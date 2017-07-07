#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <stack>
#include "Evidence.H"
#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "Rule.H"
#include "Distance.H"
#include "Potential.H"

Potential::Potential()
{
	r=gsl_rng_alloc(gsl_rng_default);
	gcnt=0;
}

Potential::~Potential()
{
	gsl_rng_free(r);
}

int 
Potential::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
Potential::setMinLeafSize(int aSize)
{
	minLeafSize=aSize;
	return 0;
}

int 
Potential::setAssocVariable(Variable* var,Potential::VariableRole vRole)
{
	varSet[var->getID()]=var;
	switch(vRole)
	{
		case Potential::FACTOR:
		{
			factorVariables[var->getID()]=0;
			break;
		}
		case Potential::MARKOV_BNKT:
		{
			markovBlnktVariables[var->getID()]=0;
			break;
		}
	}
	return 0;
}

VSET& 
Potential::getAssocVariables()
{
	return varSet;
}

int 
Potential::potZeroInit()
{
	return 0;
}


int
Potential::calculateConditionalEntropy()
{
	//We are going to use the formular \integral P(X,Mx)logP(X|Mx). 
	//However, we do not know how to estimate P(X,Mx) exactly. P(X|Mx) is given by
	//the regression tree structure. We will approximate P(X,Mx) using two strategies:
	//we will look at the joint configuration of (X,Mx) and treat it as unique and use
	//1/M, where M is the total number of data points in our dataset.
	//In this case the conditional entropy is equal to the sum of the marginal loglikelihood
	//of the child variable weighted by the n_i/M, where n_i is the number of datapoints in 
	//leaf node i. 
	//In the other strategy, we will assume that each configuration Mx which corresponds
	//to the same path in the tree is the same. That is we are partitioning our parent set
	//into P partitions where each partition is a path form the root to a leaf node.
	//We then write P(X,Mx) as P(X|Mx)P(Mx). P(Mx)=P(P_i), where P_i is the path consistent
	//with Mx. This boils down to summing over all the marginal entropies in the leaf nodes
	//weighted by n_i^2/(\sum_j _j).
	//We will begin with this first
	vector<RegressionTree*> allLeafNodes;
	//dtree->getLeafNodes(allLeafNodes);
	conditionalEntropy=0;
	jointEntropy=0;
	double totalPaths=evMgr->getNumberOfEvidences();
	for(int i=0;i<allLeafNodes.size();i++)
	{
		RegressionTree* rTree=allLeafNodes[i];
		double aMarginalEntropy=rTree->getMarginalEntropy();
		double pathCnt=(double) rTree->getDataSubset().size();
		double pathProb=pathCnt/totalPaths;
		//conditionalEntropy=aMarginalEntropy* (pathCnt*pathCnt/totalPaths);
		conditionalEntropy=aMarginalEntropy* (pathCnt/totalPaths);
		//jointEntropy=jointEntropy+(pathCnt*pathProb*log(pathProb));
		jointEntropy=jointEntropy+(pathProb*log(pathProb));
	}
	jointEntropy=jointEntropy+conditionalEntropy;
	return 0;
}


double 
Potential::getConditionalEntropy()
{
	return conditionalEntropy;
}

double 
Potential::getJointEntropy()
{
	return jointEntropy;
}

int
Potential::populateMeFromFile(const char* inPrefix,int treeCnt)
{
	dtreeSet.clear();
	for (int i=0;i<treeCnt;i++)
	{
		RegressionTree* t = new RegressionTree;
		char inname[1024];
		sprintf(inname,"%s_%d.txt",inPrefix,i);
		ifstream iFile(inname);
		t->deserialize(iFile,varSet,NULL);
		iFile.close();
		dtreeSet.push_back(t);
	}
	return 0;
}

//To turn this code into a forest, we will create multiple trees and let them work with different  datasets. In addition, instead of searching over all variables, we would randomly sub-sample a set of features (supposed to be a third) and split on that.
int
Potential::populateMe(double regval,int treeCnt)
{
	lambda=regval;
	vector<int> dataSamples;
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	for(int t=0;t<treeCnt;t++)
	{
		RegressionTree* dtree=new RegressionTree;
		dtree->setRNG(r);
		dtreeSet.push_back(dtree);
		dtree->setType(RegressionTree::LEAF);
		dtree->setCodingLength(1.0);
		int classVarID=factorVariables.begin()->first;
		dtree->setOutputVariable(classVarID);
		//Each tree now needs to maintain its own cache and leafset
		generateSamplesForTree_Bootstrap(dataSamples);			
		//generateSamplesForTree_Stability(dataSamples);			
		for(int i=0;i<dataSamples.size();i++)
		{
			dtree->setDataID(dataSamples[i]);
		}
		dataSamples.clear();	
		for(INTINTMAP_ITER mbIter=markovBlnktVariables.begin();mbIter!=markovBlnktVariables.end();mbIter++)
		{
			//dtree->setSubtreeVariable(mbIter->first);
			dtree->setSubtreeVariable(mbIter->first,varSet[mbIter->first]->getName());
		}
		dtree->setEvidenceManager(evMgr);
		int randvarCnt=markovBlnktVariables.size()/3;
		if(markovBlnktVariables.size()<=3)
		{
			randvarCnt=markovBlnktVariables.size();
		}
		struct timeval begintime2;
		struct timeval endtime2;
		gettimeofday(&begintime2,NULL);
		dtree->learn(lambda,minLeafSize,randvarCnt);
		//dtree->learn_Pairwise(lambda,minLeafSize,randvarCnt);
		gettimeofday(&endtime2,NULL);
		cout << "Time to learn  tree " << t << " " << endtime2.tv_sec-begintime2.tv_sec<< " seconds and " << endtime2.tv_usec-begintime2.tv_usec << " micro secs" << endl;
	}
	gettimeofday(&endtime,NULL);
	cout << "Time to learn trees " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
	return 0;
}

/*
int 
Potential::populateMe(double regval)
{
	//cout << "HERE!" << regval << endl;
	lambda=regval;
	dtree=new RegressionTree;
	dtree->setType(RegressionTree::LEAF);
	dtree->setCodingLength(1.0);
	currentLeafNodes.push(dtree);
	for(INTINTMAP_ITER mbIter=markovBlnktVariables.begin();mbIter!=markovBlnktVariables.end();mbIter++)
	{
		dtree->setSubtreeVariable(mbIter->first);
	}
	
	int evidCnt=evMgr->getNumberOfEvidences();
	for(int i=0;i<evidCnt;i++)
	{
		dtree->setDataID(i);
	}
	estimateMarginal(dtree);


	double currentScore=dtree->getMarginalEntropy();
	int splitCnt=1;
	map<int,int> nodeLevel;
	while(!currentLeafNodes.empty())
	{

		PARTITION* p=NULL;
		INTINTMAP* ss1=NULL;
		INTINTMAP* ss2=NULL;
		int ss1s=0;
		int ss2s=0;
		//cout << "IN THE LOOP!" << endl;
		double maxGain=0;
		int testVarID=-1;
		double testValue=0;

		RegressionTree* currNode=currentLeafNodes.front();
		INTINTMAP& currSubset=currNode->getDataSubset();
		currentLeafNodes.pop();
		int pNode;
		int branch;
		currNode->getParentInfo(pNode,branch);
		//double penalty=lambda*log(nodeLevel.size()+1);
		//double penalty=lambda*log(splitCnt);
		

		//cout <<"Splitting child  node of " << pNode <<" branch " << branch << endl;
		INTINTMAP& subtreeVars=currNode->getSubtreeVariables();
		for(INTINTMAP_ITER mbIter=subtreeVars.begin();mbIter!=subtreeVars.end();mbIter++)
		{
			//Partition the training set into sets based on the values of mbIter->first
			//but because this is real, we consider two splits
			//0 corresponds to the partition using testVar<threshold and 1 corresponds to
			//the partition using testVar>threshold
			double gain=-1;
			double splitValue;
			int varId=mbIter->first;
			
			// will return -1 if split not possible? not currently implemented.
			// constraints on leaf size?? 
			int canSplit=getPartitions_Cached(mbIter->first,splitValue,currNode->getMarginalEntropy(),currSubset,gain);
			//getPartitions(mbIter->first,splitValue,currNode->getMarginalEntropy(),currSubset,gain);
			if(canSplit==0 && gain> maxGain)
			{
				maxGain=gain;
				testVarID=mbIter->first;
				testValue=splitValue;
				//cout << "UPDATE VAR:" << testVarID << "," << testValue << endl;
				p=allPartitions[testVarID];
				ss1=(*p)[0];
				ss2=(*p)[1];
				ss1s=ss1->size();
				ss2s=ss2->size();
				//cout << "HERE WE GO!" << ss1->size() << "," << ss2->size() << endl;
			}
		}
		gettimeofday(&endtime,NULL);
		cout << "Time to find split " << testVarID << " "<< endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
		if(testVarID==-1)
		{
			clearCache();
			//cout << "No good split found " << testVarID << endl;
			continue;
		}
		//cout << "TestVar: " << testVarID << ", Split: " << testValue << endl;
		
		int testNodeLevel=0;
		if(pNode!=-1)
		{
			int plevel=nodeLevel[pNode];
			testNodeLevel=plevel+1;
		}
		double penalty=2*lambda*log(1+testNodeLevel);
		//PARTITION* p=allPartitions[testVarID];
		//INTINTMAP* ss1=(*p)[0];
		//INTINTMAP* ss2=(*p)[1];
		// I think this means we stop splitting because the leaf size hit the minimum
		if( (ss1s<minLeafSize) && (ss2s<minLeafSize))
		{
			clearCache();
			cout << "TOO SMALL!" << ss1s << "," << ss2s << endl;
			continue;
		}
		nodeLevel[testVarID]=testNodeLevel;
		splitCnt++;
		Variable* var=varSet[testVarID];
		//cout <<"Gain: "<<maxGain << " for test var "<< testVarID << " " << var->getName().c_str() << " at level " << nodeLevel[testVarID] << endl;
		
		// report parent
		RegressionTree* parent=currNode->getParent();
		if (parent!= NULL)
		{
			int parentVar=currNode->getParent()->getTestVariable();
			//cout << "Parent was " << varSet[parentVar]->getName().c_str() << endl;
		}
		
		currNode->setTestVariable(testVarID);
		currNode->setPenalizedScore(maxGain-penalty);
		//Make new leaf nodes using the values of testValue
		currNode->setType(RegressionTree::NONLEAF);
		currNode->split(testValue,allPartitions[testVarID]);
		currNode->setChildParams(allMeans[testVarID],allVariances[testVarID],allMarginalEntropy[testVarID]);
		computeCodingLength(currNode);
		map<int,RegressionTree*>& newLeaves=currNode->getChildren();
		for(map<int,RegressionTree*>::iterator lIter=newLeaves.begin();lIter!=newLeaves.end();lIter++)
		{
			//if((!lIter->second->isPureNode()) && (lIter->second->getDataSubset().size()>10) && (testNodeLevel<1))
			// DC asks -- should that be >= minLeafSize or > minLeafSize? (original)
			// doesn't matter if leaf size is 1 -- always pure... but otherwise??
			if((!lIter->second->isPureNode()) && (lIter->second->getDataSubset().size()>minLeafSize))
			//if(!lIter->second->isPureNode())
			{
				currentLeafNodes.push(lIter->second);
			}
			lIter->second->setCodingLength(1.0);
		}
		clearCache();
		INTINTMAP treeVar;
		dtree->getTreeVars(treeVar);
		//cout << "HERE! size: " << treeVar.size() << endl;
	}
	return 0;
}*/

//I'll try bootstrap or stability selection
int 
Potential::generateSamplesForTree_Bootstrap(vector<int>& dataSamples)
{
	INTINTMAP dist;
	int size=evMgr->getNumberOfEvidences();
	double step=1.0/(double)size;
	
	struct timeval begintime;
	gettimeofday(&begintime,NULL);
	//gsl_rng_set(r,begintime.tv_sec);
	gsl_rng_set(r,gcnt);
	gcnt++;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		if(rind==size)
		{
			cout <<"Oops corner case of exceeding limit. Fixing to max datasample size" <<endl;
			rind=size-1;
		}
		dataSamples.push_back(rind);
		if(dist.find(rind)==dist.end())
		{
			dist[rind]=1;
		}
		else
		{
			dist[rind]=dist[rind]+1;
		}
	}
	/*for(INTINTMAP_ITER rIter=dist.begin();rIter!=dist.end();rIter++)
	{
		if(rIter->second>5)	
		{
			cout <<"Found "<< rIter->second << " instances of " << rIter->first<<endl;
		}
	}*/
	dist.clear();
	return 0;
}

//I'll try bootstrap or stability selection
int 
Potential::generateSamplesForTree_Stability(vector<int>& dataSamples)
{
	INTINTMAP dist;
	int size=evMgr->getNumberOfEvidences();
	double step=1.0/(double)size;
	
	struct timeval begintime;
	gettimeofday(&begintime,NULL);
	//gsl_rng_set(r,begintime.tv_sec);
	gsl_rng_set(r,gcnt);
	gcnt++;
	int cnt=(int)(size/2);
	cout <<"Need to get " << cnt << " samples" << endl;
	for(int i=0;i<cnt;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		if(rind==size)
		{
			cout <<"Oops corner case of exceeding limit. Fixing to max datasample size" <<endl;
			rind=size-1;
		}
		while(dist.find(rind)!=dist.end())
		{
			double rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
			if(rind==size)
			{
				rind=size-1;
			}
		}
		dataSamples.push_back(rind);
		if((dataSamples.size()%1000)==0)
		{
			cout <<".";
		}
		dist[rind]=1;
	}
	/*for(INTINTMAP_ITER rIter=dist.begin();rIter!=dist.end();rIter++)
	{
		if(rIter->second>5)	
		{
			cout <<"Found "<< rIter->second << " instances of " << rIter->first<<endl;
		}
	}*/
	dist.clear();
	return 0;

}

int
Potential::clearMe()
{
	for(int i=0;i<dtreeSet.size();i++)
	{
		RegressionTree* dtree=dtreeSet[i];
		if(dtree!=NULL)
		{
			dtree->clear();
			delete dtree;
		}
	}
	return 0;
}

double
Potential::getJointPotValueFor(INTINTMAP& configMap)
{
	cout <<"Not implemented" << endl;
	return 0;
	/*RegressionTree* currNode=dtree;
	double pval=-1;
	while((currNode!=NULL) && (currNode->getType()!=RegressionTree::LEAF))
	{
		int testVar=currNode->getTestVariable();
		int varVal=configMap[testVar];
		RegressionTree* childNode=currNode->getChildAt(varVal);
		currNode=childNode;
	}
	int classVarID=factorVariables.begin()->first;
	int classVarVal=configMap[classVarID];
	if(currNode!=NULL)
	{
		pval=currNode->getMarginalPDF(classVarVal);
	}
	return pval;*/
}

double 
Potential::getJointPotValueForConf(string& varConf)
{
	cout <<"Not implemented " << endl;
	return 0;
}


int
Potential::generateSample(INTINTMAP& jointConf,int vId, gsl_rng *r)
{
	cout <<"Not implemented" << endl;
	return 0;
	int sampleVal=-1;
	/*RegressionTree* currNode=dtree;

	while((currNode!=NULL) && (currNode->getType()!=RegressionTree::LEAF))
	{
		int testVar=currNode->getTestVariable();
		int varVal=jointConf[testVar];
		RegressionTree* childNode=currNode->getChildAt(varVal);
		currNode=childNode;
	}
	if(currNode==NULL)
	{
		return sampleVal;
	}	
	double childMean=currNode->getMean();
	double childVariance=currNode->getVariance();
	double rval=gsl_ran_gaussian(r,sqrt(childVariance));
	//sampleVal=rval+childMean;
	return sampleVal;*/
}


double 
Potential::predictSample(EMAP* evidMap, int vId)
{
	double predValSum=0;
	for(int i=0;i<dtreeSet.size();i++)
	{
		RegressionTree* currNode=dtreeSet[i];
		while(currNode->getType()!=RegressionTree::LEAF)
		{
			int testVar=currNode->getTestVariable();
			double varVal=(*evidMap)[testVar]->getEvidVal();
			double testValue=currNode->getTestValue();
			RegressionTree* childNode=NULL;
			if(varVal<=testValue)
			{
				childNode=currNode->getChildAt(0);
			}
			else
			{
				childNode=currNode->getChildAt(1);
			}
			currNode=childNode;
		}
		double predVal=currNode->getMean();
		predValSum=predValSum+predVal;
	}
	predValSum=predValSum/dtreeSet.size();
	return predValSum;
}

double
Potential::getMSE(ofstream& oFile)
{
	double overallMSE1=0;
	for(int i=0;i<dtreeSet.size();i++)
	{
		RegressionTree* dtree=dtreeSet[i];
		double mse=getMSE(dtree);
		overallMSE1=overallMSE1+mse;
	}
	oFile << "Column\tTrueValue\tPredictedValue\tsquared_err" << endl;
	vector<double> truevals;
	vector<double> predvals;
	overallMSE1=overallMSE1/dtreeSet.size();
	double err2=0;
	int sId=factorVariables.begin()->first;
	int evidCnt=evMgr->getNumberOfEvidences();
	int vIdDistance=-1;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* v=vIter->second;
		if(strcmp(v->getName().c_str(),"Distance")==0)
		{
			vIdDistance=vIter->first;
		}
	}
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		double predval=predictSample(evidMap,sId);
		double trueval=(*evidMap)[sId]->getEvidVal();
		double diff=predval-trueval;
		double e=diff*diff;
		double distance=-1;
		if(vIdDistance>=0)
		{
			distance=(*evidMap)[vIdDistance]->getEvidVal();
		}
		err2=err2+e;
		truevals.push_back(trueval);
		predvals.push_back(predval);
		const char* name=evMgr->getSampleName(i);
		if(name==NULL)
		{
			oFile << "TestCol" << i << "\t" << trueval << "\t" << predval << "\t" << e <<  "\t" << distance << endl;
		}
		else
		{
			oFile << name << "\t" << trueval << "\t" << predval << "\t" << e << "\t" << distance <<  endl;
		}
	}
	err2=err2/((double) evidCnt);
	Distance d;
	double cc=d.computeCC(truevals,predvals);
	cout << "total train MSE " << err2 << " CC " << cc<< endl;

	
	return 0;
}

double
Potential::getMSE(RegressionTree* dtree)
{
	double mse=0;
	queue<RegressionTree*> nodes;
	nodes.push(dtree);
	while(!nodes.empty())
	{
		RegressionTree* cnode=nodes.front();
		nodes.pop();
		map<int,RegressionTree*>& children=cnode->getChildren();
		if(children.size()==0)
		{
			double variance=cnode->getVariance();
			variance=variance*(cnode->getDataSubset().size()-1);
			mse=mse+variance;
		}
		else
		{
			for(map<int,RegressionTree*>::iterator rIter=children.begin();rIter!=children.end();rIter++)
			{
				nodes.push(rIter->second);
			}
		}
	}

	int evidCnt=evMgr->getNumberOfEvidences();
	mse=mse/((double)evidCnt);
	return mse;
}

/*
* Assess MSE on an external test data set...
* assumes same variable structure (row names) as training data
* Writes to an output file
* DC Added
*/
double
Potential::getTestMSE(EvidenceManager* testEvMgr, ostream& oFile)
{	
	int evidCnt=testEvMgr->getNumberOfEvidences();
	double mse=0;
	vector<double> truevals;
	vector<double> predvals;
	int sId=factorVariables.begin()->first;
	int vIdDistance=-1;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		Variable* v=vIter->second;
		if(strcmp(v->getName().c_str(),"Distance")==0)
		{
			vIdDistance=vIter->first;
		}
	}
	oFile << "Column\tTrueValue\tPredictedValue\tsquared_err" << endl;
	
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=testEvMgr->getEvidenceAt(i);
		double predval=predictSample(evidMap,sId);
		double trueval=(*evidMap)[sId]->getEvidVal();	
		truevals.push_back(trueval);
		predvals.push_back(predval);
		double diff=predval-trueval;
		double e=diff*diff;
		const char* name=testEvMgr->getSampleName(i);

		double distance=-1;
		if(vIdDistance>=0)
		{
			distance=(*evidMap)[vIdDistance]->getEvidVal();
		}
		if(name==NULL)
		{
			oFile << "TestCol" << i << "\t" << trueval << "\t" << predval << "\t" << e << "\t" << distance <<  endl;
		}
		else
		{
			oFile << name << "\t" << trueval << "\t" << predval << "\t" << e << "\t" << distance << endl;
		}
		mse=mse+e;
	}
	mse=mse/((double) evidCnt);
	//cout << "Total test MSE " << mse << endl;
	Distance d;
	double cc=d.computeCC(truevals,predvals);
	cout <<"Total test MSE " << mse << " CC " << cc << endl;
	truevals.clear();
	predvals.clear();
	return mse;
}

int
Potential::dumpPotential(ostream& oFile)
{
	for(int i=0;i<dtreeSet.size();i++)
	{
		oFile <<"Tree"<<i<<endl;
	
		queue<RegressionTree*> nodeList;
		RegressionTree* dtree=dtreeSet[i];
		nodeList.push(dtree);
		while(!nodeList.empty())
		{
			RegressionTree* currNode=nodeList.front();
			nodeList.pop();
			map<int,RegressionTree*>& children=currNode->getChildren();
			for(map<int,RegressionTree*>::iterator dIter=children.begin();dIter!=children.end();dIter++)
			{
				oFile << currNode->getTestVariable()<< "\t"<<dIter->first << "\t" << dIter->second->getTestVariable() << endl;		
				nodeList.push(dIter->second);
			}
		}
	}
}

/**
* Dumps out gene expression to a tab-delim file.
* Leaves are separated by dummy columns.
* Doesn't have appropriate names for target variables yet (row names)
*/
int
Potential::dumpExpressionTab(RegressionTree* theTree, int sId, vector<string>* columnNames, ostream& oFile)
{
	vector<RegressionTree*> leaves;
	theTree->getLeafNodes(leaves);
	
	// order of data values
	vector<double> dataVect;
	
	// order of evidence IDs
	vector<int> evidVect;	
	
	
	// order of evidence IDs
	vector<int> parentVars;	
	
	int dummyCt=0;
	for(int i=0; i<leaves.size(); i++)
	{
		RegressionTree* currNode=leaves[i];
		
		// get info about parent
		RegressionTree* parent=currNode->getParent();
		
		// get data
		//INTINTMAP& dataSubset=currNode->getDataSubset();
		INTVECT& dataSubset=currNode->getDataSubset();
		
		//for(INTINTMAP_ITER dIter=dataSubset.begin();dIter!=dataSubset.end();dIter++)
		for(int d=0;d<dataSubset.size();d++)
		{
			// if only one leaf, then no parent.
			if (parent==NULL)
			{
				parentVars.push_back(-1);
			}
			else
			{
				parentVars.push_back(parent->getTestVariable());
			}
			evidVect.push_back(dataSubset[d]);
		}
		// dummy value in between leaves
		if (i<leaves.size()-1)
		{
			evidVect.push_back(-1);
			parentVars.push_back(-1);
		}
	}
	
	// print out headers
	oFile << "Gene";
	dummyCt=0;
	for(int i=0; i<evidVect.size(); i++)
	{
		int eId=evidVect[i];
		if (eId >= 0)
		{
			if (columnNames==NULL || columnNames->size() ==0)
			{
				oFile << "\tExpCol" << eId;
			} else
			{
				oFile << "\t" << (*columnNames)[eId];
			}
			
		}
		else 
		{
			oFile << "\tDummy" << dummyCt++;
		}
	}
	oFile << endl;
	
	// print out parents
	/*
	oFile << "ParentName";
	dummyCt=0;
	for(int i=0; i<parentVars.size(); i++)
	{
		int parId=parentVars[i];
		if (parId >= 0)
		{
			oFile << "\t" << varSet[parId]->getName();
		}
		else 
		{
			oFile << "\tDummy" << dummyCt++;
		}
	}
	oFile << endl;*/
	
	// go through -- add dummy value in between leaves
	oFile << "Target" << sId;
	for(int i=0; i<evidVect.size(); i++)
	{
		int eId=evidVect[i];
		
		// dummy value for space between leaves
		double dataval=-100;
		if (eId>=0)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(eId);
			dataval=(*evidMap)[sId]->getEvidVal();
			
			// turns out all of them in this factor have the same variable name
			// b/c one variable per factor
			//int varId=(*evidMap)[sId]->getAssocVariable();
			//Variable* var=varSet[varId];
			//cout << "\t" << varId << "\t" << var->getName();
		}
		oFile << "\t" << dataval;
	}
	oFile << endl;
	return 0;
}

/**
* Dumps out gene expression to a 3-column file for hmawk.
* Leaves are separated by vertical spacers.
*/
int
Potential::dumpExpressionHmawk(RegressionTree* theTree, int sId, vector<string>* columnNames, ostream& oFile)
{
	vector<RegressionTree*> leaves;
	theTree->getLeafNodes(leaves);
	
	// order of data values
	vector<double> dataVect;
	
	// order of evidence IDs
	vector<int> evidVect;	
	
	// order of parent IDs
	vector<int> parentVars;
	
	for(int i=0; i<leaves.size(); i++)
	{
		RegressionTree* currNode=leaves[i];
		
		// get info about parent
		RegressionTree* parent=currNode->getParent();
		
		//INTINTMAP& dataSubset=currNode->getDataSubset();
		INTVECT& dataSubset=currNode->getDataSubset();
		//oFile << dataSubset.size() << endl;
		//for(INTINTMAP_ITER dIter=dataSubset.begin();dIter!=dataSubset.end();dIter++)
		for(int d=0;d<dataSubset.size();d++)
		{
			parentVars.push_back(parent->getTestVariable());
			evidVect.push_back(dataSubset[d]);
		}
		// dummy value in between leaves
		if (i<leaves.size()-1)
		{
			evidVect.push_back(-1);
			parentVars.push_back(-1);
		}
	}
	
	// go through -- add vertical spacer in between leaves
	//cout << sId << " SID" << endl;
	int dummyCt=0;
	for(int i=0; i<evidVect.size(); i++)
	{
		int eId=evidVect[i];
		
		// dummy value for space between leaves
		double dataval=-100;
		if (eId>=0)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(eId);
			dataval=(*evidMap)[sId]->getEvidVal();
			
			// turns out all of them in this factor have the same variable name
			// b/c one variable per factor
			int varId=(*evidMap)[sId]->getAssocVariable();
			Variable* var=varSet[varId];
			//cout << "\t" << varId << "\t" << var->getName();

			if (columnNames==NULL || columnNames->size() == 0)
			{
				oFile << var->getName() << "||9 ExpCol" << eId << "||9 " << dataval << "|1" << endl;
			} else
			{
				oFile << var->getName() << "||9 " << (*columnNames)[eId] << "||9 " << dataval << "|1" << endl;
			}
			
			
		}
		else 
		{
			oFile << "|- VertSpacer" << dummyCt++ << "|Space|6" << endl;
		}	
	}
	return 0;
}

/**
* Dumps out node features for Cytoscape.
* Currently: just expression
*/
int
Potential::dumpNodeInfo(RegressionTree* theTree, int sId, vector<string>* columnNames, ostream& oFile)
{
	// check to see if names exist
	bool hasNames=(columnNames != NULL && columnNames->size()>0);
	
	vector<RegressionTree*> leaves;
	theTree->getLeafNodes(leaves);
	
	oFile << "Gene\tExpression" << endl; // header
	// for each leaf
	for(int i=0; i<leaves.size(); i++)
	{
		RegressionTree* currNode=leaves[i];
		
		//INTINTMAP& dataSubset=currNode->getDataSubset();
		INTVECT& dataSubset=currNode->getDataSubset();
		// for each datapoint
		//for(INTINTMAP_ITER dIter=dataSubset.begin();dIter!=dataSubset.end();dIter++)
		for(int d=0;d<dataSubset.size();d++)
		{
			int eId=dataSubset[d];
			if (eId>=0)
			{
				EMAP* evidMap=evMgr->getEvidenceAt(eId);
				double dataval=(*evidMap)[sId]->getEvidVal();
				// verify names provided; otherwise just print number of column
				if (hasNames)
				{
					oFile << (*columnNames)[eId] <<"\t" << dataval << endl;
				} else
				{
					oFile << "ExpCol" << eId << "\t" << dataval << endl;
				}
			}
		}
	}
	return 0;
}


int 
Potential::makeValidJPD()
{
	cout <<"Not implemented yet" << endl;
	return 0;
}

vector<RegressionTree*>&
Potential::getForest()
{
	return dtreeSet;
}

/*// Probably okay with duplicate feature values now? 
int 
Potential::getPartitions(int vId, double& splitValue,double marginalEntropy, INTINTMAP& dataSet, double& infoGain)
{
	//Here we want to partition dataSet into smaller partitions corresponding to the different values of
	//vId. We also want to compute the entropy of the class variable using the partitions. So we store
	//for each value of variable vId, the distribution of the class variable.
	int classVarID=factorVariables.begin()->first;
	vector<int> sortedInd;
	vector<double> sortedValues;
	for(INTINTMAP_ITER dIter=dataSet.begin();dIter!=dataSet.end();dIter++)
	{
		EMAP* evSet=evMgr->getEvidenceAt(dIter->first);
		if(evSet->find(vId)==evSet->end())
		{
			cout <<"No evidence value for " << vId << endl;
			return -1;
		}
		Evidence* evid=(*evSet)[vId];
		double attrVal=evid->getEvidVal();
		sortedInd.push_back(dIter->first);
		sortedValues.push_back(attrVal);
	}
	sortAttrVals(sortedValues,sortedInd);
	int partId=1; // should we start this at minLeafSize?
	//partId=minLeafSize+1; // DC ADD
	double maxGain=0;
	int splitId=0;
	double bestEntropyLeft=0;
	double bestMeanLeft=0;
	double bestVarianceLeft=0;
	double bestEntropyRight=0;
	double bestMeanRight=0;
	double bestVarianceRight=0;
	while(partId<sortedValues.size()-1)
	{
		double attrVal=sortedValues[partId];
		// this part is crucial to handling duplicate values
		int leftId=partId+1;
		while(leftId<sortedValues.size() && sortedValues[leftId]==attrVal)
		{
			leftId++;
		}

		partId=leftId-1;
		leftId--;
		
		//Now all datapoints from 0 to partID+leftID are less than or equal to attrVal
		double classEntrLeft=0;
		double meanLeft=0;
		double varianceLeft=0;
		getSubEntropy(sortedInd,0,partId,meanLeft,varianceLeft,classEntrLeft);
		double classEntrRight=0;
		double meanRight=0;
		double varianceRight=0;
		getSubEntropy(sortedInd,partId+1,sortedValues.size()-1,meanRight,varianceRight,classEntrRight);
		//Now compute the new gain
		double weightedEntropy=0;
		double total1=((double)(partId+1));
		double frac1=total1/((double)dataSet.size());
		double frac2=1-frac1;
		weightedEntropy=(frac1*classEntrLeft)+ (frac2*classEntrRight);
		double currInfoGain=marginalEntropy-weightedEntropy;
		
		if(currInfoGain>maxGain)
		{
			maxGain=currInfoGain;
			splitId=partId;
			splitValue=sortedValues[partId];
			bestEntropyLeft=classEntrLeft;
			bestMeanLeft=meanLeft;
			bestVarianceLeft=varianceLeft;
			bestEntropyRight=classEntrRight;
			bestMeanRight=meanRight;
			bestVarianceRight=varianceRight;
		}
		partId=partId+1;
	}
	infoGain=maxGain;
	PARTITION* partition=new PARTITION;
	INTDBLMAP* mean=new INTDBLMAP; 
	INTDBLMAP* variance=new INTDBLMAP; 
	INTDBLMAP* entropy=new INTDBLMAP;

	for(int i=0;i<sortedInd.size();i++)
	{
		int pseudoVal=-1;
		if(i<=splitId)
		{
			pseudoVal=0;
		}
		else
		{
			pseudoVal=1;
		}
		INTINTMAP* apart=NULL;
		if(partition->find(pseudoVal)==partition->end())
		{
			apart=new INTINTMAP;
			(*partition)[pseudoVal]=apart;
		}
		else
		{
			apart=(*partition)[pseudoVal];
		}
		(*apart)[sortedInd[i]]=0;
	}
	(*mean)[0]=bestMeanLeft;
	(*mean)[1]=bestMeanRight;
	(*variance)[0]=bestVarianceLeft;
	(*variance)[1]=bestVarianceRight;
	(*entropy)[0]=bestEntropyLeft;
	(*entropy)[1]=bestEntropyRight;
	allPartitions[vId]=partition;
	allMeans[vId]=mean;
	allVariances[vId]=variance;
	allMarginalEntropy[vId]=entropy;
	return 0;

}*/





/*
* This can be improved -- use quicksort instead of bubble sort
*/
/*int
Potential::sortAttrVals(vector<double>& sortedValues,vector<int>& sortedInd)
{
	for(int i=0;i<sortedValues.size();i++)
	{
		for(int j=i+1;j<sortedValues.size();j++)
		{
			double aVal=sortedValues[i];
			double bVal=sortedValues[j];
			if(aVal<bVal)
			{
				continue;
			}
			sortedValues[i]=bVal;
			sortedValues[j]=aVal;
			int tempRank=sortedInd[i];
			sortedInd[i]=sortedInd[j];
			sortedInd[j]=tempRank;

		}
	}
	//check for sort
	return 0;
}*/



int
Potential::showTree()
{
	for(int i=0;i<dtreeSet.size();i++)
	{
		showTree(dtreeSet[i]);
	}
	return 0;
}

int
Potential::showTree(RegressionTree* dtree)
{
	stack<RegressionTree*> nodeList;
	nodeList.push(dtree);
	int currLevel=0;
	int nodeID=0;
	int factorVariable=factorVariables.begin()->first;
	while(!nodeList.empty())
	{
		RegressionTree* currNode=nodeList.top();
		nodeList.pop();
		map<int,RegressionTree*>& children=currNode->getChildren();
		int vId=currNode->getTestVariable();
		if(vId==-1)
		{
			vId=factorVariable;
		}
		double testValue=currNode->getTestValue();
		VSET_ITER vIter=varSet.find(vId);
		if(vIter==varSet.end())
		{
			cout <<"No variable with id " << vId << endl;
			return -1;
		}
		Variable* v=vIter->second;
		if(vId!=factorVariable)
		{
			cout << v->getName().c_str() <<" <= " << testValue << endl;
		}
		else
		{
			cout << v->getName().c_str() << endl;
		}
		for(map<int,RegressionTree*>::iterator dIter=children.begin();dIter!=children.end();dIter++)
		{
			RegressionTree* cNode=dIter->second;
			nodeList.push(cNode);
		}
	}
	return 0;
}

/**
* Dumps out tree to Cytoscape-friendly network file.
* Each split in the tree becomes 2 edges. We annotate each edge with split value.
* Leaves will have many edges off them (one per target gene).
*
* Edge was:   Target GATA3 <=0 CEBPB
*		
* Edge now: nodeA	etype	nodeB	splitInfo	Target
*			GATA3_n 	split	CEBPB_n	<=0			Target
*			blah-left 	terminal	childName	""			Target
*
* _n uses to differentiate multiple instances of same regulator in the tree
*/
int
Potential::dumpTreeToNetwork(RegressionTree* currNode, ostream& oFile, const char* targetName, map<int,int> &timesSeen, int sId, vector<string>* columnNames)
{
	// if no parents, print out the header row
	if (currNode->getParent() == NULL)
	{
		oFile << "source\tinteraction\ttarget\tsplitCondition\ttreeName\tnumDataBelow" << endl;
	}
	
	int testVarID=currNode->getTestVariable();
	double testValue=currNode->getTestValue();
	
	if(testVarID==-1)
	{
		return 0;
	}

	// get the name/instance of the test variable -- applies to both left and right
	char testVarName[1024];
	sprintf(testVarName, "%s_%d", varSet[testVarID]->getName().c_str(), timesSeen[testVarID]);
	
	string splitCond;
	string leafType;
	map<int,RegressionTree*>& children=currNode->getChildren();
	for(map<int,RegressionTree*>::iterator cIter=children.begin();cIter!=children.end();cIter++)
	{
		RegressionTree* childNode=cIter->second;
		int dataBelow;
		
		// is potential leaf node left or right child?
		if(cIter->first==0)	
		{
			splitCond="<="; // first child, left
			leafType="-left"; // first child
		} else
		{
			splitCond=">";	// second child, right
			leafType="-right"; // second child
		}

		// tab-delim fields: parent "split" child splitCondition targetName
		
		
		if(childNode->getTestVariable()!=-1)
		{
			// parent "split" child
			// parent "split"
			oFile << testVarName << "\tsplit";
		
			char childName[1024];
			int cv=childNode->getTestVariable();
			
			// increment child HERE
			timesSeen[cv]++;
			
			sprintf(childName, "%s_%d", varSet[cv]->getName().c_str(), timesSeen[cv]);
			oFile  << "\t" << childName;
			
			// weird if child name is same as parent
			if (strcmp(testVarName, childName)==0)
			{
				cout << "Splitting on same variable a second time" << endl;
			}
			
			// how many nodes below the child? 
			dataBelow=childNode->getDataSubset().size();
			// parent "split" child [<=/<]testValue targetName
			oFile << "\t" << splitCond << testValue << "\t" << targetName << "\t" << dataBelow << endl;
		
		}
		else
		{
			// parent "split" parent-leafType
			// parent "split"
			oFile << testVarName << "\tsplit";
			
			// how many nodes in this leaf?
			dataBelow=childNode->getDataSubset().size();
			oFile << "\t" << testVarName << leafType;
			// parent "split" child [<=/<]testValue targetName
			oFile << "\t" << splitCond << testValue << "\t" << targetName << "\t" << dataBelow << endl;
		
			////SR is commenting this out to avoid excessive I/O
			// order of evidence IDs
			/*vector<int> evidVect;	
			//getDataAtLeaf(childNode, evidVect);
			//INTINTMAP& dataSubset=childNode->getDataSubset();
			INTVECT& dataSubset=childNode->getDataSubset();
			//for(INTINTMAP_ITER dIter=dataSubset.begin();dIter!=dataSubset.end();dIter++)
			for(int d=0;d<dataSubset.size();d++)
			{
				evidVect.push_back(dataSubset[d]);
			}
			
			// for each child...
			// print out headers
			for(int i=0; i<evidVect.size(); i++)
			{
				char childName[1024];
				double dataval=-100;
				
				int eId=evidVect[i];
				if (eId >= 0)
				{
					if (columnNames->size() ==0)
					{
						sprintf(childName, "ExprCol_%d", eId);
					} else
					{
						sprintf(childName, "%s", (*columnNames)[eId].c_str());
					}
					EMAP* evidMap=evMgr->getEvidenceAt(eId);
					dataval=(*evidMap)[sId]->getEvidVal();	
					
					// print edge with gene value
					oFile << testVarName << leafType << "\tterminal\t" << childName << "\tterminal\t" << targetName << "\t" << 1 << endl;
				}
				else
				{
					cout << "problem in potential:: dump network ??????" << endl;
				}
			}*/			
		}
		
		
		// keep going down the tree...
		if(childNode->getTestVariable()!=-1)
		{
			dumpTreeToNetwork(childNode, oFile,targetName, timesSeen, sId, columnNames);
		}
	}
		
	return 0;
}


/*
int
Potential::getDataAtLeaf(RegressionTree* currNode, vector<int>* evidVect)
{
	INTINTMAP& dataSubset=currNode->getDataSubset();
	for(INTINTMAP_ITER dIter=dataSubset.begin();dIter!=dataSubset.end();dIter++)
	{
		evidVect->push_back(dIter->first);
	}
}*/

int
Potential::prune()
{
	cout <<"Not implemented" <<endl;

/*	dtree->generateRuleSet();
	vector<Rule*>& ruleSet=dtree->getRuleSet();

	int evidCnt=evMgr->getNumberOfEvidences();
	for(int i=0;i<evidCnt;i++)
	{
		pruneDataSet[i]=0;
	}

	for(int i=0;i<ruleSet.size();i++)
	{
		if(pruneDataSet.size()==0)
		{
			continue;
		}
		Rule* arule=ruleSet[i];
		if(arule==NULL)
		{
			continue;
		}
		//Keep trying to delete conditions from this rule until the entropy does not
		//decrease
		map<int,RegressionTree*>& conditionSet=arule->getAllConditions();
		bool foundDelCondition=true;
		while(foundDelCondition)
		{
			//double currEntropy=arule->getMarginalEntropy();
			//int currCoverage=arule->getCoverage();
			int currCoverage=0;
			double currEntropy=computeEntropyIfDeleted(arule,-1,currCoverage);
			cout <<"Pruning rule " << i << " with coverage: " << currCoverage << endl;
			arule->showRule();
			int oldCoverage=currCoverage;
			double oldComplexity=arule->getRuleComplexity(-1);
			double currComplexity=oldComplexity;
			if(currCoverage==0)
			{
				foundDelCondition=false;
				delete arule;
				ruleSet[i]=NULL;
				arule=NULL;
				continue;
			}
			int toDel=-1;
			for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
			{
				if(lIter->second->getChildren().size()==0)
				{
					continue;
				}
				int newCoverage;
				double newEntropy=computeEntropyIfDeleted(arule,lIter->first,newCoverage);
				if(newCoverage<oldCoverage)
				{
					cout <<"Wierd rule. Deleting condition does not increase coverage!"<<endl;
				}
				//We are guaranteed that ruleLen>1
				double ruleLen=(double)conditionSet.size();
				double newComplexity=arule->getRuleComplexity(lIter->first);
				if(((newCoverage*newEntropy)+newComplexity)<=((currCoverage*currEntropy)+currComplexity))
				{
					currEntropy=newEntropy;
					currCoverage=newCoverage;
					toDel=lIter->first;
					currComplexity=newComplexity;
				}
			}
			if(toDel==-1)
			{
				foundDelCondition=false;
				continue;
			}
			map<int,RegressionTree*>::iterator delIter=conditionSet.find(toDel);
			conditionSet.erase(delIter);
			//If the condition set has only one node, it means it is the leaf node and can be deleted.
			if(conditionSet.size()==1)
			{
				foundDelCondition=false;
				delete arule;
				ruleSet[i]=NULL;
				arule=NULL;
				continue;
			}
			arule->setMarginalEntropy(currEntropy);
			arule->setCoverage(currCoverage);
		}
		if(arule!=NULL)
		{
			prunedRuleSet.push_back(arule);
			reducePruneSet(arule);
			//If I am using the Nguyen approach then I need to delete the datapoints
			//already covered by this rule
		}
	}
	cout <<"Number of rules after pruning" << prunedRuleSet.size() << endl;
	removeDuplicateCoverage();*/
	return 0;
}

int
Potential::getAssocVariables_PostPruning(INTINTMAP& newVarSet)
{
	for(int r=0;r<prunedRuleSet.size();r++)
	{
		Rule* arule=prunedRuleSet[r];
		if(arule==NULL)
		{
			continue;
		}
		map<int,RegressionTree*>& conditionSet=arule->getAllConditions();
		for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
		{
			RegressionTree* node=lIter->second;
			if(node->getChildren().size()==0)
			{
				continue;
			}
			int testVarID=node->getTestVariable();
			newVarSet[testVarID]=0;
		}
	}
	return 0;
}


double
Potential::computeEntropyIfDeleted(Rule* arule,int delMe,int& acoverage)
{
	cout <<"Not implemented" <<endl;
	/*vector<double> filteredVals;
	int evidCnt=evMgr->getNumberOfEvidences();
	double mean=0;
	double var=0;
	int classVarID=factorVariables.begin()->first;
	map<int,RegressionTree*>& conditionSet=arule->getAllConditions();

	//for(int e=0;e<evidCnt;e++)
	for(INTINTMAP_ITER aIter=pruneDataSet.begin();aIter!=pruneDataSet.end();aIter++)
	{
		int e=aIter->first;
		EMAP* evMap=evMgr->getEvidenceAt(e);
		double classVarVal=(*evMap)[classVarID]->getEvidVal();
		bool testTrue=true;
		for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
		{
			RegressionTree* rnode=lIter->second;
			if(lIter->first==delMe)
			{
				continue;
			}
			if(rnode->getChildren().size()==0)
			{
				continue;
			}
			int parentBranch=arule->getBranch(lIter->first);
			if(parentBranch==-1)
			{
				cout <<"Invalid branch value " << endl;
				exit(0);
			}
			int testVar=rnode->getTestVariable();
			double testValue=rnode->getTestValue();
			double varVal=(*evMap)[testVar]->getEvidVal();
			if(parentBranch==0)
			{
				if(varVal>testValue)
				{
					testTrue=false;
				}
			}
			else if(parentBranch==1)
			{
				if(varVal<=testValue)
				{
					testTrue=false;
				}
			}
			if(!testTrue)
			{
				break;
			}
		}
		if(testTrue)
		{
			filteredVals.push_back(classVarVal);
			mean=mean+classVarVal;
		}
	}
	//cout <<"Datapoints covered " << filteredVals.size() << endl;
	mean=mean/((double)filteredVals.size());
	for(int i=0;i<filteredVals.size();i++)
	{
		double diff=mean-filteredVals[i];
		var=var+(diff*diff);
	}
	var=var/((double)filteredVals.size()-1);
	double newEntropy=-1;
	RegressionTree::getSubEntropy(mean,var,newEntropy);
	acoverage=filteredVals.size();
	filteredVals.clear();
	return newEntropy;*/
}

int
Potential::reducePruneSet(Rule* r)
{
	INTINTMAP coveredSet;
	map<int,RegressionTree*>& conditionSet=r->getAllConditions();
	for(INTINTMAP_ITER aIter=pruneDataSet.begin();aIter!=pruneDataSet.end();aIter++)
	{
		EMAP* evMap=evMgr->getEvidenceAt(aIter->first);
		bool testTrue=true;
		for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
		{
			RegressionTree* rnode=lIter->second;
			if(rnode->getChildren().size()==0)
			{
				continue;
			}
			int parentBranch=r->getBranch(lIter->first);
			if(parentBranch==-1)
			{
				cout <<"Invalid branch value " << endl;
				exit(0);
			}
			int testVar=rnode->getTestVariable();
			double testValue=rnode->getTestValue();
			double varVal=(*evMap)[testVar]->getEvidVal();
			if(parentBranch==0)
			{
				if(varVal>testValue)
				{
					testTrue=false;
				}
			}
			else if(parentBranch==1)
			{
				if(varVal<=testValue)
				{
					testTrue=false;
				}
			}
			if(!testTrue)
			{
				break;
			}
		}
		if(testTrue)
		{
			coveredSet[aIter->first]=0;
			INTINTMAP* coveredByRules=NULL;
			if(dataRuleCoverage.find(aIter->first)==dataRuleCoverage.end())
			{
				coveredByRules=new INTINTMAP;
				dataRuleCoverage[aIter->first]=coveredByRules;
			}
			else
			{
				coveredByRules=dataRuleCoverage[aIter->first];
			}
			(*coveredByRules)[prunedRuleSet.size()-1]=0;
		}
	}

	for(INTINTMAP_ITER dIter=coveredSet.begin();dIter!=coveredSet.end();dIter++)
	{
		INTINTMAP_ITER eIter=pruneDataSet.find(dIter->first);
		if(eIter==pruneDataSet.end())
		{
			cout <<"No datapoint found " << dIter->first <<endl;
			exit(0);
		}
		//pruneDataSet.erase(eIter);
	}
	coveredSet.clear();
	return 0;
}




int 
Potential::removeDuplicateCoverage()
{
	int overlappedData=0;
	for(map<int,INTINTMAP*>::iterator dIter=dataRuleCoverage.begin();dIter!=dataRuleCoverage.end();dIter++)
	{
		INTINTMAP* covRuleSet=dIter->second;
		if(covRuleSet->size()==1)
		{
			continue;
		}
		overlappedData++;
		//If there are multiple rules, then keep the one with the least MDL
		double minEntropy=0;
		double minComplexity=0;
		double coverage=0;
		int bestRule=-1;
		for(INTINTMAP_ITER rIter=covRuleSet->begin();rIter!=covRuleSet->end();rIter++)
		{
			Rule* rule=prunedRuleSet[rIter->first];
			if(rIter==covRuleSet->begin())
			{
				minEntropy=rule->getMarginalEntropy();
				minComplexity=rule->getRuleComplexity(-1);
				coverage=rule->getCoverage()-1;
				bestRule=rIter->first;
			}
			else
			{
				int newCoverage=rule->getCoverage()-1;
				double newComplexity=rule->getRuleComplexity(-1);
				double newEntropy=rule->getMarginalEntropy();
				if((newComplexity+(newCoverage*newEntropy))
				  <(minComplexity+(coverage*minEntropy)))
				{
					minEntropy=newEntropy;
					coverage=newCoverage;
					minComplexity=newComplexity;
					bestRule=rIter->first;
				}
			}
		}
		for(INTINTMAP_ITER rIter=covRuleSet->begin();rIter!=covRuleSet->end();rIter++)
		{
			if(rIter->first==bestRule)
			{
				continue;
			}
			Rule* arule=prunedRuleSet[rIter->first];
			int newCov=arule->getCoverage()-1;
			arule->setCoverage(newCov);
		}
		covRuleSet->clear();
		(*covRuleSet)[bestRule]=0;
	}
	cout <<"Number of datapoints overlapped "<<  overlappedData << endl;
	//Now finally get rid of any rules that do not have any coverage;
	int blankrule=0;
	for(int r=0;r<prunedRuleSet.size();r++)
	{
		if(prunedRuleSet[r]->getCoverage()>0)
		{
			continue;
		}
		blankrule++;
		delete prunedRuleSet[r];
		prunedRuleSet[r]=NULL;
	}
	cout <<"Blank rules " << blankrule << endl;
	return 0;
}


