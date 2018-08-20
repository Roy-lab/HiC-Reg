#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <string.h>
#include <time.h>
#include "Error.H"
#include "Variable.H"
#include "SlimFactor.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "VariableManager.H"

#include "LatticeStructure.H"

#include "PotentialManager.H"
#include "FactorGraph.H"

#include "Vertex.H"
#include "Graph.H"
#include "FactorManager.H"
#include "Move.H"
#include "FGEditor.H"

#include "GenanatomyConverter.H"
#include "FGMaximizer.H"

FGMaximizer::FGMaximizer()
{
	columnNames=NULL; // assume no column names by default
	testEvMgr=NULL; // assume no test data by default
}

FGMaximizer::~FGMaximizer()
{
}

//How many mbs to use to do the greedy, approximate search of MB vars 
int 
FGMaximizer::setMBCntForApproxGreedy(int aCnt)
{
	mbCntForApproxGreedy=aCnt;
	return 0;
}

//How many top graphs we should return
int 
FGMaximizer::setBeamSize(int aSize)
{
	beamSize=aSize;
	return 0;
}

int 
FGMaximizer::setFactorManager(FactorManager* aPtr)
{
	fMgr=aPtr;
	return 0;
}

int 
FGMaximizer::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int
FGMaximizer::setPotentialManager(PotentialManager* aPtr)
{
	potMgr=aPtr;
	return 0;
}

int
FGMaximizer::setTrainEvidenceManager(EvidenceManager* aPtr)
{
	trainEvMgr=aPtr;
	return 0;
}

int
FGMaximizer::setTestEvidenceManager(EvidenceManager* aPtr)
{
	testEvMgr=aPtr;
	return 0;
}

int
FGMaximizer::setTreeCnt(int cnt)
{
	treeCnt=cnt;
	return 0;
}

int
FGMaximizer::setTreeLocation(const char* tl)
{
	strcpy(treeloc,tl);
	return 0;
}

int 
FGMaximizer::setOutputDir(const char* aPath) 
{
	strcpy(outputDir,aPath);
	return 0;
}

int
FGMaximizer::setConvergenceThreshold(double t)
{
	convThreshold=t;
	return 0;
}

int
FGMaximizer::setMinLeafSize(int aSize)
{
	minLeafSize=aSize;
	return 0;
}

int
FGMaximizer::setPriorGraph(const char* aFName)
{
	priorGraph.setDirectionality(true);
	priorGraph.makeGraph(aFName);
	return 0;
}

// DC added
int
FGMaximizer::setColumnNames(vector<string>* colNames) 
{
	columnNames=colNames;
	return 0;
}

int 
FGMaximizer::findBestGraphs()
{
	FactorGraph* currFg=fMgr->createInitialFactorGraph();
	maxMBSize=fMgr->getMaxFactorSize();
	maxMBSize_Approx=fMgr->getMaxFactorSize_Approx();
	maxMBSize--;
	maxMBSize_Approx--;

	fgeditor.setFactorManager(fMgr);
	fgeditor.setMBCntForApproxGreedy(mbCntForApproxGreedy);
	fgeditor.setMBSize_Exact(maxMBSize);
	int currK=1;
	double currScore=getScore(currFg);
	cout<<"Current score: " << currScore << endl;
	while(currK<=maxMBSize_Approx)
	{
		bool notConverged=true;
		int iterCnt=0;
		while((notConverged)&&(iterCnt<10))
		{
			fgeditor.setFactorGraph(currFg);
			fgeditor.setMaxMBSize(currK);
			if(currK<=maxMBSize)
			{
				fgeditor.setMoveType(Move::EXACT);
			}
			else
			{
				generateNextLevelClusters(currK,currFg);
				fgeditor.setMoveType(Move::GREEDY);
			}
			fgeditor.generateMoves();
			vector<Move*>& moveSet=fgeditor.getMoves();
			sortMoves(moveSet);
			//Now decide what move to make. A move is the update to the MB of the srcVertex. 
			//As a result of this one move, MBs of several variables will get updated and hence 
			//some moves will get invalidated
			//So making the move would
			//change the MBs of some variables. We store these variables in our affectedVars list
			//along with the variables that must be included in its MB.
			map<int,int> srcVertex;
			map<int,int> targetVertex;
			//double currLL=fMgr->getLikelihood_ChainRule(currFg);
			//cout <<"Current likelihood_chainrule: "<< currLL<< endl;
			//double currLL_Can=fMgr->getLikelihood(currFg);
			//cout <<"Current("<< currK <<") likelihood_canonical: "<< currLL_Can<< endl;
			//double currLL_MCMC=fMgr->getLikelihood_MCMC(currFg);
			//cout <<"Current likelihood_mcmc: "<< currLL_MCMC<< endl;
			int successMoves=0;
			int failedMoves=0;
			for(int m=0;m<moveSet.size();m++)
			{
				Move* aMove=moveSet[m];
				if(aMove->getSrcVertex()==8)
				{
					//cout <<"Stop here" << endl;
				}
				SlimFactor* sFactor=currFg->getFactorAt(aMove->getSrcVertex());
				if(sFactor->mergedMB.size()==currK)
				{
					//cout << endl <<"Nothing to do for " << sFactor->fId << endl;
				}	
				else 
				{
					int moveStat=attemptMove(currFg,m,moveSet,srcVertex,targetVertex);
					if(moveStat==-1)
					{
						return -1;
					}
					else if(moveStat==0)
					{
						successMoves++;
					}
					else
					{
						failedMoves++;
					}
				}
			}
			cout <<"Total success moves " << successMoves << " failed moves " << failedMoves << endl;
			if(!currFg->isConsistent())
			{
				cout <<"Consistent check failure " << endl;
				return -1;
			}
			double newScore=getScore(currFg);
			cout<<"New score: " << newScore<< endl;
			if((currScore-newScore) <=convThreshold)
			{
				notConverged=false;		
			}
			currScore=newScore;
			fgeditor.clearOldMoves();
			iterCnt++;
		}
		currFg->dumpVarMB_PairwiseFormat(outputDir,currK,vMgr->getVariableSet(),fMgr);
		if(maxMBSize==1)
		{
			currFg->dumpCandidateVarMB_PairwiseFormat(outputDir,mbCntForApproxGreedy,vMgr->getVariableSet(),fMgr);
		}
		currK++;
	}
	fg=currFg;
	double finalScore=getScore(currFg);
	cout<<"Final score: " << finalScore << endl;
	//double finalLL=fMgr->getLikelihood_ChainRule(currFg);
	//cout <<"final likelihood_chainrule: "<< finalLL<< endl;
	//double finalLL_Can=fMgr->getLikelihood(currFg);
	//cout <<"Final likelihood_canonical: "<< finalLL_Can<< endl;
	//double currLL_MCMC=fMgr->getLikelihood_MCMC(currFg);
	//cout <<"Current likelihood_mcmc: "<< currLL_MCMC<< endl;
	return 0;
}

int
FGMaximizer::findBestGraphs_TopDown(double lambda,const char* projectFName)
{
	//The main idea here is to take each variable and learn a regression tree using all other variables.
	//Those variables that will not contribute significantly will not end up in the tree
	FactorGraph* currFg=fMgr->createInitialFactorGraph();
	double currScore=getScore(currFg);
	cout<<"Current score: " << currScore << endl;
	VSET& varSet=vMgr->getVariableSet();
	map<int,double > varErr;
	map<int,double > varInitErr;
	char regProgFName[1024];
	sprintf(regProgFName,"%s/regpro.txt",outputDir);
	ofstream oFile(regProgFName);

	map<int,RegressionTree*> regTreeModel;
	
	for(int i=0;i<currFg->getFactorCnt();i++)
	{
		//cout << "Analyzing factor " << i << endl;
		SlimFactor* sFactor=currFg->getFactorAt(i);
		Potential* aPotFunc=new Potential;
		map<string,int> regset;
		for(int j=0;j<currFg->getFactorCnt();j++)
		{
			SlimFactor* mFactor=currFg->getFactorAt(j);
			Variable* aVar=varSet[mFactor->vIds[0]];
			if(i==j)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				Vertex* aVertex=priorGraph.getVertex(aVar->getName().c_str());
				if(aVertex==NULL)
				{
					cout << "No vertex for " << aVar->getName() << endl;
					continue;
				}
				if(aVertex->getOutDegree()>0)
				{
					aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
					regset[aVar->getName()]=0;
				}
			}
		}
		if(i==0)
		{
			cout <<"Candidate regulators  " << regset.size() << endl;
		}
		aPotFunc->potZeroInit();
		aPotFunc->setMinLeafSize(minLeafSize);
		struct timeval begintime;
		struct timeval endtime;
		gettimeofday(&begintime,NULL);
	
		//potMgr->populatePotential(aPotFunc,false,lambda);
		aPotFunc->setEvidenceManager(trainEvMgr);
		gettimeofday(&endtime,NULL);
		cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
		//aPotFunc->showTree();

		aPotFunc->calculateConditionalEntropy();
		sFactor->conditionalEntropy=aPotFunc->getConditionalEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		//The Markov blanket variables are all the nonleaf nodes of the decision tree
		/////WARNING!! THIS IS A HACK FIX. THIS CODE IS NOT COMPATIBLE WITH THE POTENTIAL FUNCTION
		vector<RegressionTree*>& dtreeSet=aPotFunc->getForest();
		RegressionTree* rtree=dtreeSet[0];
		string starter;
	//	rtree->showMe(starter,varSet);
		INTINTMAP treeVar;
		//aPotFunc->prune();
		//aPotFunc->getAssocVariables_PostPruning(treeVar);
		rtree->getTreeVars(treeVar);
		Variable* sVar=varSet[sFactor->fId];
		rtree->dumpTree(oFile,varSet,sVar->getName().c_str());
		regTreeModel[sFactor->fId]=rtree;
		for(INTINTMAP_ITER vIter=treeVar.begin();vIter!=treeVar.end();vIter++)
		{
			if(vIter->first!=sFactor->vIds[0])
			{
				sFactor->mergedMB[vIter->first];
			}
		}
		if(sFactor->mergedMB.size()==0)
		{
			sFactor->mbScore=sFactor->conditionalEntropy;
		}
		else
		{
			sFactor->mbScore=sFactor->conditionalEntropy+(lambda*log(sFactor->mergedMB.size()));
		}
		//double err=aPotFunc->getMSE();
		
		
		//varErr[sFactor->fId]=err;
		varInitErr[sFactor->fId]=rtree->getVariance();
	}
	oFile.close();
	char errFName[1024];
	sprintf(errFName,"%s/knownregulator_err.txt",outputDir);
	ofstream errFile(errFName);
	for(map<int,double>::iterator vIter=varErr.begin();vIter!=varErr.end();vIter++)
	{
		SlimFactor* sFactor=currFg->getFactorAt(vIter->first);
		errFile << varSet[vIter->first]->getName() << "\t" <<varInitErr[vIter->first] << "\t" << vIter->second << "\t" << sFactor->mergedMB.size();
		for(INTINTMAP_ITER uIter=sFactor->mergedMB.begin();uIter!=sFactor->mergedMB.end();uIter++)
		{
			errFile <<"\t" << varSet[uIter->first]->getName();
		}
		errFile << endl;
	}
	errFile.close();
	currFg->dumpVarMB_PairwiseFormat(outputDir,currFg->getFactorCnt(),vMgr->getVariableSet(),fMgr);
	double finalScore=getScore(currFg);
	cout<<"Final score: " << finalScore << endl;
	
	dispRegTree_Genanatomy(regTreeModel,projectFName);
	return 0;
}


int
FGMaximizer::findBestGraphs_PriorGraph(double lambda,const char* projectFName)
{
	//The main idea here is to take each variable and learn a regression tree using all other variables.
	//Those variables that will not contribute significantly will not end up in the tree
	FactorGraph* currFg=fMgr->createInitialFactorGraph();
	double currScore=getScore(currFg);
	cout<<"Current score: " << currScore << endl;
	
	VSET& varSet=vMgr->getVariableSet();
	map<int,double > varErr;
	map<int,double > varInitErr;
	
	char regProgFName[1024];
	char netFName[1024];
	char expFName[1024];
	char hmFName[1024];
	char nodeFName[1024];
	char testerrFName[1024];
	char trainerrFName[1024];
	
	// regulatory program file output
	sprintf(regProgFName,"%s/regpro.txt",outputDir);
	sprintf(netFName,"%s/regtree_network.tab",outputDir);
	sprintf(expFName,"%s/regtree_expression.tab",outputDir);
	sprintf(hmFName,"%s/regtree_expression_hmawk.txt",outputDir);
	sprintf(nodeFName,"%s/node_attributes.txt",outputDir);
	sprintf(testerrFName,"%s/testset_error.txt",outputDir);
	sprintf(trainerrFName,"%s/trainset_error.txt",outputDir);
	
	ofstream oFile(regProgFName);
	ofstream netFile(netFName);
	ofstream expFile(expFName);
	ofstream hmFile(hmFName);
	ofstream nodeFile(nodeFName);
	
	//if (testEvMgr!=NULL)
	//{
	ofstream testerrFile(testerrFName);
	ofstream trainerrFile(trainerrFName);
	//}
	
	// header row for network file (moved to actual function)
	//netFile << "source\tinteraction\ttarget\tsplitCondition\ttreeName\tnumDataBelow" << endl;
	
	map<int,RegressionTree*> regTreeModel;
	for(int i=0;i<currFg->getFactorCnt();i++)
	{

		SlimFactor* sFactor=currFg->getFactorAt(i);
		Variable* sVar=varSet[sFactor->fId];
		Vertex* sVertex=priorGraph.getVertex(sVar->getName().c_str());
		cout << "Analyzing factor " << i << " " << sVar->getName() << endl;
		if(sVertex==NULL)
		{
			continue;
		}
		NINFO_MAP& potentialParents=sVertex->getInNeighbours();
		if(potentialParents.size()==0)
		{
			continue;
		}
		Potential* aPotFunc=new Potential;
		aPotFunc->setAssocVariable(sVar,Potential::FACTOR);
		for(NINFO_MAP_ITER nIter=potentialParents.begin();nIter!=potentialParents.end();nIter++)
		{
			int varID=vMgr->getVarID(nIter->first.c_str());
			if(varID==-1)
			{
				continue;
			}
			//cout << "\tPotential Parent " << varSet[varID]->getName() << " for " << sVar->getName() << endl; 
			SlimFactor* mFactor=currFg->getFactorAt(varID);
			if(mFactor->fId==sFactor->fId)
			{
				continue;
			}
			Variable* aVar=varSet[varID];
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		aPotFunc->potZeroInit();
		aPotFunc->setMinLeafSize(minLeafSize);

	
		struct timeval begintime;
		struct timeval endtime;
		gettimeofday(&begintime,NULL);
	
		//potMgr->populatePotential(aPotFunc,false,lambda);
		//potMgr->populatePotential(aPotFunc,false,lambda);
		aPotFunc->setEvidenceManager(trainEvMgr);
		aPotFunc->populateMe(lambda,treeCnt);
		gettimeofday(&endtime,NULL);
		cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
		//aPotFunc->showTree();

		aPotFunc->calculateConditionalEntropy();
		sFactor->conditionalEntropy=aPotFunc->getConditionalEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		
		
		//The Markov blanket variables are all the nonleaf nodes of the decision tree
		vector<RegressionTree*>& rtreeSet=aPotFunc->getForest();
		for(int t=0;t<rtreeSet.size();t++)
		{
			netFile <<"Tree " << t << endl;
			RegressionTree* rtree=rtreeSet[t];

			char regSerFName[1024];
			sprintf(regSerFName,"%s/regtree_node_%d.txt",outputDir,t);
			ofstream sFile(regSerFName);
			rtree->serialize(sFile,varSet);
			sFile.close();

			// dump expression data
			//aPotFunc->dumpExpressionTab(rtree, sFactor->fId, columnNames, expFile);
			//aPotFunc->dumpExpressionHmawk(rtree, sFactor->fId, columnNames, hmFile);
		
			string starter;
			//rtree->showMe(starter,varSet);
			INTINTMAP treeVar;
			//aPotFunc->prune();
			//aPotFunc->getAssocVariables_PostPruning(treeVar);
			rtree->getTreeVars(treeVar);
		
			// dump tree to original file format
			rtree->dumpTree(oFile,varSet,sVar->getName().c_str());
		
			// dump tree to network file
			map<int,int> timesSeen; // empty map
			//rtree->dumpTreeToNetwork(netFile, varSet,sVar->getName().c_str(), timesSeen);
			aPotFunc->dumpTreeToNetwork(rtree, netFile, sVar->getName().c_str(), timesSeen, sFactor->fId, columnNames);
		
			//SR is commenting this out
			//dump node info to Cytoscape file
			//aPotFunc->dumpNodeInfo(rtree, sFactor->fId, columnNames, nodeFile);
		
			//regTreeModel[sFactor->fId]=rtree;
			for(INTINTMAP_ITER vIter=treeVar.begin();vIter!=treeVar.end();vIter++)
			{
				if(vIter->first!=sFactor->vIds[0])
				{
					sFactor->mergedMB[vIter->first];
				}
			}
			if(sFactor->mergedMB.size()==0)
			{
				sFactor->mbScore=sFactor->conditionalEntropy;
			}
			else
			{
				sFactor->mbScore=sFactor->conditionalEntropy+(lambda*log(sFactor->mergedMB.size()));
			}
		}
		double err=aPotFunc->getMSE(trainerrFile);
		
		// test data set
		if (testEvMgr!= NULL)
		{
			double testErr=aPotFunc->getTestMSE(testEvMgr, testerrFile);
			cout << "Test Set MSE " << testErr << endl;
		}
		
		varErr[sFactor->fId]=err;
		//varInitErr[sFactor->fId]=rtree->getVariance();
		
	}
	oFile.close(); // close tree program file
	netFile.close(); // close network file
	expFile.close(); // close expression tab file
	hmFile.close(); // close expression hmawk file
	nodeFile.close(); // close node data file
	
	//if (testEvMgr!= NULL)
	//{
	testerrFile.close(); // close error file
	//}
	trainerrFile.close();
	// Known regulator error file
	char errFName[1024];
	sprintf(errFName,"%s/knownregulator_err.txt",outputDir);
	ofstream errFile(errFName);
	for(map<int,double>::iterator vIter=varErr.begin();vIter!=varErr.end();vIter++)
	{
		SlimFactor* sFactor=currFg->getFactorAt(vIter->first);
		errFile << varSet[vIter->first]->getName() << "\t" <<varInitErr[vIter->first] << "\t" << vIter->second << "\t" << sFactor->mergedMB.size();
		for(INTINTMAP_ITER uIter=sFactor->mergedMB.begin();uIter!=sFactor->mergedMB.end();uIter++)
		{
			errFile <<"\t" << varSet[uIter->first]->getName();
		}
		errFile << endl;
	}
	errFile.close();

	currFg->dumpVarMB_PairwiseFormat(outputDir,currFg->getFactorCnt(),vMgr->getVariableSet(),fMgr);
	double finalScore=getScore(currFg);
	cout<<"Final score: " << finalScore << endl;
	
	// seg fault in here
	//dispRegTree_Genanatomy(regTreeModel,projectFName);
	
	return 0;
}

int
FGMaximizer::readBestGraphs_PriorGraph(double lambda,const char* projectFName)
{
	//The main idea here is to take each variable and learn a regression tree using all other variables.
	//Those variables that will not contribute significantly will not end up in the tree
	FactorGraph* currFg=fMgr->createInitialFactorGraph();
	double currScore=getScore(currFg);
	cout<<"Current score: " << currScore << endl;
	
	VSET& varSet=vMgr->getVariableSet();
	map<int,double > varErr;
	map<int,double > varInitErr;
	
	
	map<int,RegressionTree*> regTreeModel;
	for(int i=0;i<currFg->getFactorCnt();i++)
	{
		SlimFactor* sFactor=currFg->getFactorAt(i);
		Variable* sVar=varSet[sFactor->fId];
		Vertex* sVertex=priorGraph.getVertex(sVar->getName().c_str());
		cout << "Analyzing factor " << i << " " << sVar->getName() << endl;
		if(sVertex==NULL)
		{
			continue;
		}
		NINFO_MAP& potentialParents=sVertex->getInNeighbours();
		if(potentialParents.size()==0)
		{
			continue;
		}
		Potential* aPotFunc=new Potential;
		
		aPotFunc->setAssocVariable(sVar,Potential::FACTOR);
		for(NINFO_MAP_ITER nIter=potentialParents.begin();nIter!=potentialParents.end();nIter++)
		{
			int varID=vMgr->getVarID(nIter->first.c_str());
			if(varID==-1)
			{
				continue;
			}
			//cout << "\tPotential Parent " << varSet[varID]->getName() << " for " << sVar->getName() << endl; 
			SlimFactor* mFactor=currFg->getFactorAt(varID);
			if(mFactor->fId==sFactor->fId)
			{
				continue;
			}
			Variable* aVar=varSet[varID];
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		aPotFunc->potZeroInit();
		aPotFunc->setMinLeafSize(minLeafSize);
		//aPotFunc->setEvidenceManager(trainEvMgr);

		//aPotFunc->populateMe(lambda,treeCnt);
		aPotFunc->populateMeFromFile(treeloc,treeCnt);
		//aPotFunc->showTree();

		//Do we need this?
		//aPotFunc->calculateConditionalEntropy();
		//sFactor->conditionalEntropy=aPotFunc->getConditionalEntropy();
		//sFactor->jointEntropy=aPotFunc->getJointEntropy();
		
		
		// Report testing error per example
		// test data set
		if (testEvMgr!= NULL)
		{
			char testerrFName[1024];
			sprintf(testerrFName,"%s/testset_error.txt",outputDir);
			ofstream testerrFile(testerrFName);
			// switch on reg mode
			// Testing error
			double testErr=aPotFunc->getTestMSE(testEvMgr, testerrFile);
			testerrFile.close(); // close error file 
		}
	}
	
	return 0;
}

int
FGMaximizer::showOutput()
{
	VSET& variableSet=vMgr->getVariableSet();
	fg->dumpVarMB_PairwiseFormat(outputDir,maxMBSize_Approx,variableSet,fMgr);
	return 0;
}

double
FGMaximizer::getScore(FactorGraph* aFg)
{
	double aScore=0;
	for(int i=0;i<aFg->getFactorCnt();i++)
	{
		SlimFactor* sFactor=aFg->getFactorAt(i);
		aScore=aScore+sFactor->mbScore;
	}
	return aScore;
}

int
FGMaximizer::dispRegTree_Genanatomy(map<int,RegressionTree*>& regtreeSet,const char* projectFName)
{
	char aFName[1024];
	sprintf(aFName,"mkdir %s/genatomyfiles",outputDir);
	if(system(aFName)!=0)
	{
		cout <<"Could not execute " << aFName << endl;
		return 0;
	}
	//Need to make five files
	//Project file
	GenanatomyConverter gc;
	gc.setTemplateProject(projectFName);
	sprintf(aFName,"%s/genatomyfiles/projectfile.gpf",outputDir);
	//Expression file
	char expressionFName[1024];
	sprintf(expressionFName,"%s/genatomyfiles/expr.tab",outputDir);
	char modulenetFName[1024];
	sprintf(modulenetFName,"%s/genatomyfiles/modnet.xml",outputDir);
	ofstream oFile(modulenetFName);
	oFile <<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	oFile <<"<root>" << endl;
	VSET& varSet=vMgr->getVariableSet();

	for(map<int,RegressionTree*>::iterator rIter=regtreeSet.begin();rIter!=regtreeSet.end();rIter++)
	{
		RegressionTree* rtree=rIter->second;
		oFile <<"<Module Name=\"Mod" << rIter->first <<"\">";
		rtree->genGenanatomy(oFile,varSet);
		oFile <<"<Expr>" << endl;
		int evidCnt=rtree->getDataSubset().size();
		for(int e=0;e<evidCnt;e++)
		{
			if(e==0)
			{
				oFile << "Exp" << e;
			}
			else
			{
				oFile <<" Exp" << e;
			}
		}	
		oFile <<"</Expr>" << endl;
		oFile <<"<Genes>" << endl;
		oFile << varSet[rIter->first]->getName() ;
		INTINTMAP treeVars;
		rtree->getTreeVars(treeVars);
		for(INTINTMAP_ITER aIter=treeVars.begin();aIter!=treeVars.end();aIter++)
		{
			oFile << " " << varSet[aIter->first]->getName();
		}
		oFile << endl;
		oFile <<"</Genes>" << endl;
		oFile <<"</Module>"<< endl;
	}
	oFile <<"</root>" << endl;
	oFile.close();

	char mainDir[1024];
	sprintf(mainDir,"%s/genatomyfiles",outputDir);
	gc.genTemplateFile(aFName,expressionFName,modulenetFName,mainDir);
	
	return 0;
}


int
FGMaximizer::sortMoves(vector<Move*>& moveSet)
{
	for(int i=0;i<moveSet.size();i++)
	{
		for(int j=i+1;j<moveSet.size();j++)
		{
			double delta1=moveSet[i]->getScoreImprovement();
			double delta2=moveSet[j]->getScoreImprovement();
			if(delta1<delta2)
			{
				Move* tempMove=moveSet[i];
				moveSet[i]=moveSet[j];
				moveSet[j]=tempMove;
			}
		}
	}
	return 0;
}


int 
FGMaximizer::updateMoveOrder(int currOrder, Move* newMove,vector<Move*>& moveSet,int& newOrder)
{
	newOrder=currOrder;
	double newScoreImp=newMove->getScoreImprovement();
	while(newOrder<moveSet.size())
	{
		if(newScoreImp>moveSet[newOrder]->getScoreImprovement())
		{
			break;
		}
		newOrder++;
	}
	delete moveSet[currOrder];
	if(newOrder==moveSet.size())
	{
		newOrder--;
	}
	//Everything from currOrder to newOrder-1 is shifted one place up
	for(int i=currOrder;i<newOrder;i++)
	{
		moveSet[i]=moveSet[i+1];
	}
	moveSet[newOrder]=newMove;
	return 0;
}


int 
FGMaximizer::attemptMove(FactorGraph* currFg,int currMoveOrder, vector<Move*>& moveSet,map<int,int>& srcVertex, map<int,int>& targetVertex)
{
	//Can this move fail even if srcVertex has not been affected yet.
	//Yes, if any variable in its MB cannot be updated because it exceeds the current k
	Move* aMove=moveSet[currMoveOrder];
	//cout <<endl << "Attempting move for: " << aMove->getSrcVertex() << " MBVars: ";
	INTINTMAP& targets=aMove->getTargetSet();
	for(INTINTMAP_ITER vIter=targets.begin();vIter!=targets.end();vIter++)
	{
	//	cout << " " << vIter->first;
	}
	//cout << endl;
	int moveStatus=fgeditor.makeMove(currFg,aMove,srcVertex,targetVertex);

	if(moveStatus==0)
	{
	//	cout <<"Move success for: " << aMove->getSrcVertex() << endl;
		for(INTINTMAP_ITER tIter=targets.begin();tIter!=targets.end();tIter++)
		{
			targetVertex[tIter->first]=0;
		}
		srcVertex[aMove->getSrcVertex()]=0;
	}
	else if (moveStatus==-1)
	{
	//	cout <<"Move failure for: " <<aMove->getSrcVertex() << endl;
		/*cout <<"Identifying new move for " << aMove->getSrcVertex()<< endl;
		Move* newMove=fgeditor.getNewMove(currFg,aMove);
		//If the move is null it means that we cannot find a new move which satisfies the affectingVars
		if(newMove!=NULL)
		{
			int moveOrder=-1;
			updateMoveOrder(currMoveOrder,newMove,moveSet,moveOrder);
			if(moveOrder==currMoveOrder)
			{
				attemptMove(currFg,moveOrder,moveSet,srcVertex,targetVertex);
			}
			else
			{
				cout <<"Update move order "<< currMoveOrder <<" -> "<< moveOrder << endl;
				attemptMove(currFg,currMoveOrder,moveSet,srcVertex,targetVertex);
			}
		}
		else
		{
			cout << "No new move for " << aMove->getSrcVertex() << endl;*/
			SlimFactor* sFactor=currFg->getFactorAt(aMove->getSrcVertex());
			for(INTINTMAP_ITER tIter=sFactor->mergedMB.begin();tIter!=sFactor->mergedMB.end();tIter++)
			{
				targetVertex[tIter->first]=0;
			}
			srcVertex[aMove->getSrcVertex()]=0;
		//}
		return 1;
	}
	else
	{
		return -1;
	}

	return 0;
}

//We need to make sure that the factorManager has all the one-variable extensions
//that may be need for this level of search
int
FGMaximizer::generateNextLevelClusters(int currK,FactorGraph* currFg)
{
	for(int i=0;i<currFg->getFactorCnt();i++)
	{
		SlimFactor* sFactor=currFg->getFactorAt(i);
		if((sFactor->mergedMB.size()>=maxMBSize)&&(sFactor->mergedMB.size()<maxMBSize_Approx))
		{	
			SlimFactor* mbFactor=fMgr->getFactorAt(sFactor->goodMBIDs.begin()->first);
			if(nextLevelFlag.find(mbFactor->fId)==nextLevelFlag.end())
			{
				fMgr->generateNextLevelClusters(mbFactor);
				nextLevelFlag[mbFactor->fId]=mbFactor->vCnt;
			}
			else
			{
				int currLevel=nextLevelFlag[mbFactor->fId];
				if(currLevel<currK)
				{
					fMgr->generateNextLevelClusters(mbFactor);
					nextLevelFlag[mbFactor->fId]=mbFactor->vCnt;
				}
			}
		}
	}
	return 0;
}
