#include <iostream>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
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


FactorManager::FactorManager()
{
	globalFactorID=0;
	maxFactorSize_Approx=-1;
	mbPenalty=0;
}

FactorManager::~FactorManager()
{
}


int 
FactorManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int 
FactorManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
FactorManager::setPotentialManager(PotentialManager* aPtr)
{
	potMgr=aPtr;
	return 0;
}

int 
FactorManager::setBaseInstantiation()
{
	VSET& variableSet=vMgr->getVariableSet();
	/*for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		Variable* aVar=vIter->second;
		INTVECT& varVals=aVar->getValues();
		defaultInstantiation[aVar->getID()]=varVals[0];
	}*/
	int evidCnt=evMgr->getNumberOfEvidences();
	map<string,int> evidDist;
	map<string,int> evidStrIDMap;
	char confStr[CONSTR_LEN];
	for(int e=0;e<evidCnt;e++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(e);
		string aConfStr;
		
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
			aConfStr.append(confStr);
		}
		if(evidDist.find(aConfStr)==evidDist.end())
		{
			evidDist[aConfStr]=1;
		}
		else
		{
			evidDist[aConfStr]=evidDist[aConfStr]+1;
		}
		evidStrIDMap[aConfStr]=e;
	}
	int maxCnt=0;
	int minCnt=evidCnt;
	int maxevid=-1;
	int minevid=-1;
	for(map<string,int>::iterator aIter=evidDist.begin();aIter!=evidDist.end();aIter++)
	{
		if(aIter->second>maxCnt)
		{
			maxCnt=aIter->second;
			maxevid=evidStrIDMap[aIter->first];
		}
		if(aIter->second<minCnt)
		{
			minCnt=aIter->second;
			minevid=evidStrIDMap[aIter->first];
		}
	}
	//maxevid=0;
	EMAP* evidMap=evMgr->getEvidenceAt(maxevid);
//	EMAP* evidMap=evMgr->getEvidenceAt(minevid);
	string defConfStr;
	for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
	{
		defaultInstantiation[eIter->first]=eIter->second->getHardEvidVal();
		sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
		defConfStr.append(confStr);
		defInstMap[eIter->first]=eIter->second->getHardEvidVal();
	}
	defProb=(double)evidDist[defConfStr]/(double)evidCnt;
	defInstID=maxevid;
	cout << "Default Inst: ID " << maxevid << " Freq: " << evidDist[defConfStr] << " Def prob: " << defProb << endl;
	cout << "Default Inst: " << defConfStr.c_str() << endl;
	return 0;
}


int 
FactorManager::setBaseInstantiation_Variable()
{
	VSET& variableSet=vMgr->getVariableSet();
	int evidCnt=evMgr->getNumberOfEvidences();
	map<string,int> evidStrIDMap;
	char confStr[CONSTR_LEN];
	INTINTMAP defInstMap;
	string defInstStr;
	defProb=0;
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		map<int,double> vvFreq;
		double total=0;
		for(int e=0;e<evidCnt;e++)
		{
			EMAP* evidMap=evMgr->getEvidenceAt(e);
			EMAP_ITER eIter=evidMap->find(vIter->first);
			int vval=eIter->second->getHardEvidVal();
			if(vvFreq.find(vval)==vvFreq.end())
			{
				vvFreq[vval]=1;
			}
			else
			{
				vvFreq[vval]=vvFreq[vval]+1;
			}
			total=total+1;
		}
		double maxFreq=0;
		int maxVal=0;
		for(map<int,double>::iterator fIter=vvFreq.begin();fIter!=vvFreq.end();fIter++)
		{
			if(fIter->second>maxFreq)
			{
				maxFreq=fIter->second;
				maxVal=fIter->first;
			}
		}
		defInstMap[vIter->first]=maxVal;
		sprintf(confStr,"-%d=%d-",vIter->first,maxVal);
		defInstStr.append(confStr);
		double pval=maxFreq/total;
		//defProb=defProb+log(pval);
		if(pval>defProb)
		{
			defProb=pval;
		}
	}
	//defProb=exp(defProb);
	cout << "Default Inst prob: "<<  defProb << endl;
	cout << "Default Inst: " << defInstStr.c_str() << endl;
	return 0;
}


int 
FactorManager::setMaxFactorSize(int size)
{
	maxFactorSize=size;
	if(maxFactorSize_Approx<maxFactorSize)
	{
		maxFactorSize_Approx=size;
	}
	return 0;
}


int 
FactorManager::setRandMISdCnt(double sdCnt) 
{
	misdCnt=sdCnt;
	return 0;
}


int 
FactorManager::setMaxFactorSize_Approx(int size)
{
	if(size>maxFactorSize)
	{
		maxFactorSize_Approx=size;
	}
	return 0;
}


int 
FactorManager::getMaxFactorSize()
{
	return maxFactorSize;
}

int 
FactorManager::getMaxFactorSize_Approx()
{
	return maxFactorSize_Approx;
}

int
FactorManager::setOutputDir(const char* aPtr)
{
	strcpy(outputDir,aPtr);
	return 0;
}

int
FactorManager::setModelName(const char* aPtr)
{
	strcpy(modelName,aPtr);
	return 0;
}

int
FactorManager::setGraph(Graph* gPtr)
{
	graph=gPtr;
	return 0;
}

int
FactorManager::setBeamSize(int aSize)
{
	beamSize=aSize;
	return 0;
}

int
FactorManager::setPenalty(double p)
{
	mbPenalty=p;
	return 0;
}

int
FactorManager::allocateFactorSpace()
{
	initFactorSet();
	populateFactorSet();
	return 0;
}

int
FactorManager::allocateFactorSpace_Graph()
{
	generateSingletonFactors();
	generateCanonicalFactors();
	return 0;
}


FactorGraph* 
FactorManager::createInitialFactorGraph()
{
	FactorGraph* fg=new FactorGraph;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* factor=fIter->second;
		if(factor->vCnt>1)
		{
			break;
		}
		SlimFactor* newFactor=new SlimFactor;
		newFactor->vCnt=1;
		newFactor->vIds=new int[1];
		newFactor->vIds[0]=factor->vIds[0];
		newFactor->fId=factor->fId;
		newFactor->mutualInfo=factor->mutualInfo;
		newFactor->conditionalEntropy=factor->conditionalEntropy;
		newFactor->marginalEntropy=factor->conditionalEntropy;
		newFactor->mbScore=factor->mbScore;
		newFactor->moveScore=factor->moveScore;
		fg->setFactor(newFactor);
	}
	return fg;
}

//This function is called when we know the topology and only want to estimate parameters
int 
FactorManager::paramEstimation(const char* estMethod)
{
	setBaseInstantiation();
	VSET& variableSet=vMgr->getVariableSet();
	if(allocateFactorSpace_Graph()==-1)
	{
		return -1;
	}
	estimateCanonicalParameters(estMethod);
	double dll=getLikelihood();
	cout <<"Data likelihood canonical " << dll << endl;
	double dll_exact=getLikelihood_ChainRule();
	cout <<"Data likelihood exact " << dll_exact << endl;
	evaluateMarkovBlanket(-1);
	return 0;	
}

int
FactorManager::generateSingletonFactors()
{
	VSET& variableSet=vMgr->getVariableSet();
	//Create the factors and the Markov blanket variables using the neighbours of each variable
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		Vertex* v=graph->getVertex(vIter->second->getName().c_str());
		if(v==NULL)
		{
			cout <<"No vertex corresponding to " << vIter->second->getName() << endl;
			return -1;
		}
		SlimFactor* sFactor=new SlimFactor;
		canonicalFactorSet[globalFactorID]=sFactor;
		sFactor->vIds=new int[1];
		sFactor->vIds[0]=vIter->first;
		sFactor->vCnt=1;
		sFactor->secondPId=-1;
		sFactor->mutualInfo=0;
		sFactor->conditionalEntropy=0;
		sFactor->fId=globalFactorID;
		string key;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		factorNameToIDMap[key]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=key;
		globalFactorID++;

		NINFO_MAP& neighbours=v->getImmediateNeighbours();
		for(NINFO_MAP_ITER nmIter=neighbours.begin();nmIter!=neighbours.end();nmIter++)
		//for(VSET_ITER nmIter=variableSet.begin();nmIter!=variableSet.end();nmIter++)
		{
			int vId=vMgr->getVarID(nmIter->first.c_str());
			//int vId=nmIter->first;
			if(vId==-1)
			{
				cout <<"Did not find variable id in variableManager for " << nmIter->first.c_str() << endl;
				return -1;
			}
			//Take care of the self-loops
			if(vId==sFactor->vIds[0])
			{
				continue;
			}
			sFactor->mergedMB[vId]=0;
		}
	}
	cout <<"Populated " << canonicalFactorSet.size() << " factors and their MB from graph " << endl;
	return 0;
}

int
FactorManager::generateCanonicalFactors()
{
	VSET& variableSet=vMgr->getVariableSet();
	int** subsetSpace=new int* [variableSet.size()];
	for(int i=0;i<variableSet.size();i++)
	{
		subsetSpace[i]=new int[variableSet.size()-1];
	}

	INTINTMAP parentIDs;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactorSet.begin();aIter!=canonicalFactorSet.end();aIter++)
	{
		parentIDs[aIter->first]=0;
	}
	int fSize=2;
	while(parentIDs.size()>0)
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		INTINTMAP tempParentIDs;
		for(INTINTMAP_ITER pIter=parentIDs.begin();pIter!=parentIDs.end();pIter++)
		{
			SlimFactor* pFactor=canonicalFactorSet[pIter->first];
			for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
			{
				int newVId=vIter->first;
				if(newVId<pFactor->vIds[pFactor->vCnt-1])
				{
					continue;
				}
				//Now add newVId to if newVId is in the Markov blanket of all variables in pFactor
				bool clqChk=true;
				int vid=0;
				while((vid<pFactor->vCnt) && clqChk)
				{
					SlimFactor* qFactor=canonicalFactorSet[pFactor->vIds[vid]];
					if(qFactor->mergedMB.find(newVId)==qFactor->mergedMB.end())
					{
						clqChk=false;
					}
					vid++;
				}
				if(!clqChk)
				{
					continue;
				}
				SlimFactor* sFactor=new SlimFactor;
				sFactor->vIds=new int[pFactor->vCnt+1];
				sFactor->vCnt=pFactor->vCnt+1;
				int fIter=0;
				int dIter=0;
				while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				sFactor->vIds[dIter]=newVId;
				dIter++;
				while(fIter<pFactor->vCnt) 
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				sFactor->fId=globalFactorID;
				string fKey;
				//Get the markov blanket variables 
				for(int i=0;i<sFactor->vCnt;i++)
				{
					SlimFactor* qFactor=canonicalFactorSet[sFactor->vIds[i]];
					for(INTINTMAP_ITER mIter=qFactor->mergedMB.begin();mIter!=qFactor->mergedMB.end();mIter++)
					{
						if(mIter->first!=sFactor->vIds[0])
						{
							sFactor->mergedMB[mIter->first]=mIter->second;
						}
					}
					if(i>0)
					{
						sFactor->mergedMB[qFactor->fId]=0;
					}
				}
				getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
				factorNameToIDMap[fKey]=globalFactorID;
				factorIDToNameMap[globalFactorID]=fKey;
				canonicalFactorSet[globalFactorID]=sFactor;
				globalFactorID++;
				//Here we need to add the super-set and sub-set relationships
				//Specifically, sFactor is a super-set of pFactor
				//pFactor is a sub-set of sFactor
				sFactor->secondPId=pFactor->fId;
				addToLattice(sFactor,subsetSpace);
				tempParentIDs[sFactor->fId]=0;
			}
		}
		cout <<"Added new " << tempParentIDs.size() << " factors of size "<< fSize << endl;
		fSize++;
		parentIDs.clear();
		for(INTINTMAP_ITER pIter=tempParentIDs.begin();pIter!=tempParentIDs.end();pIter++)
		{
			parentIDs[pIter->first]=0;
		}
		tempParentIDs.clear();
	}
	for(int i=0;i<variableSet.size();i++)
	{
		delete [] subsetSpace[i];
	}
	delete[] subsetSpace;
	cout << "Global factor id " << globalFactorID << endl;
	return 0;
}


int
FactorManager::generateCanonicalFactors(FactorGraph* fg,map<int,SlimFactor*>& canonicalFactors)
{
	VSET& variableSet=vMgr->getVariableSet();
	int** subsetSpace=new int* [variableSet.size()];
	for(int i=0;i<variableSet.size();i++)
	{
		subsetSpace[i]=new int[variableSet.size()-1];
	}

	INTINTMAP parentIDs;
	map<int,SlimFactor*>& allFactors=fg->getAllFactors();
	for(map<int,SlimFactor*>::iterator aIter=allFactors.begin();aIter!=allFactors.end();aIter++)
	{
		SlimFactor* sFactor=slimFactorSet[aIter->first];
		canonicalFactors[aIter->first]=sFactor;
		sFactor->mergedMB.clear();
		SlimFactor* rFactor=aIter->second;
		for(INTINTMAP_ITER vIter=rFactor->mergedMB.begin();vIter!=rFactor->mergedMB.end();vIter++)
		{
			sFactor->mergedMB[vIter->first]=vIter->second;
		}
		parentIDs[aIter->first]=0;
	}
	cout <<"Before making cliques" << endl;
	for(map<int,SlimFactor*>::iterator sIter=canonicalFactors.begin();sIter!=canonicalFactors.end();sIter++)
	{
		string key;
		SlimFactor* sFactor=sIter->second;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		cout<< "Can. factor: " << key.c_str() <<" MB: ";
		for(INTINTMAP_ITER mvIter=sFactor->mergedMB.begin();mvIter!=sFactor->mergedMB.end();mvIter++)
		{
			cout <<" "<< mvIter->first;
		}
		cout << endl;
	}
	//Now we need to make sure that the MBs of all the variables also form cliques
	int moreChange=1;
	int maxMBSize=maxFactorSize_Approx;
	if(maxMBSize>8)
	{
		maxMBSize=8;
	}
	//Hack to get out of loop
	map<int,int> alreadyDeletedFrom;
	while(moreChange)
	{
		moreChange=0;
		for(map<int,SlimFactor*>::iterator cIter=canonicalFactors.begin();cIter!=canonicalFactors.end();cIter++)
		{
			SlimFactor* cFactor=cIter->second;
			INTINTMAP toErase;
			for(INTINTMAP_ITER uIter=cFactor->mergedMB.begin();uIter!=cFactor->mergedMB.end();uIter++)
			{
				if(toErase.find(uIter->first)!=toErase.end())
				{
					continue;
				}
				SlimFactor* uFactor=canonicalFactors[uIter->first];	
				INTINTMAP_ITER vIter=uIter;
				vIter++;
				for(;vIter!=cFactor->mergedMB.end();vIter++)
				{
					SlimFactor* vFactor=canonicalFactors[vIter->first];
					//Now we need to make sure that u is in v's MB and v is in u's MB
					if(uFactor->mergedMB.find(vIter->first)==uFactor->mergedMB.end())
					{
						if((uFactor->mergedMB.size()<maxMBSize) && (vFactor->mergedMB.size()<maxMBSize))
						{
							if((alreadyDeletedFrom.find(vIter->first)==alreadyDeletedFrom.end()) && 
							   (alreadyDeletedFrom.find(uIter->first)==alreadyDeletedFrom.end()))
						   	{
								uFactor->mergedMB[vIter->first]=vIter->second;
								vFactor->mergedMB[uIter->first]=uIter->second;
							}
						}
						else
						{
							if(uFactor->mergedMB.size()>=maxMBSize)
							{
								toErase[uIter->first]=0;
								alreadyDeletedFrom[uIter->first]=0;
							}
							else
							{
								toErase[vIter->first]=0;
								alreadyDeletedFrom[vIter->first]=0;
							}
							moreChange++;
						}
					}
				}
			}
			for(INTINTMAP_ITER eIter=toErase.begin();eIter!=toErase.end();eIter++)
			{
				//If we are eliminating eIter->first from the MB of dFactor
				//we need to make sure that every variable other than the one
				//we are elminating must also not have an edge to the variable
				//we are eliminating
				SlimFactor* dFactor=canonicalFactors[eIter->first];
				INTINTMAP_ITER dIter=cFactor->mergedMB.find(eIter->first);
				INTINTMAP_ITER fIter=dFactor->mergedMB.find(cFactor->fId);
				if(dIter!=cFactor->mergedMB.end())
				{
					cFactor->mergedMB.erase(dIter);
				}
				if(fIter!=dFactor->mergedMB.end())
				{
					dFactor->mergedMB.erase(fIter);
				}
				for(INTINTMAP_ITER uIter=cFactor->mergedMB.begin();uIter!=cFactor->mergedMB.end();uIter++)
				{
					SlimFactor* eFactor=canonicalFactors[uIter->first];
					INTINTMAP_ITER dIter=eFactor->mergedMB.find(eIter->first);
					INTINTMAP_ITER fIter=dFactor->mergedMB.find(uIter->first);
					if(dIter!=eFactor->mergedMB.end())
					{
						eFactor->mergedMB.erase(dIter);
					}
					if(fIter!=dFactor->mergedMB.end())
					{
						dFactor->mergedMB.erase(fIter);
					}

				}
			}
		}

		cout <<"After making " << moreChange<< " changes " << endl;
		for(map<int,SlimFactor*>::iterator sIter=canonicalFactors.begin();sIter!=canonicalFactors.end();sIter++)
		{
			string key;
			SlimFactor* sFactor=sIter->second;
			getFactorKey(sFactor->vIds,sFactor->vCnt,key);
			cout<< "Can. factor: " << key.c_str() <<" MB: ";
			for(INTINTMAP_ITER mvIter=sFactor->mergedMB.begin();mvIter!=sFactor->mergedMB.end();mvIter++)
			{
				cout <<" "<< mvIter->first;
			}
			cout << endl;
		}
	}
	int fSize=2;
	while(parentIDs.size()>0)
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		INTINTMAP tempParentIDs;
		for(INTINTMAP_ITER pIter=parentIDs.begin();pIter!=parentIDs.end();pIter++)
		{
			SlimFactor* pFactor=canonicalFactors[pIter->first];
			for(INTINTMAP_ITER vIter=pFactor->mergedMB.begin();vIter!=pFactor->mergedMB.end();vIter++)
			{
				int newVId=vIter->first;
				//Now add newVId to if newVId is in the Markov blanket of all variables in pFactor
				bool clqChk=true;
				int vid=0;
				while((vid<pFactor->vCnt) && clqChk)
				{
					SlimFactor* qFactor=canonicalFactors[pFactor->vIds[vid]];
					if(qFactor->mergedMB.find(newVId)==qFactor->mergedMB.end())
					{
						clqChk=false;
					}
					vid++;
				}
				if(!clqChk)
				{
					continue;
				}
				SlimFactor* sFactor=new SlimFactor;
				sFactor->vIds=new int[pFactor->vCnt+1];
				sFactor->vCnt=pFactor->vCnt+1;
				int fIter=0;
				int dIter=0;
				while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				sFactor->vIds[dIter]=newVId;
				dIter++;
				while(fIter<pFactor->vCnt) 
				{
					sFactor->vIds[dIter]=pFactor->vIds[fIter];
					fIter++;
					dIter++;
				}
				string fKey;
				getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
				if(factorNameToIDMap.find(fKey)!=factorNameToIDMap.end())
				{
					delete sFactor;
					int fId=factorNameToIDMap[fKey];
					sFactor=slimFactorSet[fId];
					sFactor->mergedMB.clear();
				}
				else
				{
					potMgr->populateFactor(slimFactorSet,variableSet,sFactor,false);
					factorNameToIDMap[fKey]=globalFactorID;
					factorIDToNameMap[globalFactorID]=fKey;
					sFactor->fId=globalFactorID;
					slimFactorSet[sFactor->fId]=sFactor;
					globalFactorID++;
				}
				canonicalFactors[sFactor->fId]=sFactor;
				//Here we need to add the super-set and sub-set relationships
				//Specifically, sFactor is a super-set of pFactor
				//Get the markov blanket variables 
				for(int i=0;i<sFactor->vCnt;i++)
				{
					SlimFactor* qFactor=canonicalFactors[sFactor->vIds[i]];
					for(INTINTMAP_ITER mIter=qFactor->mergedMB.begin();mIter!=qFactor->mergedMB.end();mIter++)
					{
						if(mIter->first!=sFactor->vIds[0])
						{
							sFactor->mergedMB[mIter->first]=mIter->second;
						}
					}
					if(i>0)
					{
						sFactor->mergedMB[qFactor->fId]=0;
					}
				}
				addToLattice(sFactor,subsetSpace);
				tempParentIDs[sFactor->fId]=0;
			}
		}
		cout <<"Added new " << tempParentIDs.size() << " factors of size "<< fSize << endl;
		fSize++;
		parentIDs.clear();
		for(INTINTMAP_ITER pIter=tempParentIDs.begin();pIter!=tempParentIDs.end();pIter++)
		{
			parentIDs[pIter->first]=0;
		}
		tempParentIDs.clear();
	}
	for(int i=0;i<variableSet.size();i++)
	{
		delete [] subsetSpace[i];
	}
	/*for(map<int,SlimFactor*>::iterator sIter=canonicalFactors.begin();sIter!=canonicalFactors.end();sIter++)
	{
		string key;
		SlimFactor* sFactor=sIter->second;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		cout<< "Can. factor: " << key.c_str() <<" MB: ";
		for(INTINTMAP_ITER mvIter=sFactor->mergedMB.begin();mvIter!=sFactor->mergedMB.end();mvIter++)
		{
			cout <<" "<< mvIter->first;
		}
		cout << endl;
	}*/
	
	delete[] subsetSpace;
	return 0;
}


//In this function we are going to associate with each canonical parameter, a joint potential table.
//Each entry is a pair corresponding to a configuration of variables of the canonical factor and its value.
//The value is estimated by taking the exponent of the sum over all subsets using the default instantiation.
int
FactorManager::estimateCanonicalParameters(const char* estMethod)
{
	VSET& variableSet=vMgr->getVariableSet();
	char aFName[256];
	sprintf(aFName,"%s/%s.txt",outputDir,estMethod);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactorSet.begin();aIter!=canonicalFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		INTINTMAP allSubsets;
		lattice.getAllSubsets(aIter->first,allSubsets);
		allSubsets[sFactor->fId]=0;
		if(strstr(estMethod,"joint")!=NULL)
		{
			potMgr->estimateCanonicalPotential_Joint(sFactor,variableSet,defaultInstantiation,allSubsets,canonicalFactorSet);
		}
		else if(strstr(estMethod,"markovblnkt")!=NULL)
		{
			potMgr->estimateCanonicalPotential(sFactor,variableSet,defaultInstantiation,allSubsets,canonicalFactorSet);
		}
		else if(strstr(estMethod,"approx")!=NULL)
		{
			potMgr->estimateCanonicalPotential_Approximate(sFactor,variableSet,defaultInstantiation,allSubsets,canonicalFactorSet);
		}
		oFile <<"Canonical potential for ID: " << aIter->first <<" ";
		sFactor->canonicalParams->dumpPotential(oFile);
	}
	return 0;
}

int 
FactorManager::learnStructure()
{
	Error::ErrorCode err=estimateClusterProperties();
	if(err!=Error::SUCCESS)
	{
		cout <<Error::getErrorString(err) << endl;
		return 0;
	}
	return 0;
}

int 
FactorManager::showStructure()
{
	VSET& variableSet=vMgr->getVariableSet();
	char aFName[1024];
	sprintf(aFName,"%s/strct_k%d.txt",outputDir,maxFactorSize);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* aFactor=fIter->second;
		oFile <<fIter->first<< " " << factorIDToNameMap[fIter->first] 
			<<" " << aFactor->mutualInfo << " " 
			<< aFactor->conditionalEntropy<< endl;
	}
	oFile.close();
	return 0;
}

int 
FactorManager::showStructure_allK()
{
	VSET& variableSet=vMgr->getVariableSet();
	char aFName[1024];
	map<int,ofstream*> filePtrs;
	for(int k=2;k<=maxFactorSize;k++)
	{
		sprintf(aFName,"%s/strct_k%d.txt",outputDir,k);
		ofstream* oFile=new ofstream(aFName);
		filePtrs[k]=oFile;
	}
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* aFactor=fIter->second;
		for(int k=2;k<=maxFactorSize;k++)
		{
			if(aFactor->vCnt<=k)
			{
				ofstream* oFile=filePtrs[k];
				(*oFile) <<fIter->first<< " " << factorIDToNameMap[fIter->first] 
				<<" " << aFactor->mutualInfo << " " 
				<< aFactor->conditionalEntropy<< endl;
			}
		}
	}
	for(int k=2;k<=maxFactorSize;k++)
	{
		ofstream* oFile=filePtrs[k];
		oFile->close();
		delete oFile;
	}
	filePtrs.clear();
	return 0;
}

int
FactorManager::readStructure()
{
	char aFName[1024];
	sprintf(aFName,"%s/strct_k%d.txt",outputDir,maxFactorSize);
	if(readClusterProperties(aFName)==-1)
	{
		return -1;
	}
	return 0;
}

//This reads a file of the format written by the above function, showStructure
int
FactorManager::readClusterProperties(const char* aFName)
{
	ifstream inFile(aFName);
	if(!inFile.good())
	{
		inFile.close();
		return -1;
	}
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int tokCnt=0;
		int fID=0;
		string fKey;
		double condEntropy=0;
		double mutualInfo=0;
		
		char* tok=strtok(buffer," ");
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					fID=atoi(tok);
					break;
				}
				case 1:
				{
					fKey.append(tok);
					break;
				}
				case 2:
				{
					mutualInfo=atof(tok);
					break;
				}
				case 3:
				{
					condEntropy=atof(tok);
					break;
				}
			}
			tok=strtok(NULL," ");
			tokCnt++;
		}
		SlimFactor* factor=slimFactorSet[fID];
		factor->mutualInfo=mutualInfo;
		factor->conditionalEntropy=condEntropy;
		if(factor->vCnt==1)
		{
			factor->marginalEntropy=condEntropy;
			factor->mbScore=condEntropy;
			factor->moveScore=condEntropy;
		}
		INTINTMAP fvarSet;
		for(int i=0;i<factor->vCnt;i++)
		{
			fvarSet[factor->vIds[i]]=0;
		}
		potMgr->setConditionalEntropy(factor->vIds[factor->vCnt-1],fvarSet,vMgr->getVariableSet(),condEntropy);
	}
	inFile.close();
	return 0;
}

//Dump all the information files in outputDir

int
FactorManager::estimateRandomInfo()
{
	VSET& variableSet=vMgr->getVariableSet();
	vector<double> randInfo;
	char fName[256];
	sprintf(fName,"%s/randmi_summary.txt",outputDir);
	ofstream misummary(fName);
	misummary << "k\tmean_rand_mi\tstd_rand_mi" << endl;
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	evMgr->randomizeEvidence(r);
	potMgr->initRandom();

	potMgr->estimateMarginalEntropies(slimFactorSet,variableSet,true);
	//potMgr->computeSummaryStats(maxFactorSize,true);
	for(int i=2;i<=maxFactorSize;i++)
	{
		int fCnt=combCnt(variableSet.size(),i);
		potMgr->estimateRandomInfo(slimFactorSet,variableSet,randInfo,i);
		double mean=0;
		sprintf(fName,"%s/randmi_k%d.txt",outputDir,i);
		ofstream oFile(fName);
		for(int j=0;j<randInfo.size();j++)
		{
			oFile <<randInfo[j] << endl;
			mean=mean+randInfo[j];
		}
		oFile.close();
		mean=mean/randInfo.size();
		double std=0;
		for(int j=0;j<fCnt;j++)
		{
			double diff=mean-randInfo[j];
			std=std+(diff*diff);
		}
		std=sqrt(std/(randInfo.size()-1));
		misummary << i << "\t" << mean << "\t" << std << endl;
		randInfo.clear();
	}
	misummary.close();
	return 0;
}

//This is similar to the function above except that for k greater than exactk it uses sampleCnt
//number of combinations
int
FactorManager::estimateRandomInfo_Approximate(int sampleCnt)
{
	VSET& variableSet=vMgr->getVariableSet();
	vector<double> randInfo;
	char fName[256];
	sprintf(fName,"%s/randmi_approx_asummary.txt",outputDir);
	ofstream misummary(fName);
	misummary << "k\tmean_rand_mi\tstd_rand_mi" << endl;
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	evMgr->randomizeEvidence(r);
	potMgr->initRandom();
	potMgr->estimateMarginalEntropies(slimFactorSet,variableSet,true);
	//potMgr->computeSummaryStats(maxFactorSize,true);
	for(int i=2;i<=maxFactorSize_Approx;i++)
	{
		cout <<"Estimating randinfo for k " << i << endl;
		double mean=0;
		sprintf(fName,"%s/randmi_k%d.txt",outputDir,i);
		ofstream oFile(fName);
		if(i<=maxFactorSize)
		{
			potMgr->estimateRandomInfo(slimFactorSet,variableSet,randInfo,i);
			for(int j=0;j<randInfo.size();j++)
			{
				oFile <<randInfo[j] << endl;
				mean=mean+randInfo[j];
			}
		}
		else
		{
			//Create a random factor with i variables
			SlimFactor* sFactor=new SlimFactor;
			sFactor->vIds=new int[i];
			sFactor->vCnt=i;
			int j=0;
			int duplicateFCnt=0;
			map<string,int> usedFactorIDs;
			double step=1.0/((double) variableSet.size());
			while(j<sampleCnt)
			{
				map<int,int> usedVarIDs;
				int vid=0;
				while(usedVarIDs.size()<i)
				{
					double rVal=gsl_ran_flat(r,0,1);
					int rind=(int)(rVal/step);
					while(usedVarIDs.find(rind)!=usedVarIDs.end())
					{
						rVal=gsl_ran_flat(r,0,1);
						rind=(int)(rVal/step);
					}
					usedVarIDs[rind]=0;
					sFactor->vIds[vid]=rind;
					vid++;
				}
				string akey;
				getFactorKey(sFactor->vIds,i,akey);
				if(usedFactorIDs.find(akey)!=usedFactorIDs.end())
				{
					duplicateFCnt++;
				}
				usedFactorIDs[akey]=0;
				//Estimate the information associated with this factor
				potMgr->populateFactor(slimFactorSet,variableSet,sFactor,true);
				double mi=sFactor->mutualInfo;
				randInfo.push_back(mi);
				oFile <<mi << endl;
				mean=mean+mi;
				j++;	
			}
			delete sFactor;
			if(duplicateFCnt>0)
			{
				cout <<"Duplicates : " << duplicateFCnt << " out of a total of "<< sampleCnt << " factors of size " << i << endl;
			}
		}
		oFile.close();
		mean=mean/randInfo.size();
		double std=0;
		for(int j=0;j<randInfo.size();j++)
		{
			double diff=mean-randInfo[j];
			std=std+(diff*diff);
		}
		std=sqrt(std/(randInfo.size()-1));
		misummary << i << "\t" << mean << "\t" << std << endl;
		randInfo.clear();
	}
	misummary.close();
	gsl_rng_free(r);
	return 0;
}


int
FactorManager::readRandomInfo()
{
	char fName[256];
	sprintf(fName,"%s/randmi_approx_asummary.txt",outputDir);
	ifstream misummary(fName);
	int lineCnt=0;
	char buffer[1024];
	while(misummary.good())
	{
		misummary.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(lineCnt==0)
		{
			lineCnt++;
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int k;
		double mi_mean;
		double mi_std;
		while(tok!=NULL)
		{
			switch(tokCnt)
			{
				case 0:
				{
					k=atoi(tok);
					break;
				}
				case 1:
				{
					mi_mean=atof(tok);
					break;
				}
				case 2:
				{
					mi_std=atof(tok);
					break;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		randMI_mean[k]=mi_mean;
		randMI_std[k]=mi_std;
	}
	
	misummary.close();
	if(randMI_mean.size()<(maxFactorSize-1))
	{
		return -1;
	}
	return 0;
}


//This takes the number of standard deviations a cluster's information must exceed the mean random mi
int
FactorManager::filterClustersWithMI(double eps)
{
	INTINTMAP badFactors;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		double randmi_m=randMI_mean[sFactor->vCnt];
		double randmi_sd=randMI_std[sFactor->vCnt];
		if(sFactor->mutualInfo < (randmi_m + (randmi_sd*eps)))
		{
			//Erase this cluster from the lattice
			deleteFromLattice(fIter->first);
			//Erase from factorNameToIDMap and from factorIDToNameMap
			map<int,string>::iterator idnameIter=factorIDToNameMap.find(fIter->first);
			map<string,int>::iterator nameidIter=factorNameToIDMap.find(idnameIter->second);
			delFactors_MI[nameidIter->first]=fIter->second->mutualInfo;
			//factorIDToNameMap.erase(idnameIter);
			//factorNameToIDMap.erase(nameidIter);
			//slimFactorSet.erase(fIter);
			if(badFactors.find(sFactor->vCnt)==badFactors.end())
			{
				badFactors[sFactor->vCnt]=1;
			}
			else
			{
				badFactors[sFactor->vCnt]=badFactors[sFactor->vCnt]+1;
			}
		}
	}
	for(map<string,double>::iterator bIter=delFactors_MI.begin();bIter!=delFactors_MI.end();bIter++)
        {
                map<string,int>::iterator nameidIter=factorNameToIDMap.find(bIter->first);
                map<int,string>::iterator idnameIter=factorIDToNameMap.find(nameidIter->second);
                map<int,SlimFactor*>::iterator gIter=slimFactorSet.find(idnameIter->first);
                if(gIter==slimFactorSet.end())
                {
                        cout <<"Problem: factor " << idnameIter->first << " not found " << endl;
                        return -1;
                }
                factorIDToNameMap.erase(idnameIter);
                factorNameToIDMap.erase(nameidIter);
                slimFactorSet.erase(gIter);
        }
	for(INTINTMAP_ITER imIter=badFactors.begin();imIter!=badFactors.end();imIter++)
	{
		cout <<"Eliminated " << imIter->second << " bad factors of size " << imIter->first << endl;
	}
	showAllFactors(eps);
	return 0;
}

int
FactorManager::applyDPICorrection(double threshold, double dpiPercent)
{
	if(maxFactorSize>2)
	{
		cout << "Cannot perform DPI correction on factors larger than 2" << endl;
		return 0;
	}
	cout <<"Applying DPI correction " << endl;

	map<int,bool> factorStatus;
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		factorStatus[aIter->first]=true;
	}

	VSET& variableSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt==1)
		{
			factorStatus[fIter->first]=false;		
			continue;
		}
		if(sFactor->mutualInfo<threshold)
		{
			factorStatus[fIter->first]=false;
			continue;
		}
		//Other wise make a triple
		int triple[3];
		triple[0]=sFactor->vIds[0];
		triple[1]=sFactor->vIds[1];
		int factorIds[3];
		for(int j=sFactor->vIds[1]+1;j<variableSet.size();j++)
		{
			triple[2]=j;
			
			int fVars[2];
			for(int k=0;k<3;k++)
			{
				if(k==0)
				{
					fVars[0]=triple[0];
					fVars[1]=triple[1];
				}
				else if(k==1)
				{
					fVars[0]=triple[0];
					fVars[1]=triple[2];
				}
				else if(k==2)
				{
					fVars[0]=triple[1];
					fVars[1]=triple[2];
				}
				int fid=getFactorIndex(fVars,2);
				factorIds[k]=fid;
			}
			for(int k=0;k<3;k++)
			{
				for(int l=k+1;l<3;l++)
				{
					if(slimFactorSet[factorIds[k]]->mutualInfo < slimFactorSet[factorIds[l]]->mutualInfo)
					{
						int tempId=factorIds[k];
						factorIds[k]=factorIds[l];
						factorIds[l]=tempId;
					}
				}
			}
			double minInfo=slimFactorSet[factorIds[2]]->mutualInfo;
			double nextToMin=slimFactorSet[factorIds[1]]->mutualInfo;
			double withDpiInfo=nextToMin-(dpiPercent*nextToMin);
			if(minInfo<=withDpiInfo)
			{
				factorStatus[factorIds[2]]=false;
			}
		
		}
	}
	char aFName[1024];
	sprintf(aFName,"%s/dpi.txt",outputDir);
	ofstream oFile(aFName);
	map<int,INTINTMAP*> varNeighbour;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		if(factorStatus[fIter->first])
		{
			SlimFactor* factor=fIter->second;
			int uid=factor->vIds[0];
			int vid=factor->vIds[1];
			oFile<<variableSet[uid]->getName()
				<<"\t" <<variableSet[vid]->getName() << endl;
			INTINTMAP* nhood1=NULL;
			INTINTMAP* nhood2=NULL;
			if(varNeighbour.find(uid)==varNeighbour.end())
			{
				nhood1=new INTINTMAP;
				varNeighbour[uid]=nhood1;
			}
			else
			{
				nhood1=varNeighbour[uid];
			}
			if(varNeighbour.find(vid)==varNeighbour.end())
			{
				nhood2=new INTINTMAP;
				varNeighbour[vid]=nhood2;
			}
			else
			{
				nhood2=varNeighbour[vid];
			}
			(*nhood1)[vid]=0;
			(*nhood2)[uid]=0;
		}
	}
	oFile.close();
	sprintf(aFName,"%s/mbscore_dpi.txt",outputDir);//,dpiPercent,threshold);
	ofstream nFile(aFName);

	nFile<< "Var\tMB\tCondEntr\tMargEntr" << endl;
	for(map<int,INTINTMAP*>::iterator fIter=varNeighbour.begin();fIter!=varNeighbour.end();fIter++)
	{
		nFile <<variableSet[fIter->first]->getName() <<"\t";
		INTINTMAP* nhood=fIter->second;
                for(INTINTMAP_ITER vIter=nhood->begin();vIter!=nhood->end();vIter++)
                {
                        if(vIter!=nhood->begin())
                        {
                                nFile <<"-";
                        }
                        nFile<< variableSet[vIter->first]->getName();
                }
		(*nhood)[fIter->first]=0;
		double condEntr=potMgr->getConditionalEntropy(fIter->first,*nhood,variableSet);
		INTINTMAP_ITER delIter=nhood->find(fIter->first);
		nhood->erase(delIter);
                nFile <<"\t" << condEntr <<"\t" << slimFactorSet[fIter->first]->conditionalEntropy << endl;
	}
	nFile.close();
	return 0;
}


int
FactorManager::findTrueHO(double sdCnt)
{
	int** subset=new int* [maxFactorSize_Approx];
	int* subsetId=new int [maxFactorSize_Approx];
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		subset[i]=new int[maxFactorSize_Approx-1];
	}
	//Count those k>2 clusters for which at the k-1 there was little information
	//but adding the kth variable increased information beyond random chance
	//Such clusters would not be detected if we used a greedy algorithm
	INTINTMAP trueHO;
	//Count up those k>2 clusters that would be detected if we used a greedy approach
	INTINTMAP greedyHO;
	filterClustersWithMI(sdCnt);
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt<3)
		{
			continue;
		}

		//Check if the mutual information
		//of this factor is greater than random
		//because of the pairwise information
		//or because of true higher order dependency
			
		//Create all subsets of variables of sFactor of size one less
		int sSize=sFactor->vCnt-1;
		//sInd will be the index that will be eliminated
		for(int sInd=0;sInd<sFactor->vCnt;sInd++)
		{
			int vInd=0;
			int currSize=0;
			while(currSize<sSize)
			{
				if(vInd!=sInd)
				{
					subset[sInd][currSize]=sFactor->vIds[vInd];
					currSize++;
				}
				vInd++;
			}
		}
		//Now that we have created sFactor->vCnt subsets of 1 less 
		//we need to check if the addition of one variable significantly increases
		//the MI.
		int foundSubsets=0;
		int ssId=0;
		while((ssId<sFactor->vCnt) && (foundSubsets==0))
		{
			int fId=getFactorIndex(subset[ssId],sSize);
			if(slimFactorSet.find(fId)!=slimFactorSet.end())
			{
				foundSubsets++;
			}
			ssId++;
		}
		//If there is atleast one subset that is found at the lower level
		//we can find the higher order dependency via a greedy approach
		if(foundSubsets)
		{
			if(greedyHO.find(sFactor->vCnt)==greedyHO.end())
			{
				greedyHO[sFactor->vCnt]=1;
			}
			else
			{
				greedyHO[sFactor->vCnt]=greedyHO[sFactor->vCnt]+1;
			}
		}
		else
		{
			if(trueHO.find(sFactor->vCnt)==trueHO.end())
			{
				trueHO[sFactor->vCnt]=1;
			}
			else
			{
				trueHO[sFactor->vCnt]=trueHO[sFactor->vCnt]+1;
			}
		}
	}
	char aFName[1024];
	sprintf(aFName,"%s/hotype.txt",outputDir);
	ofstream oFile(aFName);
	oFile <<"k\tGreedyHO\tTrueHO"<< endl;
	for(int k=3;k<=maxFactorSize_Approx;k++)
	{
		oFile << k;
		if(greedyHO.find(k)==greedyHO.end())
		{
			oFile <<"\t0";
		}
		else
		{
			oFile <<"\t" << greedyHO[k];
		}
		if(trueHO.find(k)==trueHO.end())
		{
			oFile <<"\t0";
		}
		else
		{
			oFile <<"\t" << trueHO[k];
		}
		oFile << endl;
	}
	oFile.close();
	delete [] subsetId;
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		delete [] subset[i];
	}
	delete[] subset;
	return 0;
}



//We use simple Apriori-like algorithm to find clusters. For clusters of size k<=maxFactorSize, for
//which we can compute multi-information exactly we use the random-mutual information as a threshold
//to get good clusters. For clusters of size greater than maxFactorSize, we will use confidence criteria
//similar to Apriori
int
FactorManager::generateClusters(double epsilon,bool latticeCheckLowerLevel , double reqConf, int maxClusterSize)
{
	//Start with all factors of size 2
	int currK=2;
	map<int,Variable*>& variableSet=vMgr->getVariableSet();
	int oldClusterCnt=0;
	int newClusterCnt=1;
	INTVECT parentIds;
	INTVECT newParentIds;
	INTINTMAP generatedIDs;

	//The variables of the new factor
	int* newVids=new int[maxClusterSize];
	while((newClusterCnt>0) && (currK<=maxClusterSize))//No more factors can be added
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		double randmi_mean=randMI_mean[currK];
		double randmi_sd=randMI_std[currK];
		if(currK<=maxFactorSize_Approx)
		{
			int pFidCnt=parentIds.size();
			if(pFidCnt==0)
			{
				pFidCnt=variableSet.size();
			}
			for(int p=0;p<pFidCnt;p++)
			{
				SlimFactor* pFactor;
				if(parentIds.size()==0)
				{
					pFactor=slimFactorSet[p];
				}
				else
				{
					pFactor=goodSlimFactors[parentIds[p]];
				}
				for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
				{
					int newVId=vIter->first;
					if(newVId<=pFactor->vIds[pFactor->vCnt-1])
					{
						continue;
					}
					for(int v=0;v<pFactor->vCnt;v++)
					{
						newVids[v]=pFactor->vIds[v];
					}
					newVids[pFactor->vCnt]=newVId;
					int currFid=getFactorIndex(newVids,pFactor->vCnt+1);
					SlimFactor* sFactor=slimFactorSet[currFid];

					//Check for support
					if(sFactor->mutualInfo >= (randmi_mean + (randmi_sd*epsilon)))
					{
						if((sFactor->vCnt>2) && latticeCheckLowerLevel)
						{
							INTINTMAP* subsets=lattice.getSubsets(sFactor->fId);
							//get all subsets and check for confidence
							//We need to check the last sFactor->vCnt subsets as these
							//will be the ones that correspond to the maximal subsets
							map<int,int>::reverse_iterator rIter=subsets->rbegin();
							int sscnt=0;
							int hitCnt=0;
							SlimFactor* aSubset=slimFactorSet[rIter->first];
							while((sscnt<sFactor->vCnt) && ((sFactor->vCnt-aSubset->vCnt) ==1 ))
							{
								int ssId=rIter->first;
								if(goodSlimFactors.find(ssId)!=goodSlimFactors.end())
								{
									hitCnt++;
								}
								sscnt++;
								rIter++;
								aSubset=slimFactorSet[rIter->first];
							}
							double conf=((double)hitCnt)/((double) sFactor->vCnt);
							if(conf>=reqConf)
							{
								goodSlimFactors[sFactor->fId]=sFactor;
								newParentIds.push_back(sFactor->fId);
							}
						}
						else //No need to check confidence at the lower levels where
							//we can compute multi-information correctly
						{
							goodSlimFactors[sFactor->fId]=sFactor;
							newParentIds.push_back(sFactor->fId);
						}
					}
				}
			}
			//Now update parentIds from newParentIds
			parentIds.clear();
			for(int i=0;i<newParentIds.size();i++)
			{
				parentIds.push_back(newParentIds[i]);
			}
			newParentIds.clear();
		}
		else
		{
			int pFidCnt=parentIds.size();
			for(int p=0;p<pFidCnt;p++)
			{
				SlimFactor* pFactor=goodSlimFactors[parentIds[p]];
				for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
				{
					int newVId=vIter->first;
					if(newVId<=pFactor->vIds[pFactor->vCnt-1])
					{
						continue;
					}
					//Here we have to create a new factor
					SlimFactor* sFactor=new SlimFactor;
					sFactor->vCnt=pFactor->vCnt+1;
					sFactor->vIds=new int[pFactor->vCnt+1];
					for(int j=0;j<pFactor->vCnt;j++)
					{
						sFactor->vIds[j]=pFactor->vIds[j];
					}
					sFactor->vIds[sFactor->vCnt-1]=newVId;
					sFactor->fId=getFactorIndex(sFactor->vIds,sFactor->vCnt);
					if(generatedIDs.find(sFactor->fId) !=generatedIDs.end())
					{
						delete sFactor;
						continue;
					}
					generatedIDs[sFactor->fId]=0;
					potMgr->populateFactor(slimFactorSet,variableSet,sFactor,false);
					if(sFactor->mutualInfo < (randmi_mean + (randmi_sd*epsilon)))
					{
						delete sFactor;
						continue;
					}

					//Allocate memory for subsets
					int** subsets=new int*[sFactor->vCnt];
					for(int i=0;i<sFactor->vCnt;i++)
					{
						subsets[i]=new int[sFactor->vCnt-1];
					}
					sFactor->generateMaximalSubsets(subsets);
					//Now check the confidence and in the meantime store the subset ids
					//to update the lattice structure
					int* ssIds=new int[sFactor->vCnt];
					int hitCnt=0;
					for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
					{
						int sId=getFactorIndex(subsets[sscnt],sFactor->vCnt-1);
						ssIds[sscnt]=sId;
						if(goodSlimFactors.find(sId)!=goodSlimFactors.end())
						{
							hitCnt++;
						}
					}
					double conf=((double)hitCnt)/((double)sFactor->vCnt);
					if(conf>=reqConf)
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						newParentIds.push_back(sFactor->fId);
						//Update the lattice structure
						for(int i=0;i<sFactor->vCnt;i++)
						{
							lattice.addSubset(ssIds[i],sFactor->fId);
							lattice.addSuperset(sFactor->fId,ssIds[i]);
						}
					}
					else
					{
						delete sFactor;
					}
					for(int i=0;i<sFactor->vCnt;i++)
					{
						delete[] subsets[i];
					}
					delete[] subsets;
					delete ssIds;
				}
			}
			parentIds.clear();
			for(int i=0;i<newParentIds.size();i++)
			{
				parentIds.push_back(newParentIds[i]);
			}
			newParentIds.clear();
		}
		currK++;
		newClusterCnt=goodSlimFactors.size()-oldClusterCnt;
		oldClusterCnt=goodSlimFactors.size();
		cout <<"Added " << newClusterCnt<< " new factors "<< endl; 
	}
	return 0;
}

//In this function, we go over all good slim factors and prune out the ones that have already been found at higher levels.
int
FactorManager::getMaximalClusters()
{
	//Start from the smallest clusters
	map<int,SlimFactor*>::iterator factorIter=goodSlimFactors.begin();
	while(factorIter!=goodSlimFactors.end())
	{
		int fid=factorIter->first;
		INTINTMAP* supersets=lattice.getSupersets(fid);
		if(supersets!=NULL)
		{
			//If any of the supersets of fid exist in goodSlimFactors then get rid of fid
			//This is because, the dependencies in fid are captured by the superset of fid
			INTINTMAP_ITER iIter=supersets->begin();
			bool foundSuper=false;
			while((iIter!=supersets->end()) && (!foundSuper))
			{
				if(goodSlimFactors.find(iIter->first)!=goodSlimFactors.end())
				{
					foundSuper=true;
				}
				iIter++;
			}
			if(foundSuper)
			{
				goodSlimFactors.erase(factorIter);
			}
		}
		factorIter++;
	}
	return 0;	
}

int
FactorManager::removeDupFactors()
{
	INTINTMAP factorVars;
	string factorVarKey;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		//Need to consider only those factors that may have a Markov blanket 
		if(fIter->second->vCnt==maxFactorSize_Approx)
		{
			break;
		}
		factorVars.clear();
		factorVarKey.clear();
		SlimFactor* sFactor=fIter->second;
		for(int i=0;i<sFactor->vCnt;i++)
		{
			factorVars[sFactor->vIds[i]]=0;
		}
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			factorVars[mbIter->first]=0;
		}
		char keyPart[32];
		//First append the factor size
		sprintf(keyPart,"%d:",sFactor->vCnt);
		factorVarKey.append(keyPart);
		for(INTINTMAP_ITER vIter=factorVars.begin();vIter!=factorVars.end();vIter++)
		{
			if(vIter!=factorVars.begin())
			{
				sprintf(keyPart,"-%d",vIter->first);
			}
			else
			{
				sprintf(keyPart,"%d",vIter->first);
			}
			factorVarKey.append(keyPart);
		}
		INTDBLMAP* factorSet=NULL;
		if(signFactGrpMap.find(factorVarKey)==signFactGrpMap.end())
		{
			factorSet=new INTDBLMAP;
			signFactGrpMap[factorVarKey]=factorSet;
		}
		else
		{
			factorSet=signFactGrpMap[factorVarKey];
		}
		(*factorSet)[fIter->first]=fIter->second->mbScore;
		factSignMap[fIter->first]=factorVarKey;
	}
	//Now get the non-redundant factors per signature
	char aFName[1024];
	sprintf(aFName,"%s/factsign_k%d.txt",outputDir,maxFactorSize_Approx);
	VSET& variableSet=vMgr->getVariableSet();
	ofstream oFile(aFName);
	oFile <<"GrpSign\tGrpSize\tRepFid\tFactor\tMBlanket" << endl;
	
	for(map<string,INTDBLMAP*>::iterator fSetIter=signFactGrpMap.begin();fSetIter!=signFactGrpMap.end();fSetIter++)
	{
		INTDBLMAP* factorSet=fSetIter->second;
		double minEntropy=10000;
		int minFid=-1;
		SlimFactor* minFactor=NULL;
		for(INTDBLMAP_ITER idIter=factorSet->begin();idIter!=factorSet->end();idIter++)
		{
			double fEntropy=idIter->second;
			if(fEntropy<minEntropy)
			{
				minEntropy=fEntropy;
				minFactor=slimFactorSet[idIter->first];
			}
		}
		if(minFactor==NULL)
		{
			cout <<"No factor for sign " << fSetIter->first.c_str() << endl;
			return -1;
		}
		groupFactorRepMap[fSetIter->first]=minFactor->fId;
		oFile << fSetIter->first.c_str() <<"\t" << fSetIter->second->size() << "\t" << minFactor->fId << "\t";
		minFactor->showFactor(oFile,variableSet,false);
		for(INTINTMAP_ITER mbIter=minFactor->mergedMB.begin();mbIter!=minFactor->mergedMB.end();mbIter++)
		{
			if(mbIter!=minFactor->mergedMB.begin())
			{
				oFile <<"-";
			}
			else
			{
				oFile <<"\t";
			}
			oFile << variableSet[mbIter->first]->getName();
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


//This function must be called only after we have estimated the information theoretic 
//variables: joint entropy and multiinformation
int 
FactorManager::generateClusters_Apriori(double confidence, double eps_support,int maxClusterSize, bool remDup)
{
	//Get rid of the bad clusters/factors
	filterClustersWithMI(eps_support);
	if(maxFactorSize_Approx>maxFactorSize)
	{
		generateApproximateClusters(eps_support,confidence);
	}
	if(!checkMonotonicity())
	{
		cout <<"Monotonic property of factor ids is violated. Quitting! " << endl;
		return -1;	
	}
	//Now find the best Markov blankets using the greedy approach and also merge multiple Markov blankets together
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		if(getBestMarkovBlanket(fIter->second,eps_support,confidence)==-1)
		{
			return -1;
		}
		//The moment we reach a factor which is of size maxFactorSize, there is no point
		//in continuing because we can't exactly find its Markov blanket.
		if(fIter->second->vCnt==maxFactorSize_Approx)
		{
			break;
		}
	}
	
	VSET& variableSet=vMgr->getVariableSet();
	char aFName[1024];
	sprintf(aFName,"%s/mbsize_dist_precheck%d.txt",outputDir,maxFactorSize_Approx);
	ofstream checkFile(aFName);
	INTINTMAP mbSizeDist;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt==maxFactorSize_Approx)
		{
			break;
		}
		int asize=sFactor->mergedMB.size();
		if(mbSizeDist.find(asize)==mbSizeDist.end())
		{
			mbSizeDist[asize]=1;
		}	
		else
		{
			mbSizeDist[asize]=mbSizeDist[asize]+1;
		}
		sFactor->showFactor(checkFile,variableSet,false);
		checkFile << "\t" << sFactor->vCnt <<"\t" << asize<<"\t";
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			if(mbIter!=sFactor->mergedMB.begin())
			{
				checkFile <<"-";
			}
			checkFile << variableSet[mbIter->first]->getName();
		}
		checkFile << endl;
	}
	checkFile << "MBSize\tFactorCnt"  << endl;
	for(INTINTMAP_ITER sIter= mbSizeDist.begin();sIter!=mbSizeDist.end(); sIter++)
	{
		checkFile << sIter->first <<"\t " << sIter->second << endl;
	}
	checkFile.close();
	mbSizeDist.clear();
	//Make Markov blankets consistent
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		//No need to make Markov blankets of factors of size maxFactorSize consistent.
		if(fIter->second->vCnt==(maxFactorSize_Approx-1))
		{
			break;
		}
		//This function starts from the smallest superset and makes it consistent
		//with all the children of this superset
		//makeMarkovBlanketConsistent(fIter->second);
	}
	sprintf(aFName,"%s/mbsize_dist_%d.txt",outputDir,maxFactorSize_Approx);
	ofstream oFile(aFName);
	oFile << "fId\tfVars\tsize\tmbsize\tmb" << endl;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt==maxFactorSize_Approx)
		{
			break;
		}
		int asize=sFactor->mergedMB.size();
		if(mbSizeDist.find(asize)==mbSizeDist.end())
		{
			mbSizeDist[asize]=1;
		}
		else
		{
			mbSizeDist[asize]=mbSizeDist[asize]+1;
		}
		oFile <<sFactor->fId<<"\t";
		sFactor->showFactor(oFile,variableSet,false);
		oFile << "\t" << sFactor->vCnt <<"\t" << asize <<"\t";
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			if(mbIter!=sFactor->mergedMB.begin())
			{
				oFile <<"-";
			}
			oFile << variableSet[mbIter->first]->getName();
		}
		oFile << endl;
	}
	oFile  << "MBSize\tFactorCnt"  << endl;
	for(INTINTMAP_ITER sIter= mbSizeDist.begin();sIter!=mbSizeDist.end(); sIter++)
	{
		oFile << sIter->first <<"\t " << sIter->second << endl;
	}
	oFile.close();
	if(removeDupFactors()==-1)
	{
		return -1;
	}
	if(remDup)
	{
		produceClusters_NoDup(confidence,maxClusterSize);
	}
	else
	{
		produceClusters(confidence,maxClusterSize);
	}
	return 0;
}


int 
FactorManager::learnMBStructure(double eps_support,double confidence)
{
	//Get rid of the bad clusters/factors
	filterClustersWithMI(eps_support);
	if(!checkMonotonicity())
	{
		cout <<"Monotonic property of factor ids is violated. Quitting! " << endl;
		return -1;	
	}
	//Now find the best Markov blankets using the greedy approach and also merge multiple Markov blankets together
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		if(fIter->second->vCnt>1)
		{
			break;
		}
		if(getBestMarkovBlanket(fIter->second,eps_support,confidence)==-1)
		{
			return -1;
		}
	}
	makeMBMutuallyConsistent();
	VSET& variableSet=vMgr->getVariableSet();
	INTINTMAP mbSizeDist;
	char aFName[1024];
	sprintf(aFName,"%s/mbsize_dist_%d.txt",outputDir,maxFactorSize_Approx);
	ofstream oFile(aFName);
	oFile << "fId\tfVars\tsize\tmbsize\tmb" << endl;
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		int asize=sFactor->mergedMB.size();
		if(mbSizeDist.find(asize)==mbSizeDist.end())
		{
			mbSizeDist[asize]=1;
		}
		else
		{
			mbSizeDist[asize]=mbSizeDist[asize]+1;
		}
		oFile <<sFactor->fId<<"\t";
		sFactor->showFactor(oFile,variableSet,false);
		oFile << "\t" << sFactor->vCnt <<"\t" << asize <<"\t";
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			if(mbIter!=sFactor->mergedMB.begin())
			{
				oFile <<"-";
			}
			oFile << variableSet[mbIter->first]->getName();
		}
		oFile << endl;
	}
	oFile  << "MBSize\tFactorCnt"  << endl;
	for(INTINTMAP_ITER sIter= mbSizeDist.begin();sIter!=mbSizeDist.end(); sIter++)
	{
		oFile << sIter->first <<"\t " << sIter->second << endl;
	}
	oFile.close();

	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		if(fIter->second->vCnt>1)
		{
			break;
		}
		canonicalFactorSet[fIter->first]=fIter->second;
	}
	generateCanonicalFactors();
	//getPseudoLikelihood();
	setBaseInstantiation();
	estimateCanonicalParameters("paramLearn_approx");
	double dll=getLikelihood();
	cout <<"Data likelihood canonical" << dll << endl;
	double dll_exact=getLikelihood_ChainRule();
	cout <<"Data likelihood exact " << dll_exact << endl;
	evaluateMarkovBlanket(eps_support);
	return 0;
}


double
FactorManager::getPseudoLikelihood()
{
	double pseudoLL=0;
	VSET& variableSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		if(sFactor->goodMBIDs.size()==0)
		{
			continue;
		}
		double factorLL=potMgr->getPseudoLikelihood(sFactor,variableSet);	
		pseudoLL=pseudoLL+factorLL;
	}
	cout <<"Pseudolikelihood: " << pseudoLL << endl;
	return pseudoLL;
}

double
FactorManager::getLikelihood_ChainRule()
{
	VSET& variableSet=vMgr->getVariableSet();
	double dataLL=0;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactorSet.begin();aIter!=canonicalFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		double factorLL=potMgr->getLikelihood(sFactor,variableSet);
		dataLL=factorLL+dataLL;
	}

	return dataLL;
}

double 
FactorManager::getLikelihood_ChainRule(FactorGraph* fg)
{
	VSET& variableSet=vMgr->getVariableSet();
	double dataLL=0;
	int factorCnt=fg->getFactorCnt();
	struct timeval begintime;
	gettimeofday(&begintime,NULL);
	
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r,begintime.tv_usec);
	vector<double> lls;
	double meanLL=0;
	for(int j=0;j<50;j++)
	{
		double step=1.0/(double)factorCnt;
		map<int,int> usedInit;
		vector<int> fIds;
		dataLL=0;
		for(int i=0;i<factorCnt;i++)
		{
			double rVal=gsl_ran_flat(r,0,1);
			int rind=(int)(rVal/step);
			while(usedInit.find(rind)!=usedInit.end())
			{
				rVal=gsl_ran_flat(r,0,1);
				rind=(int)(rVal/step);
			}
			usedInit[rind]=0;
			fIds.push_back(rind);
		}
		map<int,int> visitedVertices;
		for(int i=0;i<fIds.size();i++)
		{
			SlimFactor* sFactor=fg->getFactorAt(fIds[i]);
			double factorLL=potMgr->getLikelihood(sFactor,variableSet,visitedVertices);
			dataLL=factorLL+dataLL;
			visitedVertices[sFactor->fId]=0;
		}
		meanLL=meanLL+dataLL;
		lls.push_back(dataLL);
	}
	meanLL=meanLL/((double) lls.size());
	double stdLL=0;
	for(int i=0;i<lls.size();i++)
	{
		double diff=meanLL-lls[i];
		stdLL=stdLL+(diff*diff);
	}
	stdLL=stdLL/(double(lls.size()-1));
	stdLL=sqrt(stdLL);
	cout <<"Mean: " << meanLL <<" Std: " << stdLL << endl;
	dataLL=meanLL;
	gsl_rng_free(r);
	return dataLL;
	
}

double
FactorManager::getLikelihood()
{
	double dataLL=0;
	int evidCnt=evMgr->getNumberOfEvidences();
	INTDBLMAP potScore;
	VSET& varSet=vMgr->getVariableSet();
	int maxCanSize=0;
	int incorrectLL=0;
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		string aConfStr;
		char confStr[CONSTR_LEN];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
			aConfStr.append(confStr);
		}
		double dataPtll=0;
		for(map<int,SlimFactor*>::iterator fIter=canonicalFactorSet.begin();fIter!=canonicalFactorSet.end();fIter++)
		{
			SlimFactor* sFactor=fIter->second;
			if(sFactor->vCnt==1)
			{
				if(sFactor->mergedMB.size()>maxCanSize)
				{
					maxCanSize=sFactor->mergedMB.size();
				}
			}
			string confKey;
			for(int v=sFactor->vCnt-1;v>=0;v--)
			{
				int vId=sFactor->vIds[v];
				Evidence* evid=(*evidMap)[vId];
				char confStr[CONSTR_LEN];
				if(evid->getType()!=Evidence::HARD)
				{
					cout <<"Handling only hard evidence now " << endl;
					return -1;
				}
				sprintf(confStr,"-%d=%d-",vId,evid->getHardEvidVal());	
				confKey.append(confStr);
			}
			if(sFactor->canonicalParams==NULL)
			{
				cout <<"No canonical params for factor " << endl;
				sFactor->showFactor(cout,varSet);
				return -1;
			}
			double pVal=sFactor->canonicalParams->getJointPotValueForConf(confKey);
			if(pVal==-1)
			{
				cout <<"No value for " << confKey.c_str() << endl;
				return -1;
			}
			if(potScore.find(fIter->first)==potScore.end())
			{
				potScore[fIter->first]=log(pVal);
			}
			else
			{
				potScore[fIter->first]=potScore[fIter->first]+log(pVal);
			}
			dataPtll=dataPtll+log(pVal);
		}
		dataPtll=dataPtll+log(defProb);
		if(dataPtll>0)
		{
			incorrectLL++;
			cout <<"P value " << dataPtll << " for data point " << i << " is greater than 1 " << endl;
			cout << aConfStr.c_str() << endl;
		}
		dataLL=dataLL+dataPtll;
	}
	cout <<"Total Number of data points with log likelihood > 0 " << incorrectLL << endl;
	char aFName[1024];
	sprintf(aFName,"%s/canonical_score_k%d.txt",outputDir,maxCanSize);
	ofstream oFile(aFName);
	for(INTDBLMAP_ITER iIter=potScore.begin();iIter!=potScore.end();iIter++)
	{
		oFile <<iIter->first << "\t" << factorIDToNameMap[iIter->first] <<"\t" <<iIter->second << "\t";
		SlimFactor* sFactor=canonicalFactorSet[iIter->first];
		for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
		{
			if(mIter!=sFactor->mergedMB.begin())
			{
				oFile <<"-";
			}
			oFile << mIter->first;
		}
		oFile << endl;
	}

	oFile.close();

	return dataLL;
}


double
FactorManager::getLikelihood(FactorGraph* fg)
{
	double dataLL=0;
	INTDBLMAP potScore;
	VSET& varSet=vMgr->getVariableSet();
	int maxCanSize=0;
	int incorrectLL=0;
	//setBaseInstantiation_Variable();
	setBaseInstantiation();
	map<int,SlimFactor*> canonicalFactors;
	generateCanonicalFactors(fg,canonicalFactors);
	potMgr->resetPotFuncs();
	map<string,int> canFactOrder;
	int maxCanFactSize=0;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		string key;
		getFactorKey(sFactor->vIds,sFactor->vCnt,key);
		if(maxCanFactSize<sFactor->vCnt)
		{
			maxCanFactSize=sFactor->vCnt;
		}
		canFactOrder[key]=aIter->first;
	}
	double den=(double) pow(2.0,(double)maxCanFactSize+2);
	double canThreshold=0.1/den;
	int trivialFactors=0;
	//for(map<string,int>::iterator aIter=canFactOrder.begin();aIter!=canFactOrder.end();aIter++)
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		INTINTMAP allSubsets;
		lattice.getAllSubsets(aIter->first,allSubsets);
		allSubsets[sFactor->fId]=0;
		potMgr->estimateCanonicalPotential_Approximate(sFactor,varSet,defaultInstantiation,allSubsets,canonicalFactors);
		//cout <<"Canonical potential for: " << sFactor->fId << endl;
		//sFactor->canonicalParams->dumpPotential(cout);
		sFactor->thresholdToOne(canThreshold);
		if(sFactor->allEntriesInsignificant())
		{
			trivialFactors++;
		}
		//potMgr->estimateCanonicalPotential_Abbeel(sFactor,varSet,defaultInstantiation,allSubsets,canonicalFactors);
		//potMgr->estimateCanonicalPotential(sFactor,varSet,defaultInstantiation,allSubsets,canonicalFactors);
	}
	cout <<"Total number of trivial factors " << trivialFactors << endl;
	int evidCnt=evMgr->getNumberOfEvidences();
	double actualDefProb=0;
	for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		string confKey;
		char confStr[CONSTR_LEN];
		for(int v=sFactor->vCnt-1;v>=0;v--)
		{
			int vId=sFactor->vIds[v];
			int vVal=defInstMap[vId];
			sprintf(confStr,"-%d=%d-",vId,vVal);
			confKey.append(confStr);
		}
		if(sFactor->canonicalParams==NULL)
		{
			cout <<"No canonical params for factor " << endl;
			sFactor->showFactor(cout,varSet);
			return -1;
		}
		double pVal=sFactor->canonicalParams->getJointPotValueForConf(confKey);
		if(pVal==-1)
		{
			cout <<"No value for " << confKey.c_str() << endl;
			return -1;
		}
		if(pVal!=1.0)
		{
			cout <<"Bad pvalue " << pVal << " for potential " << confKey.c_str() << endl;
		}
		actualDefProb=actualDefProb+log(pVal);
	}
	double maxUnnormll=0;
	int maxevid=0;
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		string aConfStr;
		char confStr[CONSTR_LEN];
		int match=0;
		int mismatch=0;

		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sprintf(confStr,"-%d=%d-",eIter->first,eIter->second->getHardEvidVal());
			aConfStr.append(confStr);
			if(eIter->second->getHardEvidVal()==defInstMap[eIter->first])
			{
				match++;
			}
			else
			{
				mismatch++;
			}
		}
		double dataPtll=0;
		for(map<int,SlimFactor*>::iterator aIter=canonicalFactors.begin();aIter!=canonicalFactors.end();aIter++)
		{
			SlimFactor* sFactor=aIter->second;
			string confKey;
			for(int v=sFactor->vCnt-1;v>=0;v--)
			{
				int vId=sFactor->vIds[v];
				Evidence* evid=(*evidMap)[vId];
				char confStr[CONSTR_LEN];
				if(evid->getType()!=Evidence::HARD)
				{
					cout <<"Handling only hard evidence now " << endl;
					return -1;
				}
				sprintf(confStr,"-%d=%d-",vId,evid->getHardEvidVal());	
				confKey.append(confStr);
			}
			if(sFactor->canonicalParams==NULL)
			{
				cout <<"No canonical params for factor " << endl;
				sFactor->showFactor(cout,varSet);
				return -1;
			}
			double pVal=sFactor->canonicalParams->getJointPotValueForConf(confKey);
			if(pVal==-1)
			{
				cout <<"No value for " << confKey.c_str() << endl;
				return -1;
			}
			if(potScore.find(aIter->first)==potScore.end())
			{
				potScore[aIter->first]=log(pVal);
			}
			else
			{
				potScore[aIter->first]=potScore[aIter->first]+log(pVal);
			}
			dataPtll=dataPtll+log(pVal);
		}
		if(dataPtll>maxUnnormll)
		{
			maxUnnormll=dataPtll;
			maxevid=i;
		}
		dataPtll=dataPtll+log(defProb);
		if(dataPtll>0)
		{
			incorrectLL++;
		//	cout << "Pvalue>1\t"<<dataPtll <<"\t" << mismatch <<endl;
		//	cout <<"P value " << dataPtll << " for data point " << i << " is greater than 1 " << endl;
		//	cout << aConfStr.c_str() << endl;
		}
		else
		{
		//	cout << "Pvalue<=1\t"<<dataPtll <<"\t" << mismatch <<endl;
		}
		dataLL=dataLL+dataPtll;
	}
	cout <<"Max unnormalized data ll " << maxUnnormll << " for evidence " << maxevid << endl;
	cout <<"Total Number of data points with log likelihood > 0 " << incorrectLL << endl;
	//Corrected datall would be if the -maxUnnormll=log(defProb). 
	//So subtract for dataLL evidCnt*log(defProb), and then 
	//add evidCnt*(-maxUnnormll)
	if(maxUnnormll>fabs(log(defProb)))
	{
		cout <<"Applying correction" << endl;
		dataLL=dataLL-(evidCnt*log(defProb));
		dataLL=dataLL-(evidCnt*maxUnnormll);
	}

	return dataLL;
}


double 
FactorManager::getLikelihood_MCMC(FactorGraph* fg)
{
	double dataLL=0;
	potMgr->resetPotFuncs();
	double z=estimatePartitionFunction(fg);
	if(z==-1)
	{
		return -1;
	}
	int incorrectLL=0;
	VSET& varSet=vMgr->getVariableSet();
	int evidCnt=evMgr->getNumberOfEvidences();
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		INTINTMAP sample;
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sample[eIter->first]=eIter->second->getHardEvidVal();
		}
		double dataptLL=potMgr->getSampleLikelihood(fg->getAllFactors(),varSet,&sample);
		z=z+exp(dataptLL);
	}
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=evMgr->getEvidenceAt(i);
		INTINTMAP sample;
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			sample[eIter->first]=eIter->second->getHardEvidVal();
		}
		double dataptLL=potMgr->getSampleLikelihood(fg->getAllFactors(),varSet,&sample);
		dataptLL=dataptLL-log(z);
		if(dataptLL>0)
		{
			incorrectLL++;
		}
		dataLL=dataptLL+dataLL;
	}
	cout <<"Total Number of data points with log likelihood (MCMC) > 0 " << incorrectLL << endl;

	return dataLL;
}

double 
FactorManager::estimatePartitionFunction(FactorGraph* fg)
{
	double z=0;
	int burnIn=50000;
	int sampleCnt=200000;
	int sampleID=0;
	struct timeval begintime;
	gettimeofday(&begintime,NULL);
	
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r,begintime.tv_usec);

	INTINTMAP* firstSample=new INTINTMAP;
	INTINTMAP* nextSample=new INTINTMAP;
	INTINTMAP* currSample=firstSample;
	randInitSample(currSample,r);
	VSET& varSet=vMgr->getVariableSet();
	while(sampleID <= sampleCnt)
	{
		if(sampleID>burnIn)
		{
			double sampleLikelihood=potMgr->getSampleLikelihood(fg->getAllFactors(),varSet,currSample);
			if(sampleLikelihood==-1)
			{
				return -1;
			}
			z=z+exp(sampleLikelihood);
		}
		sampleID++;
		if(getNextSample(fg,currSample,nextSample,r)==-1)
		{
			return -1;
		}
		currSample=nextSample;
	}
	return z;
}

int
FactorManager::randInitSample(INTINTMAP* sample,gsl_rng* r)
{
	VSET& variableSet=vMgr->getVariableSet();
	
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		INTVECT& values=vIter->second->getValues();
		double step=1.0/(double) values.size();
		int valInd=0;
		double total=step;
		double rVal=gsl_ran_flat(r,0,1);
		while((rVal > total) && (valInd<values.size()))
		{
			total=total+step;
			valInd++;
		}
		if(valInd==values.size())
		{
			cout <<"Could not find sample value " << endl;
			return -1;
		}
		int vval=values[valInd];
		(*sample)[vIter->first]=vval;
	}
	return 0;
}


int 
FactorManager::getNextSample(FactorGraph* fg,INTINTMAP* currSample,INTINTMAP* nextSample,gsl_rng* r)
{
	map<int,SlimFactor*>& factorSet=fg->getAllFactors();
	VSET& varSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		int varSample=potMgr->getVariableSample(*currSample,varSet,fIter->first,fIter->second,r);
		if(varSample==-1)
		{
			return -1;
		}
		(*nextSample)[fIter->first]=varSample;
	}
	return 0;
}

int
FactorManager::evaluateMarkovBlanket(double eps_support)
{
	char aFName[1024];
	VSET& varSet=vMgr->getVariableSet();
	potMgr->estimateMarginalEntropies(canonicalFactorSet,varSet,false);
	if(eps_support!=-1)
	{
		sprintf(aFName,"%s/mbscore_score_n%d_x%d_e%.1f.txt",outputDir,beamSize,maxFactorSize_Approx,eps_support);
	}
	else
	{
		sprintf(aFName,"%s/mbscore_score_n%d_x%d.txt",outputDir,beamSize,maxFactorSize_Approx);
	}
	ofstream oFile(aFName);

	oFile<< "Var\tMB\tCondEntr\tMargEntr" << endl;
	for(map<int,SlimFactor*>::iterator fIter=canonicalFactorSet.begin();fIter!=canonicalFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		if(sFactor->mbScore==-1)
		{
			Potential* apot=potMgr->getPotential(sFactor->fId);
			if(apot==NULL)
			{
				cout <<"No potential estimated for factor " << sFactor->fId << endl;
				return -1;
			}
			apot->calculateConditionalEntropy();
			sFactor->mbScore=apot->getConditionalEntropy();
		}
		oFile <<varSet[sFactor->vIds[0]]->getName() <<"\t";
                for(INTINTMAP_ITER vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++)
                {
                        if(vIter!=sFactor->mergedMB.begin())
                        {
                                oFile <<"-";
                        }
                        oFile<< varSet[vIter->first]->getName();
                }
                oFile <<"\t" <<sFactor->mbScore <<"\t" <<sFactor->conditionalEntropy << endl;
	}
	/*oFile <<"Factor\tFactorStrID\tAllVars\tConditionalEntr\tJointEntr\tMutualInfo" << endl;
	for(map<int,SlimFactor*>::iterator fIter=canonicalFactorSet.begin();fIter!=canonicalFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		potMgr->populateFactor(canonicalFactorSet,varSet,sFactor,false);
		showFactor(sFactor,oFile,false);
		oFile <<"\t";
		for(int i=0;i<sFactor->vCnt;i++)
		{
			oFile <<"-"<<sFactor->vIds[i];
		}
		oFile <<"\t"<< varSet[sFactor->vIds[0]]->getName();
		for(INTINTMAP_ITER vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++)
		{
			oFile <<"-"<< varSet[vIter->first]->getName();
		}
		oFile<< "\t" << sFactor->mbScore 
			<< "\t" << sFactor->jointEntropy
			<<"\t" << sFactor->mutualInfo << endl;
	}*/
	oFile.close();
	return 0;
}

//Here we generate clusters of size greater than max factor size upto the value for which 
//we can compute random mutual information exactly

int
FactorManager::generateApproximateClusters(double eps, double confThresh)
{
	map<int,int> parentIDs;
	VSET& variableSet=vMgr->getVariableSet();
	//get all the parent factor ids.
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt<=2)
		{
			sFactor->confidence=1;
		}
		else 
		{
			INTINTMAP* subsets=lattice.getSubsets(fIter->first);
			if(subsets==NULL)
			{
				cout <<"No subsets for " << fIter->first << endl;
				continue;
			}
			if(subsets->size()==0)
			{
				continue;
			}
			//get all subsets and check for confidence
			//We need to check the last sFactor->vCnt subsets as these
			//will be the ones that correspond to the maximal subsets
			map<int,int>::reverse_iterator rIter=subsets->rbegin();
			int sscnt=0;
			double hitConf=0;
			//However check for subsets in slimFactorSet
			SlimFactor* aSubset=slimFactorSet[rIter->first];
			while((sscnt<sFactor->vCnt) && ((sFactor->vCnt-aSubset->vCnt)==1 ))
			{
				int ssId=rIter->first;
				map<int,SlimFactor*>::iterator ssIter=slimFactorSet.find(ssId);
				if(ssIter!=slimFactorSet.end())
				{
					hitConf=hitConf+ssIter->second->confidence;
				}
				sscnt++;
				rIter++;
				aSubset=slimFactorSet[rIter->first];
			}
			double conf=hitConf/((double) sFactor->vCnt);
			sFactor->confidence=conf;
		}
		if(fIter->second->vCnt==maxFactorSize)
		{
			if(sFactor->confidence >=confThresh)
			{
				parentIDs[fIter->first]=0;
			}
		}
	}
	int fSize=maxFactorSize+1;
	int oldCnt=slimFactorSet.size();
	while(fSize<=maxFactorSize_Approx)
	{
		double randmi_mean=randMI_mean[fSize];
		double randmi_sd=randMI_std[fSize];
		map<int,int> tempParentIDs;
		for(map<int,int>::iterator pIter=parentIDs.begin();pIter!=parentIDs.end();pIter++)
		{
			SlimFactor* pFactor=slimFactorSet[pIter->first];
			generateNextLevelClusters(pFactor,eps,confThresh,tempParentIDs);
			if((tempParentIDs.size() % 100)==0)
			{
				cout <<"Populated " << tempParentIDs.size() << " parents " << endl;
			}
		}
		parentIDs.clear();
		for(map<int,int>::iterator pIter=tempParentIDs.begin();pIter!=tempParentIDs.end();pIter++)
		{
			parentIDs[pIter->first]=0;
		}
		tempParentIDs.clear();
		fSize++;
		cout <<"Populated " << slimFactorSet.size()-oldCnt<< " new factors "<< endl; 
		oldCnt=slimFactorSet.size();
	}

	return 0;
}

//This function generates all k+1 possible clusters from pFactor, where k is the size of this cluster
int
FactorManager::generateNextLevelClusters(SlimFactor* pFactor,double eps_support,double confThresh,INTINTMAP& newParentIDs)
{
	int** subsetSpace=new int*[maxFactorSize_Approx];
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		subsetSpace[i]=new int[maxFactorSize_Approx-1];
	}
	int* ssIds=new int[maxFactorSize_Approx];
	VSET& variableSet=vMgr->getVariableSet();
	double randmi_mean=randMI_mean[pFactor->vCnt+1];
	double randmi_sd=randMI_std[pFactor->vCnt+1];
	double randmi_mean_prev=randMI_mean[pFactor->vCnt];
	double randmi_sd_prev=randMI_std[pFactor->vCnt];
	for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		if(pFactor->isMemberVariable(vIter->first))
		{
			continue;
		}
		int newVId=vIter->first;
		//Propose a factor
		SlimFactor* sFactor=new SlimFactor;
		sFactor->vIds=new int[pFactor->vCnt+1];
		sFactor->vCnt=pFactor->vCnt+1;
		int fIter=0;
		int dIter=0;
		while((pFactor->vIds[fIter]<newVId) && (fIter!=pFactor->vCnt))
		{
			sFactor->vIds[dIter]=pFactor->vIds[fIter];
			dIter++;
			fIter++;
		}
		sFactor->vIds[dIter]=newVId;
		dIter++;
		while(fIter<pFactor->vCnt)
		{
			sFactor->vIds[dIter]=pFactor->vIds[fIter];
			fIter++;
			dIter++;
		}
		//This factor might have been created
		string fKey;
		getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
		if(strcmp(fKey.c_str(),"-9-10-11-33-34")==0)
		{
			cout <<"Reached "<< fKey.c_str()<< endl;
		}
		if(strcmp(fKey.c_str(),"-9-11-33-34")==0)
		{
			cout <<"Reached "<< fKey.c_str()<< endl;
		}
		if(factorNameToIDMap.find(fKey)!=factorNameToIDMap.end())
		{
			delete sFactor;
			continue;
		}
		sFactor->generateMaximalSubsets(subsetSpace);
		//Now check the confidence and in the meantime store the subset ids
		//to update the lattice structure. 
		double hitConf=0;
		int foundSsets=0;
		for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
		{
			string ssKey;
			getFactorKey(subsetSpace[sscnt],sFactor->vCnt-1,ssKey);
			if(factorNameToIDMap.find(ssKey)==factorNameToIDMap.end())
			{
				continue;
			}
			int sId=factorNameToIDMap[ssKey];
			ssIds[foundSsets]=sId;
			foundSsets++;
			map<int,SlimFactor*>::iterator ssIter=slimFactorSet.find(sId);
			if(ssIter!=slimFactorSet.end())
			{
				hitConf=hitConf+ssIter->second->confidence;
			}
		}
		hitConf=hitConf/sFactor->vCnt;
		if(hitConf<confThresh)
		{
		//	delete sFactor;
		//	continue;
		}
		sFactor->confidence=hitConf;
		potMgr->populateFactor(slimFactorSet,variableSet,sFactor,false);
		if(sFactor->mutualInfo < (randmi_mean+(eps_support*randmi_sd)))
		{
			delete sFactor;
			continue;
		}
		sFactor->fId=globalFactorID;
		//Check the confidence of this factor
		factorNameToIDMap[fKey]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=fKey;
		sFactor->secondPId=pFactor->fId;
		globalFactorID++;
		if(globalFactorID<0)
		{
			cout <<"Global factor id became negative !! " << endl;
		}
		for(int i=0;i<foundSsets;i++)
		{
			int subsetId=ssIds[i];
			if(slimFactorSet.find(subsetId)!=slimFactorSet.end())
			{
				lattice.addSubset(subsetId,sFactor->fId);
				lattice.addSuperset(sFactor->fId,subsetId);
			}
		}
		slimFactorSet[sFactor->fId]=sFactor;
		newParentIDs[sFactor->fId]=0;
	}
	if(newParentIDs.size()>(variableSet.size()-pFactor->vCnt))
	{
		cout <<"Too many ("<< newParentIDs.size() <<") one-variable extensions of " << pFactor->fId << endl;
		return -1;
	}
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		delete[] subsetSpace[i];
	}
	delete[] ssIds;
	delete[] subsetSpace;
	return 0;
}


//This function generates all k+1 possible clusters from pFactor, where k is the size of this cluster
int
FactorManager::generateNextLevelClusters(SlimFactor* pFactor,double eps_support, INTDBLMAP& candidateExtensions)
{
	int** subsetSpace=new int*[maxFactorSize_Approx];
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		subsetSpace[i]=new int[maxFactorSize_Approx-1];
	}
	int* ssIds=new int[maxFactorSize_Approx];
	VSET& variableSet=vMgr->getVariableSet();
	double randmi_mean=randMI_mean[pFactor->vCnt+1];
	double randmi_sd=randMI_std[pFactor->vCnt+1];
	for(INTDBLMAP_ITER vIter=candidateExtensions.begin();vIter!=candidateExtensions.end();vIter++)
	{
		if(pFactor->isMemberVariable(vIter->first))
		{
			continue;
		}
		int newVId=vIter->first;
		//Propose a factor
		SlimFactor* sFactor=new SlimFactor;
		sFactor->vIds=new int[pFactor->vCnt+1];
		sFactor->vCnt=pFactor->vCnt+1;
		int fIter=0;
		int dIter=0;
		while((pFactor->vIds[fIter]<newVId) && (fIter!=pFactor->vCnt))
		{
			sFactor->vIds[dIter]=pFactor->vIds[fIter];
			dIter++;
			fIter++;
		}
		sFactor->vIds[dIter]=newVId;
		dIter++;
		while(fIter<pFactor->vCnt)
		{
			sFactor->vIds[dIter]=pFactor->vIds[fIter];
			fIter++;
			dIter++;
		}
		//This factor might have been created
		string fKey;
		getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
		if(strcmp(fKey.c_str(),"-9-11-33-34")==0)
		{
	                cout << " Reached " << fKey.c_str() << endl;
		}
		if(strcmp(fKey.c_str(),"-9-10-11-33-34")==0)
		{
	                cout << " Reached " << fKey.c_str() << endl;
		}
		if(factorNameToIDMap.find(fKey)!=factorNameToIDMap.end())
		{
			delete sFactor;
			continue;
		}
		sFactor->generateMaximalSubsets(subsetSpace);
		int foundSsets=0;
		for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
		{
			string ssKey;
			getFactorKey(subsetSpace[sscnt],sFactor->vCnt-1,ssKey);
			if(factorNameToIDMap.find(ssKey)==factorNameToIDMap.end())
			{
				continue;
			}
			int sId=factorNameToIDMap[ssKey];
			ssIds[foundSsets]=sId;
			foundSsets++;
		}

		potMgr->populateFactor(slimFactorSet,variableSet,sFactor,false);
		if(sFactor->mutualInfo < (randmi_mean+(eps_support*randmi_sd)))
		{
			delete sFactor;
			continue;
		}
		sFactor->fId=globalFactorID;
		//Check the confidence of this factor
		factorNameToIDMap[fKey]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=fKey;
		sFactor->secondPId=pFactor->fId;
		globalFactorID++;
		for(int i=0;i<foundSsets;i++)
		{
			int subsetId=ssIds[i];
			if(slimFactorSet.find(subsetId)!=slimFactorSet.end())
			{
				lattice.addSubset(subsetId,sFactor->fId);
				lattice.addSuperset(sFactor->fId,subsetId);
			}
		}
		slimFactorSet[sFactor->fId]=sFactor;
	}
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		delete[] subsetSpace[i];
	}
	delete[] ssIds;
	delete[] subsetSpace;
	return 0;
}

int
FactorManager::generateNextLevelClusters(SlimFactor* sFactor)
{
	INTINTMAP newExtensions;
	generateNextLevelClusters(sFactor,misdCnt,1.0,newExtensions);
	return 0;
}

int
FactorManager::generateNextLevelClusters(SlimFactor* sFactor,SlimFactor* mbFactor)
{
	INTDBLMAP candidateExtensions;
	VSET& variableSet=vMgr->getVariableSet();
	VSET_ITER vIter=variableSet.find(sFactor->fId);
	vIter++;
	for(;vIter!=variableSet.end();vIter++)
	{
		candidateExtensions[vIter->first]=0;	
	}
	generateNextLevelClusters(mbFactor,misdCnt,candidateExtensions);
	return 0;
}


int 
FactorManager::showConditionalPotentials(FactorGraph* fg)
{
	char aFName[1024];
	sprintf(aFName,"%s/condpot_k%d",outputDir,maxFactorSize_Approx);
	ofstream oFile(aFName);
	map<int,SlimFactor*>& factorSet=fg->getAllFactors();
	VSET& varSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		Potential* apot=new Potential;
		Variable* currVar=varSet[sFactor->fId];
		apot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
		for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
		{
			Variable* aVar=varSet[mIter->first];
			apot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		apot->potZeroInit();
		potMgr->populatePotential(apot,false);
		apot->makeValidJPD();
		oFile <<"Potential for " << currVar->getName() << "\t" << currVar->getID() << endl;
		apot->dumpPotential(oFile);
		delete apot;
	}
	oFile.close();
	return 0;
}

int
FactorManager::showAllFactors(double eps)
{
	char aFName[1024];
	//sprintf(aFName,"%s/mifilter_%.3ftxt",outputDir,eps);
	sprintf(aFName,"%s/mifilter_k%d.txt",outputDir,maxFactorSize);

	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		fIter->second->showFactor(oFile,vMgr->getVariableSet());
	}
	oFile.close();
	return 0;
}


bool 
FactorManager::checkMonotonicity()
{
	bool monotonic=true;
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		int predVCnt=aIter->second->vCnt;
		map<int,SlimFactor*>::iterator bIter=aIter;
		bIter++;
		for(;bIter!=slimFactorSet.end();bIter++)
		{
			int currVCnt=bIter->second->vCnt;
			if(currVCnt<predVCnt) 
			{
				monotonic=false;
				cout <<"Monotonic property violated at prev factor " << aIter->first <<
					" and succ factor " << bIter->first << endl;
			}
			break;
		}
		if(!monotonic)
		{
			break;
		}
	}
	return monotonic;
}



INTINTMAP*
FactorManager::getSupersets(int fId)
{
	return lattice.getSupersets(fId);
}


int
FactorManager::getSupersets(int fId, int level, INTINTMAP& superSetID)
{
	lattice.getSupersets(fId,superSetID,level);
	return 0;
}


double 
FactorManager::getMBScore(SlimFactor* sFactor,int supId)
{
	
	int* mbVars=new int[maxFactorSize_Approx];
	int mbvarCnt=0;
	SlimFactor* mbFactor=slimFactorSet[supId];
	//sFactor->getSetDiff(mbFactor,mbVars,mbvarCnt);
	//double mbInfo=getMIFromVars(mbVars,mbvarCnt);
	//For now we can safely assume that we are finding
	//the Markov blankets of only one variable
	mbvarCnt=mbFactor->vCnt-1;
	INTINTMAP fVars;
	for(int i=0;i<mbFactor->vCnt;i++)
	{
		fVars[mbFactor->vIds[i]]=0;
	}
	double proposedScore=potMgr->getConditionalEntropy(sFactor->fId,fVars,vMgr->getVariableSet());
	//proposedScore=proposedScore+mbInfo;
	double penalty=mbPenalty*log(double(mbvarCnt));
	proposedScore=proposedScore+penalty;
	/*if(proposedScore<0)
	{
		cout<<"Negative entropy for factor " << sFactor->fId << " mbfactor: "<< mbFactor->fId << " Setting to zero" << endl; 
		proposedScore=0;
		delete[] mbVars;
		return -1;
	}*/
	delete[] mbVars;
	return proposedScore;
}


double 
FactorManager::getMBScore(SlimFactor* sFactor)
{
	INTINTMAP fVars;
	fVars[sFactor->fId]=0;
	for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
	{
		fVars[mIter->first]=0;
	}
	double condEntropy=potMgr->getConditionalEntropy(sFactor->fId,fVars,vMgr->getVariableSet());
	double proposedScore=condEntropy;
	double penalty=mbPenalty*log(double(fVars.size()-1));
	proposedScore=proposedScore+penalty;
	return proposedScore;
}


SlimFactor*
FactorManager::getFactorAt(int fId)
{
	if(slimFactorSet.find(fId)==slimFactorSet.end())
	{
		return NULL;
	}
	return slimFactorSet[fId];
}


//Find the best Markov blankets for sFactor. 
//Populates the set of Markov blankets for sFactor
int
FactorManager::getBestMarkovBlanket(SlimFactor* sFactor,double eps_support, double confidence)
{
	VSET& variableSet=vMgr->getVariableSet();
	if(sFactor->fId==9)
	{
		cout <<"Stop here " << endl;
	}
	int currK=sFactor->vCnt+1;
	double marginalEntropy=0;
	for(int i=0;i<sFactor->vCnt;i++)
	{
		SlimFactor* vFactor=slimFactorSet[sFactor->vIds[i]];
		marginalEntropy=marginalEntropy+vFactor->marginalEntropy;
	}
	//Initialize to the sum of the marginal entropy of the variables in sFactor
	double minEntropy=marginalEntropy;
	while(currK<=maxFactorSize_Approx)
	{
		if(currK==sFactor->vCnt+1)
		{
			INTINTMAP* supersets=lattice.getSupersets(sFactor->fId);
			if(getGoodCandidateMarkovBlankets(sFactor,supersets,currK,minEntropy,marginalEntropy,sFactor->goodMBIDs)==-1)
			{
				return -1;
			}
		}
		//Implementing stage II exactly for markov blankets of size 3
		else if(currK<=maxFactorSize)
		{
			INTDBLMAP tempGoodMBIDs;
			//Get all supersets of sFactor at level 2
			INTINTMAP supersets;
			lattice.getSupersets(sFactor->fId,supersets,currK-sFactor->vCnt);
			if(getGoodCandidateMarkovBlankets(sFactor,&supersets,currK,minEntropy,marginalEntropy,tempGoodMBIDs)==-1)
			{
				return -1;
			}
			//else update sFactor's mbIds
			if(tempGoodMBIDs.size()>0)
			{
				sFactor->goodMBIDs.clear();
				for(INTDBLMAP_ITER aIter=tempGoodMBIDs.begin();aIter!=tempGoodMBIDs.end();aIter++)
				{
					if(aIter->second<=minEntropy)
					{
						sFactor->goodMBIDs[aIter->first]=minEntropy;
					}
				}
			}
		}
		else
		{
			INTDBLMAP tempGoodMBIDs;
			//This holds the id of Markov blankets that could be successfully grown
			//We must delete these ids from the set of good MBs
			INTINTMAP replaceID;
			for(INTDBLMAP_ITER mbIter=sFactor->goodMBIDs.begin();mbIter!=sFactor->goodMBIDs.end();mbIter++)
			{
				SlimFactor* currMBFactor=slimFactorSet[mbIter->first];
				if(currMBFactor->vCnt<(currK-1))
				{
					continue;
				}
				//INTINTMAP newExtensions;
				//generateNextLevelClusters(currMBFactor,eps_support,confidence,newExtensions);
				generateNextLevelClusters(currMBFactor,eps_support,sFactor->candidateNeighbours);
				INTINTMAP* supersets=lattice.getSupersets(mbIter->first);
				int mbCnt=tempGoodMBIDs.size();
				if(getGoodCandidateMarkovBlankets(sFactor,supersets,currK,minEntropy,marginalEntropy,tempGoodMBIDs)==-1)
				{
					return -1;
				}
				if(tempGoodMBIDs.size()>mbCnt)
				{
					replaceID[mbIter->first]=0;
				}
			}
			//Replace all entries from goodMBIDs which are in replaceID or are above the minEntropy threshold
			for(INTDBLMAP_ITER rIter=sFactor->goodMBIDs.begin();rIter!=sFactor->goodMBIDs.end();rIter++)
			{
				if(replaceID.find(rIter->first)!=replaceID.end())
				{
					sFactor->goodMBIDs.erase(rIter);
				}
				else if(rIter->second > minEntropy)
				{
					sFactor->goodMBIDs.erase(rIter);
				}
			}
			//Add the new MBs that have been created that have entropy equal to minEntropy
			for(INTDBLMAP_ITER aIter=tempGoodMBIDs.begin();aIter!=tempGoodMBIDs.end();aIter++)
			{
				if(aIter->second<=minEntropy)
				{
					sFactor->goodMBIDs[aIter->first]=minEntropy;
				}
			}
			tempGoodMBIDs.clear();
			replaceID.clear();
		}
		currK++;
	}
	
	//Finally merge the Markov blankets together
	SlimFactor* maxMbFactor=NULL;
	int maxMB=0;
	vector<int> mbIds;
	for(INTDBLMAP_ITER mbIter=sFactor->goodMBIDs.begin();mbIter!=sFactor->goodMBIDs.end();mbIter++)
	{
		/*SlimFactor* mbFactor=slimFactorSet[mbIter->first];
		for(int i=0;i<mbFactor->vCnt;i++)
		{
			if(!sFactor->isMemberVariable(mbFactor->vIds[i]))
			{
				sFactor->mergedMB[mbFactor->vIds[i]]=0;
			}
		}*/
		mbIds.push_back(mbIter->first);
	}
	for(int i=0;i<mbIds.size();i++)
	{
		for(int j=i+1;j<mbIds.size();j++)
		{
			if(sFactor->goodMBIDs[mbIds[i]]>sFactor->goodMBIDs[mbIds[j]])
			{
				int tempfid=mbIds[i];
				mbIds[i]=mbIds[j];
				mbIds[j]=tempfid;
			}
		}
	}
	int i=0;
	while(sFactor->mergedMB.size()<(maxFactorSize_Approx*beamSize))
	{
		if(i>=mbIds.size())
		{
			break;
		}
		SlimFactor* mbFactor=slimFactorSet[mbIds[i]];
		for(int j=0;j<mbFactor->vCnt;j++)
		{
			if(!sFactor->isMemberVariable(mbFactor->vIds[j]))
			{
				sFactor->mergedMB[mbFactor->vIds[j]]=0;
			}
			//Add variable of sFactor into the Markov blanket of mbFactor
		/*	SlimFactor* mbVarFactor=slimFactorSet[mbFactor->vIds[j]];
			for(int k=0;k<sFactor->vCnt;k++)
			{
				if(!mbVarFactor->isMemberVariable(sFactor->vIds[k]))
				{
					mbVarFactor->mergedMB[sFactor->vIds[k]]=0;
				}
			}*/
		}
		i++;
	}
	sFactor->mbScore=minEntropy;
	return 0;
}


int
FactorManager::makeMBMutuallyConsistent()
{
	VSET& varSet=vMgr->getVariableSet();
	for(map<int,SlimFactor*>::iterator aIter=slimFactorSet.begin();aIter!=slimFactorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		int totalAdd=0;
		int totalDel=0;
		if(sFactor->vCnt>1)
		{
			break;
		}
		if(sFactor->fId==4 || sFactor->fId==5)
		{
			cout << "Stop here " << endl;
		}
		INTINTMAP fVars;
		fVars[sFactor->fId]=0;
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			fVars[mbIter->first]=0;
		}
		double condEntropy=potMgr->getConditionalEntropy(sFactor->fId,fVars,varSet);
		sFactor->mbScore=condEntropy;
		INTINTMAP deleteMBVarID;
		for(INTINTMAP_ITER mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
		{
			if(mbIter->first<sFactor->fId)
			{
				continue;
			}
			SlimFactor* mbFactor=slimFactorSet[mbIter->first];
			if(mbFactor->mergedMB.find(sFactor->fId)!=mbFactor->mergedMB.end())
			{
				continue;
			}
			INTINTMAP mbVars;
			mbVars[mbFactor->fId]=0;
			for(INTINTMAP_ITER vIter=mbFactor->mergedMB.begin();vIter!=mbFactor->mergedMB.end();vIter++)
			{
				mbVars[vIter->first]=0;
			}
			double mbcondEntropy=potMgr->getConditionalEntropy(mbFactor->fId,mbVars,varSet);
			//Now we have to decide if its better to add sFactor's var to mbFactor or to delete mbFactor's 
			//var from sFactor
			INTINTMAP_ITER vIter=fVars.find(mbIter->first);
			fVars.erase(vIter->first);
			double newcondEntropy=potMgr->getConditionalEntropy(sFactor->fId,fVars,varSet);

			mbVars[sFactor->fId]=0;
			double newmbCondEntropy=potMgr->getConditionalEntropy(mbFactor->fId,mbVars,varSet);

			double gain=mbcondEntropy-newmbCondEntropy;
			double loss=newcondEntropy-condEntropy;
			double gainpenalty=mbPenalty*log(((double) mbFactor->mergedMB.size())/((double)mbFactor->mergedMB.size()+1));
			double losspenalty=mbPenalty*log(((double)sFactor->mergedMB.size())/((double)sFactor->mergedMB.size() + 1));
			gain=gain+gainpenalty;
			loss=loss+losspenalty;
			if(((gain > 0) || (loss > 0)) && (mbFactor->mergedMB.size() <maxFactorSize_Approx ))
			{
				mbFactor->mergedMB[sFactor->fId]=0;
				mbFactor->mbScore=newmbCondEntropy;
				fVars[mbFactor->fId]=0;
				totalAdd++;
			}
			else
			{
				deleteMBVarID[mbIter->first]=0;
				sFactor->mbScore=newcondEntropy;
				totalDel++;
			}
		}
		for(INTINTMAP_ITER dIter=deleteMBVarID.begin();dIter!=deleteMBVarID.end();dIter++)
		{
			INTINTMAP_ITER mIter=sFactor->mergedMB.find(dIter->first);
			sFactor->mergedMB.erase(mIter);
		}
		cout <<"Total additions " << totalAdd << " total deletions " << totalDel << endl;
	}
	return 0;
}


//This function belongs to a bottom up approach of making the Markov blankets consistent.
//We want to make sure that once we are done with the consistency check at one level l, we are guaranteed that every
//factor below l is consistent with the factors at l. 
//For example, consider the factor X and its possible 1-extensions, which can
//be given to us by the superset matrix. Let Y=X\Union Xi be one such extension. To make
//the Markov blankets of X consistent with Y, we need to check that the MBx is a subset of
//Mby \Union Y/X. We can either delete a variable from MBx or add a variable to MBy.
//While adding a variable to MBy will not requires downstream of l, deleting a variable
//from MBx will require checks because we need to make sure that does not make the 
//MBx inconsistent with the Markov blankets of the subsets of X. 
//Using a top down approach also has the same problem if we want to add variables
//to the Markov blanket, because then we have to go up the hierarchy to make sure that
//everything is consistent.
//We chose the bottom up approach and allow only the addition of random variables to the Markov
//blankets. We delete only those R.Vs. from X that are local to the Mx. If we encounter a variable
//that must be deleted which is not local to Mx, that is, was added because a subset of X
//had it in its Markov blanket, we simply add the the variable to the superset of X.
int
FactorManager::makeMarkovBlanketConsistent(SlimFactor* aFactor)
{
	INTINTMAP* superset=lattice.getSupersets(aFactor->fId);
	int delOps=0;
	int addOps=0;
	int lowmi_addOps=0;
	if(superset==NULL)
	{
		return 0;
	}
	INTINTMAP_ITER sIter=superset->begin();
	int* tempVids=new int[maxFactorSize_Approx];
	while(sIter!=superset->end())
	{
		SlimFactor* ssFactor=slimFactorSet[sIter->first];
		
		if( ( (ssFactor->vCnt-aFactor->vCnt)>1) || (ssFactor->vCnt==maxFactorSize_Approx))
		{
			break;
		}
		
		for(INTINTMAP_ITER mbvIter=aFactor->mergedMB.begin();mbvIter!=aFactor->mergedMB.end();mbvIter++)
		{
			int vid=mbvIter->first;
			if((ssFactor->mergedMB.find(vid)==ssFactor->mergedMB.end()) && (!ssFactor->isMemberVariable(vid)))	
			{
				//Now we need to decide whether to add the variable to ssFactor Markov blanket
				//or to delete it from sFactor. Let the variables of ssFactor be Y. Let the variable
				//missing from MBy be V.  Before we add to MBy we need to check that I(Y,V) is a valid factor
				int ssid=0;
				while((ssid<ssFactor->vCnt) && (ssFactor->vIds[ssid]<vid ))
				{
					tempVids[ssid]=ssFactor->vIds[ssid];
					ssid++;
				}
				tempVids[ssid]=vid;
				ssid++;
				while(ssid<ssFactor->vCnt)
				{
					tempVids[ssid]=ssFactor->vIds[ssid];
					ssid++;
				}
				int testfId=getFactorIndex(tempVids,ssFactor->vCnt+1);
				if(slimFactorSet.find(testfId)!=slimFactorSet.end())
				{
					//This just indicates that this variable in ssFactor's Markov blanket has been added
					//to ensure consistency of a Markov blanket of a subset of ssFactor
					ssFactor->mergedMB[vid]=ssFactor->vCnt-1;	
					addOps++;
				}
				else
				{
					if(aFactor->mergedMB[vid]==0)
					{
						INTINTMAP_ITER dIter=aFactor->mergedMB.find(vid);
						aFactor->mergedMB.erase(dIter);
						delOps++;
					}
					else
					{
						ssFactor->mergedMB[vid]=ssFactor->vCnt-1;	
						lowmi_addOps++;
					}
				}
			}
		}
		sIter++;
	}
	delete[] tempVids;
	//if(addOps || delOps || lowmi_addOps)
	if(lowmi_addOps)
	{
		cout<< " Factor " << aFactor->fId <<" addOps " << addOps << " delOps " << delOps << " low mutualinfo addOps " << lowmi_addOps << endl;
	}
	return 0;
}

//Here we generate the clusters using the same technique as the lattice cluster generation algorithm
//but we restrict the growth of clusters only to the Markov blankets.
//First of all we do not need to make checks of non-random mutual information for clusters of size <=maxFactorSize
//because we have already filtered for good clusters in filterClustersWithMI
int
FactorManager::produceClusters(double reqConf,int maxClusterSize)
{
	//Start with all factors of size 2
	int currK=2;
	//This plays the role of Fk in the pseudo code
	INTVECT parentIds;
	INTVECT newParentIds;
	INTVECT penultimateMaxIds;
	bool latticeCheckLowerLevel=true;

	//The variables of the new factor
	int* newVids=new int[maxClusterSize];
	//All factors of size 2 which are in slimFactorSet can be blindly added into goodSlimFactors because we have already
	//done the mutual information a.k.a. the support test and confidence does not really make sense here.
	//Since the factor ids grow monotonically with the size of the factors, we can simply iterate over slimfactor set
	//and break out when we reach the first factor whose number of variables is greater than 2.
	map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
	while(sfIter!=slimFactorSet.end())
	{
		if(sfIter->second->vCnt==2)
		{
			goodSlimFactors[sfIter->first]=sfIter->second;
			sfIter->second->confidence=1.0;
			parentIds.push_back(sfIter->first);
		}
		sfIter++;	
		if(sfIter->second->vCnt>2)
		{
			continue;
		}
	}
	if(maxFactorSize_Approx==2)
	{
		map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
		while(sfIter!=slimFactorSet.end())
		{
			if(sfIter->second->vCnt==1)
			{
				//This is needed for growing clusters using Apriori
				penultimateMaxIds.push_back(sfIter->first);
				goodSlimFactors[sfIter->first]=sfIter->second;
				sfIter->second->confidence=1.0;
			}
			sfIter++;
			if(sfIter->second->vCnt>1)
			{
				continue;
			}
		}
	}
	int oldClusterCnt=goodSlimFactors.size();
	int newClusterCnt=1;
	//Use this make map to not reanalyze factors
	INTINTMAP generatedIDs;
	while((newClusterCnt>0) && (currK<=maxClusterSize))//No more factors can be added
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		if(currK<maxFactorSize_Approx)
		{
			//These are factors that have enough mutual information but are not supported
			//at the lower level
			int failedAttempts=0;
			int totalAttempts=0;
			for(int p=0;p<parentIds.size();p++)
			{
				SlimFactor* pFactor=goodSlimFactors[parentIds[p]];
				if(pFactor->fId==160)
				{
					cout <<"Stop here " << endl;
				}
				//Use variables from mergedMB to create a cluster a of size currK+1
				for(INTINTMAP_ITER vIter=pFactor->mergedMB.begin();vIter!=pFactor->mergedMB.end();vIter++)
				{
					int newVId=vIter->first;
					//Do merge-sort of the parent variables and newVId
					int fIter=0;
					int dIter=0;
					while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					newVids[dIter]=newVId;
					dIter++;
					while(fIter<pFactor->vCnt) 
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					int currFid=getFactorIndex(newVids,pFactor->vCnt+1);
					if(slimFactorSet.find(currFid)==slimFactorSet.end())
					{
						cout << "No factor id  " << currFid << endl;
						return -1;
					}
					SlimFactor* sFactor=slimFactorSet[currFid];
					totalAttempts++;
					//If this check is turned off, we will get the same number of clusters
					//at currK as we did before executing this code.
					if(latticeCheckLowerLevel)
					{
						INTINTMAP* subsets=lattice.getSubsets(sFactor->fId);
						//get all subsets and check for confidence
						//We need to check the last sFactor->vCnt subsets as these
						//will be the ones that correspond to the maximal subsets
						map<int,int>::reverse_iterator rIter=subsets->rbegin();
						int sscnt=0;
						double hitConf=0;
						SlimFactor* aSubset=slimFactorSet[rIter->first];
						while((sscnt<sFactor->vCnt) && ((sFactor->vCnt-aSubset->vCnt) ==1 ))
						{
							int ssId=rIter->first;
							map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(ssId);
							if(ssIter!=goodSlimFactors.end())
							{
								hitConf=hitConf+ssIter->second->confidence;
							}
							sscnt++;
							rIter++;
							aSubset=slimFactorSet[rIter->first];
						}
						double conf=hitConf/((double) sFactor->vCnt);
						if(conf>=reqConf)
						{
							goodSlimFactors[sFactor->fId]=sFactor;
							sFactor->confidence=conf;
							newParentIds.push_back(sFactor->fId);
						}
						else
						{
							failedAttempts++;
						}
					}
					else //No need to check confidence at the lower levels where
						//we can compute multi-information correctly
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						newParentIds.push_back(sFactor->fId);
					}
				}//Done with one factor
			}//Done with all factors
			cout <<"Failed to add " << failedAttempts << " out of total " << totalAttempts <<" factors at level " << currK+1 << endl;
			if(currK==(maxFactorSize_Approx-1))
			{
				for(int i=0;i<parentIds.size();i++)
				{
					penultimateMaxIds.push_back(parentIds[i]);
				}
			}
			//Now update parentIds from newParentIds
			parentIds.clear();
			for(int i=0;i<newParentIds.size();i++)
			{
				parentIds.push_back(newParentIds[i]);
			}
			newParentIds.clear();
		}//Upto here we can exactly compute the Markov blankets.
		else
		{
			//The number of variables to add to the factor to create a new factor of size currK+1
			int varsToAdd=currK+1-(maxFactorSize_Approx-1);
			//penultimateMaxIds is made up of all factors of size maxFactorSize-1 that have 
			//Markov blankets associated with them.
			for(int p=0;p<penultimateMaxIds.size();p++)
			{
				SlimFactor* pFactor=goodSlimFactors[penultimateMaxIds[p]];
				//Need to make sure that pFactor's Markov blanket has currK+1-pFactor->vCnt R.V.s
				if(pFactor->mergedMB.size() < varsToAdd)
				{
					continue;
				}
				//Now we need to create subsets of the Markov blanket of size 2.
				
				pFactor->genMBSubsets(varsToAdd);
				int mbssStart=pFactor->mbSubsetStartInd[varsToAdd];
				int mbssEnd=pFactor->mbSubsetStartInd[varsToAdd+1];
				for(int mbsid=mbssStart;mbsid<mbssEnd;mbsid++)
				{
					INTINTMAP* mbset=pFactor->mbSubsets[mbsid];
					//Here we have to create a new factor
					SlimFactor* sFactor=new SlimFactor;
					sFactor->vCnt=currK+1;
					sFactor->vIds=new int[currK+1];
					//Now we need to do a merge sort into the variables of sFactor
					INTINTMAP_ITER mIter=mbset->begin();
					int fIter=0;
					int dIter=0;
					while(dIter!=sFactor->vCnt && mIter!=mbset->end() && fIter!=pFactor->vCnt)
					{
						if(mIter->first<pFactor->vIds[fIter])
						{
							sFactor->vIds[dIter]=mIter->first;
							mIter++;
						}
						else
						{
							sFactor->vIds[dIter]=pFactor->vIds[fIter];
							fIter++;
						}
						dIter++;
					}
					while(mIter!=mbset->end())
					{
						sFactor->vIds[dIter]=mIter->first;
						mIter++;
						dIter++;
					}
					while(fIter!=pFactor->vCnt)
					{
						sFactor->vIds[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					sFactor->fId=getFactorIndex(sFactor->vIds,sFactor->vCnt);
					if(generatedIDs.find(sFactor->fId) !=generatedIDs.end())
					{
						delete sFactor;
						continue;
					}
					generatedIDs[sFactor->fId]=0;
					//Allocate memory for subsets
					int** subsets=new int*[sFactor->vCnt];
					for(int i=0;i<sFactor->vCnt;i++)
					{
						subsets[i]=new int[sFactor->vCnt-1];
					}
					sFactor->generateMaximalSubsets(subsets);
					//Now check the confidence and in the meantime store the subset ids
					//to update the lattice structure
					int* ssIds=new int[sFactor->vCnt];
					double hitConf=0;
					for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
					{
						int sId=getFactorIndex(subsets[sscnt],sFactor->vCnt-1);
						ssIds[sscnt]=sId;
						map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(sId);
						if(ssIter!=goodSlimFactors.end())
						{
							hitConf=hitConf+ssIter->second->confidence;
						}
					}
					double conf=hitConf/((double)sFactor->vCnt);
					if(conf>=reqConf)
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						sFactor->confidence=conf;
						newParentIds.push_back(sFactor->fId);
						//Update the lattice structure
						for(int i=0;i<sFactor->vCnt;i++)
						{
							if(goodSlimFactors.find(ssIds[i])!=goodSlimFactors.end())
							{
								lattice.addSubset(ssIds[i],sFactor->fId);
								lattice.addSuperset(sFactor->fId,ssIds[i]);
							}
						}
					}
					else
					{
						delete sFactor;
					}
					for(int i=0;i<sFactor->vCnt;i++)
					{
						delete[] subsets[i];
					}
					delete[] subsets;
					delete ssIds;
				}
			}
		}
		currK++;
		newClusterCnt=goodSlimFactors.size()-oldClusterCnt;
		oldClusterCnt=goodSlimFactors.size();
		cout <<"Added " << newClusterCnt<< " new clusters of size " << currK<< endl; 
	}
	return 0;
}


int
FactorManager::produceClusters_NoDup(double reqConf,int maxClusterSize)
{
	//Start with all factors of size 2
	int currK=2;
	//This plays the role of Fk in the pseudo code
	INTINTMAP parentIds;
	INTINTMAP newParentIds;
	INTINTMAP penultimateMaxIds;
	bool latticeCheckLowerLevel=true;

	//The variables of the new factor
	int* newVids=new int[maxClusterSize];
	//All factors of size 2 which are in slimFactorSet can be blindly added into goodSlimFactors because we have already
	//done the mutual information a.k.a. the support test and confidence does not really make sense here.
	//Since the factor ids grow monotonically with the size of the factors, we can simply iterate over slimfactor set
	//and break out when we reach the first factor whose number of variables is greater than 2.
	map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
	while(sfIter!=slimFactorSet.end())
	{
		if(sfIter->second->vCnt==2)
		{
			//Get the group to which this factor belongs
			string& groupSign=factSignMap[sfIter->first];
			int repId=groupFactorRepMap[groupSign];
			goodSlimFactors[sfIter->first]=sfIter->second;
			sfIter->second->confidence=1.0;
			if(parentIds.find(repId)==parentIds.end())
			{
				parentIds[repId]=0;
			}
		}
		if(sfIter->second->vCnt>2)
		{
			break;
		}
		sfIter++;	
	}
	if(maxFactorSize_Approx==2)
	{
		map<int,SlimFactor*>::iterator sfIter=slimFactorSet.begin();
		while(sfIter!=slimFactorSet.end())
		{
			if(sfIter->second->vCnt==1)
			{
				//This is needed for growing clusters using Apriori
				goodSlimFactors[sfIter->first]=sfIter->second;
				sfIter->second->confidence=1.0;
				string& groupSign=factSignMap[sfIter->first];
				int repId=groupFactorRepMap[groupSign];
				if(penultimateMaxIds.find(repId)==penultimateMaxIds.end())
				{
					penultimateMaxIds[repId]=0;
				}
			}
			sfIter++;
			if(sfIter->second->vCnt>1)
			{
				break;
			}
		}
	}
	int oldClusterCnt=goodSlimFactors.size();
	int newClusterCnt=1;
	//Use this make map to not reanalyze factors
	INTINTMAP generatedIDs;
	while((newClusterCnt>0) && (currK<=maxClusterSize))//No more factors can be added
	{
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		if(currK<maxFactorSize_Approx)
		{
			//These are factors that have enough mutual information but are not supported
			//at the lower level
			int failedAttempts=0;
			int totalAttempts=0;
			for(INTINTMAP_ITER pIter=parentIds.begin();pIter!=parentIds.end();pIter++)
			{
				SlimFactor* pFactor=slimFactorSet[pIter->first];
				//Use variables from mergedMB to create a cluster a of size currK+1
				for(INTINTMAP_ITER vIter=pFactor->mergedMB.begin();vIter!=pFactor->mergedMB.end();vIter++)
				{
					int newVId=vIter->first;
					//Do merge-sort of the parent variables and newVId
					int fIter=0;
					int dIter=0;
					while((fIter<pFactor->vCnt) && (pFactor->vIds[fIter]<newVId))
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					newVids[dIter]=newVId;
					dIter++;
					while(fIter<pFactor->vCnt) 
					{
						newVids[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					int currFid=getFactorIndex(newVids,pFactor->vCnt+1);
					if(slimFactorSet.find(currFid)==slimFactorSet.end())
					{
						cout << "Factor " << currFid << " has probably been deleted " <<  endl;
						continue;
					//	return -1;
					}
					//Must select a factor for growth only from the nu dup list
					SlimFactor* sFactor=slimFactorSet[currFid];
					totalAttempts++;
					string& groupSign=factSignMap[sFactor->fId];
					int repId=groupFactorRepMap[groupSign];
					
					//If this check is turned off, we will get the same number of clusters
					//at currK as we did before executing this code.
					INTINTMAP* subsets=lattice.getSubsets(sFactor->fId);
					//get all subsets and check for confidence
					//We need to check the last sFactor->vCnt subsets as these
					//will be the ones that correspond to the maximal subsets
					map<int,int>::reverse_iterator rIter=subsets->rbegin();
					int sscnt=0;
					double hitConf=0;
					//However check for subsets in slimFactorSet
					SlimFactor* aSubset=slimFactorSet[rIter->first];
					while((sscnt<sFactor->vCnt) && ((sFactor->vCnt-aSubset->vCnt) ==1 ))
					{
						int ssId=rIter->first;
						map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(ssId);
						if(ssIter!=goodSlimFactors.end())
						{
							hitConf=hitConf+ssIter->second->confidence;
						}
						sscnt++;
						rIter++;
						aSubset=slimFactorSet[rIter->first];
					}
					double conf=hitConf/((double) sFactor->vCnt);
					sFactor->confidence=conf;
					if(latticeCheckLowerLevel)
					{
						if(conf>=reqConf)
						{
							goodSlimFactors[sFactor->fId]=sFactor;
							newParentIds[repId]=0;
						}
						else
						{
							failedAttempts++;
						}
					}
					else //No need to check confidence at the lower levels where
						//we can compute multi-information correctly
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						newParentIds[repId]=0;
					}
				}//Done with one factor
			}//Done with all factors
			cout <<"Failed to add " << failedAttempts << " out of total " << totalAttempts <<" factors at level " << currK+1 << endl;
			if(currK==(maxFactorSize_Approx-1))
			{
				for(INTINTMAP_ITER pIter=parentIds.begin();pIter!=parentIds.end();pIter++)
				{
					penultimateMaxIds[pIter->first]=0;
				}
			}
			//Now update parentIds from newParentIds
			parentIds.clear();
			for(INTINTMAP_ITER pIter=newParentIds.begin();pIter!=newParentIds.end();pIter++)
			{
				parentIds[pIter->first]=0;
			}
			newParentIds.clear();
		}//Upto here we can exactly compute the Markov blankets.
		else
		{
			//The number of variables to add to the factor to create a new factor of size currK+1
			int varsToAdd=currK+1-(maxFactorSize_Approx-1);
			int totalAttempts=0;
			int failedAttempts=0;
			//penultimateMaxIds is made up of all factors of size maxFactorSize-1 that have 
			//Markov blankets associated with them.
			for(INTINTMAP_ITER pIter=penultimateMaxIds.begin();pIter!=penultimateMaxIds.end();pIter++)
			{
				SlimFactor* pFactor=slimFactorSet[pIter->first];
				//Need to make sure that pFactor's Markov blanket has currK+1-pFactor->vCnt R.V.s
				if(pFactor->mergedMB.size() < varsToAdd)
				{
					continue;
				}
				//Now we need to create subsets of the Markov blanket of size 2.
				
				pFactor->genMBSubsets(varsToAdd);
				int mbssStart=pFactor->mbSubsetStartInd[varsToAdd];
				int mbssEnd=pFactor->mbSubsetStartInd[varsToAdd+1];
				for(int mbsid=mbssStart;mbsid<mbssEnd;mbsid++)
				{
					INTINTMAP* mbset=pFactor->mbSubsets[mbsid];
					//Here we have to create a new factor
					SlimFactor* sFactor=new SlimFactor;
					sFactor->vCnt=currK+1;
					sFactor->vIds=new int[currK+1];
					//Now we need to do a merge sort into the variables of sFactor
					INTINTMAP_ITER mIter=mbset->begin();
					int fIter=0;
					int dIter=0;
					while(dIter!=sFactor->vCnt && mIter!=mbset->end() && fIter!=pFactor->vCnt)
					{
						if(mIter->first<pFactor->vIds[fIter])
						{
							sFactor->vIds[dIter]=mIter->first;
							mIter++;
						}
						else
						{
							sFactor->vIds[dIter]=pFactor->vIds[fIter];
							fIter++;
						}
						dIter++;
					}
					while(mIter!=mbset->end())
					{
						sFactor->vIds[dIter]=mIter->first;
						mIter++;
						dIter++;
					}
					while(fIter!=pFactor->vCnt)
					{
						sFactor->vIds[dIter]=pFactor->vIds[fIter];
						fIter++;
						dIter++;
					}
					sFactor->fId=getFactorIndex(sFactor->vIds,sFactor->vCnt);
					if(generatedIDs.find(sFactor->fId) !=generatedIDs.end())
					{
						delete sFactor;
						continue;
					}
					totalAttempts++;
					generatedIDs[sFactor->fId]=0;
					//Allocate memory for subsets
					int** subsets=new int*[sFactor->vCnt];
					for(int i=0;i<sFactor->vCnt;i++)
					{
						subsets[i]=new int[sFactor->vCnt-1];
					}
					sFactor->generateMaximalSubsets(subsets);
					//Now check the confidence and in the meantime store the subset ids
					//to update the lattice structure
					int* ssIds=new int[sFactor->vCnt];
					double hitConf=0;
					for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
					{
						int sId=getFactorIndex(subsets[sscnt],sFactor->vCnt-1);
						ssIds[sscnt]=sId;
						map<int,SlimFactor*>::iterator ssIter=goodSlimFactors.find(sId);
						if(ssIter!=goodSlimFactors.end())
						{
							hitConf=hitConf+ssIter->second->confidence;
						}
					}
					double conf=hitConf/((double)sFactor->vCnt);
					if(conf>=reqConf)
					{
						goodSlimFactors[sFactor->fId]=sFactor;
						sFactor->confidence=conf;
						//Update the lattice structure
						for(int i=0;i<sFactor->vCnt;i++)
						{
							if(goodSlimFactors.find(ssIds[i])!=goodSlimFactors.end())
							{
								lattice.addSubset(ssIds[i],sFactor->fId);
								lattice.addSuperset(sFactor->fId,ssIds[i]);
							}
						}
					}
					else
					{
						failedAttempts++;
						delete sFactor;
					}
					for(int i=0;i<sFactor->vCnt;i++)
					{
						delete[] subsets[i];
					}
					delete[] subsets;
					delete ssIds;
				}
			}
			cout <<"Failed to add " << failedAttempts << " out of total " << totalAttempts <<" factors at level " << currK+1 << endl;
		}
		currK++;
		newClusterCnt=goodSlimFactors.size()-oldClusterCnt;
		oldClusterCnt=goodSlimFactors.size();
		cout <<"Added " << newClusterCnt<< " new clusters of size " << currK<< endl; 
	}
	return 0;
}



int
FactorManager::getGoodCandidateMarkovBlankets(SlimFactor* sFactor,INTINTMAP* supersets,int currK, double& minEntropy,double marginalEntropy,INTDBLMAP& goodMBs)
{
	if(supersets==NULL)
	{
		return 0;
	}
	if(supersets->size()==0)
	{
		return 0;
	}
	INTINTMAP_ITER fIter=supersets->begin();
	bool doneFlag=false;
	INTDBLMAP entropyMap;
	vector<int> mbIds;
	double currMinEntropy=minEntropy;
	//This is a temporarary array to store the variables of the Markov blanket whenever need arises
	int* tempMBVars=new int[maxFactorSize_Approx];
	//The supersets are ordered such that the smaller ones are before the larger ones
	//This is because we assign ids to the clusters based on the order in which they are created
	while((!doneFlag) &&(fIter!=supersets->end()))
	{
		int supId=fIter->first;
		SlimFactor* supFactor=slimFactorSet[supId];
		if(supFactor->vCnt>currK)
		{
			doneFlag=true;
			break;
		}//Do stage I
		double proposedEntropy=marginalEntropy-supFactor->mutualInfo;
		double penalty=0;
		if(supFactor->vCnt-sFactor->vCnt>1)
		{
			//Here we need to get the multi-information of the variables only in the Markov blanket
			int mbvarCnt=0;
			sFactor->getSetDiff(supFactor,tempMBVars,mbvarCnt);
			
			double mbInfo=getMIFromVars(tempMBVars,mbvarCnt);
			if(mbInfo<0)
			{
				cout << "Error occured while getting mutual information of sub factors" << endl;
				return -1;
			}
			proposedEntropy=proposedEntropy+mbInfo;
			penalty=mbPenalty*log(double(currK-sFactor->vCnt));
		}
		if(proposedEntropy*(1+penalty)<=currMinEntropy)
		{
			if(proposedEntropy<0)
			{
			//	cout<<"Negative entropy for factor " << sFactor->fId << " mbfactor: "<< supFactor->fId << " setting to 0 "<< endl; 
			//	proposedEntropy=0;
				//return -1;
			}
			currMinEntropy=proposedEntropy;
		}
		entropyMap[supFactor->fId]=proposedEntropy;
		//if(mbIds.size()<beamSize)
		if(mbIds.size()<10)
		{
			mbIds.push_back(supFactor->fId);
			//if(mbIds.size()==beamSize)
			if(mbIds.size()==10)
			{
				for(int i=0;i<mbIds.size();i++)
				{
					for(int j=i+1;j<mbIds.size();j++)
					{
						if(entropyMap[mbIds[i]]>entropyMap[mbIds[j]])
						{
							int tempfid=mbIds[i];
							mbIds[i]=mbIds[j];
							mbIds[j]=tempfid;
						}
					}
				}
			}
		}
		else 
		{
			//Compare with the last element
			int i=mbIds.size()-1;
			bool insertFlag=false;
			while(proposedEntropy<entropyMap[mbIds[i]] && i>=0)
			{
				i--;
				insertFlag=true;
			}

			if(insertFlag)
			{
				i++;
				if(mbIds.size()>1)
				{
					//drop off the last element
					mbIds.pop_back();
					//make a copy of the last element
					mbIds.push_back(mbIds[mbIds.size()-1]);
					for(int j=i;j<mbIds.size()-2;j++)
					{
						mbIds[j]=mbIds[j+1];
					}
				}
				mbIds[i]=supFactor->fId;
			}

		}
		fIter++;
	}
	if(currMinEntropy>=minEntropy)
	{
		return 0;
	}
	double newMinEntropy=currMinEntropy;
	/*for(INTDBLMAP_ITER idIter=entropyMap.begin();idIter!=entropyMap.end();idIter++)
	{
		if(idIter->second<=currMinEntropy)
		{
			goodMBs[idIter->first]=idIter->second;
		}
	}*/
	int i=0; 
	while(i<beamSize) 
	{
		if(i>=mbIds.size())
		{
			break;
		}
		double anentropy=entropyMap[mbIds[i]];
		if((newMinEntropy <= anentropy) && (anentropy<minEntropy))
		{
			newMinEntropy=anentropy;
		}
		if(newMinEntropy>=anentropy)
		{
			goodMBs[mbIds[i]]=anentropy;
		}
		i++;
	}
	for(int m=0;m<mbIds.size();m++)
	{
		SlimFactor* supFactor=slimFactorSet[mbIds[m]];
		for(int v=0;v<supFactor->vCnt;v++)
		{
			int nId=supFactor->vIds[v];
			if(!sFactor->isMemberVariable(nId))
			{
				sFactor->candidateNeighbours[nId]=0;
			}
		}
	}
	minEntropy=newMinEntropy;
	return 0;
}

int 
FactorManager::getFactorIndex(int* vIds, int vCnt)
{
	string aKey;
	getFactorKey(vIds,vCnt,aKey);
	int fId=-1;
	if(factorNameToIDMap.find(aKey)!=factorNameToIDMap.end())
	{
		fId=factorNameToIDMap[aKey];
	}
	return fId;
}

int
FactorManager::getFactorKey(int* vIds, int vCnt, string& key)
{
	for(int v=0;v<vCnt;v++)
	{
		char aBuff[56];
		sprintf(aBuff,"-%d",vIds[v]);
		key.append(aBuff);
	}
	return 0;
}


int 
FactorManager::getMBFactorIndex(SlimFactor* sFactor)
{
	INTINTMAP fVars;
	fVars[sFactor->fId]=0;
	for(INTINTMAP_ITER vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++)
	{
		fVars[vIter->first]=0;
	}
	string key;
	for(INTINTMAP_ITER vIter=fVars.begin();vIter!=fVars.end();vIter++)
	{
		char aBuff[56];
		sprintf(aBuff,"-%d",vIter->first);
		key.append(aBuff);
	}
	int fId=-1;
	if(factorNameToIDMap.find(key)!=factorNameToIDMap.end())
	{
		fId=factorNameToIDMap[key];
	}
	else
	{
		fId=addNewMBFactor(sFactor);
	}

	return fId;
}

int
FactorManager::addNewMBFactor(SlimFactor* sFactor)
{
	int newVId=sFactor->fId;
	SlimFactor* mbFactor=new SlimFactor;
	mbFactor->vIds=new int[sFactor->mergedMB.size()+1];
	mbFactor->vCnt=sFactor->mergedMB.size()+1;
	int dIter=0;
	INTINTMAP_ITER fIter=sFactor->mergedMB.begin();
	while(fIter!=sFactor->mergedMB.end())
	{
		if(fIter->first>newVId)
		{
			break;
		}
		mbFactor->vIds[dIter]=fIter->first;
		dIter++;
		fIter++;
	}
	mbFactor->vIds[dIter]=newVId;
	dIter++;
	while(fIter!=sFactor->mergedMB.end())
	{
		mbFactor->vIds[dIter]=fIter->first;
		fIter++;
		dIter++;
	}
	//This factor might have been created
	string fKey;
	getFactorKey(mbFactor->vIds,mbFactor->vCnt,fKey);
	int** subsetSpace=new int*[maxFactorSize_Approx];
	for(int i=0;i<maxFactorSize_Approx;i++)
	{
		subsetSpace[i]=new int[maxFactorSize_Approx-1];
	}
	int* ssIds=new int[maxFactorSize_Approx];
	mbFactor->generateMaximalSubsets(subsetSpace);
	int foundSsets=0;
	for(int sscnt=0;sscnt<sFactor->vCnt;sscnt++)
	{
		string ssKey;
		getFactorKey(subsetSpace[sscnt],sFactor->vCnt-1,ssKey);
		if(factorNameToIDMap.find(ssKey)==factorNameToIDMap.end())
		{
			continue;
		}
		int sId=factorNameToIDMap[ssKey];
		ssIds[foundSsets]=sId;
		foundSsets++;
	}
	potMgr->populateFactor(slimFactorSet,vMgr->getVariableSet(),mbFactor,false);
	double randmi_mean=randMI_mean[mbFactor->vCnt];
	double randmi_sd=randMI_std[mbFactor->vCnt];
	if(mbFactor->mutualInfo < randmi_mean+(randmi_sd*misdCnt))
	{
		cout <<"MB factor " << mbFactor->fId <<" has below random mi " << mbFactor->mutualInfo << endl;
	}
	mbFactor->fId=globalFactorID;
	factorNameToIDMap[fKey]=mbFactor->fId;
	factorIDToNameMap[mbFactor->fId]=fKey;
	globalFactorID++;
	if(globalFactorID<0)
	{
		cout <<"Global factor ID became negative" << endl;
	}
	for(int i=0;i<foundSsets;i++)
	{
		int subsetId=ssIds[i];
		if(slimFactorSet.find(subsetId)!=slimFactorSet.end())
		{
			lattice.addSubset(subsetId,sFactor->fId);
			lattice.addSuperset(sFactor->fId,subsetId);
		}
	}
	slimFactorSet[mbFactor->fId]=mbFactor;

	return mbFactor->fId;
}


SlimFactor*
FactorManager::getFactorFromVars(int* vIds, int vCnt)
{
	int fid=getFactorIndex(vIds,vCnt);
	if(slimFactorSet.find(fid)==slimFactorSet.end())
	{
		cout <<"Warning! Accessing null factor " << endl;
		return NULL;
	}
	return slimFactorSet[fid];
}

double
FactorManager::getMIFromVars(int* vIds,int vCnt)
{
	string factorKey;
	getFactorKey(vIds,vCnt,factorKey);
	double mi=-1;
	map<string,int>::iterator aIter=factorNameToIDMap.find(factorKey);
	if(aIter!=factorNameToIDMap.end())
	{
		mi=slimFactorSet[aIter->second]->mutualInfo;
	}
	else if(delFactors_MI.find(factorKey)!=delFactors_MI.end())
	{
		mi=delFactors_MI[factorKey];
	}
	else
	{
		if(mbSpecific_MI.find(factorKey)!=mbSpecific_MI.end())
		{
			mi=mbSpecific_MI[factorKey];
		}
		else
		{
			VSET& varSet=vMgr->getVariableSet();
			Potential* apot=new Potential;
			for(int j=0;j<vCnt;j++)
			{
				Variable* aVar=varSet[vIds[j]];
				if(j==vCnt-1)
				{
					apot->setAssocVariable(aVar,Potential::FACTOR);
				}
				else
				{
					apot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
				}
			}
			apot->potZeroInit();
			potMgr->populatePotential(apot,false);
			apot->calculateConditionalEntropy();
			mi=-1*apot->getJointEntropy();
			for(int j=0;j<vCnt;j++)
			{
				SlimFactor* sFactor=slimFactorSet[vIds[j]];
				mi=mi+sFactor->marginalEntropy;
			}
			mbSpecific_MI[factorKey]=mi;
		}
	}
	return mi;
}

//vId is the variable id which is used to create the clusters of size cSize, 
//from a total of N variables, starting with vId
int
FactorManager::getClusterCnt(int vId, int cSize, int N)
{
	int clusterCnt=0;
	if(cSize==1)
	{
		clusterCnt=1;
	}
	else if(cSize==2)
	{
		clusterCnt=N-vId-1;
	}
	else if(cSize==3)
	{
		int n=N-cSize+1-vId;
		clusterCnt=(n*(n+1))/2;
	}
	else if(cSize>3)
	{
		//This is the sum of all clusters of size cSize-1 starting
		//from vId+1 to N
		for(int newvId=vId+1;newvId<N;newvId++)
		{
			clusterCnt=clusterCnt+getClusterCnt(newvId,cSize-1,N);
		}
	}
	return clusterCnt;
}

int
FactorManager::initFactorSet()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	
	int factorCnt=0;
	int fSize=1;
	int vCnt=vMgr->getVariableSet().size();
	while(fSize<=maxFactorSize)
	{
		int fCnt=combCnt(vCnt,fSize);
		factorCnt=factorCnt+fCnt;
		fSize++;
	}
	//Calculate the number of factors
	cout <<" Number of factors "  << factorCnt << endl;
	fSize=1;
	int fCnt=combCnt(vCnt,fSize);
	for(int i=0;i<factorCnt;i++)
	{
		SlimFactor* sFactor=new SlimFactor;
		slimFactorSet[i]=sFactor;

		if(i==fCnt)
		{
			fSize++;
			fCnt=fCnt+combCnt(vCnt,fSize);
		}
		sFactor->vIds=new int[fSize];
		sFactor->vCnt=fSize;
		sFactor->secondPId=-1;
		sFactor->mutualInfo=0;
		sFactor->conditionalEntropy=0;
		sFactor->fId=globalFactorID;
		globalFactorID++;
	}
	cout << "Global factor id " << globalFactorID << endl;
	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;

	return 0;
}

//Here we simply write to the memory that we have allocated
int
FactorManager::populateFactorSet()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	int fSize=1;
	//These indices are simply used to create factors of size k+1
	//from factors of size k
	int startFid=0;
	int endFid=0;
	int currFid=0;
	map<int,Variable*>& variableSet=vMgr->getVariableSet();
	for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		SlimFactor* sFactor=slimFactorSet[currFid];
		currFid++;
		sFactor->vIds[sFactor->vCnt-1]=vIter->first;
		string fKey;
		getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
		factorNameToIDMap[fKey]=sFactor->fId;
		factorIDToNameMap[sFactor->fId]=fKey;
	}
	endFid=currFid;
	fSize++;
	
	//This int matrix is allocated before hand and reused
	//This is done to avoid the allocation and deallocation repeatedly
	//The maximum of subsets is equal to the maxFactorSize
	//Each subset can be at most maxFactorSize-1
	int** subsetSpace=new int*[maxFactorSize];
	for(int i=0;i<maxFactorSize;i++)
	{
		subsetSpace[i]=new int[maxFactorSize-1];
	}
	while(fSize<=maxFactorSize)
	{
		int pFidCnt=endFid-startFid;
		//The number of newFids added
		//For each parent factor, iterate over the set of variables
		//and add the variable in the parent to the new factor 
		for(int i=0;i<pFidCnt;i++)
		{
			SlimFactor* pFactor=slimFactorSet[startFid+i];
			for(map<int,Variable*>::iterator vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
			{
				int newVId=vIter->first;
				if(newVId<=pFactor->vIds[pFactor->vCnt-1])
				{
					continue;
				}
				SlimFactor* sFactor=slimFactorSet[currFid];
				//Here we need to add the super-set and sub-set relationships
				//Specifically, sFactor is a super-set of pFactor
				//pFactor is a sub-set of sFactor
				for(int j=0;j<pFactor->vCnt;j++)
				{
					sFactor->vIds[j]=pFactor->vIds[j];
				}
				sFactor->vIds[sFactor->vCnt-1]=newVId;
				string fKey;
				getFactorKey(sFactor->vIds,sFactor->vCnt,fKey);
				factorNameToIDMap[fKey]=sFactor->fId;
				factorIDToNameMap[sFactor->fId]=fKey;
				sFactor->secondPId=startFid+i;
				currFid++;
				addToLattice(sFactor,subsetSpace);
			}
		}
		fSize++;
		startFid=endFid;
		endFid=currFid;
		cout <<"Populated " << endFid-startFid << " new factors "<< endl; 
	}
	for(int i=0;i<maxFactorSize;i++)
	{
		delete [] subsetSpace[i];
	}
	delete[] subsetSpace;
	cout << "Global factor id " << globalFactorID << endl;
	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;
	return 0;
}

int
FactorManager::checkFactorIds()
{
	/*for(int i=0;i<factorCnt;i++)
	{
		SlimFactor* sFactor=slimFactorSet[i];
		int checkFid=getFactorIndex(sFactor->vIds,sFactor->vCnt);
		if(checkFid!=i)
		{
			cout <<"Index not calculated properly for ";
			showFactor(sFactor,cout);
			cout <<"Actual index "<< i << " Calculated index " << checkFid << endl;
		}
	}*/
	return 0;
}

//Computes n choose k, where k is very small
int
FactorManager::combCnt(int n, int k)
{
	int fCnt=1;
	int start=n;
	int i=0;
	while(i<k)
	{
		fCnt=fCnt*start;
		start--;
		i++;
	}
	start=k;
	i=1;
	while(i<=k)
	{
		fCnt=fCnt/i;
		i++;
	}
	return fCnt;
}


Error::ErrorCode
FactorManager::estimateClusterProperties()
{
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	
	//Error::ErrorCode err=potMgr->populatePotentialsSlimFactors_Summary(slimFactorSet,vMgr->getVariableSet(),maxFactorSize);
	Error::ErrorCode err=potMgr->populatePotentialsSlimFactors(slimFactorSet,vMgr->getVariableSet());

	gettimeofday(&endtime,NULL);
	cout << "Time elapsed " << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;

	if(err!=Error::SUCCESS)
	{
		return err;
	}
	for(map<int,SlimFactor*>::iterator fIter=slimFactorSet.begin();fIter!=slimFactorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		if(sFactor->vCnt>1)
		{
			//The multiinfo is just the difference of the joint entropy and the sum of the marginal
			//entropies of the variables
			double mi=(-1)*sFactor->jointEntropy;
			for(int j=0;j<sFactor->vCnt;j++)
			{
				SlimFactor* subFactor=getFactorFromVars(sFactor->vIds+j,1);
				mi=mi+subFactor->conditionalEntropy;
			}
			sFactor->mutualInfo=mi;
			if(mi<0)
			{
				cout <<"Negative mi for "<< fIter->first << "setting to 0" << endl;
				mi=0;
			}
		}
		else
		{
			sFactor->mbScore=sFactor->conditionalEntropy;
		}
	}
	return Error::SUCCESS;
}



//This function generates all maximal subsets of this factor
//Then it adds the subset and superset relationships
int
FactorManager::addToLattice(SlimFactor* aFactor,int** subsetSpace)
{
	aFactor->generateMaximalSubsets(subsetSpace);
	for(int i=0;i<aFactor->vCnt;i++)
	{
		int subsetId=getFactorIndex(subsetSpace[i],aFactor->vCnt-1);
		lattice.addSubset(subsetId,aFactor->fId);
		lattice.addSuperset(aFactor->fId,subsetId);
	}
	return 0;
}

int
FactorManager::deleteFromLattice(int factorId)
{
	lattice.deleteFromLattice(factorId);
	return 0;
}

//totalVars is the number of variables in the union of clusterA and clusterB
double 
FactorManager::getOverlap(SlimFactor* clusterA,SlimFactor* clusterB,int& totalVars)
{
	SlimFactor* smallCluster=clusterB;
	SlimFactor* bigCluster=clusterA;
	if(clusterA->vCnt<clusterB->vCnt)
	{
		smallCluster=clusterA;
		bigCluster=clusterB;
	}

	int commonElements=0;
	for(int i=0;i<smallCluster->vCnt;i++)
	{
		int vInd=0;
		bool found=false;
		while((vInd<bigCluster->vCnt) && (!found))
		{
			if(smallCluster->vIds[i]==bigCluster->vIds[vInd])
			{
				found=true;
			}
			vInd++;
		}
		if(found)
		{
			commonElements++;
		}
	}
	totalVars=smallCluster->vCnt+bigCluster->vCnt-commonElements;
	double overlap=((double) commonElements)/ ((double) smallCluster->vCnt);
	return overlap;
}

//Sort all factors in factorIndSet. These are all the factors of a particular
//size that satisfy some criteria
//This must be done in decreasing order of mutual information
int
FactorManager::qsort(int* factorIndSet,int startind, int endind)
{
	if(endind-startind <=1)
	{
		return 0;
	}
	//pivot is the first element
	int pivot=startind;
	for(int i=startind+1;i<endind;i++)
	{
		//If the element at i is less than element at pivot
		//place element before pivot
		if(slimFactorSet[factorIndSet[i]]->mutualInfo < slimFactorSet[factorIndSet[pivot]]->mutualInfo)
		{
			int temp=factorIndSet[i];
			factorIndSet[i]=factorIndSet[pivot];
			factorIndSet[pivot]=temp;
			pivot++;
		}
	}
	//Now make sure everything after the pivot is greater or equal to the pivot
	qsort(factorIndSet,startind,pivot-1);
	qsort(factorIndSet,pivot+1,endind);
	return 0;
}
