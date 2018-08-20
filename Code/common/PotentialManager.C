#include <iostream>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
}

PotentialManager::~PotentialManager()
{
}

int 
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int
PotentialManager::init()
{
	estimateAllMeanCov(false,globalMean,globalCovar);
	return 0;
}

int
PotentialManager::initRandom()
{
	estimateAllMeanCov(true,globalMean_Rand,globalCovar_Rand);
	return 0;
}

int
PotentialManager::estimateAllMeanCov(bool random, INTDBLMAP& gMean, map<int,INTDBLMAP*>& gCovar)
{
	int evidCnt=evMgr->getNumberOfEvidences();
	//First get the mean and then the variance
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(i);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(i);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			if(gMean.find(vId)==gMean.end())
			{
				gMean[vId]=val;
			}
			else
			{
				gMean[vId]=gMean[vId]+val;
			}
		}

	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=gMean.begin();idIter!=gMean.end();idIter++)
	{
		idIter->second=idIter->second/(double) evidCnt;
	}
	//Now the variance
	for(int i=0;i<evidCnt;i++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(i);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(i);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double vval=evid->getEvidVal();
			double vmean=gMean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(EMAP_ITER uIter=vIter;uIter!=evidMap->end();uIter++)
			{
				int uId=uIter->first;
				Evidence* evid1=uIter->second;
				double uval=evid1->getEvidVal();
				double umean=gMean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					(*vcov)[uId]=diffprod;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						(*ucov)[vId]=diffprod;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=2*vIter->second/((double)(gCovar.size()-1));
				vIter->second=2*vIter->second/((double)(evidCnt-1));
			}
			else
			{
				vIter->second=vIter->second/((double)(evidCnt-1));
				//vIter->second=vIter->second/((double)(gCovar.size()-1));
				//vIter->second=0;
			}
		}
	}

	return 0;
}

Error::ErrorCode
PotentialManager::populatePotentialsSlimFactors(map<int,SlimFactor*>& factorSet,VSET& varSet)
{
	//The set of flags to keep status of the potentials that have been calculated
	map<int,bool> doneFlag;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		doneFlag[fIter->first]=false;
	}
	int popFId=0;
	//for(map<int,SlimFactor*>::reverse_iterator rIter=factorSet.rbegin();rIter!=factorSet.rend();rIter++)
	for(map<int,SlimFactor*>::iterator rIter=factorSet.begin();rIter!=factorSet.end();rIter++)
	{
		//If we have computed the potential for this flag move one
		if(doneFlag[rIter->first])
		{
			popFId++;
			continue;
		}
		SlimFactor* sFactor=rIter->second;
		//Otherwise create the potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,false);
		aPotFunc->calculateConditionalEntropy();
		sFactor->conditionalEntropy=aPotFunc->getConditionalEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		if(sFactor->conditionalEntropy<0)
		{
		//	sFactor->jointEntropy=0;
		//	cout <<"Negative entropy for " << sFactor->fId << endl;
		}
		doneFlag[rIter->first]=true;
		delete aPotFunc;
		if(popFId%100000==0)
		{
			cout <<"Done with " << factorSet.size()-popFId << " factors " << endl;
		}
		popFId++;
	}
	return Error::SUCCESS;
}


int 
PotentialManager::estimateMarginalEntropies(map<int,SlimFactor*>& slimFactors,VSET& varSet,bool random)
{	
	for(map<int,SlimFactor*>::iterator aIter=slimFactors.begin();aIter!=slimFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		Potential* aPotFunc=new Potential;
		Variable* aVar=varSet[sFactor->vIds[0]];
		aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,random);
		aPotFunc->calculateConditionalEntropy();
		sFactor->conditionalEntropy=aPotFunc->getConditionalEntropy();
		sFactor->jointEntropy=aPotFunc->getConditionalEntropy();
		delete aPotFunc;
	}
	return 0;
}

//Estimate the random information for all factors of a particular size. Assume that
//randInfo is already allocated
Error::ErrorCode
PotentialManager::estimateRandomInfo(map<int,SlimFactor*>& factorSet, 
		VSET& varSet, vector<double>& randInfo, int fSize)
{
	int rInd=0;
	for(map<int,SlimFactor*>::iterator rIter=factorSet.begin();rIter!=factorSet.end();rIter++)
	{
		SlimFactor* sFactor=rIter->second;
		if(sFactor->vCnt<fSize)
		{
			continue;
		}
		else if(sFactor->vCnt>fSize)
		{
			break;
		}
		//Otherwise create a potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		populatePotential(aPotFunc,true);
		aPotFunc->calculateConditionalEntropy();
		double rInfo=(-1)*aPotFunc->getJointEntropy();
		for(int j=0;j<sFactor->vCnt;j++)
		{
			SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
			rInfo=rInfo+subFactor->conditionalEntropy;
		}
		randInfo.push_back(rInfo);
		rInd++;
		delete aPotFunc;
	}
	return Error::SUCCESS;

}


int 
PotentialManager::populateFactor(map<int,SlimFactor*>& factorSet,VSET& varSet,SlimFactor* sFactor,bool random)
{
	Potential* aPotFunc=new Potential;
	string fullConfStr;
	char confStr[CONSTR_LEN];
	for(int j=0;j<sFactor->vCnt;j++)
	{
		Variable* aVar=varSet[sFactor->vIds[j]];
		if(j==sFactor->vCnt-1)
		{
			
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		sprintf(confStr,"-%d",sFactor->vIds[j]);
		fullConfStr.append(confStr);
	}
	aPotFunc->potZeroInit();
	populatePotential(aPotFunc,random);
	aPotFunc->calculateConditionalEntropy();
	conditionalEntropies[fullConfStr]=aPotFunc->getConditionalEntropy();
	double rInfo=(-1)*aPotFunc->getJointEntropy();
	for(int j=0;j<sFactor->vCnt;j++)
	{
		SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
		rInfo=rInfo+subFactor->conditionalEntropy;
	}
	sFactor->mutualInfo=rInfo;
	sFactor->conditionalEntropy=aPotFunc->getConditionalEntropy();
	sFactor->jointEntropy=aPotFunc->getJointEntropy();
	delete aPotFunc;
	return 0;
}

int
PotentialManager::populatePotential(Potential* aPot, bool random)
{
	aPot->setEvidenceManager(evMgr);
	aPot->populateMe(0,1);
	return 0;
}

double
PotentialManager::getPseudoLikelihood(SlimFactor* sFactor,VSET& varSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


int
PotentialManager::populatePotential(Potential* aPot, bool random,double lambda)
{
	aPot->setEvidenceManager(evMgr);
	//We will bypass populatePotential
	aPot->populateMe(lambda,1);
	return 0;
}

double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet,map<int,int>& visitedVertices )
{
	double dll=0;
	cout <<"Not implemented" << endl;
	return dll;
}


int
PotentialManager::estimateConditionalPotential(SlimFactor* sFactor,VSET& varSet,Potential** pot, STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::populatePotential(Potential* pot,STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


int 
PotentialManager::estimateCanonicalPotential_Abbeel(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}



int 
PotentialManager::estimateCanonicalPotential_Approximate(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::resetPotFuncs()
{
	for(map<int,Potential*>::iterator pIter=potFuncs.begin();pIter!=potFuncs.end();pIter++)
	{
		delete pIter->second;
	}
	potFuncs.clear();
	return 0;
}



int 
PotentialManager::estimateCanonicalPotential_Joint(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

Potential*
PotentialManager::getPotential(int fId)
{
	if(potFuncs.find(fId)==potFuncs.end())
	{
		return NULL;
	}
	return potFuncs[fId];
}


double 
PotentialManager::getConditionalEntropy(int vId,INTINTMAP& fVars,VSET& varSet)
{
	double condEntropy=0;
	string fullConfStr;
	string partConfStr;
	char confStr[CONSTR_LEN];
	int varCnt=0;
	sprintf(confStr,"-%d",vId);
	fullConfStr.append(confStr);
	for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
	{
		if(aIter->first==vId)
		{
			continue;
		}
		sprintf(confStr,"-%d",aIter->first);
		fullConfStr.append(confStr);
	}
	if(conditionalEntropies.find(fullConfStr)!=conditionalEntropies.end())
	{
		condEntropy=conditionalEntropies[fullConfStr];
	}
	else
	{
		Potential* potFunc=new Potential;
		for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
		{
			if(aIter->first==vId)
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::FACTOR);
			}
			else
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::MARKOV_BNKT);
			}
		}
		potFunc->potZeroInit();
		populatePotential(potFunc,false);
		potFunc->calculateConditionalEntropy();
		condEntropy=potFunc->getConditionalEntropy();
		conditionalEntropies[fullConfStr]=condEntropy;
		delete potFunc;
	}
	return condEntropy;
}

int
PotentialManager::setConditionalEntropy(int vId,INTINTMAP& fVars,VSET& varSet, double condEntropy)
{
	string fullConfStr;
	char confStr[CONSTR_LEN];
	int varCnt=0;
	sprintf(confStr,"-%d",vId);
	fullConfStr.append(confStr);
	for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
	{
		if(aIter->first==vId)
		{
			continue;
		}
		sprintf(confStr,"-%d",aIter->first);
		fullConfStr.append(confStr);
	}
	conditionalEntropies[fullConfStr]=condEntropy;
	return 0;
}




double 
PotentialManager::getSampleLikelihood(map<int,SlimFactor*>& factorSet, VSET& varSet, INTINTMAP* sample)
{
	double sampleLL=0;
	cout <<"Not implemented " <<endl;
	return sampleLL;
}


int 
PotentialManager::getVariableSample(INTINTMAP& jointConf,VSET& varSet,int vId,SlimFactor* sFactor, gsl_rng* r)
{
	Potential* pot=NULL;
	if(potFuncs.find(vId)==potFuncs.end())
	{
		pot=new Potential;
		STRDBLMAP counts;
		estimateConditionalPotential(sFactor,varSet,&pot,counts);
		potFuncs[vId]=pot;
	}
	else
	{
		pot=potFuncs[vId];
	}
	int sample=pot->generateSample(jointConf,vId,r);
	return sample;
}

double
PotentialManager::estimateCanonicalValue(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* condPot)
{
	double pVal=0;
	cout <<"Not implemented " << endl;
	return 0;
}


double
PotentialManager::estimateCanonicalValue_Joint(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* jointPot)
{
	cout <<"Not implemented " << endl;
	return 0;
}
