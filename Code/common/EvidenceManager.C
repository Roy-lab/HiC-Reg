#include <fstream>
#include <iostream>
#include <math.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"


EvidenceManager::EvidenceManager()
{
}

EvidenceManager::~EvidenceManager()
{
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile(const char* inFName)
{
	ifstream inFile(inFName);
	char buffer[80000];
	while(inFile.good())
	{
		inFile.getline(buffer,80000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		else if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value
		while(tok!=NULL)
		{
			Evidence* evid;
			if(populateEvidence(&evid,tok)==-1)
			{
				cout <<"Error while populating evidence " << endl;
				return Error::DATAFILE_ERR;
			}
			(*evidMap)[evid->getAssocVariable()]=evid;
			tok=strtok(NULL,"\t");
		}
		//if(evidenceSet.size()<=20)
		//{
			evidenceSet.push_back(evidMap);
		//}
	}
	inFile.close();
	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;
	return Error::SUCCESS;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Continuous(const char* inFName)
{
	ifstream inFile(inFName);
	string buffstr;
	char* buffer=NULL;
	int bufflen=0;
	int lineNo=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		else if(strchr(buffstr.c_str(),'#')!=NULL)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(lineNo>=500)
		{
			lineNo++;
		//	continue;
		}

		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value
		while(tok!=NULL)
		{
			Evidence* evid;
			if(populateEvidence_Continuous(&evid,tok)==-1)
			{
				cout <<"Error while populating evidence " << endl;
				return Error::DATAFILE_ERR;
			}
			(*evidMap)[evid->getAssocVariable()]=evid;
			tok=strtok(NULL,"\t");
		}
		evidenceSet.push_back(evidMap);
	}
	if(buffer!=NULL)
	{
		delete[] buffer;
	}
	inFile.close();
	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;
	return Error::SUCCESS;
}


Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Simple(const char* inFName)
{
	ifstream inFile(inFName);
	string buffstr;
	char* buffer=NULL;
	int bufflen=0;
	int lineNo=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		else if(strchr(buffstr.c_str(),'#')!=NULL)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		if(lineNo==0)
		{
			//Reading in the header names. Will also populate the varManager here
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt>0)
				{
					Error::ErrorCode eCode=vMgr->addVariable(tok,tokCnt-1);		
					if(eCode!=Error::SUCCESS)
					{
						cout <<"Unable to add variable! "<< endl;
						return eCode;
					}
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
		}
		else
		{
			//All the evidences for each variable are stored in a map, indexed by the varId
			EMAP* evidMap=new EMAP;
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			string sName;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					sName.append(tok);
				}
				else
				{
					Evidence* evid=new Evidence;
					Variable* var=vMgr->getVariableAt(tokCnt-1);
					evid->assocVariable(var->getID());
					double varVal=atof(tok);
					evid->setEvidVal(varVal);
					(*evidMap)[evid->getAssocVariable()]=evid;
				}
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
			sampleNames[evidenceSet.size()]=sName;
			evidenceSet.push_back(evidMap);
		}
		lineNo++;
	}
	if(buffer!=NULL)
	{
		delete[] buffer;
	}
	inFile.close();
	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;
	return Error::SUCCESS;
}


//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	VSET& variableSet=vMgr->getVariableSet();
	int* randInds=new int[evidenceSet.size()];
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,randEvidenceSet.size());	
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[randInds[i]];
			EMAP* randEvidMap=randEvidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			(*randEvidMap)[vIter->first]=evid;
		}
	}
	return 0;
}

int 
EvidenceManager::getNumberOfEvidences()
{
	return evidenceSet.size();
}

EMAP* 
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

EMAP* 
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}

int
EvidenceManager::addToEvidence(int eSetID, int vId, INTDBLMAP& evidData)
{
	EMAP* emap=evidenceSet[eSetID];
	Evidence* evid=new Evidence;
	evid->setData(evidData);
	(*emap)[vId]=evid;
	return 0;
}

int
EvidenceManager::dumpEvidenceSet(ostream& oFile)
{
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			Evidence* evid=eIter->second;
			if(eIter!=evidMap->begin())
			{
				oFile<<"\t";
			}
			evid->dumpEvidence(oFile);
		}
		oFile << endl;
	}
	return 0;
}

int
EvidenceManager::getMLSettings(ostream& oFile)
{
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		if(i==0)
		{
			for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
			{
				if(eIter!=evidMap->begin())
				{
					oFile<<"\t";
				}
				oFile<< eIter->first;
			}
			oFile << endl;
		}

		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			if(eIter!=evidMap->begin())
			{
				oFile<<"\t";
			}
			Evidence* evid=eIter->second;
			oFile << evid->getMLVal();
		}
		oFile << endl;
	}
	return 0;
}

int
EvidenceManager::dumpEvidenceTab(const char* aFName)
{
	ofstream dataFile(aFName);
	int evidCnt=evidenceSet.size();
	for(int e=0;e<evidCnt;e++)
	{
		if(e==0)
		{
			dataFile <<"Exp" << e;
		}
		else 
		{
			dataFile <<" Exp" << e;
		}
	}
	dataFile << endl;
	VSET& varSet=vMgr->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		dataFile << vIter->second->getName();
		for(int e=0;e<evidCnt;e++)
		{
			EMAP* evidSet=getEvidenceAt(evidCnt);
			Evidence* evid=(*evidSet)[vIter->first];
			dataFile << "\t"<< evid->getEvidVal();
		}
		dataFile << endl;
	}
	dataFile.close();
	return 0;
}

int 
EvidenceManager::populateEvidence(Evidence** evid,const char* evidStr)
{
	//first check for validity of evidStr
	if(strchr(evidStr,'=')==NULL)
	{
		return -1;
	}
	*evid=new Evidence;
	
	INTDBLMAP evidData;
	int currInd=0;
	int ttInd=0;
	int tokId=0;
	char tempTok[256];
	while(evidStr[currInd]!='\0')
	{
		if((evidStr[currInd]=='=') || 
		   (evidStr[currInd]==']') ||
		   (evidStr[currInd]==',')
		  )
		{
			tempTok[ttInd]='\0';
			ttInd=0;
			if(tokId==0)
			{
				//This is the variable
				int vId=atoi(tempTok);
				Variable* var=vMgr->getVariableAt(vId);
				var->initEvidence(evidData);
				(*evid)->assocVariable(vId);
			}
			else
			{
				char* pos=strchr(tempTok,'|');
				//Hard evidence
				if(pos==NULL)
				{
					int varVal=atoi(tempTok);
					if(evidData.find(varVal)==evidData.end())
					{
						cout <<"No value "<< varVal << " in the domain of  a variable" << endl;
						return -1;
					}
					evidData[varVal]=1.0;
					(*evid)->setType(Evidence::HARD);
				}
				else
				{
					*pos='\0';
					int varVal=atoi(tempTok);
					double varValProb=atof(pos+1);
					if(evidData.find(varVal)==evidData.end())
					{
						cout <<"No value "<< varVal << " in the domain of  a variable" << endl;
						return -1;
					}
					evidData[varVal]=varValProb;
					//Will be setting it multiple times but its ok for now.
					(*evid)->setType(Evidence::SOFT);
				}
			}
			tokId++;
		}
		else if(evidStr[currInd]!='[')
		{
			tempTok[ttInd]=evidStr[currInd];
			ttInd++;
		}
		currInd++;
	}
	(*evid)->setData(evidData);
	return 0;
}


int 
EvidenceManager::populateEvidence_Continuous(Evidence** evid,const char* evidStr)
{
	//first check for validity of evidStr
	if(strchr(evidStr,'=')==NULL)
	{
		return -1;
	}
	*evid=new Evidence;
	
	int currInd=0;
	int ttInd=0;
	int tokId=0;
	char tempTok[256];
	while(evidStr[currInd]!='\0')
	{
		if((evidStr[currInd]=='=') || 
		   (evidStr[currInd]==']') ||
		   (evidStr[currInd]==',')
		  )
		{
			tempTok[ttInd]='\0';
			ttInd=0;
			if(tokId==0)
			{
				//This is the variable
				int vId=atoi(tempTok);
				Variable* var=vMgr->getVariableAt(vId);
				(*evid)->assocVariable(vId);
			}
			else
			{
				char* pos=strchr(tempTok,'|');
				//Hard evidence
				if(pos==NULL)
				{
					//double varVal=log(atof(tempTok));
					double varVal=atof(tempTok);
					(*evid)->setEvidVal(varVal);
				}
			}
			tokId++;
		}
		else if(evidStr[currInd]!='[')
		{
			tempTok[ttInd]=evidStr[currInd];
			ttInd++;
		}
		currInd++;
	}
	return 0;
}

int
EvidenceManager::dumpSummaryStat(ostream& oFile)
{
	//Need an evidence like object but to store the frequency over all evidences
	//rather than a single evidence
	map<int,INTDBLMAP*> summary;
	map<int,double> normFactors;
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			INTDBLMAP* evCnt=NULL;
			if(summary.find(eIter->first)==summary.end())
			{
				evCnt=new INTDBLMAP;
				summary[eIter->first]=evCnt;
			}
			else
			{
				evCnt=summary[eIter->first];
			}
			//Get data and add to evCnt
			INTDBLMAP& data=eIter->second->getData();
			for(INTDBLMAP_ITER idIter=data.begin();idIter!=data.end();idIter++)
			{
				if(evCnt->find(idIter->first)==evCnt->end())
				{
					(*evCnt)[idIter->first]=idIter->second;
				}
				else
				{
					(*evCnt)[idIter->first]=(*evCnt)[idIter->first]+idIter->second;
				}
				//Add the normalization factor for all the freq or exp. freq cnts
				if(normFactors.find(eIter->first)==normFactors.end())
				{
					normFactors[eIter->first]=idIter->second;
				}
				else
				{
					normFactors[eIter->first]=normFactors[eIter->first]+idIter->second;
				}
			}
		}
	}

	//Now iterate over the evidence summary, normalize and display values
	for(map<int,INTDBLMAP*>::iterator aIter=summary.begin();aIter!=summary.end();aIter++)
	{
		double normConst=normFactors[aIter->first];
		INTDBLMAP* evCnt=aIter->second;
		oFile <<"Distribution of "<< aIter->first;
		for(INTDBLMAP_ITER idIter=evCnt->begin();idIter!=evCnt->end();idIter++)
		{
			oFile << " " << idIter->first<<"=" << idIter->second/normConst;
		}
		oFile << endl;
	}
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}

const char*
EvidenceManager::getSampleName(int evidID)
{
	if(sampleNames.find(evidID)==sampleNames.end())
	{	
		return NULL;
	}
	const char* name=sampleNames[evidID].c_str();
	return name;
	return 0;
}
