#include <fstream>
#include <iostream>
#include <cstring>
#include <string.h>
#include <unistd.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "LatticeStructure.H"

#include "Vertex.H"
#include "Graph.H"

#include "FactorGraph.H"
#include "FactorManager.H"
#include "PotentialManager.H"
#include "Move.H"
#include "FGEditor.H"
#include "FGMaximizer.H"

#include "Framework.H"

Framework::Framework()
{
	epsThreshold=-1;
	constrained=true;
	predictionMode=false;
}

Framework::~Framework()
{
}

//We will use getopt here
//The options are 
//-t training_data 
//-o outputdir
//-e epsilon to control the number of standard deviations above random
//-k maxfactorsize
//-s number of samples for approximate information estimation
//-x k for which we approximate information
//-n number of trees in forest
// DC added:
// -c columnNames file that gives names of data columns (lost in fgconverting)
// -d test data

Error::ErrorCode
Framework::init(int argc, char** argv)
{
	evManager.setVariableManager(&varManager);
	potManager.setEvidenceManager(&evManager);
	fMgr.setVariableManager(&varManager);
	fgMax.setVariableManager(&varManager);
	fgMax.setPotentialManager(&potManager);
	fgMax.setTrainEvidenceManager(&evManager);
	fMgr.setEvidenceManager(&evManager);
	fMgr.setPotentialManager(&potManager);
	
	// check if minleaf size is compatible with data
	int leafcheck=0;
	int datasize=0;
	bool doTest=false; // test set provided
	
	int optret='-';
	opterr=1;
	int oldoptind=optind;
	while(optret=getopt(argc,argv,"c:t:o:k:e:x:p:l:b:u:d:n:s:")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return Error::UNKNOWN;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 't':
			{
				char fName[256];
				sprintf(fName,"%s",my_optarg);
				/*sprintf(fName,"%s.model",my_optarg);
				Error::ErrorCode eCode=varManager.readVariables(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				sprintf(fName,"%s.data",my_optarg);
				// evidence stored here
				eCode=evManager.loadEvidenceFromFile_Continuous(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}*/
				Error::ErrorCode eCode=evManager.loadEvidenceFromFile_Simple(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				
				datasize=evManager.getNumberOfEvidences();
				break;
			}
			// held-aside test data
			case 'd':
			{
				testEvManager.setVariableManager(&testVarManager);	
				fgMax.setTestEvidenceManager(&testEvManager);
				
				//potManager.setTestEvidenceManager(&testEvManager);
				// who owns the testEvManager?

				char fName[256];
				sprintf(fName,"%s",my_optarg);
				/*sprintf(fName,"%s.model",my_optarg);
				Error::ErrorCode eCode=testVarManager.readVariables(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				sprintf(fName,"%s.data",my_optarg);
				// evidence stored here
				eCode=testEvManager.loadEvidenceFromFile_Continuous(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}*/
				Error::ErrorCode eCode=testEvManager.loadEvidenceFromFile_Simple(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				
				cout << "Read test data from " << fName << endl;
				doTest=true;
				break;
			}
			case 'o':
			{
				fMgr.setOutputDir(my_optarg);
				fgMax.setOutputDir(my_optarg);
				break;
			}
			case 'p':
			{
				penalty=atof(my_optarg);
				break;
			}
			case 'e':
			{
				epsThreshold=atof(my_optarg);
				fMgr.setRandMISdCnt(epsThreshold);
				break;
			}
			case 'k':
			{
				int aSize=atoi(my_optarg);
				fMgr.setMaxFactorSize(aSize);
				break;
			}
			case 'x':
			{
				int aSize=atoi(my_optarg);
				fMgr.setMaxFactorSize_Approx(aSize);
				break;
			}
			case 'l':
			{
				int leafSize=atoi(my_optarg);
				fgMax.setMinLeafSize(leafSize);
				leafcheck=leafSize;
				break;
			}
			case 'b':
			{
				fgMax.setPriorGraph(my_optarg);
				break;
			}
			case 'u':
			{
				if(strcmp(my_optarg,"yes")==0)
				{
					constrained=true;
				}
				else
				{
					constrained=false;
				}
				break;
			}
			case 'c':
			{	
				// DC added - read column names
				readColumnNames(my_optarg);
				break;
			}
			case 'n':
			{
				fgMax.setTreeCnt(atoi(my_optarg));	
				break;
			}
			case 's':
			{
				fgMax.setTreeLocation(my_optarg);
				predictionMode=true;
				break;
			}
			default:
			{
				cout <<"Unhandled option " << c  << endl;
				return Error::UNKNOWN;
			}
		}
		oldoptind=optind;
	}
	// DC ADDED
	// Check if minleafsize is small enough for data
	if (leafcheck>=(datasize-2))
	{
		cerr << "Please specify minleafsize < (number of datapoints-2)." << endl;
		cerr << "Provided minleafsize=" << leafcheck <<", data size=" << datasize << endl;
		cerr << "Exiting." << endl;
		return Error::DATAFILE_ERR;
	}
	
	// DC added
	// if test variable manager available, make sure is compatible with data's var manager
	// this means checking to make sure they have the same variable names.
	if (doTest)
	{
		bool compatible=varManager.hasSameVars(&testVarManager);
		//cout << compatible << endl;
		if (!compatible)
		{
			cerr << "Test data file has incompatible variables with train data file. Exiting." << endl;
			return Error::DATAFILE_ERR;
		}
	}
	
	return Error::SUCCESS;
}

int 
Framework::start()
{
	// give fgmaximizer the column names
	fgMax.setColumnNames(&columnNames);
	
	fMgr.allocateFactorSpace();
/*	if(fMgr.readStructure()==-1)
	{
		fMgr.learnStructure();
		fMgr.showStructure_allK();
	}
	if(fMgr.readRandomInfo()==-1)
	{
		cout <<"Did not find random mutual informations. Calling estimateRandomInfo " << endl;
		fMgr.estimateRandomInfo_Approximate(SAMPLE_CNT);
	}*/
	fgMax.setFactorManager(&fMgr);
	if(constrained)
	{
		if(predictionMode)
		{
			if(fgMax.readBestGraphs_PriorGraph(penalty,projectFName)==-1)
			{
				cout <<"Could not find best graphs " << endl;
				return 0;
			}
		}
		else
		{
			if(fgMax.findBestGraphs_PriorGraph(penalty,projectFName)==-1)
			{
				cout <<"Could not find best graphs " << endl;
				return 0;
			}
		}
	}
	else
	{
		if(fgMax.findBestGraphs_TopDown(penalty,projectFName)==-1)
		{
			cout <<"Could not find best graphs " << endl;
			return 0;
		}
	}
        return 0;
}

// DC added this -- reads in column names from file
int Framework::readColumnNames(const char* fName)
{
	ifstream inFile(fName);
	char* buffer=NULL;
	int bufflen=0;
	string strbuff;
	int lines=0;
	bool success=true; // false if time order empty or doesn't match matrix
	
	while(inFile.good())
	{
		getline(inFile,strbuff);
		if(strbuff.length()<=0)
		{
			continue;
		}
		if(bufflen<=strbuff.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=strbuff.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,strbuff.c_str());
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;			
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				string name(tok);
				// store in member var
				columnNames.push_back(name);		
				//cout << name << endl;		
			}			
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if (!success)
		{
			break;
		}
		lines++;		
	}	
	inFile.close();
	
	if (success && lines==0) 
	{
		cerr << "Read zero lines from column file: " << fName << endl;
		return 1;
	}	
	else if (success)	
	{
		cout << "Successfully read " << lines << " column names." << endl;
		return 0;
	}
	else 
	{
		cout << "Some error reading column names from: " << fName << endl;
		return 1;
	}
}


int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"factorGraphInf " <<  endl
			<<"-t training " << endl
			<< "-o outputdir " << endl
			<< "-e epsilon [>=0]" << endl 
			 << "-p penalty"<< endl
			<< "-k maxfactorsize " << endl
			 << "-x maxfactorsize_approx" << endl
			 << "-t convergence_threshold" << endl
			 << "-l leafsize"<< endl
			 << "-b priornet"<< endl
			 << "-u [yes|no]"<< endl
			 << "-n treecnt"<< endl
			 // DC added
			 << "-s treeFilePrefix" << endl
			 << "-c dataColumnNameFile" << endl
			 << "-d testdata " << endl;

		return 0;
	}
	Framework fw;
	if(fw.init(argc,argv)!=Error::SUCCESS)
	{
		return 0;
	}
	fw.start();
	return 0;
	
}

