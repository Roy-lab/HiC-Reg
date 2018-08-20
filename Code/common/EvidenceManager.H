#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <map>
#include <fstream>

#include "gsl/gsl_randist.h"
using namespace std;

class Error;
class Evidence;

//EMAP stores the map of variable ids and their evidence. This corresponds
//to one line in the datafile
typedef map<int,Evidence*> EMAP;
typedef map<int,Evidence*>::iterator EMAP_ITER;

typedef vector<EMAP*> EVIDENCE_SET;

class VariableManager;
class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();
		//This associates the variable information with the variable id
		int setVariableManager(VariableManager*);
		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);
		Error::ErrorCode loadEvidenceFromFile_Continuous(const char*);
		//This file will just load a matrix file with columns corresponding to variables, first col is the 
		//name of the sample. 
		Error::ErrorCode loadEvidenceFromFile_Simple(const char*);
		int randomizeEvidence(gsl_rng*);
		int getNumberOfEvidences();
		//Get the evidence set at the specified id
		EMAP* getEvidenceAt(int);
		EMAP* getRandomEvidenceAt(int);
		//This is used during EM, updating missing/hidden values using 
		//expected counts. eSetId corresponds to the data entry, vId is the
		//variable whose evidence we are updating and evidData is the set
		//of values and their evidence values (soft evidence)
		int addToEvidence(int eSetId,int vId, INTDBLMAP& evidData);	
		int dumpEvidenceSet(ostream&);
		int dumpSummaryStat(ostream&);
		int getMLSettings(ostream&);
		
		int dumpEvidenceTab(const char*);
		const char* getSampleName(int);
	private:
		int populateEvidence(Evidence**,const char*);
		int populateEvidence_Continuous(Evidence**,const char*);
		int populateRandIntegers(gsl_rng*, int*,int);
		EVIDENCE_SET evidenceSet;
		EVIDENCE_SET randEvidenceSet;
		VariableManager* vMgr;
		//Suitable for the simple matrix format
		map<int,string> sampleNames;
};
#endif
