/*The class which brings it all together*/
#ifndef _FRAMEWORK
#define _FRAMEWORK

//All the managers need forward declarations
class Error;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class FactorManager;
class FGMaximizer;

class Framework
{
	public:
		Framework();
		~Framework();
		Error::ErrorCode init(int, char**);
		int start();

	private:
		VariableManager varManager;
		EvidenceManager evManager;
		
		VariableManager testVarManager; // DC added - test data
		EvidenceManager testEvManager; // DC added - test data
		
		PotentialManager potManager;
		FactorManager fMgr;
		FGMaximizer fgMax;
		double epsThreshold;
		char actionStr[256];
		double penalty;
		bool constrained;
		bool predictionMode;
		char projectFName[1024];
		
		int readColumnNames(const char*);
		
		// DC added - data column names
		vector<string> columnNames;
};
#endif
