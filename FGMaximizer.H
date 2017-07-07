/*The goal of this class is to identify the best undirected graph for a set of random variables. 
 * By best, we mean the one in which each variable has its best Markov blanket and that has
 * mutually consistent Markov blankets. The best Markov blanket in turn is measured by the reduction
 * in conditional entropy per variable. This class makes use of the FGEditor class that implements 
 * the edge driven score.
 *
 *
 * */
#ifndef _FG_MAXIMIZER
#define _FG_MAXIMIZER

class FactorGraph;
class FGEditor;
class FactorManager;
class Move;
class PotentialManager;
class RegressionTree;

class FGMaximizer
{
	public:
	FGMaximizer();
	~FGMaximizer();
	//How many mbs to use to do the greedy, approximate search of MB vars 
	int setMBCntForApproxGreedy(int);
	//How many top graphs we should return
	int setBeamSize(int);
	int setFactorManager(FactorManager* fMgr);
	int setVariableManager(VariableManager*);
	int setPotentialManager(PotentialManager*);
	
	// SR added to bypass potmag
	int setTrainEvidenceManager(EvidenceManager*);
	// DC added
	int setTestEvidenceManager(EvidenceManager*);
	
	int setTreeCnt(int);
	int setOutputDir(const char*);
	int setConvergenceThreshold(double);
	int setMinLeafSize(int);
	int setPriorGraph(const char*);
	int findBestGraphs();
	int findBestGraphs_TopDown(double,const char*);
	int findBestGraphs_PriorGraph(double,const char*);
	int readBestGraphs_PriorGraph(double,const char*);
	int showOutput();
	
	// read column names from file -- DC added
	int readColumnNames(const char*);
	int setColumnNames(vector<string>* columnNames);

	int setTreeLocation(const char*);
	private:
	double getScore(FactorGraph*);
	int sortMoves(vector<Move*>&);
	int updateMoveOrder(int, Move*,vector<Move*>& moveSet,int&);
	int attemptMove(FactorGraph*,int,vector<Move*>&,map<int,int>&,map<int,int>&);
	int generateNextLevelClusters(int,FactorGraph*);

	int dispRegTree_Genanatomy(map<int,RegressionTree*>&,const char*);
	FGEditor fgeditor;
	FactorManager* fMgr;
	FactorGraph* fg;
	VariableManager* vMgr;
	PotentialManager* potMgr;
	
	EvidenceManager* trainEvMgr; 
	// DC added to hold test data
	EvidenceManager* testEvMgr; 
	
	int maxMBSize;
	int maxMBSize_Approx;
	int mbCntForApproxGreedy;
	int beamSize;
	int minLeafSize;
	char outputDir[1024];
	char treeloc[1024];
	
	// column names - ordered by input
	// DC added
	vector<string>* columnNames;
	
	double convThreshold;
	map<int,int> nextLevelFlag;
	Graph priorGraph;
	
	int treeCnt;
};
#endif
