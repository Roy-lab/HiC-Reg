/* The goal of this class is to use the edge-based score to find the best Markov blanket
 * of size k for each variable.
 * Operations from this class update the candidate MB set.
 * The score we are going to use is edge-driven.  That is, we want to compute the sum of the effect of adding
 * the edge both on the variable for which we are searching the MB as well as the variable which is the in MB.
 * This is how the search will proceed: we will begin with the current variable set and their Markov blankets.
 * Then we will propose a  move for each variable, keeping track of the delta improvement for the variable.
 * Then we will make the moves for each variable in the order of decreasing delta. This will move the non-optimal
 * moves if any to variables for which we do not have a lot of improvement anyways.*/
#ifndef _FG_EDITOR
#define _FG_EDITOR
class FactorManager;
class Move;

class FGEditor
{
	public:
	FGEditor();
	~FGEditor();
	int setFactorGraph(FactorGraph*);
	int setMBSize_Exact(int);
	//This is mbset size we should not exceed
	int setMaxMBSize(int);
	int incrMaxMBSize();
	int setMBCntForApproxGreedy(int);
	int setFactorManager(FactorManager*);
	int setMoveType(Move::MoveType);
	int generateMoves();
	int clearOldMoves();
	vector<Move*>& getMoves();
	
	int makeMove(FactorGraph*,Move*,map<int,int>&, map<int,int>&);
	Move* getNewMove(FactorGraph*,Move*);

	private:
	int identifyMove(SlimFactor*);
	Move* identifyConstrainedMove(SlimFactor*);
	double getScoreDelta_Exact(SlimFactor* sFactor, SlimFactor* supFactor, double& newCondEntr,double&,INTINTMAP&);
	double getScoreDelta_Greedy(SlimFactor* sFactor, SlimFactor* supFactor, double& newCondEntr,double&,INTINTMAP&);
	double getScoreChange_Neighbour(SlimFactor* sFactor,INTINTMAP& ssets,int&);
	double getScoreChange_NeighbourGreedy(SlimFactor* mbFactor,SlimFactor* sFactor);
	double getLossScore(SlimFactor*);
	FactorManager* fgMgr;
	FactorGraph* fg;
	int maxMBSize_Exact;
	int maxMBSize;
	int mbCntForApproxGreedy;
	Move::MoveType moveType;
	vector<Move*> moveSet;
	map<string,double> gainExactMB;
};
#endif
