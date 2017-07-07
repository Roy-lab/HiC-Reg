#ifndef _FACTOR_MANAGER
#define _FACTOR_MANAGER

#include <vector>
#include <fstream>
using namespace std;
#include "gsl/gsl_randist.h"
#include "CommonTypes.H"
#define SAMPLE_CNT 1000

class Error;
class VariableManager;
class EvidenceManager;
class PotentialManager;
class SlimFactor;
class LatticeStructure;
class Graph;
class FactorGraph;

class FactorManager
{
	public:
		FactorManager();
		~FactorManager();
		int setVariableManager(VariableManager*);
		int setEvidenceManager(EvidenceManager*);
		int setPotentialManager(PotentialManager*);
		int setBaseInstantiation();
		int setBaseInstantiation_Variable();
		int setMaxFactorSize(int);
		int setMaxFactorSize_Approx(int);
		int setRandMISdCnt(double);
		int getMaxFactorSize();
		int getMaxFactorSize_Approx();

		int setModelName(const char*);
		int setOutputDir(const char*);
		int setGraph(Graph*);
		int setBeamSize(int);
		int setPenalty(double);
		
		int allocateFactorSpace();
		int allocateFactorSpace_Graph();

		FactorGraph* createInitialFactorGraph();

		int paramEstimation(const char*);

		int generateSingletonFactors();
		int generateCanonicalFactors();
		int generateCanonicalFactors(FactorGraph*, map<int,SlimFactor*>& canonicalFactors);
		int estimateCanonicalParameters(const char*);

		int learnStructure();
		int showStructure();
		int showStructure_allK();
		
		int readStructure();
		int readClusterProperties(const char*);
		
		int estimateRandomInfo();
		int estimateRandomInfo_Approximate(int);
		int readRandomInfo();

		
		//This function is called to remove all factors that have less than the reqd. number of sds
		//random mutual information
		int filterClustersWithMI(double eps);
		
		int applyDPICorrection(double mithreshold, double dpiPercent);

		//This is used to check if the information between k+1 variables
		//is significantly greater than the mutual information between
		//any set of k variables. If so there is a higher order dependency
		//If not then there is no higher order dependency.
		int findTrueHO(double eps);

		int getFactorIndex(int*, int);
		int getFactorKey(int*, int,string&);
		int getMBFactorIndex(SlimFactor*);

		int addNewMBFactor(SlimFactor*);

		SlimFactor* getFactorFromVars(int*, int);
		double getMIFromVars(int*,int);

		int generateClusters(double epsilon, bool latticeCheck, double confidence,int);
		int getMaximalClusters();
		//Here we populate the factSignMap
		int removeDupFactors();

		int generateClusters_Apriori(double,double,int,bool);
		int learnMBStructure(double,double);
		double getPseudoLikelihood();
		double getLikelihood_ChainRule();
		double getLikelihood_ChainRule(FactorGraph*);
		double getLikelihood();
		double getLikelihood(FactorGraph*);

		double getLikelihood_MCMC(FactorGraph*);
		double estimatePartitionFunction(FactorGraph*);
		int randInitSample(INTINTMAP* sample,gsl_rng* r);
		int getNextSample(FactorGraph* fg,INTINTMAP* currSample,INTINTMAP* nextSample,gsl_rng* r);

		int evaluateMarkovBlanket(double);
		
		int showAllFactors(double);

		bool checkMonotonicity();

		INTINTMAP* getSupersets(int);
		int getSupersets(int,int,INTINTMAP&);
		double getMBScore(SlimFactor*,int supFid);
		double getMBScore(SlimFactor*);
		SlimFactor* getFactorAt(int);
		int generateNextLevelClusters(SlimFactor*,double,INTDBLMAP&);
		int generateNextLevelClusters(SlimFactor*);
		int generateNextLevelClusters(SlimFactor*,SlimFactor*);
		int showConditionalPotentials(FactorGraph*);
	private:
		int generateApproximateClusters(double,double);
		int generateNextLevelClusters(SlimFactor*,double,double,INTINTMAP&);
		int getBestMarkovBlanket(SlimFactor*,double,double);
		int makeMBMutuallyConsistent();
		int makeMarkovBlanketConsistent(SlimFactor*);
		int produceClusters(double,int);
		int produceClusters_NoDup(double,int);

		int getGoodCandidateMarkovBlankets(SlimFactor*,INTINTMAP*, int,double&,double,INTDBLMAP&);
			
		int getClusterCnt(int,int,int);
		int initFactorSet();
		int populateFactorSet();
		int checkFactorIds();
		int combCnt(int,int);
		int qsort(int*,int,int);
		//Create all subsets of the random variables
		//without taking into account the information
		//of what associations make sense
		Error::ErrorCode createAllFactors();
		Error::ErrorCode estimateClusterProperties();

		int addToLattice(SlimFactor*,int**);
		int deleteFromLattice(int);
		
		double getOverlap(SlimFactor*,SlimFactor*,int&);
		
		int globalFactorID;	
		int maxFactorSize;
		int maxFactorSize_Approx;
		VariableManager* vMgr;
		PotentialManager* potMgr;
		EvidenceManager* evMgr;

		//This is used to store the ids of factors which have only one variable
		//This will be used to compute the multi-information
		map<int,SlimFactor*> slimFactorSet;
		map<int,SlimFactor*> canonicalFactorSet;
		map<string,int> factorNameToIDMap;
		map<int,string> factorIDToNameMap;
		//Because we delete the old factors, by our pruning step, but
		//still need the multi-information associated with the factor
		//we store it here.
		map<string,double> delFactors_MI;
		
		int* pairIds;
		int* clusterIds;

		map<int,INTDBLMAP*> markovBlanketSet;
		INTINTMAP defaultInstantiation;
		LatticeStructure lattice;
		map<int,SlimFactor*> goodSlimFactors;
		
		map<int,double> randMI_mean;
		map<int,double> randMI_std;
		char outputDir[1024];
		char modelName[1024];

		//Store a map of factor signature and the set of factors that have that signature.
		//Factor signature is defined by the union of variables in the factor and the variables in its Markov
		//blanket. We must not generate clusters for every possible factor and its Markov blanket
		//but per signature. The question we have to address here is the following:
		//Suppose we have two factors F and G with Markov blankets Mf and Mg respectively such that
		//F \Union Mf= G \Union Gf.  Should we pick F \Union Mf or should we pick G \Union Gf.
		//So which of these decompositions estimate the joint F\Union Mf better. We will pick
		//the decomposition that minimizes the interaction between the variables of Mf or in other
		//words the decomposition that lowers the conditional entropy the most. So if H(G|Mg) < H(F|Mf)
		//then we pick decomposition from G or we pick decomposition from F. Unfortunately, this may 
		//increase the interaction between the Markov blanket and the outside world.
		map<string,INTDBLMAP*> signFactGrpMap;
		map<int,string> factSignMap;
		map<string,int> groupFactorRepMap;
		map<string,double> mbSpecific_MI;

		//Only used to pairwise format of factors
		map<string,int> shownPairs;
		Graph* graph;
		double defProb;
		int defInstID;
		INTINTMAP defInstMap;
		int beamSize;
		double mbPenalty;
		double misdCnt;
};
#endif
