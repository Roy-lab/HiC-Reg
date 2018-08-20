#ifndef _RULE_
#define _RULE_

#include <map>
using namespace std;

class RegressionTree;
class Rule
{
	public:
	Rule();
	~Rule();
	int addCondition(RegressionTree*,int);
	map<int,RegressionTree*>& getAllConditions();
	int setMarginalEntropy(double);
	double getMarginalEntropy();
	int setCoverage(int);
	int getCoverage();
	int showRule();
	int getBranch(int);
	double getRuleComplexity(int);
	private:
	map<int,RegressionTree*> conditionSet;
	map<int,int> conditionTest;
	double marginalEntropy;
	int coverage;
	
};
#endif
