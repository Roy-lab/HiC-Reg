#include <iostream>
using namespace std;

#include "RegressionTree.H"
#include "Rule.H"

Rule::Rule()
{
}

Rule::~Rule()
{
}

//True means gt and false means lteq
int 
Rule::addCondition(RegressionTree* aNode, int test)
{	
	int cid=conditionSet.size();
	conditionSet[cid]=aNode;
	conditionTest[cid]=test;
	return 0;
}

map<int,RegressionTree*>& 
Rule::getAllConditions()
{
	return conditionSet;
}

int
Rule::setMarginalEntropy(double anEntropy)
{
	marginalEntropy=anEntropy;
	return 0;
}

double
Rule::getMarginalEntropy()
{
	return marginalEntropy;
}

int
Rule::setCoverage(int c)
{
	coverage=c;
}

int
Rule::getCoverage()
{
	return coverage;
}

int
Rule::showRule()
{
	int parentBranch=-1;
	int shownCondition=0;
	for(map<int,RegressionTree*>::iterator lIter=conditionSet.begin();lIter!=conditionSet.end();lIter++)
	{
		RegressionTree* node=lIter->second;
		if(conditionTest[lIter->first]==-1)
		{
			continue;
		}
		double testValue=node->getTestValue();
		int testVar=node->getTestVariable();
		int gt=conditionTest[lIter->first];
		if(shownCondition>0)
		{
			cout <<" && ";
		}
		if(gt==1)
		{
			cout << "(" << testVar << " > " << testValue << ")";
		}
		else if(gt==0)
		{
			cout << "("<< testVar<< "<= " << testValue << ")";
		}
		shownCondition++;
	}
	cout << endl;
	return 0;
}

int
Rule::getBranch(int cId)
{
	int bval=-1;
	if(conditionTest.find(cId)==conditionTest.end())
	{
		return bval;
	}
	bval=conditionTest[cId];
}

double
Rule::getRuleComplexity(int delMe)
{
	double complexity=0;
	for(map<int,RegressionTree*>::iterator cIter=conditionSet.begin();cIter!=conditionSet.end();cIter++)
	{
		if(cIter->first==delMe)
		{
			continue;
		}
		complexity=complexity+cIter->second->getCodingLength();
	}
	return complexity;
}
