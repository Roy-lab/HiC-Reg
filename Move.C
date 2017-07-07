#include "Move.H"
Move::Move()
{
}

Move::~Move()
{
}

int 
Move::setScoreImprovement(double aVal)
{
	scoreDelta=aVal;
	return 0;
}

int
Move::setScore(double aScore)
{
	score=aScore;
	return 0;
}

int
Move::setMBScore(double aScore)
{
	mbscore=aScore;
	return 0;
}

int 
Move::setSrcVertex(int vid)
{
	src=vid;
	return 0;
}

int 
Move::setSrcMBID(int mbid)
{
	srcmbid=mbid;
	return 0;
}

int 
Move::setTargetVertex(int vid,int mbid)
{
	targetSet[vid]=mbid;
	return 0;
}

int 
Move::setTargetset(INTINTMAP& vSet)
{
	for(INTINTMAP_ITER vIter=vSet.begin();vIter!=vSet.end();vIter++)
	{
		targetSet[vIter->first]=vIter->second;
	}

	return 0;
}


int 
Move::getSrcVertex()
{
	return src;
}

int
Move::getSrcMBID()
{
	return srcmbid;
}

INTINTMAP& 
Move::getTargetSet()
{
	return targetSet;
}

double 
Move::getScore()
{
	return score;
}

double 
Move::getMBScore()
{
	return mbscore;
}


double 
Move::getScoreImprovement()
{
	return scoreDelta;
}
