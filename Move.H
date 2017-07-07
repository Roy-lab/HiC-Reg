/*The goal of this class is to describe a graph edit move we make to generate a new graph.
 * We will make a move on a set of variables, typically a pair of variables, but may be more
 * if we are doing exact search. We need to store the type of move and the score
 * because of this new move. The move could be an edge addition, edge deletion, edge reversal,
 * Markov blanket addition, which is equal to addition of a set of edges */
#ifndef _MOVE_
#define _MOVE_

#include "CommonTypes.H"

class Move
{
	public:
	Move();
	~Move();
	typedef enum
	{
		EXACT,
		GREEDY,
		GREEDY_APPROX
	} MoveType;
	int setScoreImprovement(double);
	int setMBScore(double);
	int setScore(double);
	int setSrcVertex(int);
	int setSrcMBID(int);
	int setTargetVertex(int,int);
	int setTargetset(INTINTMAP&);

	int getSrcVertex();
	int getSrcMBID();
	INTINTMAP& getTargetSet();
	double getScore();
	double getMBScore();
	double getScoreImprovement();
	private:
	double scoreDelta;
	double score;
	double mbscore;
	int src;
	int srcmbid;
	INTINTMAP targetSet;

};

#endif
