#include "Error.H"
#include "Variable.H"
#include "SlimFactor.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "VariableManager.H"

#include "LatticeStructure.H"

#include "PotentialManager.H"
#include "FactorGraph.H"

#include "Vertex.H"
#include "Graph.H"
#include "FactorManager.H"
#include "Move.H"
#include "FGEditor.H"


FGEditor::FGEditor()
{
}

FGEditor::~FGEditor()
{
}

int 
FGEditor::setFactorGraph(FactorGraph* aPtr)
{
	fg=aPtr;
	return 0;
}

int 
FGEditor::setMBSize_Exact(int kexact)
{
	maxMBSize_Exact=kexact;
	return 0;
}

//This is mbset size we should not exceed
int 	
FGEditor::setMaxMBSize(int k)
{
	maxMBSize=k;
	return 0;
}

int
FGEditor::incrMaxMBSize()
{
	maxMBSize++;
	return 0;
}

int
FGEditor::setMBCntForApproxGreedy(int aCnt)
{
	mbCntForApproxGreedy=aCnt;
	return 0;
}

int 
FGEditor::setFactorManager(FactorManager* aPtr)
{
	fgMgr=aPtr;
	return 0;
}

int 
FGEditor::setMoveType(Move::MoveType mtype)
{
	moveType=mtype;
	return 0;
}

int 
FGEditor::generateMoves()
{
	gainExactMB.clear();
	for(int i=0;i<fg->getFactorCnt();i++)
	{
		SlimFactor* sFactor=fg->getFactorAt(i);
		if(sFactor->mergedMB.size()>=maxMBSize)
		{
			continue;
		}
		identifyMove(sFactor);	
	}
	return 0;
}

int
FGEditor::clearOldMoves()
{
	for(int i=0;i<moveSet.size();i++)
	{
		delete moveSet[i];
	}
	moveSet.clear();
	return 0;
}

vector<Move*>&
FGEditor::getMoves()
{
	return moveSet;
}

//This function makes the move m on the factor graph fg
//What this means is that it is going to update the Markov blanket
//of the src vertex and the target vertices 
//There are three types of target vertices here
//a) Those whose MBs have been updated completely because they preceede the src vertex
//b) Those whose MBs have been affected because of some other src vertex
//c) Those that have not been affected yet.
//However vertices of type (a) will not arise here because of our check in the FGMaximizer of the
//affectedVertices. For (b) and (c) what we do depends upon if we are doing exact or greedy search.
//If we are doing greedy search we do one-variable extensions so we will simply return from
//this function if any target vertex's MB size >= maxMBSize.
//If we are doing exact search, then for (c) we must actually not consider the current MB size of the
//variable. If we do reach the point of the function of making the move, then we must first
//clear the vertex's MB and make addition of the src vertex to its MB.
//This move can fail if any of the target vertices cannot be updated 
//because it violates the current maximum markov blanket size

int 
FGEditor::makeMove(FactorGraph* newfg,Move* m, map<int,int>& srcVertex, map<int,int>& targetVertex)
{
	SlimFactor* srcFactor=newfg->getFactorAt(m->getSrcVertex());
	if(srcFactor->mergedMB.size()>=maxMBSize)
	{
		return 0;
	}
	INTINTMAP& targetVars=m->getTargetSet();
	int nonRemovals=0;
	for(INTINTMAP_ITER vIter=targetVars.begin();vIter!=targetVars.end();vIter++)
	{
		SlimFactor* sFactor=newfg->getFactorAt(vIter->first);
		if(moveType==Move::EXACT)
		{
			//This variable's MB will be erased as part of making the move so we don't care
			//about it at this point
			if((srcVertex.find(vIter->first)==srcVertex.end()) && (targetVertex.find(vIter->first)==targetVertex.end()))
			{
				continue;
			}
			else
			{
				if(sFactor->mergedMB.find(m->getSrcVertex())!=sFactor->mergedMB.end())
				{
					nonRemovals++;
					continue;
				}
				if(sFactor->mergedMB.size()>=maxMBSize)
				{
					cout << "MB size limit reached " << vIter->first << endl;
					return -1;
				}
			}
		}
		//Doing greedy search of one variable extensions to current MB
		else
		{
			if(sFactor->mergedMB.find(m->getSrcVertex())!=sFactor->mergedMB.end())
			{
				continue;
			}
			if(sFactor->mergedMB.size()>=maxMBSize)
			{
				cout << "MB size limit reached " << vIter->first << endl;
				return -1;
			}
		}
	}
	//If we are doing a greedy search we just add to the existing
	//MBs. For exact search we erase current MBs of the src vertex and target vertices NOT in the affectedVars set 
	//and then make updates
	if(moveType==Move::EXACT)
	{
		//Need to get rid of srcFactor variable from the MBs of all the old MB vars
		vector<INTINTMAP_ITER> toErase;
		for(INTINTMAP_ITER vIter=srcFactor->mergedMB.begin();vIter!=srcFactor->mergedMB.end();vIter++)
		{
			if(srcVertex.find(vIter->first)!=srcVertex.end())
			{
				continue;
			}
			SlimFactor* mbFactor=newfg->getFactorAt(vIter->first);
			INTINTMAP_ITER dIter=mbFactor->mergedMB.find(srcFactor->fId);
			if(dIter==mbFactor->mergedMB.end())
			{
				cout << " FATAL ERROR!! No variable " << srcFactor->fId << " in MB of " << mbFactor->fId << endl;
				return -2;
			}
			mbFactor->mergedMB.erase(dIter);
			mbFactor->goodMBIDs.clear();
			if(mbFactor->mergedMB.size()>0)
			{
				int newmbid=fgMgr->getMBFactorIndex(mbFactor);
				if(newmbid==-1)
				{
					cout <<" FATAL ERROR!! No mb factor " << mbFactor->fId << endl;
					return -2;
				}
				mbFactor->goodMBIDs[newmbid]=0;
				mbFactor->mbScore=fgMgr->getMBScore(mbFactor);
			}
			else
			{
				mbFactor->mbScore=mbFactor->marginalEntropy;
			}
			toErase.push_back(vIter);
		}
		for(int i=0;i<toErase.size();i++)
		{
			srcFactor->mergedMB.erase(toErase[i]);
		}
		int carryOverDegree=srcFactor->mergedMB.size()-toErase.size();
		int newDegree=carryOverDegree + targetVars.size()-nonRemovals;
		if(newDegree>maxMBSize)
		{
			cout <<"Src var MB size limited " << endl;
			return -1;
		}
	//	cout <<"MB-update-Exact: "<<srcFactor->fId << " MBVars:";
	}
	else
	{
	//	cout <<"MB-update-Greedy: " << srcFactor->fId <<" MBVars:";
		for(INTINTMAP_ITER vIter=srcFactor->mergedMB.begin();vIter!=srcFactor->mergedMB.end();vIter++)
		{
	//		cout << " " << vIter->first;
		}

	}
	for(INTINTMAP_ITER vIter=targetVars.begin();vIter!=targetVars.end();vIter++)
	{
		SlimFactor* tgtFactor=newfg->getFactorAt(vIter->first);
	//	cout << " " << tgtFactor->fId;
		/*if(moveType==Move::EXACT)
		{
			if( (srcVertex.find(vIter->first)==srcVertex.end()) && (targetVertex.find(vIter->first)==targetVertex.end()))
			{
				//This means that the tgtFactor was neither targeted nor is a source yet 
				//for this iteration. Hence we can get rid of the edges associated with this node
				for(INTINTMAP_ITER mIter=tgtFactor->mergedMB.begin();mIter!=tgtFactor->mergedMB.end();mIter++)
				{
					SlimFactor* mbFactor=newfg->getFactorAt(mIter->first);
					INTINTMAP_ITER dIter=mbFactor->mergedMB.find(vIter->first);
					if(dIter==mbFactor->mergedMB.end())
					{
						cout << "FATAL ERROR!! No variable " << srcFactor->fId << " in MB of " << mbFactor->fId << endl;
						return -2;
					}
					mbFactor->mergedMB.erase(dIter);
					mbFactor->goodMBIDs.clear();
					if(mbFactor->mergedMB.size()>0)
					{
						int newmbid=fgMgr->getMBFactorIndex(mbFactor);
						if(newmbid==-1)
						{
							cout <<"FATAL ERROR!! No mbfactor for " << mbFactor->fId << endl;
							return -2;
						}
						mbFactor->goodMBIDs[newmbid]=0;
						mbFactor->mbScore=fgMgr->getMBScore(mbFactor);
					}
					else
					{
						mbFactor->mbScore=mbFactor->marginalEntropy;
					}
				}
				tgtFactor->mergedMB.clear();
			}
		}*/
		srcFactor->mergedMB[vIter->first]=0;
		if(tgtFactor->mergedMB.find(m->getSrcVertex())!=tgtFactor->mergedMB.end())
		{
			continue;
		}
		tgtFactor->mergedMB[m->getSrcVertex()]=0;
		//update the mbscore of tgtFactor
		tgtFactor->mbScore=fgMgr->getMBScore(tgtFactor);
		//Always store the lowest possible score to prevent bad moves
		if(tgtFactor->moveScore > m->getScore())
		{
			tgtFactor->moveScore=m->getScore();
		}
		tgtFactor->goodMBIDs.clear();
		if(tgtFactor->fId==30)
		{
			cout <<"Stop here " << endl;
		}
		int newmbId=fgMgr->getMBFactorIndex(tgtFactor);
		if(newmbId==-1)
		{
			cout <<" FATAL ERROR!! No mbfactor for " << tgtFactor->fId << endl;
			return -2;
		}
		tgtFactor->goodMBIDs[newmbId]=0;
	}
	//cout << endl;
	srcFactor->mbScore=m->getMBScore();
	srcFactor->moveScore=m->getScore();
	srcFactor->goodMBIDs.clear();
	srcFactor->goodMBIDs[m->getSrcMBID()]=0;
	return 0;
}

//This function is called when the MB of the src vertex has been affected by some
//previous iterations. We now need to find a new MB for this variable constrained by the
//vars currently in its MB
Move* 
FGEditor::getNewMove(FactorGraph* fg, Move* m)
{
	SlimFactor* sFactor=fg->getFactorAt(m->getSrcVertex());
	if(sFactor->mergedMB.size()>=maxMBSize)
	{
		return NULL;
	}

	Move* newMove=identifyConstrainedMove(sFactor);
	return newMove;
}


//Here we compute the score by using the exact move formula. Let X be the variable
//who's MB we are searching. Let Mx be its current MB and Mxnew be its new MB. Let Y be a 
//variable in Mxnew. Then the score of making My the MB of X instead of Mx is:
//H(X|Mx)-H(X|Mxnew) + sum_y H(Y|My)-min (H(Y|Mynew)), where Mynew must include X.
int
FGEditor::identifyMove(SlimFactor* sFactor)
{
	INTINTMAP supersets;
	if((moveType==Move::EXACT) || (sFactor->goodMBIDs.size()==0))
	{
		int level=sFactor->mergedMB.size()+1;
		fgMgr->getSupersets(sFactor->fId,level,supersets);
		//fgMgr->getSupersets(sFactor->fId,maxMBSize,supersets);
	}
	else
	{
		int mbid=sFactor->goodMBIDs.begin()->first;
		fgMgr->getSupersets(mbid,1,supersets);
	}
	if(supersets.size()==0)
	{
		return 0;
	}
	double currScore=sFactor->moveScore;
	//double currScore=0;
	double currMBScore=sFactor->mbScore;
	double currDelta=0;
	int bestMBID=-1;
	INTINTMAP currTargetMBVar;
	for(INTINTMAP_ITER ssIter=supersets.begin();ssIter!=supersets.end();ssIter++)
	{
		SlimFactor* supFactor=fgMgr->getFactorAt(ssIter->first);
		if(maxMBSize==3)
		{
			if(sFactor->fId==15 && supFactor->fId==8851)
			{
				cout <<"Stop here" << endl;
			}
		}
		double scoreDelta=-1;
		double scoreDelta1=-1;
		double scoreDelta2=-1;
		double newScore=-1;
		double newMBScore;
		INTINTMAP targetVarMB;
		if(moveType==Move::EXACT)
		{
			scoreDelta2=getScoreDelta_Greedy(sFactor,supFactor,newMBScore,newScore,targetVarMB);
			//scoreDelta2=getScoreDelta_Exact(sFactor,supFactor,newMBScore,newScore,targetVarMB);
			/*targetVarMB.clear();
			scoreDelta1=getScoreDelta_Exact(sFactor,supFactor,newMBScore,newScore,targetVarMB);
			if(scoreDelta1<scoreDelta2)
			{
				//cout << "Greedy Putative MB " << scoreDelta2 <<" is better than Exact Putative MB " << scoreDelta1 << " using supfactor " << supFactor->fId << " factor " << sFactor->fId <<endl;
			}*/
			scoreDelta=scoreDelta2;
			double lossScore=getLossScore(sFactor);
			scoreDelta=scoreDelta-lossScore;
		}
		else
		{
			scoreDelta=getScoreDelta_Greedy(sFactor,supFactor,newMBScore,newScore,targetVarMB);
		}
		if(scoreDelta>currDelta)
		//if(newScore<currScore)
		{
			bestMBID=ssIter->first;
			currDelta=scoreDelta;
			currScore=newScore;
			currMBScore=newMBScore;
			currTargetMBVar.clear();
			for(INTINTMAP_ITER mbIter=targetVarMB.begin();mbIter!=targetVarMB.end();mbIter++)
			{
				currTargetMBVar[mbIter->first]=mbIter->second;
			}
		}
	}
	if(bestMBID!=-1)
	{
		//Generate a move
		Move* amove=new Move;
		amove->setScoreImprovement(currDelta);
		amove->setScore(currScore);
		amove->setMBScore(currMBScore);
		amove->setSrcVertex(sFactor->fId);
		amove->setSrcMBID(bestMBID);
		amove->setTargetset(currTargetMBVar);
		moveSet.push_back(amove);
	}
	return 0;
}

Move* 
FGEditor::identifyConstrainedMove(SlimFactor* sFactor)
{
	
	INTINTMAP supersets;
	if(sFactor->goodMBIDs.size()==0)
	{
		fgMgr->getSupersets(sFactor->fId,maxMBSize,supersets);
	}
	else
	{
		fgMgr->getSupersets(sFactor->goodMBIDs.begin()->first,1,supersets);
	}
	if(supersets.size()==0)
	{
		return NULL;
	}
	double currScore=sFactor->moveScore;
	double currMBScore=sFactor->mbScore;
	double currDelta=0;
	int bestMBID=-1;
	INTINTMAP currTargetMBVar;
	for(INTINTMAP_ITER ssIter=supersets.begin();ssIter!=supersets.end();ssIter++)
	{
		SlimFactor* supFactor=fgMgr->getFactorAt(ssIter->first);
		double newScore=-1;
		INTINTMAP targetVarMB;
		double scoreDelta=-1;
		double newMBScore=-1;
		if(moveType==Move::EXACT)
		{
			//scoreDelta=getScoreDelta_Exact(sFactor,supFactor,newScore,targetVarMB);
			scoreDelta=getScoreDelta_Greedy(sFactor,supFactor,newMBScore,newScore,targetVarMB);
		}
		else if(moveType==Move::GREEDY)
		{
			scoreDelta=getScoreDelta_Greedy(sFactor,supFactor,newMBScore,newScore,targetVarMB);
		}
		//if(scoreDelta>currDelta)
		if(newScore<currScore)
		{
			bestMBID=ssIter->first;
			currScore=newScore;
			currMBScore=newMBScore;
			currDelta=scoreDelta;
			currTargetMBVar.clear();
			for(INTINTMAP_ITER mbIter=targetVarMB.begin();mbIter!=targetVarMB.end();mbIter++)
			{
				currTargetMBVar[mbIter->first]=mbIter->second;
			}
		}
	}
	Move* amove=NULL;
	if(bestMBID!=-1)
	{
		//Generate a move
		amove=new Move;
		amove->setScoreImprovement(currDelta);
		amove->setScore(currScore);
		amove->setSrcVertex(sFactor->fId);
		amove->setSrcMBID(bestMBID);
		amove->setTargetset(currTargetMBVar);
	}
	return amove;
}

double
FGEditor::getScoreDelta_Exact(SlimFactor* sFactor, SlimFactor* supFactor, double& mbScore,  double& newScore,INTINTMAP& mbVarFid)
{
	double condEntropy=fgMgr->getMBScore(sFactor,supFactor->fId);
	/*if(condEntropy==-1)
	{
		return -1;
	}*/
	newScore=condEntropy;
	mbScore=newScore;
	double scoreDelta=sFactor->mbScore-newScore;
	if(maxMBSize==1)
	{
		int mbVar=supFactor->vIds[0];
		if(supFactor->vIds[0]==sFactor->vIds[0])
		{
			mbVar=supFactor->vIds[1];
		}
		SlimFactor* mbVarFactor=fg->getFactorAt(mbVar);
		if(mbVarFactor->mergedMB.size()>=maxMBSize)
		{
			return -1;
		}
		mbVarFid[mbVar]=-1;
		double mbCondEntropy=fgMgr->getMBScore(mbVarFactor,supFactor->fId);
		newScore=newScore+mbCondEntropy;
		if(sFactor->candidateNeighbours.size()==0)
		{
			sFactor->candidateNeighbours[mbVar]=condEntropy;
			sFactor->candidateNeighbours_vect.push_back(mbVar);
		}
		else 
		{
			//Compare with the last element
			bool insertFlag=false;
			INTDBLMAP& entropyMap=sFactor->candidateNeighbours;
			vector<int>& mbIds=sFactor->candidateNeighbours_vect;
			int i=mbIds.size()-1;
			if(entropyMap.find(mbIds[i])==entropyMap.end())
			{
				cout << "No mb factor id " << mbIds[i] << " for factor " << sFactor->fId << endl;
				return -1;
			}
			double currMapMin=entropyMap[mbIds[i]];
			if(condEntropy==currMapMin)
			{
				cout <<"mbVar " << mbVar << " is as good as " << mbIds[i] << " for factor " << sFactor->fId << endl;
			}
			while((condEntropy<currMapMin) && (i>=0))
			{
				i--;
				if(i>=0)
				{
					if(entropyMap.find(mbIds[i])==entropyMap.end())
					{
						cout << "No mb factor id " << mbIds[i] << " for factor " << sFactor->fId << endl;
						return -1;
					}
					currMapMin=entropyMap[mbIds[i]];
				}
				insertFlag=true;
			}

			i++;
			if(mbIds.size()<mbCntForApproxGreedy)
			{
				insertFlag=true;
			}
			if((insertFlag) && (entropyMap.find(mbVar)==entropyMap.end()))
			{
				//drop off the last element only if there are the upper limit number of elements
				if(mbIds.size()==mbCntForApproxGreedy)
				{
					int delID=mbIds[mbIds.size()-1];
					mbIds.pop_back();
					INTDBLMAP_ITER idIter=entropyMap.find(delID);
					if(idIter==entropyMap.end())
					{
						cout <<"FATAL error!! No mb factor with id " << delID << " for factor " << sFactor->fId<<endl;
						return -1;
					}
					entropyMap.erase(idIter);
				}
				if(mbIds.size()>0)
				{
					//make a copy of the last element
					mbIds.push_back(mbIds[mbIds.size()-1]);
					for(int j=mbIds.size()-2;j>=i;j--)
					{
						mbIds[j+1]=mbIds[j];
					}
				}
				if(i<mbIds.size())
				{
					mbIds[i]=mbVar;
					entropyMap[mbVar]=condEntropy;
				}
				else
				{
					mbIds.push_back(mbVar);
					entropyMap[mbVar]=condEntropy;
				}
			}
		}
		return scoreDelta;
	}
	//Need to account for the loss in score if we were to use supFactor for the MB of sFactor.
	//This is because if were to make supFactor as the MB, then we would delete some edges
	//and add new. 
	//Now consider all variables in the supFactor other than the factor variable
	for(int m=0;m<supFactor->vCnt;m++)
	{
		if(supFactor->vIds[m]==sFactor->vIds[0])
		{
			continue;
		}
		SlimFactor* mbVarFactor=fg->getFactorAt(supFactor->vIds[m]);
		if(mbVarFactor->mergedMB.find(sFactor->vIds[0])!=mbVarFactor->mergedMB.end())
		{
			//continue;
		}
		//Cannot consider a supFactor which has one variable with already assigned MB in this iteration
		if(mbVarFactor->mergedMB.size()>=maxMBSize)
		{
			return -1;
		}
		//Figure out the best Markov blanket for each MB variable under the constraint 
		//that sFactor's variable is in the MB of each MB variable 
		int fVars[2];
		if(sFactor->vIds[0]<supFactor->vIds[m])
		{
			fVars[0]=sFactor->vIds[0];
			fVars[1]=supFactor->vIds[m];
		}
		else
		{
			fVars[0]=supFactor->vIds[m];
			fVars[1]=sFactor->vIds[0];
		}
		int fid=fgMgr->getFactorIndex(fVars,2);
		double mbVarDelta=-1;
		char key[256];
		sprintf(key,"%d-%d",fid,mbVarFactor->fId);
		string keyStr(key);

		if(gainExactMB.find(keyStr)!=gainExactMB.end())
		{
			mbVarDelta=gainExactMB[keyStr];
		}
		else
		{
			INTINTMAP mbSSet;
			int level=mbVarFactor->mergedMB.size();
			if(level==0)
			{
				mbSSet[fid]=0;
			}
			else
			{
				fgMgr->getSupersets(fid,level,mbSSet);
			}
			int newMBVarMB=-1;
			double newmbVarDelta=getScoreChange_Neighbour(mbVarFactor,mbSSet,newMBVarMB);
			if((mbVarDelta!=-1) && (mbVarDelta!=newmbVarDelta))
			{
				cout << "Bad cache value for " << fid << " old val: " << mbVarDelta << " new val: " << newmbVarDelta << endl;
			}
			mbVarDelta=newmbVarDelta;
			if(mbVarDelta!=-1)
			{
				gainExactMB[keyStr]=newmbVarDelta;
			}
		}
		if(mbVarDelta==-1)
		{
			return -1;
		}
		double mbVarDelta_greedy=getScoreChange_NeighbourGreedy(mbVarFactor,sFactor);
		//Get the varDelta if neighbour was doing greedy neighbor search
		if(mbVarDelta < mbVarDelta_greedy)
		{
			cout <<"Delta with argmin " << mbVarDelta << " is worse than one variable extension " << mbVarDelta_greedy << " for MBVar "<< mbVarFactor->fId<<" and factor " << sFactor->fId << endl;
		}
		//mbVarFid[supFactor->vIds[m]]=newMBVarMB;
		//mbVarFid[supFactor->vIds[m]]=fid;
		double loss=getLossScore(mbVarFactor);
		mbVarFid[supFactor->vIds[m]]=-1;
		scoreDelta=scoreDelta+mbVarDelta-loss;
	}
	return scoreDelta;
}


double
FGEditor::getScoreDelta_Greedy(SlimFactor* sFactor, SlimFactor* supFactor,double& newCondEntr, double& newScore, INTINTMAP& mbVarFid)
{
	int* mbmbVars=new int[maxMBSize+1];
	double condEntropy=fgMgr->getMBScore(sFactor,supFactor->fId);
	/*if(condEntropy==-1)
	{
		return -1;
	}*/
	newCondEntr=condEntropy;
	newScore=newCondEntr;
	double scoreDelta=sFactor->mbScore-newCondEntr;
	if(maxMBSize==1)
	{
		int mbVar=supFactor->vIds[0];
		if(supFactor->vIds[0]==sFactor->vIds[0])
		{
			mbVar=supFactor->vIds[1];
		}
		SlimFactor* mbVarFactor=fg->getFactorAt(mbVar);
		if(mbVarFactor->mergedMB.size()>=maxMBSize)
		{
			return -1;
		}
		mbVarFid[mbVar]=-1;
		double mbCondEntropy=fgMgr->getMBScore(mbVarFactor,supFactor->fId);
		newScore=newScore+mbCondEntropy;
		if(sFactor->candidateNeighbours.size()==0)
		{
			sFactor->candidateNeighbours[mbVar]=condEntropy;
			sFactor->candidateNeighbours_vect.push_back(mbVar);
		}
		else 
		{
			//Compare with the last element
			bool insertFlag=false;
			INTDBLMAP& entropyMap=sFactor->candidateNeighbours;
			vector<int>& mbIds=sFactor->candidateNeighbours_vect;
			int i=mbIds.size()-1;
			if(entropyMap.find(mbIds[i])==entropyMap.end())
			{
				cout << "No mb factor id " << mbIds[i] << " for factor " << sFactor->fId << endl;
				return -1;
			}
			double currMapMin=entropyMap[mbIds[i]];
			if(condEntropy==currMapMin)
			{
				cout <<"mbVar " << mbVar << " is as good as " << mbIds[i] << " for factor " << sFactor->fId << endl;
			}
			while((condEntropy<currMapMin) && (i>=0))
			{
				i--;
				if(i>=0)
				{
					if(entropyMap.find(mbIds[i])==entropyMap.end())
					{
						cout << "No mb factor id " << mbIds[i] << " for factor " << sFactor->fId << endl;
						return -1;
					}
					currMapMin=entropyMap[mbIds[i]];
				}
				insertFlag=true;
			}

			i++;
			if(mbIds.size()<mbCntForApproxGreedy)
			{
				insertFlag=true;
			}
			if((insertFlag) && (entropyMap.find(mbVar)==entropyMap.end()))
			{
				//drop off the last element only if there are the upper limit number of elements
				if(mbIds.size()==mbCntForApproxGreedy)
				{
					int delID=mbIds[mbIds.size()-1];
					mbIds.pop_back();
					INTDBLMAP_ITER idIter=entropyMap.find(delID);
					if(idIter==entropyMap.end())
					{
						cout <<"FATAL error!! No mb factor with id " << delID << " for factor " << sFactor->fId<<endl;
						return -1;
					}
					if(sFactor->fId==1 && delID==106)
					{
						cout <<"Stop here " << endl;
					}
					entropyMap.erase(idIter);
				}
				if(mbIds.size()>0)
				{
					//make a copy of the last element
					mbIds.push_back(mbIds[mbIds.size()-1]);
					for(int j=mbIds.size()-2;j>=i;j--)
					{
						mbIds[j+1]=mbIds[j];
					}
				}
				if(i<mbIds.size())
				{
					mbIds[i]=mbVar;
					entropyMap[mbVar]=condEntropy;
				}
				else
				{
					mbIds.push_back(mbVar);
					entropyMap[mbVar]=condEntropy;
				}
			}
		}
		return scoreDelta;
	}
	//Now consider all variables in the supFactor other than
	for(int m=0;m<supFactor->vCnt;m++)
	{
		if(supFactor->vIds[m]==sFactor->vIds[0])
		{
			continue;
		}
		SlimFactor* mbVarFactor=fg->getFactorAt(supFactor->vIds[m]);
		if(mbVarFactor->mergedMB.find(sFactor->vIds[0])!=mbVarFactor->mergedMB.end())
		{
			if(moveType==Move::EXACT)
			{
				mbVarFid[supFactor->vIds[m]]=-1;
			}
			continue;
		}
		//Cannot consider a supFactor which has one variable with already assigned MB in this iteration
		if(mbVarFactor->mergedMB.size()>=maxMBSize)
		{
			return -1;
		}
		//Figure out the best Markov blanket for each MB variable under the constraint 
		//that sFactor's variable is in the MB of each MB variable 
		int xVar=sFactor->vIds[0];
		int uid=0;
		if(mbVarFactor->goodMBIDs.size()==0)
		{
			if(xVar<mbVarFactor->vIds[0])
			{
				mbmbVars[0]=xVar;
				mbmbVars[1]=mbVarFactor->vIds[0];
			}
			else
			{
				mbmbVars[0]=mbVarFactor->vIds[0];
				mbmbVars[1]=xVar;
			}
			uid=2;
		}
		else
		{
			int mbid=mbVarFactor->goodMBIDs.begin()->first;
			SlimFactor* mbVarcurrMBFactor=fgMgr->getFactorAt(mbid);
			//Check the effect of adding sFactor's variable to this factor
			int vid=0;
			uid=0;
			while(vid<mbVarcurrMBFactor->vCnt && mbVarcurrMBFactor->vIds[vid]<xVar)
			{
				mbmbVars[uid]=mbVarcurrMBFactor->vIds[vid];
				vid++;
				uid++;
			}
			mbmbVars[uid]=xVar;
			uid++;
			while(vid<mbVarcurrMBFactor->vCnt)
			{
				mbmbVars[uid]=mbVarcurrMBFactor->vIds[vid];
				uid++;
				vid++;
			}
		}
		int fid=fgMgr->getFactorIndex(mbmbVars,uid);
		if(fid==-1)
		{
			/*cout <<"Bad variables in mbmbVars while computing delta for factor var "<< sFactor->fId 
				<< " mbvar " << supFactor->vIds[m]<< endl;
			for(int i=0;i<uid;i++)
			{
				cout << " " << mbmbVars[i];
			}
			cout << endl;*/
			return -1;
		}
		mbVarFid[supFactor->vIds[m]]=-1;
		double mbCondEntropy=fgMgr->getMBScore(mbVarFactor,fid);
		/*if(mbCondEntropy==-1)
		{
			return -1;
		}*/
		double mbVarDelta=mbVarFactor->mbScore-mbCondEntropy;
		scoreDelta=scoreDelta+mbVarDelta;
		newScore=newScore+mbCondEntropy;
	}
	delete[] mbmbVars;
	return scoreDelta;
}

double
FGEditor::getScoreChange_Neighbour(SlimFactor* sFactor,INTINTMAP& ssets, int& mbId)
{
	if(ssets.size()==0)
	{
		return 0;
	}

	mbId=-1;
	double currScore=sFactor->mbScore;
	for(INTINTMAP_ITER aIter=ssets.begin();aIter!=ssets.end();aIter++)
	{
		double proposedScore=fgMgr->getMBScore(sFactor,aIter->first);
		/*if(proposedScore<0)
		{
			return -1;
		}*/
		if((proposedScore<=currScore) || (currScore==-1))
		{
			currScore=proposedScore;
			mbId=aIter->first;
		}
	}
	double gain=sFactor->mbScore-currScore;
	return gain;
}



double
FGEditor::getScoreChange_NeighbourGreedy(SlimFactor* mbVarFactor,SlimFactor* sFactor)
{
	int xVar=sFactor->vIds[0];
	int uid=0;
	int* mbmbVars=new int[maxMBSize+1];
	if(mbVarFactor->goodMBIDs.size()==0)
	{
		if(xVar<mbVarFactor->vIds[0])
		{
			mbmbVars[0]=xVar;
			mbmbVars[1]=mbVarFactor->vIds[0];
		}
		else
		{
			mbmbVars[0]=mbVarFactor->vIds[0];
			mbmbVars[1]=xVar;
		}
		uid=2;
	}
	else
	{
		int mbid=mbVarFactor->goodMBIDs.begin()->first;
		SlimFactor* mbVarcurrMBFactor=fgMgr->getFactorAt(mbid);
		//Check the effect of adding sFactor's variable to this factor
		int vid=0;
		uid=0;
		while(vid<mbVarcurrMBFactor->vCnt && mbVarcurrMBFactor->vIds[vid]<xVar)
		{
			mbmbVars[uid]=mbVarcurrMBFactor->vIds[vid];
			vid++;
			uid++;
		}
		mbmbVars[uid]=xVar;
		uid++;
		while(vid<mbVarcurrMBFactor->vCnt)
		{
			mbmbVars[uid]=mbVarcurrMBFactor->vIds[vid];
			uid++;
			vid++;
		}
	}
	int fid=fgMgr->getFactorIndex(mbmbVars,uid);
	if(fid==-1)
	{
		return -1;
	}
	double mbCondEntropy=fgMgr->getMBScore(mbVarFactor,fid);
	/*if(mbCondEntropy==-1)
	{
		return -1;
	}*/
	double gain=mbVarFactor->mbScore-mbCondEntropy;
	delete[] mbmbVars;
	return gain;
}

//The goal of this function is to take the current Markov blanket of sFactor and 
//compute the net loss in the score (gain in conditional entropy) if we were to delete
//sFactor's variable from the Markov blanket variables of sFactor
double
FGEditor::getLossScore(SlimFactor* sFactor)
{
	double loss=0;
	int* newMBVars=new int[maxMBSize];
	for(INTINTMAP_ITER uIter=sFactor->mergedMB.begin();uIter!=sFactor->mergedMB.end();uIter++)
	{
		SlimFactor* mbFactor=fg->getFactorAt(uIter->first);
		if(mbFactor->mergedMB.size()==1)
		{
			loss=loss+mbFactor->marginalEntropy-mbFactor->mbScore;
		}
		else if(mbFactor->mergedMB.size()>0)
		{
			double currScore=mbFactor->mbScore;
			int vId=0;
			bool varAssigned=false;
			for(INTINTMAP_ITER vIter=mbFactor->mergedMB.begin();vIter!=mbFactor->mergedMB.end();vIter++)
			{
				if(vIter->first==sFactor->fId)
				{
					continue;
				}
				if((mbFactor->fId < vIter->first) && (!varAssigned))
				{
					newMBVars[vId]=mbFactor->fId;
					vId++;
					varAssigned=true;
				}
				newMBVars[vId]=vIter->first;
				vId++;
			}
			if(!varAssigned)
			{
				newMBVars[vId]=mbFactor->fId;
				vId++;
			}
			int fid=fgMgr->getFactorIndex(newMBVars,vId);
			if(fid!=-1)
			{
				continue;
			}
			double newScore=fgMgr->getMBScore(mbFactor,fid);
		//	if(newScore!=-1)
		//	{
				loss=loss+(newScore-currScore);
		/*	}
			else
			{
				cout <<"Negative score for factor " << fid << endl;
				return -1;
			}*/
		}
	}
	delete[] newMBVars;
	return loss;
}
