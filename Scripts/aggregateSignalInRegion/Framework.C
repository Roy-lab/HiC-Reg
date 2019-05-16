#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}


int 
Framework::readTSS(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[4096];
	while(inFile.good())
	{
		inFile.getline(buffer,4095);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string chromosome;
		string genename;
		char strand;
		int position_b;
		int position_e;
		int position;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chromosome.append(tok);
			}
			else if(tokCnt==3)
			{
				position_b=atoi(tok);
				//flipping
				//position_e=atoi(tok);
			}
			else if(tokCnt==4)
			{
				position_e=atoi(tok);
				//flipping
				//position_b=atoi(tok);
			}
			else if(tokCnt==6)
			{
				strand=tok[0];
			}
			else if(tokCnt==8)
			{
				char* pos=strchr(tok,'.');
				if(pos!=NULL)
				{	
				//	*pos='\0';
				}
				pos=strchr(tok,' ');
				if(pos!=NULL)
				{
					tok=pos+1;
				}
				pos=strchr(tok,' ');
				if(pos!=NULL)
				{	
					*pos='\0';
				}
				genename.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(genename.length()<=0)
		{
			continue;
		}
		map<string,Framework::Gene*>* geneSet=NULL;
		map<int,Framework::Gene*>* tsses=NULL;
		if(tssSet.find(chromosome)==tssSet.end())
		{
			geneSet=new map<string,Framework::Gene*>;
			tssSet[chromosome]=geneSet;
			tsses=new map<int,Framework::Gene*>;
			tsscoordSet[chromosome]=tsses;
		}
		else
		{
			geneSet=tssSet[chromosome];
			tsses=tsscoordSet[chromosome];
		}
		if(strand=='+')
		{
			position=position_b;
		}
		else	
		{
			position=position_e;
		}
		Framework::Gene* gene=NULL;
		if(geneSet->find(genename)==geneSet->end())
		{
			gene=new Framework::Gene;
			(*geneSet)[genename]=gene;
			gene->genename.append(genename);
			gene->strand=strand;
			gene->tsses[position]=0;
			if(strand=='+')
			{
				gene->begin=position_b;
				gene->end=position_e;
			}
			else
			{
				gene->begin=position_e;
				gene->end=position_b;
			}
			//arbitrarily use one tss per gene in the set of tsses. Ideally this should be the left most.
			(*tsses)[position]=gene;
		}
		else
		{
			gene=(*geneSet)[genename];
			gene->tsses[position]=0;
			//(*tsses)[position]=gene;
		}
	}
	inFile.close();
	return 0;
}

int
Framework::readChromosomeSizes(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string chromosome;
		int size;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chromosome.append(tok);
			}
			else if(tokCnt==1)	
			{
				size=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		chromosomeSizes[chromosome]=size;
	}
	inFile.close();
	return 0;
}
	

int 
Framework::readInput(const char* aFName)
{
	readData(aFName,bg,inputAvg);
	return 0;
}

int 
Framework::readSignals(const char* aFName)
{
	readData(aFName,signal,signalAvg);
	return 0;
}

int 
Framework::normalize()
{
	for(map<string,vector<double>*>::iterator cIter=signal.begin();cIter!=signal.end();cIter++)
	{
		vector<double>* values_fg=signal[cIter->first];
		//map<int,double>* lr=new map<int,double>;
		//logRatios[cIter->first]=lr;
		for(int i=0;i<values_fg->size();i++)
		{
			double normVal=(*values_fg)[i]/signalAvg;
			(*values_fg)[i]=normVal;
		}
	}
	return 0;
}

int 
Framework::readData(const char* fName,map<string,vector<double>*>& data, double& globalAvg)
{
	ifstream inFile(fName);
	char buffer[1024];
	int cnt=0;
	double totalBins=0;
	globalAvg=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer," \t");
		int tokCnt=0;
		string chr;
		int begin;
		int end;
		double val;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chr.append(tok);
			}
			else if(tokCnt==1)
			{
				begin=atoi(tok);
			}
			else if(tokCnt==2)
			{
				end=atoi(tok);
			}
			else if(tokCnt==3)
			{
				val=atof(tok);
			}
			tok=strtok(NULL," \t");
			tokCnt++;
		}
		//avoid 0
		val=val+1e-3;
		vector<double>* values=NULL;
		if(data.find(chr)==data.end())
		{
			if(chromosomeSizes.find(chr)==chromosomeSizes.end())	
			{
				cout <<"Chromosome sizes are not read in or chromosome " << chr << " does not exist" <<endl;
			}
			int chrSize=chromosomeSizes[chr];
			values=new vector<double>(chrSize);
			for(int i=0;i<chrSize;i++)
			{
				(*values)[i]=0;
			}
			data[chr]=values;
			cout <<"Created chromosome " << chr<< endl;	
		}
		else
		{
			values=data[chr];
		}
		//Now store for every location
		for(int i=begin;i<end;i++)
		{	
			(*values)[i]=val;
		}
		cnt++;
		if(cnt%500000==0)
		{
			cout <<".";
			cout.flush();
		}
		globalAvg=globalAvg+val;
	}
	
	inFile.close();
	globalAvg=globalAvg/cnt;
	cout << endl <<"Read "<< cnt << " values " << " Avg: " << globalAvg  << endl;
	return 0;
}

int 
Framework::getSignalInRegion(const char* aFName)
{
	cout << "in getSignalInRegion()" << endl;
	char sigFName[1024];
//	sprintf(sigFName,"%s_profile.txt",aFName);
//	ofstream sigFile(sigFName);
	char avgFName[1024];
	sprintf(avgFName,"%s_avg.txt",aFName);
	cout << aFName << endl;	
	ofstream avgFile(aFName);
	//Go over each chromosome and output the signal and averaged values
	cout << "tssSet.size(): " << tssSet.size() << endl;
	for(map<string,map<string,Framework::Gene*>*>::iterator cIter=tssSet.begin();cIter!=tssSet.end();cIter++)
	{
		map<string,Framework::Gene*>* geneSet=cIter->second;
		if(signal.find(cIter->first)==signal.end())
		{
			cout << "Skipping chromosome " << cIter->first << endl;
			continue;
		}
		vector<double>* svals=signal[cIter->first];
		for(map<string,Framework::Gene*>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			//We will focus on the start and end of this gene/regon
			Framework::Gene* g=gIter->second;

			map<int,map<int,double>*> allValues;
			int binstart=g->begin;
			int binend=g->end;
			double mean=0;
			map<int,double>* vset=NULL;
			for(int b=binstart;b<=binend;b++)
			{
				int coord_key=b;
				/*if(svals->find(coord_key)==svals->end())
				{
					continue;
				}*/
				mean=mean+(*svals)[coord_key];
				if(vset==NULL)
				{
					vset=new map<int,double>;
				}
				(*vset)[b]=(*svals)[coord_key];
			}
			if(vset==NULL)
			{
				continue;
			}
			//mean=mean/vset->size();
			mean=mean/(binend-binstart);
			if(mean<=0)
			{
				cout <<"This region" << g->genename << " really has no signal. But still printing " << endl;
				continue;
			}
			//Now we have the best TSS
			/*double regVal=0;
			int regCnts=0;
			sigFile <<g->genename;
			for(int tval=binstart;tval<binend;tval++)
			{
				if(vset->find(tval)==vset->end())
				{
					sigFile <<"\t0";
					continue;	
				}
				regVal=regVal+(*vset)[tval];
				regCnts++;
				sigFile <<"\t" << (*vset)[tval];
			}
			sigFile << endl;*/
			//Now free the damn memory
			vset->clear();
			delete vset;
			avgFile << g->genename <<"\t" << mean << endl;
		}
	}
	//sigFile.close();
	avgFile.close();
	return 0;
}


int
main(int argc, const char** argv)
{
	if(argc!=5)
	{
		//cout <<"Usage: aggregateSignals tss signal input outfile" << endl;
		cout <<"Usage: aggregateSignals tss chromosomesizes signal outfile" << endl;
		return 0;
	}
	Framework fw;
	fw.readTSS(argv[1]);
	fw.readChromosomeSizes(argv[2]);
	fw.readSignals(argv[3]);
	//fw.readInput(argv[3]);
	//fw.normalize();
	//fw.mapSignalToGene(argv[3],atoi(argv[4]),atoi(argv[5]));
	fw.getSignalInRegion(argv[4]);
	return 0;
}
