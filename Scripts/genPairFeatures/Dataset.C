#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "Dataset.H"

Dataset::Dataset()
{
	//by default we assume are reading a chip-seq data
	type='C';
}

Dataset::~Dataset()
{
}

int 
Dataset::readDataSet(const char* aFName)
{
	type='C';
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t_");
		int tokCnt=0;
		string chrom;
		int begin;
		int end;
		double val=-1;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chrom.append(tok);
			}
			else if(tokCnt==1)
			{
				begin=atoi(tok);
			}
			else if(tokCnt==2)
			{
				end=atoi(tok);
				//temp fix to make sure things are the same as in the cnts file
				if(end%5000!=0)
				{
					end=end+1;
				}	
			}
			else if(tokCnt==3)
			{
				val=atof(tok);
				//temp fix to normalize things by 1000 to reduce variance
				val=val/1000;
			}
			tok=strtok(NULL,"\t_");
			tokCnt++;
		}	
		Dataset::Peak* r=new Dataset::Peak;	
		r->start=begin;
		r->end=end;
		r->lqval=val;
		map<int,Peak*>* rSet=NULL;
		if(peakSet.find(chrom)==peakSet.end())
		{
			rSet=new map<int,Peak*>;
			peakSet[chrom]=rSet;
		}
		else
		{
			rSet=peakSet[chrom];
		}
		if(rSet->find(begin)==rSet->end())
		{
			(*rSet)[begin]=r;
		}
		else
		{
			//If for some reason the other dataset has a hit at he same peak location, then use the greater value
			Peak* oldpeak=(*rSet)[begin];
			if(r->lqval>oldpeak->lqval)
			{
				(*rSet)[begin]=r;
			}
		}
	}
	inFile.close();
	cout <<"Done reading "<< peakSet.size() << " different chromosomes"<<endl;
	return 0;
}


int 
Dataset::readDataSetAsPIQMotif(const char* aFName, char ty)
{
	type=ty; // 'M' or 'D'
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
		string chrom;
		string strand;
		int begin;
		int end;
		double val=-1;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chrom.append(tok);
			}
			/*else if(tokCnt==2)
			{
				//stretch on both sides
				begin=atoi(tok)-10;
				end=atoi(tok)+10;
				if(begin<0)
				{
					begin=0;
				}
				//This is the start of the motif, so we 
			}*/
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
			else if(tokCnt==4)
                        {
                                strand.append(tok);
                        }
			tok=strtok(NULL,"\t");
			tokCnt++;
		//	cout << "Read PIQ:--"<< chrom <<" - " <<begin <<" - " << end<< " , "<<val<<endl;
		}
		//cout << "Read PIQ:--"<< chrom <<" - " <<begin <<" - " << end<< " , "<<val<<endl;	
		Dataset::Peak* r=new Dataset::Peak;	
		r->start=begin;
		r->end=end;
		r->lqval=val;
		map<int,Peak*>* rSet=NULL;
		if(peakSet.find(chrom)==peakSet.end())
		{
			rSet=new map<int,Peak*>;
			peakSet[chrom]=rSet;
		}
		else
		{
			rSet=peakSet[chrom];
		}
		if(rSet->find(begin)==rSet->end())
		{
			(*rSet)[begin]=r;
		}
		else
		{
			//If for some reason the other dataset has a hit at he same peak location, then use the greater value
			Peak* oldpeak=(*rSet)[begin];
			if(r->lqval>oldpeak->lqval)
			{
				(*rSet)[begin]=r;
			}
		}
	}
	inFile.close();
	cout <<"Done reading "<< peakSet.size() << " different chromosomes"<<endl;
	return 0;

}


//Return the value with the greatest overlap
//In theory we need not do an overlap search since our total space of regions are the same between the pairs and the signals.
//So we could speed this up greatly, provided the regions DO NOT overlap
double
Dataset::getFeatureVal(string& chrom, int begin, int end)
{
	if(peakSet.find(chrom)==peakSet.end())
	{
		return -1;
	}	
	map<int,Peak*>* peaks=peakSet[chrom];
	double maxVal=-1;
	double peakcnt=0;

	//added by Shilu, to count CTCF motifs:
	//double ctcfcnt=0;
	if(type=='C')
	{
		if(peaks->find(begin)!=peaks->end())
		{
			Peak* peak=(*peaks)[begin];
			maxVal=peak->lqval;
		}
	}
	else if(type=='M')
	{
		maxVal=0;
		for(map<int,Peak*>::iterator pIter=peaks->begin();pIter!=peaks->end();pIter++)
		{
			Peak* peak=pIter->second;
			int begin1=begin;
			int end1=end;
			int begin2=peak->start;
			int end2=peak->end;
			if(begin1>begin2)
			{
				begin2=begin;
				end2=end;
				begin1=peak->start;
				end1=peak->end;
			}
			if(end1<begin2)
			{
				continue;
			}
			/*if(maxVal<peak->lqval)
			{
				maxVal=peak->lqval;
				peakcnt++;
			}*/
			maxVal+=peak->lqval;
			peakcnt++;
		}
	}
	//------------------------------------------------
	// added by Shilu, to count CTCF motitfs with + or -
	else if(type=='D')
        {
                maxVal=0;
                for(map<int,Peak*>::iterator pIter=peaks->begin();pIter!=peaks->end();pIter++)
                {
                        Peak* peak=pIter->second;
                        int begin1=begin;
                        int end1=end;
                        int begin2=peak->start;
                        int end2=peak->end;
                        if(begin1>begin2)
                        {
                                begin2=begin;
                                end2=end;
                                begin1=peak->start;
                                end1=peak->end;
                        }
                        if(end1<begin2)
                        {
                                continue;
                        }
			//ctcfcnt++;
			//maxVal+=peak->lqval;
                        peakcnt++;
                }
		maxVal=peakcnt;
		cout << "CTCF motif counts: " << peakcnt << endl;
        }	
	else
	{
		cout <<"Type mismatch! Throwing up hands in the air and quitting" << endl;
		exit(0);
	}
	 
	return maxVal;
	
	//return peakcnt;
}
