#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include "gsl/gsl_randist.h"
#include "Dataset.H"
#include "Distance.H"
#include "Framework.H"
#include <algorithm>
#include <vector>

Framework::Framework()
{
  filterZeros=true;
  totalFilteredPairs=0;
  maxDist=1000000;
  preRandomize=false;
  featype="null";
  regionDistNeighborSize=200;
  initial=0;
  BinSize=5000;
}

Framework::~Framework()
{
}

int
  Framework::setMaxDist(int mDist)
  {
    maxDist=mDist;
    return 0;
  }

int
  Framework::setFeatype(string feat)
  {
    
    featype=feat;
    return 0;
  }

//Read the enhancer-promoter interactions
int 
  Framework::readPairs(const char* aFName)
  {
    ifstream inFile(aFName);
    char buffer[1024];
    int totalPairs=0;
    while(inFile.good())
    {
      inFile.getline(buffer,1023);
      if(strlen(buffer)<=0)
      {
        continue;
      }
      if(strchr(buffer,'-')==NULL && strchr(buffer,'\t')==NULL)
      {
        continue;
      }
      char* tok=strchr(buffer,'\t');
      if(tok!=NULL)
      {
        *tok='\0';
      }
      char* e=buffer;
      char* p=tok+1;
      char* tok2=strchr(p,'\t');
      double paircnt=0;
      if(tok2!=NULL)
      {
        *tok2='\0';
        char* cntval=tok2+1;
        char* end=strchr(cntval,'\t');
        if(end!=NULL)
        {
          *end='\0';
        }
        paircnt=atof(cntval);
      }
      char* primerpref=strrchr(e,'|');
      if(primerpref!=NULL)
      {
        e=primerpref+1;
      }
      
      int i=0;
      while(e[i]!='\0')
      {
        if(e[i]=='-' || e[i]==':' )
        {
          e[i]='_';
        }
        i++;
      }
      primerpref=strrchr(p,'|');
      if(primerpref!=NULL){
        p=primerpref+1;
      }
      i=0;
      while(p[i]!='\0')
      {
        if(p[i]=='-' || p[i]==':')
        {
          p[i]='_';
        }
        i++;
      }
      string eKey(e);
      if(regionSet.find(eKey)==regionSet.end())
      {
        Framework::Region* eregion=new Framework::Region;
        regionSet[eKey]=eregion;	
        string chrom;	
        int start;
        int end;
        tok=strtok(e,"_");
        int tokCnt=0;
        while(tok!=NULL)
        {
          if(tokCnt==0)
          {
            eregion->chromosome.append(tok);
          }
          else if(tokCnt==1)
          {	
            eregion->begin=atoi(tok);
          }
          else if(tokCnt==2)
          {
            eregion->end=atoi(tok);
          }
          tok=strtok(NULL,"_");
          tokCnt++;
          
        } 
        if(initial==0)
        {
          BinSize=eregion->end - eregion->begin;
          regionDistNeighborSize=maxDist/BinSize;
          cout << "regionDistNeighborSize is : " << regionDistNeighborSize << " BinSize is: " << BinSize << endl;
          initial=1;
        }
      }
      string pKey(p);
      if(regionSet.find(pKey)==regionSet.end())
      {
        Framework::Region* pregion=new Framework::Region;
        regionSet[pKey]=pregion;	
        tok=strtok(p,"_");
        int tokCnt=0;
        while(tok!=NULL)
        {
          if(tokCnt==0)
          {
            pregion->chromosome.append(tok);
          }
          else if(tokCnt==1)
          {	
            pregion->begin=atoi(tok);
          }
          else if(tokCnt==2)
          {
            pregion->end=atoi(tok);
          }
          tok=strtok(NULL,"_");
          tokCnt++;
        }
        regionSet[pKey]=pregion;
      }
      int dist=getDistance(eKey,pKey);
      if(dist>maxDist)
      {
        continue;
      }
      totalPairs++;
      map<string,double>* promSet=NULL;
      if(pairSet.find(eKey)==pairSet.end())
      {
        promSet=new map<string,double>;
        pairSet[eKey]=promSet;
      }	
      else
      {
        promSet=pairSet[eKey];
      }
      
      (*promSet)[pKey]=log(paircnt+1);
      //(*promSet)[pKey]=paircnt;
    }
    inFile.close();
    cout <<"Found " << totalPairs << " pairs at dist " << maxDist << endl;
    return 0;
  }

//Split pairs

int
  Framework::splitPairs(int nCV,bool regionWise,const char* outputDir)
  {
    if(regionWise)
    {
      if(preRandomize)
      {
        splitRegionsGenPairs(nCV,outputDir);
      }
      else
      {
        splitRegionsGenPairs_NoPrerandomize(nCV,outputDir);
      }
    }
    else
    {
      if(preRandomize)
      {
        genPairs_Prerandomize(nCV,outputDir);
      }
      else
      {
        genPairs(nCV,outputDir);
      }
    }
    return 0;
  }



int
  Framework::splitRegionsGenPairs(int foldCnt,const char* outputDir)
  {
    int totalRegions=regionFeatures.size();
    int testSetSize=(int)(totalRegions/foldCnt);
    int testBegin=0;
    int testEnd=0;
    int rID=0;
    map<string,int> regionNameIDMap;
    int* randInt=new int[regionFeatures.size()];
    for(map<string,map<string,double>*>::iterator rIter=regionFeatures.begin();rIter!=regionFeatures.end();rIter++)
    {
      regionNameIDMap[rIter->first]=rID;
      randInt[rID]=rID;
      rID++;
    }
    cout << "Total regions: " << totalRegions << endl;	
    if(preRandomize)
    {
      //This will do what we are doing in the above for loop, but also randomize
      populateRandIntegers(randInt,regionFeatures.size());
    }
    
    for(int i=0;i<foldCnt;i++)
    {
      char testFName[1024];
      char trainFName[1024];
      char test1FName[1024];
      char test2FName[1024];
      //	char test3FName[1024];
      
      sprintf(testFName,"%s/test%d.txt",outputDir,i);
      sprintf(trainFName,"%s/train%d.txt",outputDir,i);
      sprintf(test1FName,"%s/test%d_difP.txt",outputDir,i);
      sprintf(test2FName,"%s/test%d_difE.txt",outputDir,i);
      //sprintf(test3FName,"%s/test%d_3.txt",outputDir,i);
      
      ofstream teFile(testFName);
      ofstream trFile(trainFName);
      ofstream test1File(test1FName);
      ofstream test2File(test2FName);
      //ofstream test3File(test3FName);
      
      showHeader(teFile);
      showHeader(trFile);
      showHeader(test1File);
      showHeader(test2File);
      //showHeader(test3File);
      
      testEnd=(i+1)*testSetSize;
      if(testEnd>totalRegions)
      {
        testEnd=totalRegions;
      }
      map<int,int> trainID;
      map<int,int> testID;
      for(int tid=0;tid<regionNameIDMap.size();tid++)
      {
        int rID=randInt[tid];
        //	cout << "tid: " <<tid<< " rID: " << rID <<"test start end: "  <<testBegin << "\t" << testEnd<<endl;
        if(tid<testEnd && tid>=testBegin)
        {
          testID[rID]=0;
        }
        else
        {
          trainID[rID]=0;
        }
        
      }
      //We already have the pairs. We just need to figure out if they should go to the train set or training set..?
      for(map<string,map<string,map<string,double>*>*>::iterator pIter=pairSet_Filtered.begin();pIter!=pairSet_Filtered.end();pIter++)
      {
        map<string,double>* efeatures=regionFeatures[pIter->first];
        map<string,map<string,double>*>* neighbors=pIter->second;
        map<string,double>* enfeatures=regionNeighborFeatures[pIter->first];
        int eID=regionNameIDMap[pIter->first];
        for(map<string,map<string,double>*>::iterator nIter=neighbors->begin();nIter!=neighbors->end();nIter++)
        {
          map<string,double>* pfeatures=regionFeatures[nIter->first];
          int pID=regionNameIDMap[nIter->first];
          //We just use the randomized IDs to decide whether we will put a pair in the training set or the test set 
          if( (testID.find(eID)!=testID.end()) && (testID.find(pID)!=testID.end()))
          {
            teFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(teFile,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(eID)!=trainID.end()) && (trainID.find(pID)!=trainID.end()))
          {
            trFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(trFile,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(eID)!=trainID.end()) )
          {
            test1File <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(test1File,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(pID)!=trainID.end()) )
          {
            test2File <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(test2File,efeatures,pfeatures,enfeatures,nIter->second);
            
          }
        }	 
      }
      teFile.close();
      trFile.close();
      test1File.close();
      test2File.close();
      testBegin=testEnd;
      trainID.clear();
      testID.clear();
    }
    return 0;
  }





int
  Framework::splitRegionsGenPairs_NoPrerandomize(int foldCnt,const char* outputDir)
  {
    int totalRegions=regionFeatures.size();
    int testSetSize=(int)(totalRegions/foldCnt);
    int testBegin=0;
    int testEnd=0;
    int rID=0;
    map<string,int> regionNameIDMap;
    
    // added by Shilu 06/26 to sort regions by coordinates
    
    int* randInt=new int[regionFeatures.size()];
    vector <int> randIntnum;
    int binsize=0;
    for(map<string,map<string,double>*>::iterator rIter=regionFeatures.begin();rIter!=regionFeatures.end();rIter++)
    {
      char* region;
      strcpy(region, rIter->first.c_str());
      char* tok=strtok(region,"_");
      int tokCnt=0;
      int start=0;
      int end=0;
      while(tok!=NULL)
      {
        if(binsize>0)
        {
          // get the start coordinate 
          if(tokCnt==1)
          {
            rID=int (atoi(tok)/binsize);
          }
        }
        else{
          if(tokCnt==1)
          {
            start=int (atoi(tok));
          }
          else if(tokCnt==2)
          {
            end=int (atoi(tok));
          }
          binsize=end-start;
        }
        tok=strtok(NULL,"_");
        tokCnt++;
      }
      regionNameIDMap[rIter->first]=rID;        		
      //randInt[rID]=rID;
      randIntnum.push_back(rID);
      //rID++;
      //cout << "rIter->first: "<<rIter->first << " the coord is: " << rID << endl;
    }
    
    std::sort(randIntnum.begin(), randIntnum.end());    
    int id=0;
    for (std::vector<int>::iterator it=randIntnum.begin(); it!=randIntnum.end(); ++it)
    {
      //std::cout << ' ' << *it;
      randInt[id]=*it;
      id++;
    }
    cout << "Total regions: " << totalRegions << endl;	
    /*
     if(preRandomize)
     {
     //This will do what we are doing in the above for loop, but also randomize
     populateRandIntegers(randInt,regionFeatures.size());
     }*/
    
    for(int i=0;i<foldCnt;i++)
    {
      char testFName[1024];
      char trainFName[1024];
      char test1FName[1024];
      char test2FName[1024];
      //	char test3FName[1024];
      
      sprintf(testFName,"%s/test%d.txt",outputDir,i);
      sprintf(trainFName,"%s/train%d.txt",outputDir,i);
      sprintf(test1FName,"%s/test%d_difP.txt",outputDir,i);
      sprintf(test2FName,"%s/test%d_difE.txt",outputDir,i);
      //sprintf(test3FName,"%s/test%d_3.txt",outputDir,i);
      
      ofstream teFile(testFName);
      ofstream trFile(trainFName);
      ofstream test1File(test1FName);
      ofstream test2File(test2FName);
      //ofstream test3File(test3FName);
      
      showHeader(teFile);
      showHeader(trFile);
      showHeader(test1File);
      showHeader(test2File);
      //showHeader(test3File);
      
      testEnd=(i+1)*testSetSize;
      if(testEnd>totalRegions)
      {
        testEnd=totalRegions;
      }
      map<int,int> trainID;
      map<int,int> testID;
      /*
       for(int tid=0;tid<regionNameIDMap.size();tid++)
       {
       
       cout <<tid<<" testID[tid] "  <<testID[tid] << endl;
       }*/
      for(int tid=0;tid<regionNameIDMap.size();tid++)
      {
        int rID=randInt[tid];
        //cout << "tid: " <<tid<< " rID: " << rID <<"test start end: "  <<testBegin << "\t" << testEnd<<endl;
        if(tid<testEnd && tid>=testBegin)
        {
          testID[rID]=0;
        }
        else
        {
          trainID[rID]=0;
        }
        
        
      }
      cout <<"test start_end regions: "  <<testBegin << "_" << testEnd<<endl;
      
      //We already have the pairs. We just need to figure out if they should go to the train set or training set..?
      for(map<string,map<string,map<string,double>*>*>::iterator pIter=pairSet_Filtered.begin();pIter!=pairSet_Filtered.end();pIter++)
      {
        map<string,double>* efeatures=regionFeatures[pIter->first];
        map<string,map<string,double>*>* neighbors=pIter->second;
        map<string,double>* enfeatures=regionNeighborFeatures[pIter->first];
        int eID=regionNameIDMap[pIter->first];
        for(map<string,map<string,double>*>::iterator nIter=neighbors->begin();nIter!=neighbors->end();nIter++)
        {
          map<string,double>* pfeatures=regionFeatures[nIter->first];
          int pID=regionNameIDMap[nIter->first];
          //We just use the randomized IDs to decide whether we will put a pair in the training set or the test set 
          if( (testID.find(eID)!=testID.end()) && (testID.find(pID)!=testID.end()))
          {
            teFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(teFile,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(eID)!=trainID.end()) && (trainID.find(pID)!=trainID.end()))
          {
            trFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(trFile,efeatures,pfeatures,enfeatures,nIter->second);
          }
          // SZ: added by Shilu to generate two easy test pairs
          else if( (trainID.find(eID)!=trainID.end()) )
          {
            test1File <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(test1File,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(pID)!=trainID.end()) )
          {
            test2File <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(test2File,efeatures,pfeatures,enfeatures,nIter->second);
            
          }
        }	 
      }
      teFile.close();
      trFile.close();
      test1File.close();
      test2File.close();
      testBegin=testEnd;
      trainID.clear();
      testID.clear();
    }
    return 0;
  }


int
  Framework::splitRegionsGenPairs_Prerandomize(int foldCnt,const char* outputDir)
  {
    int totalRegions=regionFeatures.size();
    int testSetSize=(int)(totalRegions/foldCnt);
    int testBegin=0;
    int testEnd=0;
    int rID=0;
    map<string,int> regionNameIDMap;
    int* randInt=new int[regionFeatures.size()];
    for(map<string,map<string,double>*>::iterator rIter=regionFeatures.begin();rIter!=regionFeatures.end();rIter++)
    {
      regionNameIDMap[rIter->first]=rID;
      rID++;
    }
    
    //This will do what we are doing in the above for loop, but also randomize
    populateRandIntegers(randInt,regionFeatures.size());
    
    for(int i=0;i<foldCnt;i++)
    {
      char testFName[1024];
      char trainFName[1024];
      char test1FName[1024];
      char test2FName[1024];
      sprintf(testFName,"%s/test%d.txt",outputDir,i);
      sprintf(trainFName,"%s/train%d.txt",outputDir,i);
      sprintf(test1FName,"%s/test%d_difP.txt",outputDir,i);
      sprintf(test2FName,"%s/test%d_difE.txt",outputDir,i);
      
      ofstream teFile(testFName);
      ofstream trFile(trainFName);
      ofstream test1File(test1FName);
      ofstream test2File(test2FName);
      
      showHeader(teFile);
      showHeader(trFile);
      showHeader(test1File);
      showHeader(test2File);
      
      testEnd=(i+1)*testSetSize;
      if(testEnd>totalRegions)
      {
        testEnd=totalRegions;
      }
      map<int,int> testID;
      map<int,int> trainID;
      for(int tid=0;tid<regionNameIDMap.size();tid++)
      {
        int rID=randInt[tid];
        if(tid<testEnd && tid>=testBegin)
        {
          testID[rID]=0;
        }
        else
        {
          trainID[rID]=0;
        }
        
      }
      //We already have the pairs. We just need to figure out if they should go to the train set or training set..?
      for(map<string,map<string,map<string,double>*>*>::iterator pIter=pairSet_Filtered.begin();pIter!=pairSet_Filtered.end();pIter++)
      {
        map<string,double>* efeatures=regionFeatures[pIter->first];
        map<string,map<string,double>*>* neighbors=pIter->second;
        map<string,double>* enfeatures=regionNeighborFeatures[pIter->first];
        int eID=regionNameIDMap[pIter->first];
        for(map<string,map<string,double>*>::iterator nIter=neighbors->begin();nIter!=neighbors->end();nIter++)
        {
          map<string,double>* pfeatures=regionFeatures[nIter->first];
          int pID=regionNameIDMap[nIter->first];
          //We just use the randomized IDs to decide whether we will put a pair in the training set or the test set 
          if( (testID.find(eID)!=testID.end()) && (testID.find(pID)!=testID.end()))
          {
            teFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(teFile,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(eID)!=trainID.end()) && (trainID.find(pID)!=trainID.end()))
          {
            trFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(trFile,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(eID)!=trainID.end()) )
          {
            test1File <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(test1File,efeatures,pfeatures,enfeatures,nIter->second);
          }
          else if( (trainID.find(pID)!=trainID.end()) )
          {
            test2File <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(test2File,efeatures,pfeatures,enfeatures,nIter->second);
            
          }
          
        }	 
      }
      teFile.close();
      trFile.close();	
      test1File.close();
      test2File.close();
      
      testBegin=testEnd;
      trainID.clear();
      testID.clear();
    }
    return 0;
  }

int
  Framework::genPairs(int nCV, const char* outputDir)
  {
    int testSetSize=(int)(totalFilteredPairs/nCV);
    int testBegin=0;
    int testEnd=0;
    for(int f=0;f<nCV;f++)
    {
      char testFName[1024];
      char trainFName[1024];
      sprintf(testFName,"%s/test%d.txt",outputDir,f);
      sprintf(trainFName,"%s/train%d.txt",outputDir,f);
      ofstream teFile(testFName);
      ofstream trFile(trainFName);
      showHeader(teFile);
      showHeader(trFile);
      testEnd=testBegin+testSetSize;
      if(testEnd>totalFilteredPairs)
      {
        testEnd=totalFilteredPairs;
      }
      int testID=testBegin;
      for(map<string,map<string,map<string,double>*>*>::iterator pIter=pairSet_Filtered.begin();pIter!=pairSet_Filtered.end();pIter++)
      {
        map<string,double>* efeatures=regionFeatures[pIter->first];
        map<string,map<string,double>*>* neighbors=pIter->second;
        for(map<string,map<string,double>*>::iterator nIter=neighbors->begin();nIter!=neighbors->end();nIter++)
        {
          map<string,double>* pfeatures=regionFeatures[nIter->first];
          if(testID>=testBegin && testID<testEnd)
          {
            teFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(teFile,efeatures,pfeatures,nIter->second);
          }
          else
          {
            trFile <<pIter->first<<"-"<<nIter->first;
            showFeaturesPair(trFile,efeatures,pfeatures,nIter->second);
          }
          testID++;
        }	 
      }
      teFile.close();
      trFile.close();				
      testBegin=testEnd;
    }
    return 0;
  }


int
  Framework::genPairs_Prerandomize(int nCV, const char* outputDir)
  {
    //first make a giant vector of IDs that map to pairs for totalFilteredPairs
    vector<Pair*> mypairs;	
    for(map<string,map<string,map<string,double>*>*>::iterator pIter=pairSet_Filtered.begin();pIter!=pairSet_Filtered.end();pIter++)
    {
      map<string,map<string,double>*>* neighbors=pIter->second;
      for(map<string,map<string,double>*>::iterator nIter=neighbors->begin();nIter!=neighbors->end();nIter++)
      {
        Pair* p=new Pair;
        p->e.append(pIter->first.c_str());
        p->p.append(nIter->first.c_str());
        mypairs.push_back(p);
      }
    }
    int* randInt=new int[mypairs.size()];
    populateRandIntegers(randInt,mypairs.size());
    int testSetSize=(int)(totalFilteredPairs/nCV);
    int testBegin=0;
    int testEnd=0;
    for(int f=0;f<nCV;f++)
    {
      char testFName[1024];
      char trainFName[1024];
      sprintf(testFName,"%s/test%d.txt",outputDir,f);
      sprintf(trainFName,"%s/train%d.txt",outputDir,f);
      ofstream teFile(testFName);
      ofstream trFile(trainFName);
      showHeader(teFile);
      showHeader(trFile);
      testEnd=testBegin+testSetSize;
      if(testEnd>totalFilteredPairs)
      {
        testEnd=totalFilteredPairs;
      }
      int testID=testBegin;
      //for(map<string,map<string,map<string,double>*>*>::iterator pIter=pairSet_Filtered.begin();pIter!=pairSet_Filtered.end();pIter++)
      for(int i=0;i<mypairs.size();i++)
      {
        int rid=randInt[i];	
        Pair* pair=mypairs[rid];
        map<string,double>* efeatures=regionFeatures[pair->e];
        map<string,double>* pfeatures=regionFeatures[pair->p];
        map<string,map<string,double>*>* neighbors=NULL;
        if(pairSet_Filtered.find(pair->e)==pairSet_Filtered.end())
        {
          cout <<"Did not find enhancer " << pair->e << endl;
          exit(0);
        }
        neighbors=pairSet_Filtered[pair->e];
        map<string,double>* epfeatures=NULL;
        if(neighbors->find(pair->p)==neighbors->end())
        {
          cout <<"Did not find promoter " << pair->p << endl;
          exit(0);
        }
        epfeatures=(*neighbors)[pair->p];
        if(testID>=testBegin && testID<testEnd)
        {
          teFile <<pair->p<<"-"<<pair->e;
          showFeaturesPair(teFile,efeatures,pfeatures,epfeatures);
        }
        else
        {
          trFile <<pair->p<<"-"<<pair->e;
          showFeaturesPair(trFile,efeatures,pfeatures,epfeatures);
        }
        testID++;
      }	 
      teFile.close();
      trFile.close();				
      testBegin=testEnd;
    }
    return 0;
  }


int 
  Framework::populateRandIntegers(int* randInds,int size)
  {
    gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
    struct timeval begintime;
    gettimeofday(&begintime,NULL);
    gsl_rng_set(r,begintime.tv_sec);
    //double step=1.0/(double)size;
    //map<int,int> usedInit;
    for(int i=0;i<size;i++)
    {	
      randInds[i]=i;
    }
    gsl_ran_shuffle(r,randInds,size,sizeof(int));
    
    return 0;
  }


int
  Framework::setFilterZeros(bool flag)
  {
    filterZeros=flag;
    return 0;
  }

int
  Framework::setCorrelation(bool flag)
  {
    correlation=flag;
    return 0;
  }

int 
  Framework::readFeatures(const char* aFName)
  {
    ifstream inFile(aFName);
    char buffer[1024];
    int fID=0;
    while(inFile.good())
    {
      inFile.getline(buffer,1023);
      if(strlen(buffer)<=0)
      {
        continue;
      }	
      if(strstr(buffer,"#")!=NULL)
      {
        continue;
      }
      string featName;
      char fileName[1024];
      char type;
      char* tok=strtok(buffer,"\t");
      int tokCnt=0;
      while(tok!=NULL)
      {
        if(tokCnt==0)
        {
          featName.append(tok);
        }
        else if(tokCnt==1)
        {
          strcpy(fileName,tok);
        }
        else if(tokCnt==2)
        {
          type=tok[0];
        }
        
        tok=strtok(NULL,"\t");
        tokCnt++;
      }
      cout <<"Reading " << fID <<" "<< featName <<" ";
      if(strcmp(featName.c_str(),"Exp")==0)
      {
        readExpression(fileName);
      }
      else if(type=='C')
      {
        Dataset* d=new Dataset;
        if(datasets.find(featName)==datasets.end())
        {
          //featureOrder.push_back(featName);
          datasets[featName]=d;
          d->readDataSet(fileName);
        }
        else
        {
          d=datasets[featName];
          d->readDataSet(fileName);
          
        }
        fID++;	
      }
      else if(type=='M'|| type=='D')
      {
        Dataset* d=new Dataset;
        if(datasets.find(featName)==datasets.end())
        {
          //featureOrder.push_back(featName);
          datasets[featName]=d;
          d->readDataSetAsPIQMotif(fileName,type);
        }
        else
        {
          d=datasets[featName];
          d->readDataSetAsPIQMotif(fileName,type);
        }
        
      }
    }
    inFile.close();
    for(map<string,Dataset*>::iterator dIter=datasets.begin();dIter!=datasets.end();dIter++)
    {
      featureOrder.push_back(dIter->first);
    }
    // added by Shilu to read neighbor features:
    if(featype.compare("Window")==0)
    {
      int k=1;
      while(k<=regionDistNeighborSize)
      {
        char fNameNeighbor[1024];
        for(int f=0;f<featureOrder.size();f++)
        {
          sprintf(fNameNeighbor,"%s+%d",featureOrder[f].c_str(),k);
          string key(fNameNeighbor);
          featureOrderNeighbor.push_back(key);
          //cout << "Feature neighbors: "<<key <<" " <<fNameNeighbor<< endl;
          //sprintf(fNameNeighbor,"%s-%d",featureOrder[f].c_str(),k);
          //string key2(fNameNeighbor);
          //featureOrderNeighbor.push_back(key2);
        }
        k++;
      }
    }
    if(correlation)
    {
      pairFeatureOrder.push_back("Correlation");
    }
    if(geneexp.size()>0)
    {
      pairFeatureOrder.push_back("Exp");
    }
    pairFeatureOrder.push_back("Distance");
    pairFeatureOrder.push_back("Count");
    return 0;
  }

int
  Framework::readExpression(const char* aFName)
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
      string region;
      double exp;
      char* tok=strtok(buffer,"\t");
      int tokCnt=0;
      while(tok!=NULL)
      {
        if(tokCnt==0)
        {
          region.append(tok);
        }
        else
        {
          exp=atof(tok);
        }
        tok=strtok(NULL,"\t");
        tokCnt++;
      }
      if(geneexp.find(region)==geneexp.end())
      {
        geneexp[region]=exp;
      }
      else
      {
        double cval=geneexp[region];
        if(cval<exp)
        {
          geneexp[region]=exp;
        }
      }
    }
    inFile.close();
    return 0;
  }

double
  Framework::getExp(string& region)
  {
    double exp=0;
    if(geneexp.find(region)!=geneexp.end())
    {
      exp=geneexp[region];
    }
    return exp;
  }

// modified by Shilu to generate region features only:
int 
  Framework::generateFeatureFiles_Concat(const char* aFName)
  {
    char aFName_E[1024];
    sprintf(aFName_E,"%s/region.txt",aFName);
    ofstream eFile(aFName_E);
    eFile <<"Region";
    for(int f=0;f<featureOrder.size();f++)
    {
      eFile <<"\t" << featureOrder[f];
    }
    eFile << endl;
    map<string,int> shownStatus;
    map<string,int> featureFrequency;
    //These are the pairs from the sparse matrix
    for(map<string,map<string,double>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
    {
      string& e=(string&)eIter->first;
      map<string,double>* efeatures=NULL;
      if(regionFeatures.find(e)==regionFeatures.end())
      {
        continue;
      }
      else
      {
        efeatures=regionFeatures[e];
      }
      if(shownStatus.find(e)==shownStatus.end())	
      {
        showFeatures(e,efeatures,eFile);
        shownStatus[e]=0;
      }
      map<string,double>* pSet=eIter->second;
      for(map<string,double>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
      {
        string& p=(string&)pIter->first;
        if(strcmp(e.c_str(),p.c_str())==0)
        {
          //Don't allow self loops now
          continue;
        }
        map<string,double>* pfeatures=NULL;
        if(regionFeatures.find(p)==regionFeatures.end())
        {
          continue;
        }
        else
        {
          pfeatures=regionFeatures[p];
        }
        if(shownStatus.find(p)==shownStatus.end())	
        {
          // modified by Shilu, should be pfeatures not efeatures
          showFeatures(p,pfeatures,eFile);
          shownStatus[p]=0;
        }
        //At this state the e and p regions must have features
        //generate any pairwise feature
        map<string,map<string,double>*>* regionNeighs=NULL;
        if(pairSet_Filtered.find(e)==pairSet_Filtered.end())
        {
          regionNeighs=new map<string,map<string,double>*>;
          pairSet_Filtered[e]=regionNeighs;
        }
        else
        {
          regionNeighs=pairSet_Filtered[e];
        }
        //Now make any pairwise features. For now we only a have correlation.
        map<string,double>* pairFeatures=new map<string,double>;
        if(correlation)
        {
          double cc=getCCFeature(pfeatures,efeatures);
          (*pairFeatures)["Correlation"]=cc;
        }
        if(geneexp.size()>0)
        {
          double exp=getExp(p);
          (*pairFeatures)["Exp"]=exp;
        }
        int dist=getDistance(e,p);
        (*pairFeatures)["Distance"]=dist;
        (*pairFeatures)["Count"]=pIter->second;
        (*regionNeighs)[p]=pairFeatures;
        totalFilteredPairs++;
      }
    }
    
    return 0;
    
  }

int
  Framework::showHeader(ofstream& oFile)
  {
    oFile <<"Pair";
    for(int f=0;f<featureOrder.size();f++)
    {
      
      oFile <<"\t"<< featureOrder[f]<<"_E";
      
    }
    for(int f=0;f<featureOrder.size();f++)
    {
      
      oFile <<"\t"<< featureOrder[f]<<"_P";
      
    }
    if (featype.compare("Outerprod")==0){
      for(int f=0;f<featureOrder.size();f++)
      {	
        for(int j=0;j<featureOrder.size();j++)
        {	
          oFile <<"\t"<< featureOrder[f]<<"_E_" << featureOrder[j] << "_P";
        }
      }
    }
    else if (featype.compare("Window")==0){
      for(int f=0;f<featureOrder.size();f++)
      {
        oFile <<"\t"<< featureOrder[f]<<"_W";
      }	
    }
    for(int f=0;f<pairFeatureOrder.size();f++)
    {
      oFile <<"\t" << pairFeatureOrder[f];
    }
    oFile << endl;
    return 0;
  } 

int
  Framework::showFeaturesPair(ofstream& oFile, map<string,double>* efeatures, map<string,double>* pfeatures,map<string,double>* pairFeatures)
  {
    for(int f=0;f<featureOrder.size();f++)
    {
      double featVal_e=0.01;
      if(efeatures->find(featureOrder[f])!=efeatures->end())
      {
        featVal_e=(*efeatures)[featureOrder[f]]+0.01;
      }
      oFile <<"\t" << featVal_e;
    }
    for(int f=0;f<featureOrder.size();f++)
    {
      double featVal_p=0.01;
      if(pfeatures->find(featureOrder[f])!=pfeatures->end())
      {
        featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
      }
      oFile <<"\t" << featVal_p;
    }
    for(int f=0;f<pairFeatureOrder.size();f++)
    {
      double featVal=0;
      if(pairFeatures->find(pairFeatureOrder[f])!=pairFeatures->end())
      {
        featVal=(*pairFeatures)[pairFeatureOrder[f]];
      }
      oFile <<"\t"  <<featVal; 
    }
    oFile <<endl;
    return 0;
  }

// added by Shilu to generate CTCF directionality:
int
  Framework::showFeaturesPair(ofstream& oFile, map<string,double>* efeatures, map<string,double>* pfeatures,map<string,double>* enfeatures,map<string,double>* pairFeatures)
  {
    double Eplus=0;
    double Eminus=0;
    double Pplus=0,Pminus=0;
    double featVal=0;
    //map<string,double>* enfeatures=regionNeighborFeatures[pIter->first];
    for(int f=0;f<featureOrder.size();f++)
    {
      double featVal_e=0.01;
      if(efeatures->find(featureOrder[f])!=efeatures->end())
      {
        //Eplus=(*efeatures)[featureOrder[f]];
        featVal_e=(*efeatures)[featureOrder[f]]+0.01;
      }
      //cout << featureOrder[f] << endl;
      
      oFile <<"\t" << featVal_e;
      
    }
    
    
    
    for(int f=0;f<featureOrder.size();f++)
    {
      double featVal_p=0.01;
      if(pfeatures->find(featureOrder[f])!=pfeatures->end())
      {
        featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
      }
      oFile <<"\t" << featVal_p;
      
    }
    
    //add a wrapper to do this or not:	
    
    if(featype.compare("Outerprod")==0)
    {	
      for(int f=0;f<featureOrder.size();f++)
      {
        double featVal_e=0;
        if(efeatures->find(featureOrder[f])!=efeatures->end())
        {
          
          featVal_e=(*efeatures)[featureOrder[f]];
        }	
        for(int j=0;j<featureOrder.size();j++)
        {
          double featVal_p=0;
          if(pfeatures->find(featureOrder[j])!=pfeatures->end())
          {
            featVal_p=(*pfeatures)[featureOrder[j]];
          }
          featVal=featVal_e*featVal_p+0.01;
          oFile <<"\t"  <<featVal;
        }
        
      }
      
    }else if(featype.compare("Window")==0)
    {
      int neighborsize=(*pairFeatures)[pairFeatureOrder[0]]/BinSize;
      //cout << "Neighbor size is " << neighborsize << endl;
      //add winow feature between two regions:
      
      showFeaturesPairNeighbors(oFile, enfeatures,neighborsize);
      
    }
    
    for(int f=0;f<pairFeatureOrder.size();f++)
    {
      double featVal=0;
      if(pairFeatures->find(pairFeatureOrder[f])!=pairFeatures->end())
      {
        featVal=(*pairFeatures)[pairFeatureOrder[f]];
        
      }
      oFile <<"\t"  <<featVal;
    }
    oFile <<endl;
    return 0;
  }

//added by Shilu to show window average features:
int
  Framework::showFeaturesPairNeighbors(ofstream& oFile, map<string,double>* enfeatures,int neighborsize)
  {
    //loop for each feature: (get average signal of window feature)
    for(int i=0;i<featureOrder.size();i++)
    {
      double featVal_e=0;
      for(int f=i;f<(neighborsize*featureOrder.size());f+=featureOrder.size())
      {	
        string& fKey=featureOrderNeighbor[f];
        if(enfeatures->find(fKey)!=enfeatures->end())
        {
          featVal_e+=(*enfeatures)[fKey];
        }
        
      }
      if(neighborsize>0){
        featVal_e=featVal_e/neighborsize+0.01;
      }
      else{
        featVal_e=0.01;
      }
      oFile <<"\t" << featVal_e;
    }
    return 0;
  }


double
  Framework::getCCFeature(map<string,double>* pfeatures,map<string,double>* efeatures)
  {	
    double cc=0;
    vector<double> pvect;
    vector<double> evect;
    Distance dist;
    for(int i=0;i<featureOrder.size();i++)
    {
      if(pfeatures->find(featureOrder[i])==pfeatures->end())
      {
        pvect.push_back(0.001);
      }
      else
      {
        pvect.push_back((*pfeatures)[featureOrder[i]]);
      }
      if(efeatures->find(featureOrder[i])==efeatures->end())
      {
        evect.push_back(0.001);
      }
      else
      {
        evect.push_back((*efeatures)[featureOrder[i]]);
      }
    }
    cc=dist.computeCC(pvect,evect);
    pvect.clear();
    evect.clear();
    return cc;
  }


double
  Framework::getCCFeature_Bin(map<string,double>* pfeatures,map<string,double>* efeatures)
  {	
    double cc=0;
    double hit=0;
    Distance dist;
    for(map<string,double>::iterator pIter=pfeatures->begin();pIter!=pfeatures->end();pIter++)
    {
      if(efeatures->find(pIter->first)!=efeatures->end())
      {
        hit=hit+1;
      }
    }
    cc=hit/(((double)(pfeatures->size()+efeatures->size()))-hit);
    //cc=cc/((double)featureOrder.size());
    //cc=1-cc;
    return cc;
  }


int
  Framework::populateFeaturesForRegion(map<string,double>* populateMe, string& key)
  {
    if(regionSet.find(key)==regionSet.end())
    {
      cout <<"No region with name " << key<< endl;
      exit(0);
    }
    Region* r=regionSet[key];
    for(map<string,Dataset*>::iterator dIter=datasets.begin();dIter!=datasets.end();dIter++)
    {
      Dataset* d=dIter->second;
      double fval=d->getFeatureVal(r->chromosome,r->begin,r->end);
      if(fval<0)
      {
        continue;
      }
      (*populateMe)[dIter->first]=fval;
    }	
    return 0;
  }

// added 0.01 by Shilu to make consistent with pair features
int
  Framework::showFeatures(string& fName,map<string,double>* fVect,ofstream& oFile)
  {
    oFile<< fName;
    for(int i=0;i<featureOrder.size();i++)
    {
      if(fVect->find(featureOrder[i])==fVect->end())
      {
        oFile<<"\t-1";
      }
      else
      {
        oFile<<"\t" <<(*fVect)[featureOrder[i]]+0.01;
      }
    }
    oFile<< endl;
    return 0;
  }

int
  Framework::filterRegions()
  {
    map<string,int> featureFrequency;
    for(map<string,Region*>::iterator rIter=regionSet.begin();rIter!=regionSet.end();rIter++)
    {
      string& e=(string&)rIter->first;
      map<string,double>* efeatures=NULL;
      if(regionFeatures.find(e)==regionFeatures.end())
      {
        efeatures=new map<string,double>;
        populateFeaturesForRegion(efeatures,e);
        if(efeatures->size()>0)
        {
          regionFeatures[e]=efeatures;
        }
        else
        {
          delete efeatures;
        }
        for(map<string,double>::iterator fIter=efeatures->begin();fIter!=efeatures->end();fIter++)
        {
          if(featureFrequency.find(fIter->first)==featureFrequency.end())
          {
            featureFrequency[fIter->first]=1;
          }
          else
          {
            featureFrequency[fIter->first]=featureFrequency[fIter->first]+1;
          }
        }	
      }
    }
    
    
    cout <<"Filtered " << regionFeatures.size() << " out of total " << regionSet.size() << endl; 
    //Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
    //although the user can decide
    for(int f=0;f<featureOrder.size();f++)
    {
      if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
      {
        cout <<"Get rid of " << featureOrder[f] << endl;	
      }
    }
    return 0;
  }

//Now we add features for neighbor regions
int
  Framework::addNeighborRegions()
  {
    //int bin=5000;
    for(map<string,map<string,double>*>::iterator rIter=regionFeatures.begin();rIter!=regionFeatures.end();rIter++)
    {
      map<string,double>* rfeatures=rIter->second;
      Region* r=regionSet[rIter->first];
      //bin=r->end - r->begin;
      //cout << "Bin is: " << (r->end - r->begin) << endl;
      
      //For now assume we have a fixed 5000bp window
      //add features of the k next neighbors
      int k=1;
      
      map<string,double>* additionalFeatures=new map<string,double>;
      while(k<=regionDistNeighborSize)
      {
        int neighCoord=r->begin+(k*BinSize);
        char neighborName[1024];
        sprintf(neighborName,"%s_%d_%d",r->chromosome.c_str(),neighCoord,neighCoord+BinSize);
        string key(neighborName);
        if(regionFeatures.find(key)==regionFeatures.end())
        {
          k++;
          continue;		
        }
        map<string,double>* neighborFeatures=regionFeatures[key];
        for(map<string,double>::iterator sIter=neighborFeatures->begin();sIter!=neighborFeatures->end();sIter++)
        {
          char fKey[1024];
          sprintf(fKey,"%s+%d",sIter->first.c_str(),k);
          string fKeyS(fKey);
          (*additionalFeatures)[fKeyS]=sIter->second;
          //cout << "Neihgbor features: " << sIter->second << endl;
        }
        k++;
      }
      regionNeighborFeatures[rIter->first]=additionalFeatures;
    }
    return 0;
  }





int
  Framework::getDistance(string& ekey, string& pkey)
  {
    Region* e=regionSet[ekey];
    Region* p=regionSet[pkey];
    int dist=0;
    int first1=e->begin;
    int first2=e->end;
    int second1=p->begin;
    int second2=p->end;
    if(first1>second1)
    {
      first1=p->begin;
      first2=p->end;
      second1=e->begin;
      second2=e->end;
    }
    if(second1<first2) 
    {
      //there is overlap
      dist=0;
    }
    else
    {
      dist=second1-first2;
    }
    return dist;
  }
int
  Framework::setPreRandomize(bool flag)
  {	
    preRandomize=flag;
    return 0;
  }


int
  main(int argc, const char** argv)
  {
    // added an argument by Shilu to generate different additional features:
    if(argc!=10)
    {
      cout <<"Usage: getFeatures sparseMatrix dist foldcv splittype[regionwise|pairwise] featurefile correlation[yes|no] outputdir prerandomize_pairs[yes|no] featuretype[Window|Pconcat|Outerprod]" << endl;
      return 0;
    }
    Framework fw;
    fw.setMaxDist(atoi(argv[2]));
    fw.setFeatype(argv[9]);
    if(strcmp(argv[6],"yes")==0)
    {
      fw.setCorrelation(true);
    }
    else
    {	
      fw.setCorrelation(false);
    }
    //strcpy(featype,argv[9]);
    fw.readPairs(argv[1]);
    fw.readFeatures(argv[5]);
    fw.filterRegions();
    
    // add feature regions:
    if(strcmp(argv[9],"Window")==0){
      fw.addNeighborRegions();
    }
    
    //Generate features one time
    fw.generateFeatureFiles_Concat(argv[7]);
    if(strcmp(argv[8],"yes")==0)
    {
      fw.setPreRandomize(true);
    }
    if(strcmp(argv[4],"regionwise")==0)
    {
      fw.splitPairs(atoi(argv[3]),true,argv[7]);
    }
    else
    {
      fw.splitPairs(atoi(argv[3]),false,argv[7]);
    }
    return 0;
  }
