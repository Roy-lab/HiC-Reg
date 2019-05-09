#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
class Dataset;
class Framework
{
	public:
	Framework();
	~Framework();
	int setMaxDist(int);
	int setFeatype(string feat);
	int readPairs(const char*);
	//Must include features in required order
	int readFeatures(const char*);
	int readExpression(const char*);
	int generateFeatureFiles_Concat(const char* aFName);
	int setFilterZeros(bool);
	int setCorrelation(bool);
	double getCCFeature(map<string,double>*, map<string,double>*);
	double getCCFeature_Bin(map<string,double>*, map<string,double>*);
	double getExp(string&);
	int filterRegions();
	int addNeighborRegions();

	int splitPairs(int nCV,bool regionWise, const char*);
	int setPreRandomize(bool);
	struct Region
	{
		string chromosome;
		int begin;
		int end;	
	};
	struct Pair
	{
		string e;
		string p;
	};
	private:
	int populateRandIntegers(int* randInds,int size);
	int showFeatures(string&,map<string,double>*,ofstream&);
	int populateFeaturesForRegion(map<string,double>*,string&);
	int showHeader(ofstream& );
	int showFeaturesPair(ofstream& oFile, map<string,double>* efeatures, map<string,double>* pfeatures,map<string,double>* enfeatures,map<string,double>* pairFeatures);
	int showFeaturesPair(ofstream& oFile, map<string,double>* efeatures, map<string,double>* pfeatures,map<string,double>* pairFeatures);
	int splitRegionsGenPairs(int foldCnt,const char*);
	int splitRegionsGenPairs_Prerandomize(int foldCnt,const char*);
	int splitRegionsGenPairs_NoPrerandomize(int foldCnt,const char* outputDir);
	int genPairs(int nCV, const char* outputDir);
	int genPairs_Prerandomize(int nCV, const char* outputDir);
	int showFeaturesPairNeighbors(ofstream& oFile, map<string,double>* enfeatures,int neighborsize);

	map<string,map<string,double>*> regionFeatures;
	map<string,map<string,double>*> regionNeighborFeatures;
	map<string,Dataset*> datasets;
	vector<string> featureOrder;
	vector<string> pairFeatureOrder;
	vector<string> featureOrderNeighbor;
	map<string,map<string,double>*> pairSet;
	map<string,map<string,map<string,double>*>*> pairSet_Filtered;
	map<string,Region*> regionSet;
	int getDistance(string&,string&);
	bool filterZeros;
	map<string,double> geneexp;
	bool correlation;
	int maxDist;
	int totalFilteredPairs;
	bool preRandomize;
	string featype;
	int regionDistNeighborSize;
	int BinSize;
	int initial;
};

#endif
