#ifndef _DATASET
#define _DATASET

#include <map>
#include <string>
using namespace std;
class Dataset
{
	public:
	Dataset();
	~Dataset();
	int readDataSet(const char*);
	int readDataSetAsPIQMotif(const char*,char type);
	double getFeatureVal(string& chrom, int begin, int end);
	struct Peak
	{
		int start;
		int end;
		double lqval;
		double signal;
		int peakCnt;
	};
	
	private:
	map<string,map<int,Peak*>*> peakSet;
	char type;
};
#endif
