#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include <vector>
using namespace std;
class Framework
{
	public:
	Framework();
	~Framework();
	int readTSS(const char*);
	int readSignals(const char*);
	int readInput(const char*);
	int normalize();
	int getSignalInRegion(const char*);
	int readChromosomeSizes(const char*);
	private:
	int readData(const char* fName,map<string,vector<double>*>& data,double& );
	struct Gene
	{
		string genename;
		map<int,int> tsses;
		char strand;
		int begin;
		int end;
	};
	map<string,map<int,double>*> binInfo;
	map<string,map<string,Gene*>*> tssSet;
	map<string,map<int,Gene*>*> tsscoordSet;
	map<string,vector<double>*> logRatios;
	map<string,vector<double>*> signal;
	map<string,vector<double>*> bg;
	double signalAvg;
	double inputAvg;
	map<string,int> chromosomeSizes;
};
#endif
