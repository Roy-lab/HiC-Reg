#ifndef _GENANATOMY_CONV
#define _GENANATOMY_CONV
class GenanatomyConverter
{
	public:
		GenanatomyConverter();
		~GenanatomyConverter();
		int setTemplateProject(const char*);
		int genTemplateFile(const char* outputFName,const char* exprFName, const char* modulenetFName, const char* mainDirName);
	private:
		int readFileData(const char*);
		char* fileData;
		char templateFName[1024];
};

#endif

