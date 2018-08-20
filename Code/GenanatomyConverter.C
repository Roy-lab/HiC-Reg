#include <iostream>
#include <fstream>
#include <string.h>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
using namespace std;
#include "GenanatomyConverter.H"
GenanatomyConverter::GenanatomyConverter()
{
	fileData=NULL;
}

GenanatomyConverter::~GenanatomyConverter()
{
	if(fileData!=NULL)
	{
		delete[] fileData;
	}
}

int 
GenanatomyConverter::setTemplateProject(const char* aFName)
{
	readFileData(aFName);
	return 0;
}

int 
GenanatomyConverter::genTemplateFile(const char* outputFName,const char* exprFName, const char* modulenetFName, const char* mainDirName)
{
	char* pos=strstr(fileData,"</GUI>");
	if(pos==NULL)
	{
		cout <<"Bad project template" << endl;
		return -1;
	}
	pos=pos+strlen("</GUI>");
	*pos='\0';
	ofstream oFile(outputFName);
	oFile << fileData << endl;
	oFile <<"<DataProvider Type=\"File\">" << endl;;
	oFile <<"Inner Filename=\"" << exprFName << "\" Species=\"Saccharomyces_cerevisiea\" DataType=\"GENE\" NameColumn=\"0\" NumOfNameColumns=\"1\"/>" << endl;
	oFile <<"</DataProvider>" << endl;
	oFile <<"<AttributeManagers><X xAxis=\"true\"> <Used/><GOProperties/><Provider DisplayOnlyNonEmpty=\"false\"/></X>" <<endl; 
	oFile <<"<Y xAxis=\"false\"> <Used/><GOProperties/><Provider DisplayOnlyNonEmpty=\"false\"/></Y></AttributeManagers>" <<endl; 

	oFile <<"<FiltersManager>" << endl;
	oFile <<"<Inner><Filter Type=\"ModNetFilter\"><Inner FileName=\""<<modulenetFName << "\"/></Filter></Inner>" << endl;
	oFile <<"</FiltersManager>" << endl;

	oFile <<"<DataCenter Specie=\"Saccharomyces_cerevisiae\"> <LocalTables/><Use/></DataCeneter>" << endl;
	oFile <<"<ProjectParameters FileSeparator=\" Specie=\"Saccharomyces_cerevisiae\" Name=\"FGModules\" MainDirector=\""<< mainDirName
		<<"\" DataType=\"GENE\">" << endl;
     	oFile <<"<Keys/></ProjectParamaters><pValues/></root>" << endl;
	oFile.close();       
	return 0;
}

int 
GenanatomyConverter::readFileData(const char* aFName)
{
	struct stat fileInfo;
	if(stat(aFName,&fileInfo)==-1)
	{
		return -1;
	}
	int size=fileInfo.st_size;
	fileData=new char[size];
	char* buffPtr=fileData;
	ifstream inFile(aFName);
	int readChunk=0;
	while(inFile.good() && ((readChunk*8196)<size))
	{
		inFile.read(buffPtr,8196);
		buffPtr=buffPtr+8196;
	}
	inFile.close();
	return 0;
}
