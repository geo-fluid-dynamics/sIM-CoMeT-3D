#include "log.hpp"
#include <time.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

Log::Log()
{
	char * timestr = (char*)malloc(16);
	char * logFileName = (char*)malloc(30);

	time_t timestamp;
	timestamp = time(NULL);
	strftime(timestr, 16, "%Y%m%d_%H%M%S", localtime(&timestamp));

	/* assert(logFileName); */
	snprintf(logFileName, 30, "logs/log_%s.dat", timestr);
	free(timestr);

	struct stat st = {0};
	if (stat("logs", &st) == -1)
		mkdir("logs", 0755);
	if (stat("outputs", &st) == -1)
		mkdir("outputs", 0755);

	logPTR = fopen(logFileName, "w");
	filename = std::string(logFileName);
	free(logFileName);


}

Log::Log(std::string filepath)
{
	filename = filepath + "/log.dat";
	logPTR = fopen(filename.c_str(), "w");

	std::string inputsFileName = filepath + "/inputs.dat";
	inputsPTR = fopen(inputsFileName.c_str(), "w");

	fprintf(inputsPTR, "##### BEGIN HEADER #####\n");
	FILE *inputFile;

	/*write inputs as header in log file*/
	inputFile = fopen( "inputs.ini", "r");
	char ch;
	if(inputFile)
	{
		while((ch = fgetc(inputFile)) != EOF)
		{
			fputc(ch, inputsPTR);
		}
	}
	fclose(inputFile);
	fprintf(inputsPTR, "##### END HEADER #####\n\n");

}

Log::~Log()
{
	fclose(logPTR);
}

void Log::writeHeader()
{

	fprintf(logPTR, "##### BEGIN HEADER #####\n");
	FILE *inputFile;

	/*write inputs as header in log file*/
	inputFile = fopen( "inputs.ini", "r");
	char ch;
	if(inputFile)
	{
		while((ch = fgetc(inputFile)) != EOF)
		{
			fputc(ch, logPTR);
		}
	}
	fclose(inputFile);
	fprintf(logPTR, "##### END HEADER #####\n\n");

}

void Log::writeFooter()
{
	fprintf(logPTR, "\n***** EOF *****");
}
