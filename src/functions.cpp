#include "functions.hpp"
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>

std::vector<std::string> setup_directory ()
{
	char * timestr = (char*)malloc(16);
	time_t timestamp = time(NULL);
	strftime(timestr, 16, "%Y%m%d_%H%M%S", localtime(&timestamp));
	assert(timestr);
	std::string solveID(timestr);
	free(timestr);

	std::string directory = "outputs/" + solveID;
	std::string avgDir = directory + "/avg";
	std::string fieldsDir = directory + "/fields";

	std::vector<std::string> directories = {directory, avgDir, fieldsDir};

	struct stat st = {0};
	if (stat(directory.c_str(), &st) == -1)
		mkdir(directory.c_str(), 0755);
	if (stat(avgDir.c_str(), &st) == -1)
		mkdir(avgDir.c_str(), 0755);
	if (stat(fieldsDir.c_str(), &st) == -1)
		mkdir(fieldsDir.c_str(), 0755);

	return directories;
}

