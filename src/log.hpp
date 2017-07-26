#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>

class Log
{
	public:
		std::string filename;
		FILE * ptr;
		Log();
		~Log();
		void writeHeader();
		void writeFooter();
		/* void write(std::string str); */
};

#endif
