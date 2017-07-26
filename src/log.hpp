#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>

class Log
{
	private:
		std::string filename;

	public:
		FILE * ptr;
		Log();
		~Log();
		void writeHeader();
		void writeFooter();
		/* void write(std::string str); */
};

#endif
