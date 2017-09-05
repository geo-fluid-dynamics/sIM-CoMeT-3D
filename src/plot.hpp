#ifndef PLOT_H
#define PLOT_H

#include "field.hpp"

class Plot
{
	public:
		FILE * gnu;

		Plot();
		~Plot();

		void image(Field * field);
		void surface(Field * field);
		void pallete(Field * field);
		void save(Field * field);
		void setSurface();
		void saveImage(std::map< std::string, Field *> map, std::string directory);
		void saveSurface(std::map< std::string, Field *> map, std::string directory);

};

#endif
