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
		void setSurface();

};

#endif
