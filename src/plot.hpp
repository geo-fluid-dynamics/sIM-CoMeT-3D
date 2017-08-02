#ifndef PLOT_H
#define PLOT_H

#include "packages/gnuplot_i/gnuplot_i.h"
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
