#ifndef PLOT_H
#define PLOT_H

#include "packages/gnuplot_i/gnuplot_i.h"
#include "field.hpp"
/* #include "packages/matplotlib-cpp/matplotlibcpp.h" */

/* namespace plt = matplotlibcpp; */

class Plot
{
	public:
		gnuplot_ctrl * gnu;

		Plot();
		~Plot();

		void image(Field * field);
		void surface(Field * field);
		void pallete(Field * field);
		void setSurface();

};

#endif
