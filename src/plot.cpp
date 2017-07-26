#include "plot.hpp"

Plot::Plot()
{
	gnu = gnuplot_init();
}

Plot::~Plot()
{
	gnuplot_close(gnu);
}

void Plot::image(Field * field)
{
	gnuplot_resetplot(gnu);
	gnuplot_cmd(gnu, "plot '-' u 1:2:3 with image pixels");
	/* gnuplot_cmd(gnu, "splot '-' u 1:2:3 with pm3d"); */
	for(int i = 1; i<field->nx-1; i++)
	{
		for(int j=1; j<field->ny-1; j++)
		{
			gnuplot_cmd(gnu, "%d %d %e", i, j, field->get(i,j));
		}
		/* gnuplot_cmd(gnu, ""); */
	}
	gnuplot_cmd(gnu, "e");
}

void Plot::surface(Field * field)
{

	gnuplot_resetplot(gnu);
	gnuplot_cmd(gnu, "splot '-' u 1:2:3 with pm3d");
	for(int i = 1; i<field->nx-1; i++)
		for(int j=1; j<field->ny-1; j++)
		{
			gnuplot_cmd(gnu, "%d %d %e", i, j, field->get(i, j));
		}
	gnuplot_cmd(gnu, "e");

}

void Plot::setSurface()
{
	gnuplot_cmd(gnu, "set dgrid3d");
	gnuplot_cmd(gnu, "set pm3d interpolate 1,1");
}

void Plot::pallete(Field * field)
{
	gnuplot_resetplot(gnu);
	gnuplot_cmd(gnu, "splot '-' u 1:2:3:4 palette pt 5");
	for(int i=0; i<field->nx; i++)
	{
		for(int j=0; j<field->ny; j++)
			for(int k=0; k<field->nz; k++)
			{
				gnuplot_cmd(gnu, "%d %d %d %e", i, j, k, field->get(i,j,k));
			}
		gnuplot_cmd(gnu, "\n");
	}
	gnuplot_cmd(gnu, "e");

}
