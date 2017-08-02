#include "plot.hpp"

Plot::Plot()
{
	gnu = popen("gnuplot -p", "w");
}

Plot::~Plot()
{
	pclose(gnu);
}

void Plot::image(Field * field)
{
	fprintf(gnu, "plot '-' u 1:2:3 with image pixels\n");
	for(int i = 1; i<field->nx-1; i++)
	{
		for(int j=1; j<field->ny-1; j++)
		{
			fprintf(gnu, "%d %d %e\n", i, j, field->get(i,j));
		}
		fprintf(gnu, "\n");
	}
	fprintf(gnu, "e\n");
}

void Plot::surface(Field * field)
{
	fprintf(gnu, "splot '-' u 1:2:3 with pm3d\n");
	for(int i = 1; i<field->nx-1; i++)
		for(int j=1; j<field->ny-1; j++)
		{
			fprintf(gnu, "%d %d %e\n", i, j, field->get(i, j));
		}
	fprintf(gnu, "e\n");

}

void Plot::setSurface()
{
	fprintf(gnu, "set dgrid3d\n");
	fprintf(gnu, "set pm3d interpolate 1,1\n");
}

void Plot::pallete(Field * field)
{
	fprintf(gnu, "splot '-' u 1:2:3:4 palette pt 5\n");
	for(int i=0; i<field->nx; i++)
	{
		for(int j=0; j<field->ny; j++)
			for(int k=0; k<field->nz; k++)
			{
				fprintf(gnu, "%d %d %d %e\n", i, j, k, field->get(i,j,k));
			}
		fprintf(gnu, "\n");
	}
	fprintf(gnu, "e\n");

}
