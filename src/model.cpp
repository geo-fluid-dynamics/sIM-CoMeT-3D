/* model.cpp: Implementation of the Model class
 *
 * This is a God class that handles ALL of the tasks involved in the program.
 * The class has one public method: solve() that runs the ccmsolve algorithm
 * based on the inputs available in the public members of the class (defined in the
 * header file) and the initialization of Tw/qw;
 *
 * solve() invokes many private functions that serve as modules in the ccmsolve
 * algorithm.
 *
 */
#include "model.hpp"

#include <cmath>
#include <iostream>
#include <vector>
#include <assert.h>
#include <time.h>

#include "PDE.hpp"
#include "solver.hpp"
#include "estimator.hpp"
#include "log.hpp"
#include "plot.hpp"

#include "packages/inih/INIReader.h"

/* Model::Model() */
/* { */

/* 	Tw      = new Field(nx, ny, 1, Lx, Ly); */
/* 	qw      = new Field(nx, ny, 1, Lx, Ly); */
/* 	delta   = new Field(nx, ny, 1, Lx, Ly); */
/* 	p       = new Field(nx, ny, 1, Lx, Ly); */
/* 	qStefan = new Field(nx, ny, 1, Lx, Ly); */
/* 	qNorth  = new Field(nx, ny, 1, Lx, Ly); */
/* 	T       = new Field(nx, ny, nz, Lx, Ly); */
/* 	u       = new Field(nx, ny, nz, Lx, Ly); */
/* 	v       = new Field(nx, ny, nz, Lx, Ly); */
/* 	w       = new Field(nx, ny, nz, Lx, Ly); */

/* 	init_fields(); */

/* } */

/*
 * Constructor for model.
 * Parses input file for inputs and initializes fields.
 *
 */
Model::Model(std::string iniPath)
{

	INIReader reader(iniPath);

	if (reader.ParseError() < 0) {
		std::cout << "Can't load 'test.ini'\n";
		exit(-1);
	}

	mu     = reader.GetReal("constants", "mu"    , 0.001  );
	Tm     = reader.GetReal("constants", "Tm"    , 0      );
	Tinf   = reader.GetReal("constants", "Tinf"  , -20    );
	hm     = reader.GetReal("constants", "hm"    , 3.337e5);
	cpL    = reader.GetReal("constants", "cpL"   , 4222.22);
	cpS    = reader.GetReal("constants", "cpS"   , 2049.41);
	rhoL   = reader.GetReal("constants", "rhoL"  , 1000   );
	rhoS   = reader.GetReal("constants", "rhoS"  , 920    );
	kL     = reader.GetReal("constants", "kL"    , 0.57   );
	Lx     = reader.GetReal("constants", "Lx"    , 0.075  );
	Ly     = reader.GetReal("constants", "Ly"    , 0.075  );
	Fscrew = reader.GetReal("constants", "Fscrew", 375    );

	southBC = reader.GetBoundary("boundaryConditions", "southBC", DIRICHLET);
	sidesBC = reader.GetBoundary("boundaryConditions", "sidesBC", NEUMANN);

	TwExp   = reader.Get("boundaryConditions", "TwExp", "40+10*(x/Lx+y/Ly)");
	qwExp   = reader.Get("boundaryConditions", "qwExp", "10");

	nx = reader.GetInteger("gridSizes", "nx", 30);
	ny = reader.GetInteger("gridSizes", "ny", 30);
	nz = reader.GetInteger("gridSizes", "nz", 30);

	maxFindIter         = reader.GetInteger("parameters", "maxFindIter"     , 100);
	maxMainIter         = reader.GetInteger("parameters", "maxMainIter"     , 9999);

	MTol                = reader.GetReal("parameters", "MTol"               , 1e-10);
	FTol                = reader.GetReal("parameters", "FTol"               , 1e-6 );
	allowedMaxFluxError = reader.GetReal("parameters", "allowedMaxFluxError", 100  );
	allowedRelativeMFE  = reader.GetReal("parameters", "allowedRelativeMFE" , 1e-6 );
	deltaCoeffMin       = reader.GetReal("parameters", "deltaCoeffMin"      , 0.7  );
	deltaCoeffMax       = reader.GetReal("parameters", "deltaCoeffMax"      , 1.3  );
	deltaRelax          = reader.GetReal("parameters", "deltaRelax"         , 0.05 );

	dx = 2.0*Lx/(nx-1);
	dy = 2.0*Ly/(ny-1);
	dz = 1.0/(nz - 1);

	hmStar = hm + cpS*(Tm-Tinf);
	alpha = kL/(rhoL*cpL);

	/* theta   = reader.GetReal("boundaryConditions", "theta", theta); */
	/* radianTheta = theta*3.1415927/180; */

	Tw      = new Field(nx, ny, 1, Lx, Ly);
	qw      = new Field(nx, ny, 1, Lx, Ly);
	delta   = new Field(nx, ny, 1, Lx, Ly);
	p       = new Field(nx, ny, 1, Lx, Ly);
	qStefan = new Field(nx, ny, 1, Lx, Ly);
	qNorth  = new Field(nx, ny, 1, Lx, Ly);
	T       = new Field(nx, ny, nz, Lx, Ly);
	u       = new Field(nx, ny, nz, Lx, Ly);
	v       = new Field(nx, ny, nz, Lx, Ly);
	w       = new Field(nx, ny, nz, Lx, Ly);

	init_fields();

	Field * num = bcSouth->differentiate(CX2, Y);
	Field * den = bcSouth->differentiate(CX2, X);

	double arg = num->average() / den->average();
	radianTheta = (std::isnan(arg))? 0 : atan(arg);
	theta = radianTheta * 180 / M_PI;

	delete num;
	delete den;

}


Model::~Model()
{
	delete Tw;
	delete qw;
	delete delta;
	delete p;
	delete qStefan;
	delete qNorth;
	delete T;
	delete u;
	delete v;
	delete w;
	delete bcSouth;
}

void Model::solve()
{
	clock_t start = clock();

	Log log;
	fprintf(log.ptr, "iter\t\t\t\tU0\t\t\t\tr\t\t\t\tdelta\t\t\t\tp\t\t\t\tT\t\t\t\tu\t\t\t\tv\t\t\t\tw\t\t\t\tMFE\t\t\t\trelMFE\n");

	Plot plot;

	char * gnucmd = (char*)malloc(100);
	snprintf(gnucmd, 100, "plot '%s' using 1:10 w l\n", log.filename.c_str());

	int iter=0;

	int itsolve=0;

	while(iter < maxMainIter)
	{
		printf("Main Loop: %d\n", ++iter);


		find_U();

		find_r();

		if(std::fabs(Fscrew - F) >= FTol)
		{
			printf("Continuing!\n");
			continue;
		}

		init_uvw();

		if(iter == 1)
		{
			TSolveWrapper();
			itsolve++;

		}

		calc_maxFluxError();

		/* break based on absolute MFE */

		if(maxFluxError <= allowedMaxFluxError)
		{
			TSolveWrapper();
			itsolve++;
			calc_maxFluxError();
			if(maxFluxError <= allowedMaxFluxError)
			{
				printf("Breaking with MFE = %f\n", maxFluxError);
				break;
			}
		}

		/* if(relativeMFE <= allowedRelativeMFE) */
		/* { */
		/* 	printf("Breaking with MFE = %f\n", maxFluxError); */
		/* 	break; */
		/* } */


		adjustDelta();

		recalcDelta = 0;

		dumper();

		fprintf(log.ptr, "%4d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
				iter, U0, r, delta->average(), p->average(), T->average(),
				u->average(),v->average(), w->average(), maxFluxError, relativeMFE);

		plot.image(delta);
		/* fprintf(plot.gnu, gnucmd); */

	}

	if( iter >= maxMainIter )
	{
		printf("Reached max specified limit for main iterations\n");
	}

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);
	printf("solved temp %d times\n", itsolve);

	/* log.writeHeader(); */
	/* log.writeFooter(); */

}


double Model::xVal(int i)
{
	return (nx==1)? 0 : -Lx + 2*Lx/(nx-1)*i;
}

double Model::yVal(int j)
{
	return (ny==1)? 0 : -Ly + 2*Ly/(ny-1)*j;
}

double Model::zVal(int i, int j, int k)
{
	return (double)k/(nz-1)*delta->get(i, j, k);
}

double Model::UVal(int i, int j)
{
	double vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j);
	return U0*(1-vector/r);
}

void Model::init_fields()
{

	if(southBC == DIRICHLET)
	{
		Tw->set(TwExp);
	}
	else
	{
		Tw->setAll(Tm + 0.1);
	}

	qw->set(qwExp);

	r = 0.1;
	U0 = 1e-4;

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			delta->set(i,j, kL * (Tw->get(i,j) - Tm)/(rhoS * hmStar * UVal(i,j)));
		}

	bcSouth = (southBC == DIRICHLET) ? Tw->copy() : qw->copy()->divide(-kL);

}

void Model::update_fields()
{
	if(recalcDelta == 1)
	{
		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				delta->set(i,j, kL * (Tw->get(i,j) - Tm)/(rhoS * hmStar * UVal(i,j)));
			}
	}

	PSolveWrapper();
	/* F = p->integrateXY(); */

	Field * pv = new Field (nx, ny, 1, Lx, Ly);

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			pv->set(i, j, p->get(i,j,0)*( cos(radianTheta) * xVal(i)) + sin(radianTheta) * yVal(j) )  ;
		}

	Mtheta = pv->integrateXY();

	delete pv;


}

void Model::dumper()
{
	printf("\n<<<<<<<<<< STARTING DUMP >>>>>>>>>>\n");
	printf("r = %e\tM = %e\n", r, Mtheta);
	printf("U0 = %e\tF = %e\tFeff = %e\n", U0, F, Fscrew);
	printf("u = %e\tv=%e\nw=%e\n", u->average(), v->average(), w->average());
	printf("Temp = %e\n", T->average());
	printf("Pres = %e\n", p->average());
	printf("MFE = %e\n", maxFluxError);
	printf("stefan = %e\n", qStefan->average());
	printf("northflux= %e\n", qNorth->average());
	printf("delta = %e\n", delta->average());
	printf("<<<<<<<<<<<<< END OF DUMP >>>>>>>>>>>\n\n");

}

/*
 * prints input information to stdout
 */
void Model::print()
{
	//printf might be better for formatting after all
	using namespace std;

	cout << "Printing input variables\n";
	cout << "mu     = " << mu      << "\t\t" << "Lx                  = " << Lx                  << "\n"
		<< "Tm     = " << Tm      << "\t\t" << "Ly                  = " << Ly                  << "\n"
		<< "hm     = " << hm      << "\t\t" << "nx                  = " << nx                  << "\n"
		<< "cpL    = " << cpL     << "\t" << "ny                  = " << ny                  << "\n"
		<< "cpS    = " << cpS     << "\t" << "nz                  = " << nz                  << "\n"
		<< "rhoL   = " << rhoL    << "\t\t" << "MTol                = " << MTol                << "\n"
		<< "rhoS   = " << rhoS    << "\t\t" << "FTol                = " << FTol                << "\n"
		<< "kL     = " << kL      << "\t\t" << "allowedMaxFluxError = " << allowedMaxFluxError << "\n"
		<< "Fscrew = " << Fscrew  << "\t\t" << "allowedRelativeMFE  = " << allowedRelativeMFE  << "\n"
		<< "Tinf   = " << Tinf    << "\t\t" << "deltaCoeffMin       = " << deltaCoeffMin       << "\n"
		<< "\t\t\t" << "deltaCoeffMax       = " << deltaCoeffMax       << "\n"
		<< "\t\t\t" << "deltaRelax          = " << deltaRelax          << "\n";


	if(southBC == DIRICHLET)
		printf("southBC = DIRICHLET\n");
	else
		printf("southBC = NEUMANN\n");

	if(sidesBC == DIRICHLET)
		printf("sidesBC = DIRICHLET\n");
	else
		printf("sidesBC = NEUMANN\n");

	cout << "TwExp = " << TwExp << "\n";
	cout << "qwExp = " << qwExp << "\n";
	cout << "theta = " << theta << "\n";

}
