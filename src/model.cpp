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
#include "PDE.hpp"
#include "solver.hpp"
#include "estimator.hpp"
#include <assert.h>
#include <time.h>
#include "packages/muparser/muParser.h"
#include "packages/inih/INIReader.h"

#include "log.hpp"

Model::Model()
{

	Tw      = new Field(nx, ny, 1, Lx, Ly, 0);
	qw      = new Field(nx, ny, 1, Lx, Ly, 0);
	delta   = new Field(nx, ny, 1, Lx, Ly, 0);
	p       = new Field(nx, ny, 1, Lx, Ly, 0);
	qStefan = new Field(nx, ny, 1, Lx, Ly, 0);

	T = new Field(nx, ny, nz, Lx, Ly, 0);
	u = new Field(nx, ny, nz, Lx, Ly, 0);
	v = new Field(nx, ny, nz, Lx, Ly, 0);
	w = new Field(nx, ny, nz, Lx, Ly, 0);

	init_fields();

}

/*
 * Constructor for model.
 * Parses input file for inputs and initializes fields.
 *
 */
Model::Model(std::string iniPath)
{

	/* INIReader reader("../examples/test.ini"); */
	INIReader reader(iniPath);

	if (reader.ParseError() < 0) {
		std::cout << "Can't load 'test.ini'\n";
		exit(-1);
	}

	mu     = reader.GetReal("constants", "mu", mu);
	Tm     = reader.GetReal("constants", "Tm", Tm);
	Tinf   = reader.GetReal("constants", "Tinf", Tinf);
	hm     = reader.GetReal("constants", "hm", hm);
	cpL    = reader.GetReal("constants", "cpL", cpL);
	cpS    = reader.GetReal("constants", "cpS" ,cpS);
	rhoL   = reader.GetReal("constants", "rhoL", rhoL);
	rhoS   = reader.GetReal("constants", "rhoS" ,rhoS);
	kL     = reader.GetReal("constants", "kL", kL);
	Lx     = reader.GetReal("constants", "Lx", Lx);
	Ly     = reader.GetReal("constants", "Ly", Ly);
	Fscrew = reader.GetReal("constants", "Fscrew", Fscrew);

	southBC = reader.GetBoundary("boundaryConditions", "southBC", southBC);
	sidesBC = reader.GetBoundary("boundaryConditions", "sidesBC", sidesBC);

	TwExp   = reader.Get("boundaryConditions", "TwExp", TwExp);
	qwExp   = reader.Get("boundaryConditions", "qwExp", qwExp);


	nx = reader.GetInteger("gridSizes", "nx", nx);
	ny = reader.GetInteger("gridSizes", "ny", ny);
	nz = reader.GetInteger("gridSizes", "nz", nz);

	/* theta   = reader.GetReal("boundaryConditions", "theta", theta); */



	MTol                = reader.GetReal("parameters", "MTol", MTol);
	FTol                = reader.GetReal("parameters", "FTol", FTol);
	allowedMaxFluxError = reader.GetReal("parameters", "allowedMaxFluxError", allowedMaxFluxError);
	allowedRelativeMFE  = reader.GetReal("parameters", "allowedRelativeMFE", allowedRelativeMFE);
	deltaCoeffMin       = reader.GetReal("parameters", "deltaCoeffMin", deltaCoeffMin);
	deltaCoeffMax       = reader.GetReal("parameters", "deltaCoeffMax", deltaCoeffMax);
	deltaRelax          = reader.GetReal("parameters", "deltaRelax", deltaRelax);

	dx = 2.0*Lx/(nx-1);
	dy = 2.0*Ly/(ny-1);
	dz = 1.0/(nz - 1);

	hmStar = hm + cpS*(Tm-Tinf);
	alpha = kL/(rhoL*cpL);
	/* radianTheta = theta*3.1415927/180; */

	Tw      = new Field(nx, ny, 1, Lx, Ly, 0);
	qw      = new Field(nx, ny, 1, Lx, Ly, 0);
	delta   = new Field(nx, ny, 1, Lx, Ly, 0);
	p       = new Field(nx, ny, 1, Lx, Ly, 0);
	qStefan = new Field(nx, ny, 1, Lx, Ly, 0);
	T       = new Field(nx, ny, nz, Lx, Ly, 0);
	u       = new Field(nx, ny, nz, Lx, Ly, 0);
	v       = new Field(nx, ny, nz, Lx, Ly, 0);
	w       = new Field(nx, ny, nz, Lx, Ly, 0);

	init_fields();

	double arg = bcSouth->differentiate(CX2, Y)->average() / bcSouth->differentiate(CX2, X)->average();
	radianTheta = (std::isnan(arg))? 0 : atan(arg);
	theta = radianTheta * 180 / M_PI;

	/* std::cout <<bcSouth->differentiate(CX2, Y)->average() << bcSouth->differentiate(CX2, X)->average() << std::endl; */
	/* std::cout << radianTheta << "\t" << theta << std::endl; */
	/* exit(-1); */
}


Model::~Model()
{
	delete Tw;
	delete qw;
	delete delta;
	delete p;
	delete qStefan;

	delete T;
	delete u;
	delete v;
	delete w;
}

void Model::legacySolve()
{
	clock_t start = clock();

	int iter=0;

	while(1)
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

		TSolveWrapper();

		calc_maxFluxError();

		/* break based on absolute MFE */

		if(maxFluxError <= allowedMaxFluxError)
		{
			printf("Breaking with MFE = %f\n", maxFluxError);
			break;
		}

		/* if(relativeMFE <= allowedRelativeMFE) */
		/* { */
		/* 	printf("Breaking with MFE = %f\n", maxFluxError); */
		/* 	break; */
		/* } */

		adjustDelta();

		recalcDelta = 0;

		dumper();
	}

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);

}

void Model::solve()
{
	clock_t start = clock();

	/* init_fields(); */
	int iter=0;

	/* log(); */
	Log log;
	log.writeHeader();
	fprintf(log.ptr, "U0\t\t\t\tr\t\t\t\tdelta\t\t\t\tp\t\t\t\tT\t\t\t\tu\t\t\t\tv\t\t\t\tw\t\t\t\tMFE\t\t\t\trelMFE\n");

	while(1)
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
			TSolveWrapper();

		calc_maxFluxError();

		/* break based on absolute MFE */

		if(maxFluxError <= allowedMaxFluxError)
		{
			TSolveWrapper();
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
		fprintf(log.ptr, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", U0, r, delta->average(), p->average(), T->average(), u->average(),
				v->average(), w->average(), maxFluxError, relativeMFE);

	}

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);

	log.writeFooter();

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
	double x;
	double y;

	mu::Parser p;
	p.DefineVar("x", &x);
	p.DefineVar("y", &y);
	p.DefineVar("Lx", &Lx);
	p.DefineVar("Ly", &Ly);

	if(southBC == DIRICHLET)
	{
		for(int i = 0; i < nx; i++)
			for(int j = 0; j < ny; j++)
			{

				x = xVal(i);
				y = yVal(j);
				p.SetExpr(TwExp);
				Tw->set(i,j,0, p.Eval());

				/* Tw->set(i, j, 0, 40+10*(xVal(i)/Lx + yVal(j)/Ly)); */
				/* Tw->set(i, j, 0, 40+10*(xVal(i)/Lx)); */

			}
	}
	else
	{
		Tw->setAll(0.1);
	}

	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
		{
			x = xVal(i);
			y = yVal(j);
			p.SetExpr(qwExp);
			qw->set(i, j, 0, p.Eval());
		}

	r = 0.1;
	U0 = 1e-4;
	/* delta->setAll(1e-2); */

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			delta->set(i,j, kL * (Tw->get(i,j) - Tm)/(rhoS * hmStar * UVal(i,j)));
		}

	bcSouth = (southBC == DIRICHLET) ? Tw : qw->copy()->divide(-kL);

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

	Field * pv = new Field (nx, ny, 1, Lx, Ly, 0);

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
	printf("r = %e\tMx = %e\t My = %e\n", r, Mx, My);
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

/* void Model::log() */
/* { */
/* 	char * timestr = (char*)malloc(16); */
/* 	char * logFileName = (char*)malloc(30); */

/* 	time_t timestamp; */
/* 	timestamp = time(NULL); */
/* 	strftime(timestr, 16, "%Y%m%d_%H%M%S", localtime(&timestamp)); */

/* 	assert(logFileName); */
/* 	snprintf(logFileName, 30, "logs/log_%s.dat", timestr); */

/* 	free(timestr); */

/* 	struct stat st = {0}; */
/* 	if (stat("logs", &st) == -1) */
/* 		mkdir("logs", 0755); */
/* 	if (stat("outputs", &st) == -1) */
/* 		mkdir("outputs", 0755); */

/* 	logFilePtr = fopen(logFileName, "w"); */
/* 	fprintf(logFilePtr, "##### BEGIN HEADER #####\n"); */
/* 	FILE *inputFile; */

/* 	inputFile = fopen( "inputs.ini", "r"); */
/* 	char ch; */
/* 	if(inputFile) */
/* 	{ */
/* 		while((ch = fgetc(inputFile)) != EOF) */
/* 		{ */
/* 			fputc(ch, logFilePtr); */
/* 		} */
/* 	} */
/* 	fclose(inputFile); */
/* 	fprintf(logFilePtr, "##### END HEADER #####\n\n"); */
/* 	free(logFileName); */
/* } */