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
#include <stdlib.h>

#include "PDE.hpp"
#include "solver.hpp"
#include "estimator.hpp"
#include "log.hpp"
#include "plot.hpp"
#include "functions.hpp"

#include "packages/inih/INIReader.h"
#include "operators.hpp"

#include <sys/types.h>
#include <sys/stat.h>

/**
 * \brief Constructor for the Model class
 * \param[in] iniPath String with path to the ini file with inputs
 *
 * \details
 * The constructor first reads the input file. If not possible, directly exits the program.
 * Once the file is read, it extracts input variable values from the specified sections in the input file.
 * If the parser isn't able to find the input variable, a default value (hardcoded) is used in the program.
 *
 * After parsing the ini file and initializing input variables, it calculates the derived constants, and
 * initializes the essential Field members of the class.
 */
Model::Model(std::string iniPath)
{
	parseINI(iniPath);
}

void Model::init()
{
	if(!curvilinearMelting)
		r = 1.0/0.0;

	dx = 2.0*Lx/(nx-1);
	dy = 2.0*Ly/(ny-1);
	dz = 1.0/(nz - 1);

	hmStar = hm + cpS*(Tm-Tinf);
	alpha = kL/(rhoL*cpL);

	Tw      = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	qw      = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	p       = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	qStefan = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	qNorth  = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	T       = std::make_unique<Field>(nx, ny, nz, Lx, Ly);
	u       = std::make_unique<Field>(nx, ny, nz, Lx, Ly);
	v       = std::make_unique<Field>(nx, ny, nz, Lx, Ly);
	w       = std::make_unique<Field>(nx, ny, nz, Lx, Ly);

	/* q      = std::make_unique<Field>(nx, ny, 1, Lx, Ly); */

	variables["pi"] = M_PI;

	if(southBC == DIRICHLET)
		Tw->set(TwExp, variables);
	else
		Tw->setAll(Tm + 0.1);

	qw->set(qwExp, variables);

	bcSouth = (southBC == DIRICHLET) ? Tw->copy() : *qw / (-kL);

	auto num = bcSouth->differentiate(CX2, Y);
	auto den = bcSouth->differentiate(CX2, X);

	double arg = num->average() / den->average();
	radianTheta = (std::isnan(arg))? 0 : atan(arg);
	theta = radianTheta * 180 / M_PI;

	radianThetaField = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	vecField = std::make_unique<Field>(nx, ny, 1, Lx, Ly);

	auto argField = *num/ *den;
	double val;
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			val = argField->get(i,j);
			radianThetaField->set(i,j, std::isnan(val)? 0 : atan(val));
		}

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			val = cos(radianThetaField->get(i,j))*xVal(i) + sin(radianThetaField->get(i,j))*yVal(j);
			vecField->set(i,j, val );
 		}

	vec       = std::make_unique<Field>(nx, ny, 1, Lx, Ly);
	variables["radianTheta"] = radianTheta;
	vec->set("cos(radianTheta) * x + sin(radianTheta) * y", variables);


	U = *(*(*vec/(-r)) + 1.0) * U0;
	/* U = *(*(*vecField/(-r)) + 1.0) * U0; */
	delta = *( *(*Tw - Tm)/(*U) ) * (kL/(rhoS*hmStar));

}

/**
 * \brief Destructor for the Model class
 *
 * \details
 * Deletes the allocated Field members of the class
 */

Model::~Model()
{
}

/**
 *
 * \brief Uses the CCMSOLVE algorithm to solve the PCM melting problem.
 *
 * \details
 * The code solves for the final values of melting velocity (U0), curve radius (r),
 * melt film interface (delta) and other variables based on the Navier Stokes Eqn with the
 * appropriate constraints and boundary conditions.
 *
 * The program logs the progress of the variables(averaged fields) in a .dat file with a timestamped name.
 * The same file is read by GNUPLOT while displaying the realtime averaged plots.
 *
 * Note that the main loop will only run for maxMainIter times. Ensure that the value of the variable is
 * large enough for the given case.
 */
void Model::solve()
{
	clock_t start = clock();

	char * timestr = (char*)malloc(16);
	time_t timestamp = time(NULL);
	strftime(timestr, 16, "%Y%m%d_%H%M%S", localtime(&timestamp));
	assert(timestr);
	std::string solveID(timestr);
	free(timestr);

	std::string directory = "outputs/" + solveID;
	std::string avgDir = directory + "/avg";
	std::string fieldsDir = directory + "/fields";

	struct stat st = {0};
	if (stat(directory.c_str(), &st) == -1)
		mkdir(directory.c_str(), 0755);
	if (stat(avgDir.c_str(), &st) == -1)
		mkdir(avgDir.c_str(), 0755);
	if (stat(fieldsDir.c_str(), &st) == -1)
		mkdir(fieldsDir.c_str(), 0755);

	std::vector<std::string> logColumns = {"NULL", "iter", "U0", "r", "delta", "p", "T", "u", "v", "w", "qStefan", "qNorth","MFE", "RMFE", "RU0"};

	/* Log log; */
	Log log(directory);
	/* fprintf(log.logPTR, "iter\t\t\t\tU0\t\t\t\tr\t\t\t\tdelta\t\t\t\tp\t\t\t\tT\t\t\t\tu\t\t\t\tv\t\t\t\tw\t\t\t\tMFE\t\t\t\trelMFE\n"); */

	Plot plot;
	Plot plot2;
	Plot realTime;

	std::string gnucmd = "plot '" + log.filename +"' using 1:" + std::to_string(index_var) + " w l\n";

	int iter=0;
	int itsolve=0;


	/* double * allowedExitVarPtr = (exitVarFlag == 0)? &allowedMaxFluxError : &allowedRelativeMFE; */
	/* double * exitVarPtr = (exitVarFlag == 0)? &maxFluxError : &relativeMFE; */

	double * allowedExitVarPtr = (exitVarFlag == 0)? &allowedMaxFluxError : &allowedMaxFluxError;
	double * exitVarPtr = (exitVarFlag == 0)? &maxFluxError : &intFE;

	/* fprintf(plot.gnu, "set terminal gif animate delay 100\n"); */
	/* fprintf(plot.gnu, "set output 'test.gif'\n"); */

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


		//break loop
		if(*exitVarPtr<= *allowedExitVarPtr)
		{
			if(itsolve == maxTempSolve)
				break;
			TSolveWrapper();
			itsolve++;
			calc_maxFluxError();
			if(*exitVarPtr<= *allowedExitVarPtr)
			{
				printOutputs();
				printf("Breaking with MFE = %f\n", maxFluxError);
				break;
			}
		}

		/* double oldDeltaAvg = delta->average(); */
		adjustDelta();
		recalcDelta = 0;

		printOutputs();

		fprintf(log.logPTR, "%4d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
				iter, U0, r, delta->average(), p->average(), T->average(),
				u->average(),v->average(), w->average(), qStefan->average(), qNorth->average(), maxFluxError, relativeMFE, relativeU0);

		plot.image(qStefan.get());
		plot2.image(qNorth.get());
		/* plot.setSurface(); */
		/* plot.surface(delta.get()); */
		fprintf(realTime.gnu, gnucmd.c_str());

	}

	auto DT_Dz = T->differentiate(FX2, Z) ;
	q = DT_Dz->getSubfield(0,nx-1, 0,ny-1, 0,0) ;
	q = *(*q * (-kL) ) / *delta;
	Q = q->integrateXY();

	//render average plots
	fflush(log.logPTR);
	fprintf(plot.gnu, "unset output\n");
	fprintf(plot.gnu, "set terminal png\n");
	for(int i=2; i<=14; i++)
	{
		fprintf(plot.gnu, "set output '%s/plot%02d.png'\n", avgDir.c_str(), i);
		fprintf(plot.gnu, "set title '%s'\n", logColumns[i].c_str());
		fprintf(plot.gnu, "plot '%s' using 1:%d w l\n", log.filename.c_str(), i );
		fflush(plot.gnu);
	}

	updateMap();
	plot.saveImage(fieldMap, fieldsDir);
	/* plot.saveSurface(fieldMap, fieldsDir); */


	if( iter >= maxMainIter )
	{
		printf("Reached max specified limit for main iterations\n");
	}

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);
	printf("solved temp %d times\n", itsolve);

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

void Model::update_fields()
{
	if(recalcDelta == 1)
	{
		U = *(*(*vec/(-r)) + 1.0) * U0;
		/* U = *(*(*vecField/(-r)) + 1.0) * U0; */
		delta = *( *(*Tw - Tm)/(*U) ) * (kL/(rhoS*hmStar));
	}

	PSolveWrapper();

	auto pv = *p * *vec;
	/* auto pv = *p * *vecField; */

	Mtheta = pv->integrateXY();
	/* F = p->integrateXY(); */

}

void Model::printOutputs()
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
	printf("theta = %.2f\n", theta);
	printf("<<<<<<<<<<<<< END OF DUMP >>>>>>>>>>>\n\n");

}

void Model::printInputs()
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

void Model::parseINI(std::string iniPath)
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

	r = reader.GetReal("initialValues", "r", 0.1);
	U0 = reader.GetReal("initialValues", "U0", 1e-4);

	maxFindIter         = reader.GetInteger("parameters", "maxFindIter"     , 100);
	maxMainIter         = reader.GetInteger("parameters", "maxMainIter"     , 9999);
	maxTempSolve		= reader.GetInteger("parameters", "maxTempSolve"	, 10   );

	MTol                = reader.GetReal("parameters", "MTol"               , 1e-10);
	FTol                = reader.GetReal("parameters", "FTol"               , 1e-6 );
	allowedMaxFluxError = reader.GetReal("parameters", "allowedMaxFluxError", 100  );
	allowedRelativeMFE  = reader.GetReal("parameters", "allowedRelativeMFE" , 1e-6 );
	deltaCoeffMin       = reader.GetReal("parameters", "deltaCoeffMin"      , 0.7  );
	deltaCoeffMax       = reader.GetReal("parameters", "deltaCoeffMax"      , 1.3  );
	deltaRelax          = reader.GetReal("parameters", "deltaRelax"         , 0.05 );

	exitVarFlag         = reader.GetInteger("flags", "exitVarFlag"     , 0);
	index_var			= reader.GetInteger("flags", "index_var"     , 12);

	printOutputsToScreen = reader.GetBoolean("flags", "printOutputsToScreen", true);
	curvilinearMelting = reader.GetBoolean("flags", "curvilinearMelting", true);

}

void Model::updateMap()
{

	fieldMap = {
		{ "Tw", Tw.get()},
		{ "qw", qw.get()},
		{ "bcSouth", bcSouth.get()},
		{ "p", p.get()},
		{ "qStefan", qStefan.get()},
		{ "qNorth", qNorth.get()},
		{ "delta", delta.get()}
	};

}

void Model::solve2()
{
	clock_t start = clock();

	std::vector<std::string> dirs = setup_directory();

	std::vector<std::string> logColumns = {"NULL", "iter", "U0", "r", "delta",
		"p", "T", "u", "v", "w", "qStefan", "qNorth","MFE", "RMFE", "RU0", "intFR"};

	Log log(dirs[0]);

	Plot plot, plot2;
	Plot realTime;
	updateMap();
	plot.saveImage(fieldMap, dirs[2]);
	plot2.image(bcSouth.get());

	std::string gnucmd = "plot '" + log.filename +"' using 1:" + std::to_string(index_var) + " w l\n";

	int iter=0;
	int itsolve=0;

	double * allowedExitVarPtr = (exitVarFlag == 0)? &allowedMaxFluxError : &allowedMaxFluxError;
	double * exitVarPtr = (exitVarFlag == 0)? &maxFluxError : &avgFluxError;

	find_U();
	if(curvilinearMelting)
		find_r();
	init_uvw();
	recalcDelta=0;

	double oldTMFE = 0;
	double TMFE = maxFluxError;
	double relTMFE = 0;

	do
	{
		TSolveWrapper();
		itsolve++;
		oldTMFE = TMFE;
		calc_maxFluxError();
		TMFE = maxFluxError;
		relTMFE = fabs(TMFE - oldTMFE) / TMFE;

		while(*exitVarPtr > *allowedExitVarPtr)
		{
			++iter;
			find_U();

			if(curvilinearMelting)
				find_r();

			if(std::fabs(Fscrew - F) >= FTol)
			{
				printf("Continuing!\n");
				continue;
			}
			init_uvw();
			calc_maxFluxError();
			adjustDelta();

			if(printOutputsToScreen)
			{
				printf("Main Loop: %d\n", iter);
				printOutputs();
			}

			fprintf(log.logPTR, "%4d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
					iter, U0, r, delta->average(), p->average(), T->average(),
					u->average(),v->average(), w->average(), qStefan->average(),
					qNorth->average(), maxFluxError, relativeMFE, relativeU0, intFR);

			fprintf(realTime.gnu, gnucmd.c_str());

		}

		/* exit(0); */

	}
	while( (relTMFE > 1e-2) && (TMFE > allowedMaxFluxError) && (itsolve < maxTempSolve));

	/* auto f = *(delta->differentiate(CX2, X)) + *(delta->differentiate(CX2, Y)); */
	/* auto bcd = *(qNorth->differentiate(CX2, X)) + *(qNorth->differentiate(CX2, Y)); */
	/* auto g = *(*f / f->average()) - *(*bcd / bcd->average()); */
	/* plot2.image(g.get()); */
	/* g->print(); */

	fflush(log.logPTR);
	fprintf(plot.gnu, "unset output\n");
	fprintf(plot.gnu, "set terminal png\n");
	for(int i=2; i<=15; i++)
	{
		fprintf(plot.gnu, "set output '%s/%s.png'\n", dirs[1].c_str(), logColumns[i].c_str());
		fprintf(plot.gnu, "set title '%s'\n", logColumns[i].c_str());
		fprintf(plot.gnu, "plot '%s' using 1:%d w l\n", log.filename.c_str(), i );
		fflush(plot.gnu);
	}

	updateMap();
	plot.saveImage(fieldMap, dirs[2]);

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);
	printf("solved temp %d times\n", itsolve);
	printf("relTMFE = %e\n", relTMFE);

}

void Model::solve3()
{
	clock_t start = clock();

	std::vector<std::string> dirs = setup_directory();

	std::vector<std::string> logColumns = {"NULL", "iter", "U0", "r", "delta",
		"p", "T", "u", "v", "w", "qStefan", "qNorth","MFE", "RMFE", "RU0", "intFR"};

	Log log(dirs[0]);

	Plot plot, plot2;
	Plot realTime;
	updateMap();
	plot.saveImage(fieldMap, dirs[2]);

	std::string gnucmd = "plot '" + log.filename +"' using 1:" + std::to_string(index_var) + " w l\n";

	int iter=0;
	int itsolve=0;

	double * allowedExitVarPtr = (exitVarFlag == 0)? &allowedMaxFluxError : &allowedMaxFluxError;
	double * exitVarPtr = (exitVarFlag == 0)? &maxFluxError : &avgFluxError;

	find_U();
	if(curvilinearMelting)
		find_r();
	init_uvw();
	TSolveWrapper();
	recalcDelta=0;

	double oldTMFE = 0;
	double TMFE = maxFluxError;
	double relTMFE = 0;

	do
	{

		while(*exitVarPtr > *allowedExitVarPtr)
		{
			++iter;
			/* find_U(); */

			/* if(curvilinearMelting) */
			/* 	find_r(); */

			/* if(std::fabs(Fscrew - F) >= FTol) */
			/* { */
			/* 	printf("Continuing!\n"); */
			/* 	continue; */
			/* } */
			update_fields();
			F = p->integrateXY();
			init_uvw();
			calc_maxFluxError();
			adjustDelta();

			if(printOutputsToScreen)
			{
				printf("Main Loop: %d\n", iter);
				printOutputs();
			}

			fprintf(log.logPTR, "%4d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
					iter, U0, r, delta->average(), p->average(), T->average(),
					u->average(),v->average(), w->average(), qStefan->average(),
					qNorth->average(), maxFluxError, relativeMFE, relativeU0, intFR);


			fprintf(realTime.gnu, gnucmd.c_str());

		}

		find_U();

		if(curvilinearMelting)
			find_r();

		if(std::fabs(Fscrew - F) >= FTol)
		{
			printf("Continuing!\n");
			continue;
		}

		init_uvw();

		/* TSolveWrapper(); */
		/* itsolve++; */
		/* oldTMFE = TMFE; */
		calc_maxFluxError();
		/* TMFE = maxFluxError; */
		/* relTMFE = fabs(TMFE - oldTMFE) / TMFE; */
		/* exit(0); */


	}
	while( *exitVarPtr > *allowedExitVarPtr);

	/* auto f = *(delta->differentiate(CX2, X)) + *(delta->differentiate(CX2, Y)); */
	/* auto bcd = *(qNorth->differentiate(CX2, X)) + *(qNorth->differentiate(CX2, Y)); */
	/* auto g = *(*f / f->average()) - *(*bcd / bcd->average()); */
	/* plot2.image(g.get()); */
	/* g->print(); */

	fflush(log.logPTR);
	fprintf(plot.gnu, "unset output\n");
	fprintf(plot.gnu, "set terminal png\n");
	for(int i=2; i<=15; i++)
	{
		fprintf(plot.gnu, "set output '%s/%s.png'\n", dirs[1].c_str(), logColumns[i].c_str());
		fprintf(plot.gnu, "set title '%s'\n", logColumns[i].c_str());
		fprintf(plot.gnu, "plot '%s' using 1:%d w l\n", log.filename.c_str(), i );
		fflush(plot.gnu);
	}

	updateMap();
	plot.saveImage(fieldMap, dirs[2]);

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);
	printf("solved temp %d times\n", itsolve);
	printf("relTMFE = %e\n", relTMFE);

}
