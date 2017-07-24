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
#include "muparser/include/muParser.h"
#include "INIReader.h"

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
	}

	double duration = (double) (clock() - start)/CLOCKS_PER_SEC;
	printf("Time of Exec = %.2fs\n", duration);

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

/* void Model::find_U() */
/* { */
/* 	int iter = 0; */
/* 	double tolerance = FTol; */
/* 	double U0new; */

/* 	std::vector<double> U0s; */
/* 	std::vector<double> Fs; */

/* 	while(1) */
/* 	{ */
/* 		iter++; */
/* 		update_fields(); */
/* 		F = p->integrateXY(); */

/* 		/1* printf("U[%d]:\t%e\t%e\n", iter, U0, F); *1/ */

/* 		if(F<0) */
/* 		{ */
/* 			std::cout << "Negative Force\n"; */
/* 			printf("AVERAGE p = %e\n", p->average()); */
/* 			exit(-1); */
/* 		} */

/* 		U0s.push_back(U0); */
/* 		Fs.push_back(F); */

/* 		if(std::fabs(F-Fscrew) < tolerance) */
/* 			break; */
/* 		else */
/* 		{ */
/* 			if(iter == 1) */
/* 			{ */
/* 				if( F > Fscrew) */
/* 					U0 *= 0.5; */
/* 				else */
/* 					U0 *= 1.5; */
/* 			} */
/* 			else */
/* 			{ */
/* 				Estimator est; */
/* 				U0new = est.loglinear(Fs, U0s, Fscrew); */
/* 				if(!std::isfinite(U0new)) */
/* 					tolerance *= 10; */
/* 				else if(U0new > 0) */
/* 					U0 = U0new; */
/* 				else */
/* 					U0 *= 0.5; */

/* 			} */
/* 		} */

/* 	} */

/* } */

/* void Model::init_uvw() */
/* { */
/* 	Field * Dp_Dx     = p->differentiate(CX2, X); */
/* 	Field * Dp_Dy     = p->differentiate(CX2, Y); */
/* 	Field * D2p_Dx2   = p->differentiate(CXX2, X); */
/* 	Field * D2p_Dy2   = p->differentiate(CXX2, Y); */
/* 	Field * Ddelta_Dx = delta->differentiate(CX2, X); */
/* 	Field * Ddelta_Dy = delta->differentiate(CX2, Y); */
/* 	Field * lapl_p = D2p_Dx2->copy()->add(D2p_Dy2); */

/* 	double z; */
/* 	double value; */
/* 	double fo_term; */


/* 	if(nx!=1) */
/* 	{ */
/* 		for(int i = 0; i < nx; i++) */
/* 			for(int j = 0; j < ny; j++) */
/* 				for(int k = 0; k < nz; k++) */
/* 				{ */
/* 					z = (double)k * delta->get(i, j) / (nz-1); */
/* 					value = Dp_Dx->get(i, j) * z * (z - delta->get(i,j)) / (2*mu); */
/* 					u->set(i, j, k, value); */
/* 				} */
/* 		assert(u->isFinite()); */
/* 	} */

/* 	if(ny!=1) */
/* 	{ */
/* 		for(int i = 0; i < nx; i++) */
/* 			for(int j = 0; j < ny; j++) */
/* 				for(int k = 0; k < nz; k++) */
/* 				{ */
/* 					z =(double) k * delta->get(i, j)/(nz-1); */
/* 					value = Dp_Dy->get(i, j) * z * (z - delta->get(i,j)) / (2*mu); */
/* 					v->set(i, j, k, value); */
/* 				} */
/* 		assert(v->isFinite()); */

/* 	} */

/* 	for(int i = 0; i < nx; i++) */
/* 		for(int j = 0; j < ny; j++) */
/* 		{ */
/* 			fo_term = Dp_Dx->get(i,j) * Ddelta_Dx->get(i,j) + Dp_Dy->get(i,j) * Ddelta_Dy->get(i,j); */
/* 			for(int k = 0; k < nz; k++) */
/* 			{ */
/* 				z =(double) k * delta->get(i, j) / (nz-1); */
/* 				value = - (lapl_p->get(i, j) * (pow(z, 3) / 3 - delta->get(i,j) * pow(z, 2) / 2) - pow(z, 2)/2 * fo_term)  / (2*mu); */
/* 				w->set(i,j,k, value); */
/* 			} */

/* 		} */

/* 	delete Dp_Dx; */
/* 	delete Dp_Dy; */
/* 	delete D2p_Dx2; */
/* 	delete D2p_Dy2; */
/* 	delete Ddelta_Dx; */
/* 	delete Ddelta_Dy; */
/* 	delete lapl_p; */



/* } */

/* void Model::calc_maxFluxError() */
/* { */
/* 	Field * DT_Dz = T->differentiate(BX2, Z); */
/* 	qNorth = DT_Dz->getSubfield(0, nx-1, 0, ny-1, nz-1, nz-1)->divide(delta)->multiply(-1); */

/* 	double value; */

/* 	for(int i = 0; i < nx; i++) */
/* 		for( int j = 0; j < ny; j++) */
/* 		{ */
/* 			value = - hmStar * rhoL * w->get(i,j, nz-1)/kL; */
/* 			qStefan->set(i,j, value); */
/* 		} */

/* 	maxFluxError = 0; */
/* 	double fluxError; */

/* 	int sx = (ny==1)? 0 : 1; */
/* 	int ex = (nx==1)? 1 : nx-1; */
/* 	int sy = (ny==1)? 0 : 1; */
/* 	int ey = (ny==1)? 1 : ny-1; */

/* 	for(int i = sx; i < ex; i++) */
/* 		for( int j = sy; j < ey; j++) */
/* 		{ */
/* 			fluxError = std::fabs(qStefan->get(i,j) - qNorth->get(i,j)); */
/* 			/1* printf("%e\n", fluxError ); *1/ */
/* 			if(maxFluxError < fluxError) */
/* 				maxFluxError = fluxError; */
/* 		} */

/* 	relativeMFE = maxFluxError/qStefan->average(); */

/* 	/1* qStefan->print(); *1/ */
/* 	/1* qNorth->print(); *1/ */

/* 	delete DT_Dz; */

/* } */

/* void Model::adjustDelta() */
/* { */
/* 	double coeff, tempvar, value; */

/* 	for(int i = 0; i < nx; i++) */
/* 		for(int j = 0; j < ny; j++) */
/* 		{ */
/* 			coeff = qNorth->get(i,j) / qStefan->get(i,j); */
/* 			tempvar = coeff < deltaCoeffMax? coeff : deltaCoeffMax; */
/* 			coeff = tempvar > deltaCoeffMin? tempvar : deltaCoeffMin; */

/* 			value = (1 + deltaRelax * (coeff-1)) * delta->get(i,j); */
/* 			delta->set(i,j, value); */
/* 		} */
/* } */

/* void Model::find_r() */
/* { */
/* 	int iter = 0; */
/* 	double tolerance = MTol; */
/* 	double reciprocal; */

/* 	/1* std::vector<double> Mxs; *1/ */
/* 	/1* std::vector<double> Mys; *1/ */
/* 	std::vector<double> Mthetas; */
/* 	std::vector<double> rs; */

/* 	while(1) */
/* 	{ */
/* 		iter++; */
/* 		update_fields(); */
/* 		/1* F = p->integrateXY(); *1/ */


/* 		/1* printf("r[%d]:\t%e\t%e\t%e\t%e\n", iter, r, Mx, My, Mtheta); *1/ */

/* 		/1* Mxs.push_back(Mx); *1/ */
/* 		/1* Mys.push_back(My); *1/ */
/* 		Mthetas.push_back(Mtheta); */
/* 		rs.push_back(-1.0/r); */

/* 		if(std::fabs(Mtheta) < tolerance) */
/* 			break; */
/* 		else */
/* 		{ */
/* 			if(iter == 1) */
/* 			{ */
/* 				r = - Mtheta/fabs(Mtheta) * 1.1 * std::sqrt(Lx*Lx + Ly*Ly); */
/* 			} */
/* 			else */
/* 			{ */
/* 				Estimator est; */
/* 				reciprocal = est.linear(Mthetas, rs, 0.0); */

/* 				if(!std::isfinite(1/reciprocal) || !std::isfinite(reciprocal) ) */
/* 					tolerance *= 10; */
/* 				else */
/* 					r = -1.0/reciprocal; */

/* 			} */
/* 		} */

/* 	} */

/* } */

/* void Model::TSolveWrapper() */
/* { */
/* 	PDE TEqn(T); */

/* 	Field * TBC = new Field(nx, ny, nz, Lx, Ly, 0); */
/* 	TBC->setAll(0); */
/* 	Field * TBCFlag = new Field(nx, ny, nz, Lx, Ly, 0); */
/* 	TBCFlag->setAll(0); */
/* 	Field * zero = new Field(nx, ny, nz, Lx, Ly, 0); */
/* 	zero->setAll(0); */

/* 	for(int i=0; i<TBC->nx; i++) */
/* 		for(int j=0; j<TBC->ny; j++) */
/* 		{ */
/* 			TBC->set(i,j,0,bcSouth->get(i,j,0)); */
/* 		} */

/* 	for(int i=0; i<TBCFlag->nx; i++) */
/* 		for(int j=0; j<TBCFlag->ny; j++) */
/* 			for(int k=1; k<TBCFlag->nz-1; k++) */
/* 			{ */
/* 				if(TBCFlag->getSide(i,j,k) != INTERIOR) */
/* 					TBCFlag->set(i,j,k, (int)NEUMANN ); */
/* 			} */

/* 	if(southBC == NEUMANN) */
/* 	{ */
/* 		for(int i=0; i<TBCFlag->nx; i++) */
/* 			for(int j=0; j<TBCFlag->ny; j++) */
/* 				TBCFlag->set(i,j,0,(int)NEUMANN); */
/* 	} */

/* 	Field * alphaField = new Field(T); */
/* 	alphaField->setAll(alpha); */

/* 	Field * zeta = new Field(T); */
/* 	for(int i=0; i<zeta->nx; i++) */
/* 		for(int j=0; j<zeta->ny; j++) */
/* 			for(int k=0; k<zeta->nz; k++) */
/* 			{ */
/* 				zeta->set(i,j,k, (double)k/(zeta->nz-1)); */
/* 			} */

/* 	Field * delta_3D = delta->replicateZ(nz); */
/* 	Field * Ddelta_Dx = delta_3D->differentiate(CX2, X); */
/* 	Field * Ddelta_Dy = delta_3D->differentiate(CX2, Y); */
/* 	Field * temp = Ddelta_Dy->copy()->pow(2); */
/* 	Field * ssDdelta = Ddelta_Dx->copy()->pow(2)->add(temp); */
/* 	Field * temp2 = delta_3D->differentiate(CXX2, Y); */
/* 	Field * sD2delta = delta_3D->differentiate(CXX2, X)->add(temp2); */



/* 	Field * z4 = Ddelta_Dy->copy()->multiply(v); */
/* 	Field * z1 = ssDdelta->copy()->multiply(2)->divide(delta_3D)->subtract(sD2delta)->multiply(alphaField)->multiply(zeta)->divide(delta_3D); */
/* 	Field * z2 = Ddelta_Dx->copy()->multiply(u)->add(z4)->multiply(zeta)->divide(delta_3D); */
/* 	Field * z3 = w->copy()->multiply(-1)->divide(delta_3D); */
/* 	Field * term = alphaField->copy()->multiply(zeta)->multiply(-2)->divide(delta_3D); */

/* 	Field * zz1 = delta_3D->copy()->pow(2); */
/* 	Field * zz2 = zeta->copy()->pow(2); */
/* 	Field * zz3 = ssDdelta->copy()->multiply(zz2)->add(1.0); */

/* 	TEqn.xx = alphaField; */
/* 	TEqn.yy = alphaField; */
/* 	TEqn.zz = alphaField->copy()->multiply(zz3)->divide(zz1); */

/* 	TEqn.x  = u->copy()->multiply(-1); */
/* 	TEqn.y  = v->copy()->multiply(-1); */
/* 	TEqn.z  = z1->copy()->add(z2)->add(z3); */

/* 	TEqn.xy = zero; */
/* 	TEqn.yz = term->copy()->multiply(Ddelta_Dy); */
/* 	TEqn.zx = term->copy()->multiply(Ddelta_Dx); */

/* 	TEqn.rhs = zero; */

/* 	Solver TSolver(TEqn, TBC, TBCFlag); */
/* 	TSolver.init(); */
/* 	TSolver.solve(T); */

/* 	delete delta_3D; */
/* 	delete Ddelta_Dx; */
/* 	delete Ddelta_Dy; */
/* 	delete ssDdelta; */
/* 	delete sD2delta; */
/* 	delete z1; */
/* 	delete z2; */
/* 	delete z3; */
/* 	delete term; */
/* 	delete zeta; */
/* 	delete alphaField; */

/* 	delete TEqn.zz; */
/* 	delete TEqn.x; */
/* 	delete TEqn.y; */
/* 	delete TEqn.z; */
/* 	delete TEqn.yz; */
/* 	delete TEqn.zx; */

/* 	delete temp; */
/* 	delete temp2; */

/* 	delete zero; */
/* 	delete TBC; */
/* 	delete TBCFlag; */

/* } */

/* void Model::PSolveWrapper() */
/* { */
/* 	/1* Field pBC(nx, ny, 1, Lx, Ly, 0); *1/ */
/* 	/1* Field zero(nx, ny, 1, Lx, Ly, 0); *1/ */

/* 	Field * zero = new Field(nx, ny, 1, Lx, Ly, 0); */
/* 	zero->setAll(0); */

/* 	/1* Field * pBCFlag = new Field(nx, ny, 1, Lx, Ly, 0); *1/ */
/* 	/1* pBCFlag.setAll(0); *1/ */


/* 	PDE pEqn(p); */

/* 	Field * lapl_coeff = (delta->copy()->pow(3)->multiply(1.0/6.0)); */
/* 	Field * term = delta->copy()->pow(2)->multiply(0.5); */
/* 	Field * x_coeff = (delta->differentiate(CX2, X))->multiply(term); */
/* 	Field * y_coeff = (delta->differentiate(CX2, Y))->multiply(term); */


/* 	Field * rhs = new Field(p); */
/* 	for(int i=0; i < nx; i++) */
/* 		for(int j=0; j<ny; j++) */
/* 		{ */
/* 			rhs->set(i,j, -2*mu*rhoS/rhoL*UVal(i,j)); */
/* 		} */

/* 	pEqn.xx = lapl_coeff; */
/* 	pEqn.yy = lapl_coeff; */
/* 	pEqn.x  = x_coeff; */
/* 	pEqn.y  = y_coeff; */
/* 	pEqn.rhs = rhs; */

/* 	pEqn.zz = zero; */
/* 	pEqn.z  = zero; */
/* 	pEqn.xy = zero; */
/* 	pEqn.yz = zero; */
/* 	pEqn.zx = zero; */

/* 	Solver pSolver(pEqn, zero, zero); */
/* 	pSolver.init(); */
/* 	pSolver.solve(p); */

/* 	/1* p->print(); *1/ */
/* 	/1* exit(-1); *1/ */

/* 	delete lapl_coeff; */
/* 	delete term; */
/* 	delete x_coeff; */
/* 	delete y_coeff; */
/* 	delete rhs; */
/* 	delete zero; */

/* } */

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

