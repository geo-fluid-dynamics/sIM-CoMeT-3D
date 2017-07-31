#include "model.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include "PDE.hpp"
#include "solver.hpp"
#include "estimator.hpp"
#include <assert.h>

void Model::find_U()
{
	int iter = 0;
	double tolerance = FTol;
	double U0new;

	std::vector<double> U0s;
	std::vector<double> Fs;

	while(1)
	{
		iter++;
		update_fields();
		F = p->integrateXY();

		/* printf("U[%d]:\t%e\t%e\n", iter, U0, F); */

		if(F<0)
		{
			std::cout << "Negative Force\n";
			printf("AVERAGE p = %e\n", p->average());
			exit(-1);
		}

		U0s.push_back(U0);
		Fs.push_back(F);

		if(std::fabs(F-Fscrew) < tolerance)
			break;
		else
		{
			if(iter == 1)
			{
				if( F > Fscrew)
					U0 *= 0.5;
				else
					U0 *= 1.5;
			}
			else
			{
				Estimator est;
				U0new = est.loglinear(Fs, U0s, Fscrew);
				if(!std::isfinite(U0new))
					tolerance *= 10;
				else if(U0new > 0)
					U0 = U0new;
				else
					U0 *= 0.5;

			}
		}

	}

}

void Model::init_uvw()
{
	Field * Dp_Dx     = p->differentiate(CX2, X);
	Field * Dp_Dy     = p->differentiate(CX2, Y);
	Field * D2p_Dx2   = p->differentiate(CXX2, X);
	Field * D2p_Dy2   = p->differentiate(CXX2, Y);
	Field * Ddelta_Dx = delta->differentiate(CX2, X);
	Field * Ddelta_Dy = delta->differentiate(CX2, Y);
	Field * lapl_p = D2p_Dx2->copy()->add(D2p_Dy2);

	double z;
	double value;
	double fo_term;


	if(nx!=1)
	{
		for(int i = 0; i < nx; i++)
			for(int j = 0; j < ny; j++)
				for(int k = 0; k < nz; k++)
				{
					z = (double)k * delta->get(i, j) / (nz-1);
					value = Dp_Dx->get(i, j) * z * (z - delta->get(i,j)) / (2*mu);
					u->set(i, j, k, value);
				}
		assert(u->isFinite());
	}

	if(ny!=1)
	{
		for(int i = 0; i < nx; i++)
			for(int j = 0; j < ny; j++)
				for(int k = 0; k < nz; k++)
				{
					z =(double) k * delta->get(i, j)/(nz-1);
					value = Dp_Dy->get(i, j) * z * (z - delta->get(i,j)) / (2*mu);
					v->set(i, j, k, value);
				}
		assert(v->isFinite());

	}

	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
		{
			fo_term = Dp_Dx->get(i,j) * Ddelta_Dx->get(i,j) + Dp_Dy->get(i,j) * Ddelta_Dy->get(i,j);
			for(int k = 0; k < nz; k++)
			{
				z =(double) k * delta->get(i, j) / (nz-1);
				value = - (lapl_p->get(i, j) * (pow(z, 3) / 3 - delta->get(i,j) * pow(z, 2) / 2) - pow(z, 2)/2 * fo_term)  / (2*mu);
				w->set(i,j,k, value);
			}

		}

	delete Dp_Dx;
	delete Dp_Dy;
	delete D2p_Dx2;
	delete D2p_Dy2;
	delete Ddelta_Dx;
	delete Ddelta_Dy;
	delete lapl_p;



}

void Model::calc_maxFluxError()
{
	Field * DT_Dz = T->differentiate(BX2, Z);
	delete qNorth;
	qNorth = DT_Dz->getSubfield(0, nx-1, 0, ny-1, nz-1, nz-1)->divide(delta)->multiply(-1);

	double value;

	for(int i = 0; i < nx; i++)
		for( int j = 0; j < ny; j++)
		{
			value = - hmStar * rhoL * w->get(i,j, nz-1)/kL;
			qStefan->set(i,j, value);
		}

	maxFluxError = 0;
	double fluxError;

	int sx = (nx==1)? 0 : 1;
	int ex = (nx==1)? 1 : nx-1;
	int sy = (ny==1)? 0 : 1;
	int ey = (ny==1)? 1 : ny-1;

	for(int i = sx; i < ex; i++)
		for( int j = sy; j < ey; j++)
		{
			fluxError = std::fabs(qStefan->get(i,j) - qNorth->get(i,j));
			/* printf("%e\n", fluxError ); */
			if(maxFluxError < fluxError)
				maxFluxError = fluxError;
		}

	relativeMFE = maxFluxError/qStefan->average();



	delete DT_Dz;

}

void Model::adjustDelta()
{
	double coeff, tempvar, value;

	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
		{
			coeff = qNorth->get(i,j) / qStefan->get(i,j);
			tempvar = coeff < deltaCoeffMax? coeff : deltaCoeffMax;
			coeff = tempvar > deltaCoeffMin? tempvar : deltaCoeffMin;

			value = (1 + deltaRelax * (coeff-1)) * delta->get(i,j);
			delta->set(i,j, value);
		}
}

void Model::find_r()
{
	int iter = 0;
	double tolerance = MTol;
	double reciprocal;

	/* std::vector<double> Mxs; */
	/* std::vector<double> Mys; */
	std::vector<double> Mthetas;
	std::vector<double> rs;

	while(1)
	{
		iter++;
		update_fields();
		/* F = p->integrateXY(); */


		/* printf("r[%d]:\t%e\t%e\t%e\t%e\n", iter, r, Mx, My, Mtheta); */

		/* Mxs.push_back(Mx); */
		/* Mys.push_back(My); */
		Mthetas.push_back(Mtheta);
		rs.push_back(-1.0/r);

		if(std::fabs(Mtheta) < tolerance)
			break;
		else
		{
			if(iter == 1)
			{
				r = - Mtheta/fabs(Mtheta) * 1.1 * std::sqrt(Lx*Lx + Ly*Ly);
			}
			else
			{
				Estimator est;
				reciprocal = est.linear(Mthetas, rs, 0.0);

				if(!std::isfinite(1/reciprocal) || !std::isfinite(reciprocal) )
					tolerance *= 10;
				else
					r = -1.0/reciprocal;

			}
		}

	}

}

void Model::TSolveWrapper()
{
	PDE TEqn(T);

	//Boundary condition setting
	Field * TBC = new Field(nx, ny, nz, Lx, Ly, 0);
	Field * TBCFlag = new Field(nx, ny, nz, Lx, Ly, 0);
	Field * zero = new Field(nx, ny, nz, Lx, Ly, 0);

	TBC->setSubfield(0,nx-1, 0,ny-1, 0,0, bcSouth);
	TBCFlag->set(FRONT, (int)sidesBC);
	TBCFlag->set(BACK , (int)sidesBC);
	TBCFlag->set(LEFT , (int)sidesBC);
	TBCFlag->set(RIGHT, (int)sidesBC);
	TBCFlag->set(SOUTH, (int)southBC);

	//PDE set up
	Field * alphaField = new Field(T);
	alphaField->setAll(alpha);

	Field * zeta = new Field(T);
	for(int i=0; i<zeta->nx; i++)
		for(int j=0; j<zeta->ny; j++)
			for(int k=0; k<zeta->nz; k++)
			{
				zeta->set(i,j,k, (double)k/(zeta->nz-1));
			}

	Field * delta_3D = delta->replicateZ(nz);
	Field * Ddelta_Dx = delta_3D->differentiate(CX2, X);
	Field * Ddelta_Dy = delta_3D->differentiate(CX2, Y);
	Field * temp = Ddelta_Dy->copy()->pow(2);
	Field * ssDdelta = Ddelta_Dx->copy()->pow(2)->add(temp);
	/* Field * temp2 = delta_3D->differentiate(CXX2, Y); */
	/* Field * sD2delta = delta_3D->differentiate(CXX2, X)->add(temp2); */
	Field * temp2 = Ddelta_Dy->copy();
	Field * sD2delta = Ddelta_Dx->copy()->add(temp2);

	Field * z4 = Ddelta_Dy->copy()->multiply(v);
	Field * z1 = ssDdelta->copy()->multiply(2)->divide(delta_3D)->subtract(sD2delta)->multiply(alphaField)->multiply(zeta)->divide(delta_3D);
	Field * z2 = Ddelta_Dx->copy()->multiply(u)->add(z4)->multiply(zeta)->divide(delta_3D);
	Field * z3 = w->copy()->multiply(-1)->divide(delta_3D);
	Field * term = alphaField->copy()->multiply(zeta)->multiply(-2)->divide(delta_3D);

	Field * zz1 = delta_3D->copy()->pow(2);
	Field * zz2 = zeta->copy()->pow(2);
	Field * zz3 = ssDdelta->copy()->multiply(zz2)->add(1.0);

	TEqn.xx = alphaField;
	TEqn.yy = alphaField;
	TEqn.zz = alphaField->copy()->multiply(zz3)->divide(zz1);

	TEqn.x  = u->copy()->multiply(-1);
	TEqn.y  = v->copy()->multiply(-1);
	TEqn.z  = z1->copy()->add(z2)->add(z3);

	TEqn.xy = zero;
	TEqn.yz = term->copy()->multiply(Ddelta_Dy);
	TEqn.zx = term->copy()->multiply(Ddelta_Dx);

	TEqn.rhs = zero;

	Solver TSolver(&TEqn, TBC, TBCFlag);
	TSolver.init();
	TSolver.solve(T);

	delete delta_3D;
	delete Ddelta_Dx;
	delete Ddelta_Dy;
	delete ssDdelta;
	delete sD2delta;
	delete z1;
	delete z2;
	delete z3;
	delete z4;
	delete term;
	delete zeta;
	delete alphaField;

	delete zz1;
	delete zz2;
	delete zz3;

	delete TEqn.zz;
	delete TEqn.x;
	delete TEqn.y;
	delete TEqn.z;
	delete TEqn.yz;
	delete TEqn.zx;

	delete temp;
	delete temp2;

	delete zero;
	delete TBC;
	delete TBCFlag;

}

void Model::PSolveWrapper()
{
	/* Field pBC(nx, ny, 1, Lx, Ly, 0); */
	/* Field zero(nx, ny, 1, Lx, Ly, 0); */

	Field * zero = new Field(nx, ny, 1, Lx, Ly, 0);
	zero->setAll(0);

	/* Field * pBCFlag = new Field(nx, ny, 1, Lx, Ly, 0); */
	/* pBCFlag.setAll(0); */


	PDE pEqn(p);

	Field * lapl_coeff = delta->copy()->pow(3)->multiply(1.0/6.0);
	Field * term = delta->copy()->pow(2)->multiply(0.5);
	Field * Ddelta_Dx = delta->differentiate(CX2, X);
	Field * Ddelta_Dy = delta->differentiate(CX2, Y);
	/* Field * x_coeff = Ddelta_Dx->copy()->multiply(term); */
	/* Field * y_coeff = Ddelta_Dy->copy()->multiply(term); */
	Field * x_coeff = delta->differentiate(CX2, X)->multiply(term);
	Field * y_coeff = delta->differentiate(CX2, Y)->multiply(term);


	Field * rhs = new Field(p);
	for(int i=0; i < nx; i++)
		for(int j=0; j<ny; j++)
		{
			rhs->set(i,j, -2*mu*rhoS/rhoL*UVal(i,j));
		}

	pEqn.xx = lapl_coeff;
	pEqn.yy = lapl_coeff;
	pEqn.x  = x_coeff;
	pEqn.y  = y_coeff;
	pEqn.rhs = rhs;

	pEqn.zz = zero;
	pEqn.z  = zero;
	pEqn.xy = zero;
	pEqn.yz = zero;
	pEqn.zx = zero;

	Solver pSolver(&pEqn, zero, zero);
	pSolver.init();
	pSolver.solve(p);

	/* p->print(); */
	/* exit(-1); */

	delete lapl_coeff;
	delete term;
	delete x_coeff;
	delete y_coeff;
	delete Ddelta_Dx;
	delete Ddelta_Dy;
	delete rhs;
	delete zero;

}

void Model::combinedUpdate()
{
	Field * zero = new Field(nx, ny, 1, Lx, Ly, 0);
	zero->setAll(0);

	PDE pEqn(p);

	double rOld, U0Old;
	int iter =0;

	while(1)
	{
		printf("\n\nNR iter : %d\n", ++iter);

		if(recalcDelta == 1)
		{
			for(int i=0; i<nx; i++)
				for(int j=0; j<ny; j++)
				{
					delta->set(i,j, kL * (Tw->get(i,j) - Tm)/(rhoS * hmStar * UVal(i,j)));
				}
		}

		std::cout << r << "\t" << U0 << "\n";

		Field * lapl_coeff = delta->copy()->pow(3)->divide(6.0);
		Field * term = delta->copy()->pow(2)->multiply(0.5);
		Field * Ddelta_Dx = delta->differentiate(CX2, X);
		Field * Ddelta_Dy = delta->differentiate(CX2, Y);
		Field * x_coeff = Ddelta_Dx->multiply(term);
		Field * y_coeff = Ddelta_Dy->multiply(term);

		Field * rhs = new Field(p);
		for(int i=0; i < nx; i++)
			for(int j=0; j<ny; j++)
			{
				rhs->set(i,j, -2*mu*rhoS/rhoL*UVal(i,j));
			}

		pEqn.xx = lapl_coeff;
		pEqn.yy = lapl_coeff;
		pEqn.x  = x_coeff;
		pEqn.y  = y_coeff;
		pEqn.rhs = rhs;

		pEqn.zz = zero;
		pEqn.z  = zero;
		pEqn.xy = zero;
		pEqn.yz = zero;
		pEqn.zx = zero;


		Solver pSolver(&pEqn, zero, zero);
		pSolver.init();
		pSolver.solve(p);

		/* std::cout << "p :\t" << p->average() << "\n"; */

		Field * D2p_Dx2 = p->differentiate(CXX2, X);
		Field * D2p_Dy2 = p->differentiate(CXX2, Y);
		Field * Dp_Dx = p->differentiate(CX2, X);
		Field * Dp_Dy = p->differentiate(CX2, Y);

		Field * term1 = delta->differentiate(CX2, Y)->multiply(Dp_Dy);
		Field * term2 = delta->differentiate(CX2, X)->multiply(Dp_Dx)->add(term1)->multiply(delta);
		Field * rhsDdelta = D2p_Dx2->copy()->add(D2p_Dy2)->multiply(term)->add(term2)->multiply(-1);
		/* Field * newterm2 = Dp_Dy->copy()->divide(dy); */
		/* Field * newterm = Dp_Dx->copy()->divide(dx)->add(newterm2)->multiply(term); */
		/* Field * rhsDdelta = D2p_Dx2->copy()->add(D2p_Dy2)->multiply(term)->add(term2)->add(newterm)->multiply(-1); */


		Field * Dp_Ddelta = new Field(p);
		Field * Dp_DU0 = new Field(p);
		Field * Dp_Dr = new Field(p);

		pEqn.rhs = rhsDdelta;
		pSolver.buildb();
		pSolver.solve(Dp_Ddelta);

		/* std::cout << "rhs :\t" << rhs->average() << "\n"; */
		/* std::cout << "rhsDdelta :\t" << rhsDdelta->average() << "\n"; */
		/* std::cout << "Dp_Ddelta :\t" << Dp_Ddelta->average() << "\n"; */

		Field * rhsDU0 = rhs->copy()->divide(U0);
		pEqn.rhs = rhsDU0;
		pSolver.buildb();
		pSolver.solve(Dp_DU0);

		/* std::cout << "Dp_DU0 :\t" << Dp_DU0->average() << "\n"; */

		Field * rhsDr = new Field(p);
		double vector;
		for(int i=0; i < nx; i++)
			for(int j=0; j<ny; j++)
			{
				vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j);
				rhsDr->set(i,j, -2*mu*rhoS/rhoL*U0*(vector/(r*r)));
			}

		/* printf("rhsDr = %e\n", rhsDr->average()); */

		pEqn.rhs = rhsDr;
		pSolver.buildb();
		pSolver.solve(Dp_Dr);

		/* std::cout << "Dp_Dr :\t" << Dp_Dr->average() << "\n"; */

		Field * Ddelta_DU0 = delta->copy()->divide(-U0);
		Field * Ddelta_Dr = delta->copy();
		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j);
				Ddelta_Dr->set(i,j, Ddelta_Dr->get(i,j) * (-vector)/(r * (r - vector)));
			}


		Field * temp = Dp_Ddelta->copy()->multiply(Ddelta_DU0);
		Dp_DU0->add(temp);
		delete temp;
		temp = Dp_Ddelta->copy()->multiply(Ddelta_Dr);
		Dp_Dr->add(temp);
		delete temp;

		F = p->integrateXY();

		Field * pv = new Field (nx, ny, 1, Lx, Ly, 0);
		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				pv->set(i, j, p->get(i,j,0)*( cos(radianTheta) * xVal(i)) + sin(radianTheta) * yVal(j) )  ;
			}
		Mtheta = pv->integrateXY();

		double DF_DU0 = Dp_DU0->integrateXY();
		double DF_Dr = Dp_Dr->integrateXY();

		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				pv->set(i, j, Dp_DU0->get(i,j,0)*( cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j) ))  ;
			}
		double DM_DU0 = pv->integrateXY();

		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				pv->set(i, j, Dp_Dr->get(i,j,0)*( cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j) ))  ;
			}
		double DM_Dr = pv->integrateXY();

		std::cout << "===== Jac =====\n" << DM_Dr << "\t" << DM_DU0 << "\n" << DF_Dr << "\t" << DF_DU0 << "\n";

		double onebydet = 1.0/(DM_Dr * DF_DU0 - DM_DU0 * DF_Dr);
		std::cout << "onebydet = " << onebydet << "\n";

		double rhat = onebydet * (DF_DU0 * (-Mtheta) - DM_DU0 * (-(F-Fscrew)) );
		double U0hat = onebydet * (-DF_Dr * (-Mtheta) + DM_Dr * (-(F-Fscrew)));


		rOld = r;
		r = r + rhat;
		U0Old = U0;
		U0 = U0 + U0hat;

		std::cout << "r = " << r << "\t" << "U0 = " << U0 << "\n";

		std::cout << Dp_Ddelta->average() << "\t" << Dp_DU0->average() << "\t" << Dp_Dr->average() << "\n";
		std::cout << rhsDdelta->average() << "\t" << rhsDU0->average() << "\t" << rhsDr->average() << "\n";

		if( (fabs(U0-U0Old) < 1e-6) && (fabs(r - rOld) < 1e-6) )
			break;



		delete pv;
		delete D2p_Dx2;
		delete D2p_Dy2;
		delete Dp_Dx;
		delete Dp_Dy;
		delete term1;
		delete term2;
		delete rhsDdelta;
		delete rhsDU0;
		delete rhsDr;

		delete Dp_Ddelta;
		delete Dp_DU0;
		delete Dp_Dr;
		delete Ddelta_DU0;
		delete Ddelta_Dr;

		delete lapl_coeff;
		delete term;
		delete x_coeff;
		delete y_coeff;
		delete rhs;

		exit(-1);


	}


	delete zero;



}

void Model::combinedUpdate2()
{
	Field * zero = new Field(nx, ny, 1, Lx, Ly, 0);
	zero->setAll(0);

	PDE pEqn(p);

	double rOld, U0Old;
	int iter =0;

	r = 0.1;

	while(1)
	{
		printf("\n\nNR iter : %d\n", ++iter);

		if(recalcDelta == 1)
		{
			for(int i=0; i<nx; i++)
				for(int j=0; j<ny; j++)
				{
					delta->set(i,j, kL * (Tw->get(i,j) - Tm)/(rhoS * hmStar * UVal(i,j)));
				}
		}



		double vector;
		double value;

		/* double dr = 1e-6; */
		/* double dU0 = 1e-6; */
		/* Field * Ddelta_DU0 = new Field(delta); */
		/* Field * Ddelta_Dr = new Field(delta); */
		/* for(int i=0; i<nx; i++) */
		/* 	for(int j=0; j<ny; j++) */
		/* 	{ */
		/* 		vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j); */
		/* 		value = (kL*(Tw->get(i,j) - Tm)/(rhoS * hmStar * (U0+dU0) * (1-vector/r) ) - kL*(Tw->get(i,j) - Tm)/(rhoS * hmStar * (U0-dU0) * (1-vector/r) ) ) / (2*dU0); */
		/* 		Ddelta_DU0->set(i,j, value ); */
		/* 	} */

		/* for(int i=0; i<nx; i++) */
		/* 	for(int j=0; j<ny; j++) */
		/* 	{ */
		/* 		vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j); */
		/* 		value =( kL*(Tw->get(i,j) - Tm)/(rhoS * hmStar * (U0) * (1-vector/(r+dr) )) - kL*(Tw->get(i,j) - Tm)/(rhoS * hmStar * (U0) * (1-vector/(r-dr) )) )/ (2*dr); */
		/* 		Ddelta_DU0->set(i,j, value ); */
		/* 	} */

		std::cout << r << "\t" << U0 << "\n";

		Field * lapl_coeff = delta->copy()->pow(3)->divide(6.0);
		Field * term = delta->copy()->pow(2)->multiply(0.5);
		Field * x_coeff = (delta->differentiate(CX2, X))->multiply(term);
		Field * y_coeff = (delta->differentiate(CX2, Y))->multiply(term);

		Field * rhs = new Field(p);
		for(int i=0; i < nx; i++)
			for(int j=0; j<ny; j++)
			{
				rhs->set(i,j, -2*mu*rhoS/rhoL*UVal(i,j));
			}

		pEqn.xx = lapl_coeff;
		pEqn.yy = lapl_coeff;
		pEqn.x  = x_coeff;
		pEqn.y  = y_coeff;
		pEqn.rhs = rhs;

		pEqn.zz = zero;
		pEqn.z  = zero;
		pEqn.xy = zero;
		pEqn.yz = zero;
		pEqn.zx = zero;


		Solver pSolver(&pEqn, zero, zero);
		pSolver.init();
		pSolver.solve(p);

		Field * Ddelta_DU0 = delta->copy()->divide(-U0);
		Field * Ddelta_Dr = delta->copy();
		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j);
				Ddelta_Dr->set(i,j, Ddelta_Dr->get(i,j) * (-vector)/(r * (r - vector)));
			}

		Field * Dp_Ddelta = new Field(p);
		Field * Dp_DU0 = new Field(p);
		Field * Dp_Dr = new Field(p);

		Field * D2p_Dx2 = p->differentiate(CXX2, X);
		Field * D2p_Dy2 = p->differentiate(CXX2, Y);
		Field * Dp_Dx = p->differentiate(CX2, X);
		Field * Dp_Dy = p->differentiate(CX2, Y);
		Field * lapl_p = D2p_Dx2->copy()->add(D2p_Dy2);
		Field * Ddelta_Dx = delta->differentiate(CX2, X);
		Field * Ddelta_Dy = delta->differentiate(CX2, Y);

		Field * temp = Dp_Dy->copy()->multiply(Ddelta_Dy);
		Field * cross_term_delta = Dp_Dx->copy()->multiply(Ddelta_Dx)->add(temp)->multiply(delta);
		Field * cross_term_delta_DU0 = cross_term_delta->copy()->multiply(Ddelta_DU0);

		delete temp;
		Field * Ddelta_DU0_Dx = Ddelta_DU0->differentiate(CX2, X);
		Field * Ddelta_DU0_Dy = Ddelta_DU0->differentiate(CX2, Y);
		temp = Ddelta_DU0_Dx->multiply(Dp_Dx)->add(Ddelta_DU0_Dy->multiply(Dp_Dy))->multiply(term);
		Field * add_DU0 = lapl_p->copy()->multiply(Ddelta_DU0)->multiply(term)->add(cross_term_delta_DU0)->add(temp);
		delete temp;

		Field * rhsDU0 = rhs->copy()->divide(U0)->subtract(add_DU0);
		pEqn.rhs = rhsDU0;
		pSolver.buildb();
		pSolver.solve(Dp_DU0);

		Field * rhsDr = new Field(p);
		for(int i=0; i < nx; i++)
			for(int j=0; j<ny; j++)
			{
				vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j);
				rhsDr->set(i,j, -2*mu*rhoS/rhoL*U0*(vector/(r*r)));
			}

		Field * Ddelta_Dr_Dx = Ddelta_Dr->differentiate(CX2, X);
		Field * Ddelta_Dr_Dy = Ddelta_Dr->differentiate(CX2, Y);
		temp = Ddelta_Dr_Dx->multiply(Dp_Dx)->add(Ddelta_Dr_Dy->multiply(Dp_Dy))->multiply(term);
		Field * cross_term_delta_Dr = cross_term_delta->copy()->multiply(Ddelta_Dr);
		Field * add_Dr = lapl_p->copy()->multiply(Ddelta_Dr)->multiply(term)->add(cross_term_delta_Dr)->add(temp);
		delete temp;

		rhsDr->subtract(add_Dr);
		pEqn.rhs = rhsDr;
		pSolver.buildb();
		pSolver.solve(Dp_Dr);

		F = p->integrateXY();

		Field * pv = new Field (nx, ny, 1, Lx, Ly, 0);
		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				pv->set(i, j, p->get(i,j,0)*( cos(radianTheta) * xVal(i)) + sin(radianTheta) * yVal(j) )  ;
			}
		Mtheta = pv->integrateXY();

		double DF_DU0 = Dp_DU0->integrateXY();
		double DF_Dr = Dp_Dr->integrateXY();

		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				pv->set(i, j, Dp_DU0->get(i,j,0)*( cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j) ))  ;
			}
		double DM_DU0 = pv->integrateXY();

		for(int i=0; i<nx; i++)
			for(int j=0; j<ny; j++)
			{
				pv->set(i, j, Dp_Dr->get(i,j,0)*( cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j) ))  ;
			}
		double DM_Dr = pv->integrateXY();

		std::cout << "===== Jac =====\n" << DM_Dr << "\t" << DM_DU0 << "\n" << DF_Dr << "\t" << DF_DU0 << "\n";

		double onebydet = 1.0/(DM_Dr * DF_DU0 - DM_DU0 * DF_Dr);
		std::cout << "onebydet = " << onebydet << "\n";

		double rhat = onebydet * (DF_DU0 * (-Mtheta) - DM_DU0 * (-(F-Fscrew)) );
		double U0hat = onebydet * (-DF_Dr * (-Mtheta) + DM_Dr * (-(F-Fscrew)));


		rOld = r;
		r = r + rhat;
		U0Old = U0;
		U0 = U0 + U0hat;

		std::cout << "r = " << r << "\t" << "U0 = " << U0 << "\n";

		if( (fabs(U0-U0Old) < 1e-6) && (fabs(r - rOld) < 1e-6) )
			break;

		delete pv;
		delete D2p_Dx2;
		delete D2p_Dy2;
		delete Dp_Dx;
		delete Dp_Dy;
		delete rhsDU0;
		delete rhsDr;

		delete lapl_p;
		delete cross_term_delta;
		delete cross_term_delta_DU0;
		delete cross_term_delta_Dr;
		delete add_DU0;
		delete add_Dr;
		delete Ddelta_DU0_Dx;
		delete Ddelta_DU0_Dy;

		delete Dp_Ddelta;
		delete Dp_DU0;
		delete Dp_Dr;
		delete Ddelta_DU0;
		delete Ddelta_Dr;

		delete lapl_coeff;
		delete term;
		delete x_coeff;
		delete y_coeff;
		delete rhs;

		/* exit(-1); */


	}


	delete zero;



}
