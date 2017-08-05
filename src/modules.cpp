#include "model.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include "PDE.hpp"
#include "solver.hpp"
#include "estimator.hpp"
#include <assert.h>

#include "operators.hpp"

void Model::find_U()
{
	int iter = 0;
	double tolerance = FTol;
	double U0new;

	std::vector<double> U0s;
	std::vector<double> Fs;

	while(iter < maxFindIter)
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

	if( iter >= maxFindIter)
	{
		printf("U0 loop did not converge in %d iterations\nExiting\n", maxFindIter);
		exit(-1);
	}

}

void Model::init_uvw()
{

	auto Dp_Dx     = p->differentiate(CX2, X);
	auto Dp_Dy     = p->differentiate(CX2, Y);
	auto D2p_Dx2   = p->differentiate(CXX2, X);
	auto D2p_Dy2   = p->differentiate(CXX2, Y);
	auto Ddelta_Dx = delta->differentiate(CX2, X);
	auto Ddelta_Dy = delta->differentiate(CX2, Y);
	auto lapl_p = (*D2p_Dx2) + (*D2p_Dy2);

	/* auto Dp_Dx     = p->differentiate(CX2, X)->replicateZ(nz); */
	/* auto Dp_Dy     = p->differentiate(CX2, Y)->replicateZ(nz); */
	/* auto D2p_Dx2   = p->differentiate(CXX2, X)->replicateZ(nz); */
	/* auto D2p_Dy2   = p->differentiate(CXX2, Y)->replicateZ(nz); */
	/* auto Ddelta_Dx = delta->differentiate(CX2, X)->replicateZ(nz); */
	/* auto Ddelta_Dy = delta->differentiate(CX2, Y)->replicateZ(nz); */
	/* auto lapl_p = ((*D2p_Dx2) + (*D2p_Dy2))->replicateZ(nz); */

	/* auto z3d = std::make_unique<Field>(*u); */
	/* z3d->set("z"); */
	/* z3d = (*z3d) * (*delta); */

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

		/* u = *(  *( ( *Dp_Dx ) * (*z3d) ) * *(*z3d - *(delta->replicateZ(nz)) ) ) / (2*mu); */
		/* u->print(); */
		/* std::cout << z3d->dz << "\n"; */
		/* assert(u->isFinite()); */

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


}

void Model::calc_maxFluxError()
{
	auto DT_Dz = T->differentiate(BX2, Z);
	qNorth = *( *(DT_Dz->getSubfield(0, nx-1, 0, ny-1, nz-1, nz-1)) / *delta ) * -1;

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

	std::vector<double> Mthetas;
	std::vector<double> rs;

	while(iter < maxFindIter)
	{
		iter++;
		update_fields();


		/* printf("r[%d]:\t%e\t%e\t%e\t%e\n", iter, r, Mx, My, Mtheta); */

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

	if( iter >= maxFindIter)
	{
		printf("r loop did not converge in %d iterations\nExiting\n", maxFindIter);
		exit(-1);
	}

}

void Model::TSolveWrapper()
{
	PDE TEqn(*T);

	//Boundary condition setting

	Field * TBC     = new Field(*T);
	Field * TBCFlag = new Field(*T);
	Field * zero    = new Field(*T);

	TBC->setSubfield(0,nx-1, 0,ny-1, 0,0, *bcSouth);
	TBCFlag->set(FRONT, (int)sidesBC);
	TBCFlag->set(BACK , (int)sidesBC);
	TBCFlag->set(LEFT , (int)sidesBC);
	TBCFlag->set(RIGHT, (int)sidesBC);
	TBCFlag->set(SOUTH, (int)southBC);

	//PDE set up
	auto alphaField = std::make_unique<Field>(T.get());
	alphaField->setAll(alpha);

	auto zeta = std::make_unique<Field>(T.get());
	for(int i=0; i<zeta->nx; i++)
		for(int j=0; j<zeta->ny; j++)
			for(int k=0; k<zeta->nz; k++)
			{
				zeta->set(i,j,k, (double)k/(zeta->nz-1));
			}

	auto delta_3D = delta->replicateZ(nz);
	auto Ddelta_Dx = delta_3D->differentiate(CX2, X);
	auto Ddelta_Dy = delta_3D->differentiate(CX2, Y);
	auto D2delta_Dx2 = delta_3D->differentiate(CXX2, X);
	auto D2delta_Dy2 = delta_3D->differentiate(CXX2, Y);

	auto sD2delta = *D2delta_Dx2 + (*D2delta_Dy2);
	auto ssDdelta = *(*Ddelta_Dx ^(2)) + *(*Ddelta_Dy ^ 2);

	auto z1 = *( *(*(*ssDdelta * 2) / *delta_3D) - *sD2delta ) * *( *(*alphaField * *zeta) / *delta_3D);
	auto z2 = *(*( *Ddelta_Dx * *u ) + *(*Ddelta_Dy * *v) ) * *(*zeta / *delta_3D);
	auto z3 = *(*w * -1) / *delta_3D;

	auto zz3 = *(*ssDdelta * *(*zeta^(2)) ) + 1;
	auto zz = *(*alphaField * *zz3) / *(*delta_3D ^ 2);

	auto x = *u * -1.0;
	auto y  = *v * -1.0;
	auto z  = *(*z1 + *z2) + *z3;


	auto term = *(*(*alphaField * *zeta) * (-2) ) / (*delta_3D);
	auto yz = (*term) * (*Ddelta_Dy);
	auto zx = (*term) * (*Ddelta_Dx);

	TEqn.xx  = alphaField.get();
	TEqn.yy  = alphaField.get();
	TEqn.zz  = zz.get();
	TEqn.x   = x.get();
	TEqn.y   = y.get();
	TEqn.z   = z.get();
	TEqn.yz  = yz.get();
	TEqn.zx  = zx.get();
	TEqn.rhs = zero;
	TEqn.xy  = zero;

	Solver TSolver(&TEqn, TBC, TBCFlag);
	TSolver.init();
	TSolver.solve(T.get());

	delete TBC;
	delete TBCFlag;
	delete zero;


}

void Model::PSolveWrapper()
{

	auto zero = std::make_unique<Field>(*p);

	PDE pEqn(*p);

	auto Ddelta_Dx  = delta->differentiate(CX2, X);
	auto Ddelta_Dy  = delta->differentiate(CX2, Y);
	auto lapl_coeff = *(*delta^3)/6;
	auto term       = *(*delta^2)/2;
	auto x_coeff    = *(Ddelta_Dx) * (*term);
	auto y_coeff    = *(Ddelta_Dy) * (*term);

	auto rhs = std::make_unique<Field>(*p);
	for(int i=0; i < nx; i++)
		for(int j=0; j<ny; j++)
		{
			rhs->set(i,j, -2*mu*rhoS/rhoL*UVal(i,j));
		}

	pEqn.xx  = lapl_coeff.get();
	pEqn.yy  = lapl_coeff.get();
	pEqn.x   = x_coeff.get();
	pEqn.y   = y_coeff.get();
	pEqn.rhs = rhs.get();

	pEqn.zz = zero.get();
	pEqn.z  = zero.get();
	pEqn.xy = zero.get();
	pEqn.yz = zero.get();
	pEqn.zx = zero.get();

	Solver pSolver(&pEqn, zero.get(), zero.get());
	pSolver.init();
	pSolver.solve(p.get());

}

