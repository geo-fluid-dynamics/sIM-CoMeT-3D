#include "solver.hpp"
#include "stencil.hpp"
#include <assert.h>
#include <iostream>

#if __has_include(<umfpack.h>)
#include <umfpack.h>
#elif __has_include(<suitesparse/umfpack.h>)
#include <suitesparse/umfpack.h>
#endif

Solver::Solver(PDE * ipde, Field * bc, Field * bcflag)
{
	pde = ipde;
	field = bc;
	BCFlag = bcflag;

	lda = field->nx * field->ny * field->nz;

	A = new Sparse();
	x = new std::vector<double>(lda, 0);
	b = new std::vector<double>(lda, 0);

}

Solver::~Solver()
{
	delete A;
	delete x;
	delete b;
}

void Solver::init()
{

	Stencil stencil(field);
	/* int lda = field.nx * field.ny * field.nz; */
	int row;

	//interior
	int sx = (field->nx == 1)? 0 : 1;
	int sy = (field->ny == 1)? 0 : 1;
	int sz = (field->nz == 1)? 0 : 1;
	int ex = (field->nx == 1)? 1 : field->nx -1 ;
	int ey = (field->ny == 1)? 1 : field->ny -1 ;
	int ez = (field->nz == 1)? 1 : field->nz -1 ;

	/* assert(pde.isFinite()); */

	for(int i=sx; i<ex; i++)
		for(int j=sy; j<ey; j++)
			for(int k=sz; k<ez; k++)
			{
				stencil.reset();

				stencil.discretize(pde->xx->get(i,j,k), CXX2, X);
				stencil.discretize(pde->yy->get(i,j,k), CXX2, Y);
				stencil.discretize(pde->xy->get(i,j,k), CXX2, XY);
				stencil.discretize(pde->x->get(i,j,k), CX2, X);
				stencil.discretize(pde->y->get(i,j,k), CX2, Y);
				stencil.discretize(pde->zz->get(i,j,k), CXX2, Z);
				stencil.discretize(pde->yz->get(i,j,k), CXX2, YZ);
				stencil.discretize(pde->zx->get(i,j,k), CXX2, ZX);
				stencil.discretize(pde->z->get(i,j,k), CX2, Z);

				assert(stencil.isFinite());

				/* row = i + field->nx * j + field->nx * field->ny * k; */
				row = field->index(i,j,k);
				A->pushStencil(stencil, row);

				(*b)[row] = pde->rhs->get(i,j,k);

			}




	//boundaries
	//set BC from field.
	for(int i=0; i<field->nx; i++)
		for(int j=0; j<field->ny; j++)
			for(int k=0; k<field->nz; k++)
			{
				Side pointSide = field->getSide(i,j,k);
				if(pointSide != INTERIOR)
				{
					/* int row = i + field->nx * j + field->nx * field->ny * k; */
					row = field->index(i,j,k);
					stencil.reset();

					//a bithack to shorten the code for boundary values in constructing A
					//involves a specifically numbered enum set and an edge case if statement in stencil.discretize()
					//aside from the bithack in the following line.
					stencil.discretize(1, static_cast<dMode>((pointSide%2+1) & (int)(BCFlag->get(i,j,k))), static_cast<dDir>(pointSide/2));
					/* printf("%d\n", pointSide/2); */

					A->pushStencil(stencil, row);
					(*b)[row] = field->get(i,j,k);
				}
			}

	assert(A->isFinite());

	/* for(int i=0; i<field->nx; i++) */
	/* 	for(int j=0; j<field->ny; j++) */
	/* 		for(int k=0; k<field->nz; k++) */
	/* 		{ */
	/* 			/1* f->set(i,j,k, x->at(i + f->nx*j + f->nx * f->ny * k)); *1/ */
	/* 			printf("%d\t%d\t%d\t%e\n", i, j, k, b->at(i + field->nx*j + field->nx * field->ny * k)); */
	/* 		} */


}

//check umfpackcpp ex.

void Solver::solve(Field * f)
{
	int size = A->values.size();

	/* for(std::vector<double>::iterator it=x->begin(); it!=x->end(); ++it) */
	/* { */
	/* 	*it = 0; */
	/* } */

	/* double avg=0; */
	/* for(std::vector<double>::iterator it=b->begin(); it!=b->end(); ++it) */
	/* { */
	/* 	avg += *it; */
	/* } */
	/* avg /= b->size(); */
	/* printf("\t\t\t\tb avg = %e\n", avg); */

	std::vector<int> Ap(lda+1);
	std::vector<int> Ai(size);
	std::vector<double> Av(size);

	int status;
	status = umfpack_di_triplet_to_col(lda, lda, size, A->Ai.data(), A->Aj.data(), A->values.data(), Ap.data(), Ai.data(), Av.data(), NULL);

	void *Symbolic, *Numeric;

	status = umfpack_di_symbolic(lda, lda, Ap.data(), Ai.data(), Av.data(), &Symbolic, NULL, NULL);
	assert(status == UMFPACK_OK);

	status = umfpack_di_numeric(Ap.data(), Ai.data(), Av.data(), Symbolic, &Numeric, NULL, NULL);
	assert(status == UMFPACK_OK);

	umfpack_di_free_symbolic(&Symbolic);

	status = umfpack_di_solve(UMFPACK_A, Ap.data(), Ai.data(), Av.data(), x->data(), b->data(), Numeric, NULL, NULL);
	assert(status == UMFPACK_OK);

	umfpack_di_free_numeric(&Numeric);

	for(int i=0; i<f->nx; i++)
		for(int j=0; j<f->ny; j++)
			for(int k=0; k<f->nz; k++)
			{
				/* f->set(i,j,k, x->at(i + f->nx*j + f->nx * f->ny * k)); */
				f->set(i,j,k, x->at(f->index(i,j,k)));
				/* printf("%d\t%d\t%d\t%e\n", i, j, k, x->at(i + f->nx*j + f->nx * f->ny * k)); */
			}


}

void Solver::setA(Sparse iA)
{

}

void Solver::buildb()
{
	int row;

	//interior

	int sx = (field->nx == 1)? 0 : 1;
	int sy = (field->ny == 1)? 0 : 1;
	int sz = (field->nz == 1)? 0 : 1;
	int ex = (field->nx == 1)? 1 : field->nx -1 ;
	int ey = (field->ny == 1)? 1 : field->ny -1 ;
	int ez = (field->nz == 1)? 1 : field->nz -1 ;

	for(int i=sx; i<ex; i++)
		for(int j=sy; j<ey; j++)
			for(int k=sz; k<ez; k++)
			{
				/* row = i + field->nx * j + field->nx * field->ny * k; */
				row = field->index(i,j,k);
				(*b)[row] = pde->rhs->get(i,j,k);
			}



	//boundaries
	//set BC from field.

	for(int i=0; i<field->nx; i++)
		for(int j=0; j<field->ny; j++)
			for(int k=0; k<field->nz; k++)
			{
				Side pointSide = field->getSide(i,j,k);
				if(pointSide != INTERIOR)
				{
					/* int row = i + field->nx * j + field->nx * field->ny * k; */
					row = field->index(i,j,k);
					(*b)[row] = field->get(i,j,k);
				}
			}


	/* for(int i=0; i<field->nx; i++) */
	/* 	for(int j=0; j<field->ny; j++) */
	/* 		for(int k=0; k<field->nz; k++) */
	/* 		{ */
	/* 			row = i + field->nx * j + field->nx * field->ny * k; */
	/* 			(*b)[row] = pde.rhs->get(i,j,k); */
	/* 		} */


}
