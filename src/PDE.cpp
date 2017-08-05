#include "PDE.hpp"

PDE::PDE(Field & field)
{
	dx = field.dx;
	dy = field.dy;
	dz = field.dz;

	nx = field.nx;
	ny = field.ny;
	nz = field.nz;
}

/* PDE::PDE() */
/* { */
/* 	dx = 0; */
/* 	dy = 0; */
/* 	dz = 0; */

/* 	nx = 0; */
/* 	ny = 0; */
/* 	nz = 0; */
/* } */

/* PDE::PDE(int nx, int ny, int nz, double dx, double dy, double dz) */
/* { */
/* 	this->nx = nx; */
/* 	this->ny = ny; */
/* 	this->nz = nz; */
/* 	this->dx = dx; */
/* 	this->dy = dy; */
/* 	this->dz = dz; */
/* } */

bool PDE::isFinite()
{
	if(!this->xx->isFinite())
		return false;
	if(!this->yy->isFinite())
		return false;
	if(!this->zz->isFinite())
		return false;
	if(!this->xy->isFinite())
		return false;
	if(!this->yz->isFinite())
		return false;
	if(!this->zx->isFinite())
		return false;
	if(!this->x->isFinite())
		return false;
	if(!this->y->isFinite())
		return false;
	if(!this->z->isFinite())
		return false;

	return true;
}

