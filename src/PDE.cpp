#include "PDE.hpp"

PDE::PDE(Field & field)
{
	dx = field.dx;
	dy = field.dy;
	dz = field.dz;

	nx = field.nx;
	ny = field.ny;
	nz = field.nz;

	zero = std::make_unique<Field>(field);

	xx = zero.get();
	yy = zero.get();
	zz = zero.get();
	x = zero.get();
	y = zero.get();
	z = zero.get();
	xy = zero.get();
	yz = zero.get();
	zx = zero.get();
	rhs = zero.get();
}


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

