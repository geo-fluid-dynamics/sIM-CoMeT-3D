#ifndef PDE_H
#define PDE_H
#include "field.hpp"

class PDE {
	public:
		Field *xx;
		Field *yy;
		Field *zz;
		Field *xy;
		Field *yz;
		Field *zx;
		Field *x;
		Field *y;
		Field *z;
		Field *rhs;

		double nx;
		double ny;
		double nz;

		double dx;
		double dy;
		double dz;

		PDE(Field * field);
		PDE();
		PDE(int nx, int ny, int nz, double dx, double dy, double dz);
		bool isFinite();
};
#endif
