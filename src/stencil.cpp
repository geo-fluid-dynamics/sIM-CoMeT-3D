#include "stencil.hpp"
#include <vector>
#include <cmath>

Stencil::Stencil(double idx, double idy, double idz)
{
	dx = idx;
	dy = idy;
	dz = idz;

	this->reset();

}

Stencil::Stencil(Field * field)
{
	dx = field->dx;
	dy = field->dy;
	dz = field->dz;

	nx = field->nx;
	ny = field->ny;
	nz = field->nz;

	stride_x = field->index(1,0,0) - field->index(0,0,0);
	stride_y = field->index(0,1,0) - field->index(0,0,0);
	stride_z = field->index(0,0,1) - field->index(0,0,0);

	this->reset();
}

void Stencil::reset()
{
	for(int i = 0; i <5; i++)
		for(int j = 0; j <5; j++)
			for(int k = 0; k <5; k++)
			{
				_stencil[i][j][k] = 0;
			}
}

void Stencil::discretize(double coeff, dMode mode, dDir dir)
{
	if(mode == ONE)
	{
		_stencil[2][2][2] += 1;
		return;
	}

	std::vector<double> vec = coeffs[mode];
	double h = (mode == CXX2 || mode == FXX1 || mode == BXX1) ? 2 : 1;


	switch(dir)
	{
		case X:
			if(std::isfinite(dx))
				for(int i=0; i<5; i++)
					_stencil[i][2][2] += coeff*vec[i]/pow(dx, h);
			break;
		case Y:
			if(std::isfinite(dy))
				for(int j=0; j<5; j++)
					_stencil[2][j][2] += coeff*vec[j]/pow(dy, h);
			break;
		case Z:
			if(std::isfinite(dz))
				for(int k=0; k<5; k++)
					_stencil[2][2][k] += coeff*vec[k]/pow(dz, h);
			break;
		case XY:
			if(std::isfinite(dy) && std::isfinite(dx))
			{
				_stencil[3][3][2] += coeff/(4*dx*dy);
				_stencil[3][1][2] += coeff/(4*dx*dy);
				_stencil[1][3][2] -= coeff/(4*dx*dy);
				_stencil[1][1][2] -= coeff/(4*dx*dy);
			}
			break;
		case YZ:
			if(std::isfinite(dy) && std::isfinite(dz))
			{
				_stencil[2][3][3] += coeff/(4*dz*dy);
				_stencil[2][3][1] += coeff/(4*dz*dy);
				_stencil[2][1][3] -= coeff/(4*dz*dy);
				_stencil[2][1][1] -= coeff/(4*dz*dy);
			}
			break;
		case ZX:
			if(std::isfinite(dz) && std::isfinite(dx))
			{
				_stencil[3][2][3] += coeff/(4*dx*dz);
				_stencil[3][2][1] += coeff/(4*dx*dz);
				_stencil[1][2][3] -= coeff/(4*dx*dz);
				_stencil[1][2][1] -= coeff/(4*dx*dz);
			}
			break;

	}

}

bool Stencil::isFinite()
{
	for(int i = 0; i <5; i++)
		for(int j = 0; j <5; j++)
			for(int k = 0; k <5; k++)
			{
				if(!std::isfinite(_stencil[i][j][k]))
					return false;
			}

	return true;

}

void Stencil::print()
{
	for(int i = 0; i <5; i++)
		for(int j = 0; j <5; j++)
			for(int k = 0; k <5; k++)
			{
				printf("(%d\t%d\t%d):\t%e\n", i, j, k, _stencil[i][j][k]);
			}
}
