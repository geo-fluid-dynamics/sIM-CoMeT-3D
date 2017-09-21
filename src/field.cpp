#include "field.hpp"
#include <cmath>
#include "packages/lepton/Lepton.h"
#include <assert.h>

#include <iostream>

Field::Field(int inx, int iny, int inz, double iLx, double iLy)
{
	nx = inx;
	ny = iny;
	nz = inz;

	Lx = iLx;
	Ly = iLy;

	dx = 2*Lx/(nx-1);
	dy = 2*Ly/(ny-1);
	dz = 1.0/(nz-1);

	values.reserve(nx*ny*nz);
	this->setAll(0);
}

Field::Field()
{
	nx = 0;
	ny = 0;
	nz = 0;
	Lx = 0;
	Ly = 0;
	dx = 0;
	dy = 0;
	dz = 0;

	values.reserve(0);
}

Field::Field(Field * field)
{
	nx = field->nx;
	ny = field->ny;
	nz = field->nz;
	Lx = field->Lx;
	Ly = field->Ly;

	dx = 2*Lx/(nx-1);
	dy = 2*Ly/(ny-1);
	dz = 1.0/(nz-1);

	values.reserve(nx*ny*nz);
	this->setAll(0);

	/* values = field.values; */
}

Field::Field(Field & field)
{
	nx = field.nx;
	ny = field.ny;
	nz = field.nz;
	Lx = field.Lx;
	Ly = field.Ly;

	dx = field.dx;
	dy = field.dy;
	dz = field.dz;

	/* dx = 2*Lx/(nx-1); */
	/* dy = 2*Ly/(ny-1); */
	/* dz = 1.0/(nz-1); */

	values.reserve(nx*ny*nz);
	this->setAll(0);

}

void Field::init(int inx, int iny, int inz, double iLx, double iLy)
{
	nx = inx;
	ny = iny;
	nz = inz;

	Lx = iLx;
	Ly = iLy;

	dx = 2*Lx/(nx-1);
	dy = 2*Ly/(ny-1);
	dz = 1.0/(nz-1);

	values.reserve(nx*ny*nz);
	std::fill(values.begin(), values.end(), 0);
}

void Field::set(int i, int j, int k, double value)
{
	/* values[i + nx*j + nx*ny*k] = value; */
	values[index(i,j,k)] = value;
}

void Field::set(Side side, double value)
{
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				if(this->getSide(i,j,k) == side)
				{
					this->set(i,j,k, value);
				}
			}
}

void Field::set(int i, int j, double value)
{
	/* values[i + nx*j] = value; */
	values[index(i,j,0)] = value;
}

void Field::setAll(double value)
{

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				/* values[i + nx*j + nx*ny*k] = value; */
				values[index(i,j,k)] = value;
			}


}

double Field::get(int i, int j, int k)
{
	/* return values[i + nx*j + nx*ny*k]; */
	return values[index(i,j,k)];
}

double Field::get(int i, int j)
{
	/* return values[i + nx*j]; */
	return values[index(i,j,0)];
}

std::unique_ptr<Field> Field::getSubfield(int i1, int i2, int j1, int j2, int k1, int k2)
{
	/* Field * subfield = new Field(i2-i1+1, j2-j1+1, k2-k1+1, 0, 0); */
	auto subfield = std::make_unique<Field>(i2-i1+1, j2-j1+1, k2-k1+1, Lx, Ly);
	subfield->dx = dx;
	subfield->dy = dy;
	subfield->dz = dz;

	for(int i = i1; i <= i2; i++)
		for(int j = j1; j <= j2; j++)
			for(int k = k1; k <= k2; k++)
			{
				subfield->set(i-i1, j-j1, k-k1, this->get(i,j,k));
			}

	return subfield;
}

void Field::setSubfield(int i1, int i2, int j1, int j2, int k1, int k2, Field & field)
{
	for(int i = i1; i <= i2; i++)
		for(int j = j1; j <= j2; j++)
			for(int k = k1; k <= k2; k++)
			{
				this->set(i, j, k, field.get(i,j,k));
			}

}

bool Field::isFinite()
{
	//doesn't work for some reason.
	/* for(std::vector<double>::iterator i = this->values.begin(); i != this->values.end(); ++i) */
	/* { */
	/* 	printf("checking\n"); */
	/* 	if(!std::isfinite(*i)) */
	/* 	{ */
	/* 		printf("it's infinite\n"); */
	/* 		return false; */

	/* 	} */
	/* } */

	for(int i=0; i<this->nx; i++)
		for(int j=0; j<this->ny; j++)
			for(int k=0; k<this->nz; k++)
			{
				if(!std::isfinite(this->get(i,j,k)))
					return false;
			}


	return true;
}

std::unique_ptr<Field> Field::fabs()
{
	auto newField = std::make_unique<Field>(this);

	for(int i=0; i<this->nx; i++)
		for(int j=0; j<this->ny; j++)
			for(int k=0; k<this->nz; k++)
			{
				newField->set(i,j,k, std::fabs(this->get(i,j,k)));
			}

	return  newField;

}

/*
 * Differentiator for field datatypes.
 *
 * Currently always applies forward and backward FD derivatives at respective boundaries.
 *
 * Interior region is handled based on requested scheme.
 *
 * However, if the requested scheme has a larger stencil than available grid points, the calculation proceeds as normal
 * with the unavailable grid points not considered.
 *
 *
 */
std::unique_ptr<Field> Field::differentiate(dMode mode, dDir dir)
{
	std::vector<double> vec = coeffs[mode];
	int h = (mode == CXX2 || mode == FXX1 || mode == BXX1) ? 2 : 1;

	std::vector<double> vecLeft = (h==1)? coeffs[FX1] : coeffs[FXX1];
	std::vector<double> vecRight= (h==1)? coeffs[BX1] : coeffs[BXX1];

	double value;

	/* Field * dfield = new Field(this->nx, this->ny, this->nz, this->Lx, this->Ly); */
	/* dfield->setAll(0); */
	auto dfield = std::make_unique<Field>(this->nx, this->ny, this->nz, this->Lx, this->Ly);

	if( (ny == 1 && dir == Y) || (nz == 1 && dir == Z) || (nx == 1 && dir == X) )
		return dfield;

	switch(dir)
	{
		case X:
			for(int i=1; i<this->nx-1; i++)
				for(int j=0; j<this->ny; j++)
					for(int k=0; k<this->nz; k++)
					{
						value = 0;
						for(int p=-2; p<=2; p++)
							if( i+p >= 0 && i+p < this->nx)
								value += this->get(i+p, j, k)*vec[p+2]/std::pow(dx, h);
						dfield->set(i,j,k,value);
					}
			for(int j=0; j<this->ny; j++)
				for(int k=0; k<this->nz; k++)
				{
					int i = 0;
					value = 0;
					for(int p=0; p<=2; p++)
						value += this->get(i+p, j, k)*vecLeft[p+2]/std::pow(dx, h);
					dfield->set(i,j,k,value);

					i = this->nx-1;
					value = 0;
					for(int p=-2; p<=0; p++)
						value += this->get(i+p, j, k)*vecRight[p+2]/std::pow(dx, h);
					dfield->set(i,j,k,value);
				}
			break;
		case Y:
			for(int i=0; i<this->nx; i++)
				for(int j=1; j<this->ny-1; j++)
					for(int k=0; k<this->nz; k++)
					{
						value = 0;
						for(int p=-2; p<=2; p++)
							if( j+p >= 0 && j+p < this->ny)
								value += this->get(i, j+p, k)*vec[p+2]/std::pow(dy, h);
						dfield->set(i,j,k,value);
					}
			for(int i=0; i<this->nx; i++)
				for(int k=0; k<this->nz; k++)
				{
					int j = 0;
					value = 0;
					for(int p=0; p<=2; p++)
						value += this->get(i, j+p, k)*vecLeft[p+2]/std::pow(dy, h);
					dfield->set(i,j,k,value);

					j = this->ny-1;
					value = 0;
					for(int p=-2; p<=0; p++)
						value += this->get(i, j+p, k)*vecRight[p+2]/std::pow(dy, h);
					dfield->set(i,j,k,value);
				}
			break;
		case Z:
			for(int i=0; i<this->nx; i++)
				for(int j=0; j<this->ny; j++)
					for(int k=1; k<this->nz-1; k++)
					{
						value = 0;
						for(int p=-2; p<=2; p++)
							if( k+p >= 0 && k+p < this->nz)
								value += this->get(i, j, k+p)*vec[p+2]/std::pow(dz, h);
						dfield->set(i,j,k,value);
					}
			for(int i=0; i<this->nx; i++)
				for(int j=0; j<this->ny; j++)
				{
					int k = 0;
					value = 0;
					for(int p=0; p<=2; p++)
						value += this->get(i, j, k+p)*vecLeft[p+2]/std::pow(dz, h);
					dfield->set(i,j,k,value);

					k = this->nz-1;
					value = 0;
					for(int p=-2; p<=0; p++)
						value += this->get(i, j, k+p)*vecRight[p+2]/std::pow(dz, h);
					dfield->set(i,j,k,value);
				}
			break;
		case XY:
		case YZ:
		case ZX:
			printf("incorrect differentiation direction for this->n");
			exit(-1);
			break;
		default:
			printf("incorrect switch, breaking \n");

	}
	return dfield;
}

/*
 * Integrates fields in X/Y/XY based on the field itself
 */
double Field::integrateXY()
{
	double integ=0;
	if(ny != 1 && nx != 1)
	{
		integ += dx * dy * ( this->get(0, 0) + this->get(0, ny-1) + this->get(nx-1, 0) + this->get( nx-1, ny-1) ) / 4;

		for(int i = 1; i < this->nx-1; i++)
		{
			integ += dx * dy * this->get(i, 0) / 2;
			integ += dx * dy * this->get(i, ny-1) / 2;
		}

		for(int j = 1; j < this->ny-1; j++)
		{
			integ += dx * dy * this->get(0, j)/2;
			integ += dx * dy * this->get(nx-1, j)/2;

		}

		for(int i = 1; i < this->nx-1; i++)
			for(int j = 1; j < this->ny-1; j++)
			{
				integ += dx * dy * this->get(i, j);
			}
	}
	else if(ny == 1)
	{
		integ += this->get(0,0) + this->get(nx-1, 0);

		for(int i=1; i<nx-1; i++)
			integ += 2*this->get(i,0);

		integ = (double)(2*Lx)*integ/(2*nx);

	}
	else if(nx == 1)
	{
		integ += this->get(0,0) + this->get(0, ny-1);

		for(int j=1; j<ny-1; j++)
			integ += 2*this->get(0,j);

		integ = (double)(2*Ly)*integ/(2*ny);

	}

	return integ;

}

std::unique_ptr<Field> Field::copy()
{
	/* Field * newField = new Field(this); */
	auto newField = std::make_unique<Field>(this);

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k) );
			}

	return newField;
}

std::unique_ptr<Field> Field::add(double value)
{
	auto newField = std::make_unique<Field>(this);

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k) + value);
			}

	return newField;

}

std::unique_ptr<Field> Field::add(Field & field)
{
	auto newField = std::make_unique<Field>(this);

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k) + field.get(i,j,k));
			}

	return newField;

}

std::unique_ptr<Field> Field::subtract(Field & field)
{
	auto newField = std::make_unique<Field>(this);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k) - field.get(i,j,k));
			}

	return newField;

}

std::unique_ptr<Field> Field::multiply(double factor)
{

	auto newField = std::make_unique<Field>(this);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k) * factor);
			}

	return newField;

}

std::unique_ptr<Field> Field::multiply(Field & field)
{

	auto newField = std::make_unique<Field>(this);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k) * field.get(i,j,k));
			}

	return newField;

}

std::unique_ptr<Field> Field::divide(Field & field)
{
	auto newField = std::make_unique<Field>(this);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k)/field.get(i,j,k));
			}

	return newField;
}

std::unique_ptr<Field> Field::divide(double value)
{
	auto newField = std::make_unique<Field>(this);
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				newField->set(i,j,k, this->get(i,j,k)/value);
			}

	return newField;
}

std::unique_ptr<Field> Field::pow(double n)
{

	auto newField = std::make_unique<Field>(this);
	for(int i=0; i<nx*ny*nz; i++)
		newField->values[i] = std::pow(values[i], n);

	return newField;
}



std::unique_ptr<Field> Field::replicateZ(int inz)
{
	/* Field * newField = new Field(nx, ny, inz, Lx, Ly); */
	auto newField = std::make_unique<Field>(nx, ny, inz, Lx, Ly);

	for(int k=0; k<inz; k++)
		for(int i = 0; i<nx; i++)
			for(int j = 0; j<ny; j++)
			{
				newField->set(i,j,k, this->get(i,j));
			}

	return newField;
}


double Field::average()
{
	double avg=0;

	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				avg += this->get(i,j,k) ;
			}

	avg /= nx*ny*nz;

	return avg;

}

void Field::print()
{
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
			for(int k=0; k<nz; k++)
			{
				printf("%d\t%d\t%d\t%e\n", i, j, k, this->get(i,j,k));
			}
}

/*
 * Field.getSide() : returns the location of the given point,
 * i.e., whether the point is INTERIOR/LEFT/RIGHT/...
 *
 */
Side Field::getSide(int i, int j, int k)
{
	if((this->nz != 1) && (k==0))
	{
		return SOUTH;
	}
	else if((this->nz != 1) && (k==this->nz-1))
	{
		return NORTH;
	}
	else if((this->ny != 1) && j == 0)
	{
		return FRONT;
	}
	else if((this->ny != 1) && j == this->ny-1)
	{
		return BACK;
	}
	else if(i == 0)
	{
		return LEFT;
	}
	else if(i == this->nx-1)
	{
		return RIGHT;
	}
	else return INTERIOR;

}

double Field::xVal(int i)
{
	return (nx==1)? 0 : -Lx + 2*Lx/(nx-1)*i;
}

double Field::yVal(int j)
{
	return (ny==1)? 0 : -Ly + 2*Ly/(ny-1)*j;
}

void Field::set(std::string expression_string, std::map<std::string, double> variables)
/* void Field::set(std::string expression_string) */
{

	const std::string e = expression_string;

	/* std::map<std::string, double> variables; */
	variables["Lx"]     = this->Lx;
	variables["Ly"]     = this->Ly;

	/* double radianTheta = variables["radianTheta"]; */
	/* double vector; */

	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
			for(int k = 0; k < nz; k++)
			{
				/* vector = cos(radianTheta) * xVal(i) + sin(radianTheta) * yVal(j); */
				/* variables["U"] = variables["U0"]*(1-vector/variables["r"]); */
				variables["x"] = xVal(i);
				variables["y"] = yVal(j);
				variables["z"] = (double)k * dz;
				this->set(i,j,k, Lepton::Parser::parse(e).evaluate(variables));

			}

}

int Field::index(int i, int j, int k)
{
	/* return (k + nz*j + nz*ny*i); */
	return (i + nx*j + nx*ny*k);
}

/*******************************************************************************************/
std::unique_ptr<Field> Field::operator=(Field & field)
{
	/* Field * newField = new Field(this); */
	auto newField = std::make_unique<Field>(field);

	for(int i=0; i<field.nx; i++)
		for(int j=0; j<field.ny; j++)
			for(int k=0; k<field.nz; k++)
			{
				newField->set(i,j,k, field.get(i,j,k) );
			}

	return newField;
}
