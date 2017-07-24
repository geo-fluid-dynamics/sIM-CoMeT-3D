#ifndef STENCIL_H
#define STENCIL_H

#include <vector>
#include <map>
#include "enum.hpp"
#include "field.hpp"

class Stencil {
	private:
		std::map < dMode , std::vector<double> > coeffs	= {
			{ CX2 , {0.0, -0.5  , 0.0   , 0.5  , 0.0} }   ,
			{ CXX2, {0.0, 1.0   , -2.0  , 1.0  , 0.0} }   ,
			{ FX1, {0.0, 0.0   , -1.0  , 1.0  , 0.0} }   ,
			{ FX2, {0.0, 0.0   , -1.5  , 2.0  , -0.5} }  ,
			{ FXX1, {0.0, 0.0   , 1.0   , -2.0 , 1.0} }   ,
			{ BXX1, {1.0, -2.0  , 1.0   , 0.0  , 0.0} }   ,
			{ BX1, {0.0, -1.0  , 1.0   , 0.0  , 0.0} }   ,
			{ BX2, {0.5, -2.0  , 1.5   , 0.0  , 0.0} } ,
			{ ONE, {0.0, 0.0  , 1   , 0.0  , 0.0} }
		};


	public:
		double _stencil[5][5][5];
		double dx; double dy; double dz;
		double nx; double ny; double nz;
		Stencil(double idx, double idy, double idz);
		Stencil(Field * field);
		void reset();
		void discretize(double coeff, dMode mode, dDir dir);
		bool isFinite();
		void print();


};
#endif
