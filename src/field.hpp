#ifndef FIELD_H
#define FIELD_H

#include <map>
#include <vector>
#include <cmath>
#include "enum.hpp"

class Field {

	private:
		std::vector<double> values;
		std::map < dMode , std::vector<double> > coeffs = {
			{ CX2 , {0.0, -0.5  , 0.0   , 0.5  , 0.0} }   ,
			{ CXX2, {0.0, 1.0   , -2.0  , 1.0  , 0.0} }   ,
			{ FX1, {0.0, 0.0   , -1.0  , 1.0  , 0.0} }   ,
			{ FX2, {0.0, 0.0   , -1.5  , 2.0  , -0.5} }  ,
			{ FXX1, {0.0, 0.0   , 1.0   , -2.0 , 1.0} }   ,
			{ BXX1, {1.0, -2.0  , 1.0   , 0.0  , 0.0} }   ,
			{ BX1, {0.0, -1.0  , 1.0   , 0.0  , 0.0} }   ,
			{ BX2, {0.5, -2.0  , 1.5   , 0.0  , 0.0} }
		};

	public:
		int nx;
		int ny;
		int nz;

		double Lx;
		double Ly;
		double Lz;

		double dx;
		double dy;
		double dz;

		/* Field(Model&  model); */
		Field();
		Field(Field * field);
		Field(int inx, int iny, int inz, double iLx, double iLy, double iLz);
		void init(int inx, int iny, int inz, double iLx, double iLy, double iLz);
		void set(int i, int j, int k, double value);
		double get(int i, int j, int k);
		double get(int i, int j);
		void set(int i, int j, double value);
		void setAll(double value);
		bool isFinite();
		double integrateXY();
		double average();
		Field * getSubfield(int i1, int i2, int j1, int j2, int k1, int k2);
		Field * differentiate(dMode mode, dDir dir);
		Field * add(Field * field);
		Field * subtract(Field * field);
		Field * add(double value);
		Field * multiply(double factor);
		Field * multiply(Field * field);
		Field * copy();
		Field * replicateZ(int inz);
		Field * pow(double n);
		Field * divide(Field * field);
		Field * divide(double value);
		void print();
		Side getSide(int i, int j, int k);

};

#endif
