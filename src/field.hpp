#ifndef FIELD_H
#define FIELD_H

#include <map>
#include <vector>
#include <cmath>
#include <string>
#include <memory>
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

		double dx;
		double dy;
		double dz;

		/* Field(Model&  model); */
		Field();
		Field(Field * field);
		Field(Field & field);
		Field(int inx, int iny, int inz, double iLx, double iLy);
		void init(int inx, int iny, int inz, double iLx, double iLy);

		void set(int i, int j, int k, double value);
		void set(Side side, double value);
		void set(int i, int j, double value);
		/* void set(std::string expression_string); */
		void set(std::string expression_string, std::map<std::string, double> variables);
		void setSubfield(int i1, int i2, int j1, int j2, int k1, int k2, Field & field);
		void setAll(double value);
		void print();

		double get(int i, int j, int k);
		double get(int i, int j);
		double integrateXY();
		double average();

		bool isFinite();

		std::unique_ptr<Field> fabs();
		std::unique_ptr<Field> getSubfield(int i1, int i2, int j1, int j2, int k1, int k2);
		std::unique_ptr<Field> differentiate(dMode mode, dDir dir);
		std::unique_ptr<Field> add(Field & field);
		std::unique_ptr<Field> subtract(Field & field);
		std::unique_ptr<Field> add(double value);
		std::unique_ptr<Field> multiply(double factor);
		std::unique_ptr<Field> multiply(Field & field);
		std::unique_ptr<Field> copy();
		std::unique_ptr<Field> replicateZ(int inz);
		std::unique_ptr<Field> pow(double n);
		std::unique_ptr<Field> divide(Field & field);
		std::unique_ptr<Field> divide(double value);

		std::unique_ptr<Field> operator=(Field & field);
		/* std::unique_ptr<Field> Field::operator%(dMode mode, dDir dir); */

		Side getSide(int i, int j, int k);

		double xVal(int i);
		double yVal(int j);
		int index(int i, int j, int k);

};

#endif
