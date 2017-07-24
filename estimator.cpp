#include "estimator.hpp"
#include <cmath>
#include <iostream>

/* Estimator::Estimator(std::vector<double>& x, std::vector<double>& y) { */
/* 	this->x = x; */
/* 	this->y = y; */
/* } */

double Estimator::loglinear(std::vector<double>& x, std::vector<double>& y, double x0)
{
	std::vector<double> xlog(x);
	std::vector<double> ylog(y);

	double x0log = std::log10(x0);

	for(unsigned int i=0; i<x.size(); i++)
	{
		xlog[i] = std::log10(x[i]);
		ylog[i] = std::log10(y[i]);
	}

	double value = this->linear(xlog, ylog, x0log);
	value = std::pow(10, value);

	return value;
}

double Estimator::linear(std::vector<double>& x, std::vector<double>& y, double x0)
{
	double m, c, value;

	std::vector<double>::iterator xi = x.end();
	std::vector<double>::iterator yi = y.end();

	m = (*(yi-1) - *(yi-2))/(*(xi-1) - *(xi-2));
	c = *(yi-1) - (*(xi-1))*m;
	value = m*(x0) + c;


	return value;

}
