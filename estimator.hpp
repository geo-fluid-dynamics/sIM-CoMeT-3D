#ifndef ESTIMATOR_H
#define ESTIMATOR_H
#include <vector>

class Estimator {
	public:
		/* Estimator(); */
		/* Estimator(std::vector<double>& x, std::vector<double>& y); */
		double linear(std::vector<double>& x, std::vector<double>& y, double x0);
		double loglinear(std::vector<double>& x, std::vector<double>& y, double x0);

};
#endif
