#include "model.hpp"
#include <iostream>

#include "optimization.hpp"

/* std::string circular(double r, double offsetTheta, int num); */
void calc_factor(double num_sensors, double * halfwidth, double * factor);

int main(int argc, char * argv[])
{

	/* Model model("inputs.ini"); */
	/* model.init(); */
	/* model.solve3(); */

	/* optimize_xy(); */
	optimize_r();
	/* optiSolve(); */
	/* optiSolveDelta(); */



}

void calc_factor(double num_sensors, double * halfwidth, double * factor)
{
	/* *halfwidth = sqrt(2.25e-3/(num_sensors*M_PI)); */
	*halfwidth = 0.01;
	double sigma = *halfwidth/sqrt(2*log(2));
	*factor = 1/(2*sigma*sigma);
}
