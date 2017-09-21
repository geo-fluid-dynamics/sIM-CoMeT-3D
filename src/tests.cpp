#include<iostream>
#include "field.hpp"
#include "solver.hpp"
#include "operators.hpp"
#include "plot.hpp"
#include <assert.h>

void ssade2d()
{
	std::unique_ptr<Field> T = std::make_unique<Field>(50, 50, 1, 10, 10);
	std::unique_ptr<Field> u = std::make_unique<Field>(50, 50, 1, 10, 10);
	std::unique_ptr<Field> alpha = std::make_unique<Field>(50, 50, 1, 10, 10);
	std::unique_ptr<Field> erf = std::make_unique<Field>(50, 50, 1, 10, 10);
	alpha->setAll(0.57);
	u->setAll(1);
	std::map<std::string, double> variables;
	erf->set("10*(1-erf(y/(2*sqrt(0.57*(x+10)))))", variables);

	auto zero = new Field(*T);
	auto TBC = new Field(*T);
	auto TBCFlag = new Field (*T);

	auto eqn_x = *u * (-1);
	PDE ssade2d_eqn(*T);
	ssade2d_eqn.x = eqn_x.get();
	ssade2d_eqn.yy = alpha.get();

	/* ssade2d_eqn.xx = zero; */
	/* ssade2d_eqn.zz = zero; */
	/* ssade2d_eqn.y  = zero; */
	/* ssade2d_eqn.z  = zero; */
	/* ssade2d_eqn.xy = zero; */
	/* ssade2d_eqn.yz = zero; */
	/* ssade2d_eqn.zx = zero; */
	/* ssade2d_eqn.rhs = zero; */

	TBCFlag->set(FRONT, 0);
	TBCFlag->set(BACK, 0);
	TBCFlag->set(LEFT, 0);
	TBCFlag->set(RIGHT, 3);

	TBC->set(FRONT, 20);
	TBC->set(BACK, 0);
	for(int i=0; i<25; i++)
		TBC->set(0, i, 20);

	Solver ssade2d_solve(&ssade2d_eqn, TBC, TBCFlag);
	ssade2d_solve.init();
	ssade2d_solve.solve(T.get());

	Plot plot, erfplot;
	plot.image(T.get());
	erfplot.image(erf.get());

	auto err = *T - *erf;
	Plot errplot;
	errplot.image(err.get());

	assert( fabs( err->average() ) < 1e-10 );

}


int main()
{
	ssade2d();
}
