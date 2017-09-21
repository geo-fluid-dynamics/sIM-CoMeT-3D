#include "optimization.hpp"
#include "model.hpp"
#include "operators.hpp"
#include "plot.hpp"

#include <nlopt.hpp>
#include <iostream>


double obj_func_xy(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	int * count = (int *)func_data;
	*count += 1;

	Model model("inputs.ini");
	model.southBC = DIRICHLET;
	model.TwExp = "40 + ampl*exp(-factor*((x-(" + std::to_string(x[0]) + "))^2 + (y-(" + std::to_string(x[1]) + "))^2))";

	double factor, r_sensor, num_sensors = 1;
	r_sensor = 0.01;
	double sigma = r_sensor/sqrt(2*log(2));
	factor = 1/(2*sigma*sigma);

	model.variables["ampl"] = 30;
	model.variables["factor"] = factor;

	model.init();
	std::cout << ">>>>> Call " << *count << " x = " << x[0] << " y = " << x[1] <<"\n";
	model.solve2();

	return -model.U0;

}

double c1_xy(const std::vector<double> &x)
{
	return x[0]*x[0] + x[1]*x[1] - 0.06*0.06;
}

void optimize_xy()
{
	nlopt::opt opt(nlopt::LN_BOBYQA, 2);
	std::vector<double> lb(2), ub(2);
	lb[0] = -0.05; lb[1] = -0.05;
	ub[0] = 0.05; ub[1] = 0.05;

	int count=0;

	opt.set_min_objective(obj_func_xy, (void*) &count);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.set_xtol_rel(1e-4);

	std::vector<double> x(2);
	x[0] = 0.01;  x[1] = 0.01;
	double minf;
	nlopt::result result = opt.optimize(x, minf);

	std::cout << "RESULT: " << x[0] << "\t" << x[1] << "\t" << minf <<"\n";
	std::cout << "#eval = " << count << "\n";

}

void optimize_r()
{

	nlopt::opt opt(nlopt::LN_COBYLA, 3);
	std::vector<double> lb(3), ub(3);
	lb[0] = 0.1; ub[0] = 0.3;
	lb[1] = 0.05; ub[1] = 0.15;
	lb[2] = 0; ub[2] = 90;

	int count=0;

	opt.set_max_objective(obj_func_r, (void*) &count);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.set_xtol_rel(1e-4);

	std::vector<double> x(3);
	x[0] = 0.1; x[1] = 0.1; x[2] = 20;
	double minf;
	nlopt::result result = opt.optimize(x, minf);

	std::cout << "RESULT: " << x[0] << "\t" << x[1] <<"\t" << x[2] << "\t"<< minf <<"\n";
	std::cout << "#eval = " << count << "\n";
}


double obj_func_r(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	int * count = (int *)func_data;
	*count += 1;

	Model model("inputs.ini");
	model.southBC = DIRICHLET;

	double factor, r_sensor, num_sensors = 4;
	r_sensor = 0.05;
	double sigma = r_sensor/sqrt(2*log(2));
	factor = 1/(2*sigma*sigma);

	model.TwExp = "40" + circular(x[0], 0, num_sensors) + circular(x[0]+x[1], x[2], num_sensors);

	model.variables["ampl"] = 30;
	model.variables["factor"] = factor;

	model.init();
	std::cout << ">>>>> Call " << *count << " r = " << x[0] << " delta r = " << x[1] << " phase = " << x[2]<< "\n";
	model.solve2();

	return model.U0;

}

std::string circular(double r, double offsetTheta, int num)
{
	double offsetRadian = M_PI * offsetTheta / 180;
	double radian = M_PI * 2 / (num);

	std::string expr = "";

	double x, y;

	for(int i=0; i<num; i++)
	{
		x = r*cos(offsetRadian + i*radian);
		y = r*sin(offsetRadian + i*radian);
		expr = expr + "- ampl*exp(-factor*((x-(" + std::to_string(x) + "))^2 + (y-(" + std::to_string(y) + "))^2))";
	}

	return expr;

}

double F_constraint(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	Model * model = (Model *) func_data;

	model->U0 = x[0];
	model->r = x[1];

	/* model->U0 = std::exp(x[0]); */
	/* model->r = -1/x[1]; */

	auto qsba = *(*(model->qNorth->differentiate(CX2, X)) + *(model->qNorth->differentiate(CX2, Y))) / model->qNorth->average();
	auto slope = *qsba * x[2];
	std::unique_ptr<Field> xy = std::make_unique<Field>(*model->delta);
	xy->set("x+y", model->variables);
	model->delta = *(*slope * *xy) + x[2];

	model->update_fields();
	model->F = model->p->integrateXY();
	model->init_uvw();

	/* return (model->F - model->Fscrew); */
	return (log(model->F/model->Fscrew));
}

double M_constraint(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	Model * model = (Model *) func_data;

	model->U0 = x[0];
	model->r = x[1];

	/* model->U0 = std::exp(x[0]); */
	/* model->r = -1/x[1]; */

	auto qsba = *(*(model->qNorth->differentiate(CX2, X)) + *(model->qNorth->differentiate(CX2, Y))) / model->qNorth->average();
	auto slope = *qsba * x[2];
	std::unique_ptr<Field> xy = std::make_unique<Field>(*model->delta);
	xy->set("x+y", model->variables);
	model->delta = *(*slope * *xy) + x[2];

	model->update_fields();
	model->F = model->p->integrateXY();
	model->init_uvw();
	return (model->Mtheta);
}

double MF_constraint(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	Model * model = (Model *) func_data;

	model->U0 = x[0];
	model->r = x[1];

	auto qsba = *(*(model->qNorth->differentiate(CX2, X)) + *(model->qNorth->differentiate(CX2, Y))) / model->qNorth->average();
	auto slope = *qsba * x[2];
	std::unique_ptr<Field> xy = std::make_unique<Field>(*model->delta);
	xy->set("x+y", model->variables);
	model->delta = *(*slope * *xy) + x[2];

	model->update_fields();
	model->F = model->p->integrateXY();
	model->init_uvw();
	return std::pow((model->F-model->Fscrew),2) + std::pow((model->Mtheta),2);

}

double F_pos(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	Model * model = (Model *) func_data;

	model->U0 = x[0];
	model->r = x[1];

	/* model->U0 = std::exp(x[0]); */
	/* model->r = -1/x[1]; */

	auto qsba = *(*(model->qNorth->differentiate(CX2, X)) + *(model->qNorth->differentiate(CX2, Y))) / model->qNorth->average();
	auto slope = *qsba * x[2];
	std::unique_ptr<Field> xy = std::make_unique<Field>(*model->delta);
	xy->set("x+y", model->variables);
	model->delta = *(*slope * *xy) + x[2];

	model->update_fields();
	model->F = model->p->integrateXY();
	model->init_uvw();

	/* return (model->F - model->Fscrew); */
	return (-model->F);

}

int count =0;
double optiSolveObjFunc(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	Model * model = (Model *) func_data;

	++count;

	model->U0 = x[0];
	model->r = x[1];

	/* model->U0 = std::exp(x[0]); */
	/* model->r = -1/x[1]; */

	auto qsba = *(*(model->qNorth->differentiate(CX2, X)) + *(model->qNorth->differentiate(CX2, Y))) / model->qNorth->average();
	auto slope = *qsba * x[2];
	std::unique_ptr<Field> xy = std::make_unique<Field>(*model->delta);
	xy->set("x+y", model->variables);
	model->delta = *(*slope * *xy) + x[2];


	model->update_fields();
	model->init_uvw();
	/* model->TSolveWrapper(); */
	model->calc_maxFluxError();

	/* model->printOutputs(); */
	std::cout << model->U0 << "\t"<< model->r << "\n";

	/* return std::pow((model->F-model->Fscrew),2) + std::pow((model->Mtheta),2); */
	return model->maxFluxError;

}

void optiSolve()
{
	nlopt::opt opt(nlopt::LN_COBYLA, 3);

	std::vector<double> lb(3), ub(3);
	lb[0] = 1e-10; ub[0] = 1e-2;
	lb[1] = -1; ub[1] = 1;
	lb[2] = 1e-10; ub[2] = 1e-3;

	Model model("inputs.ini");
	model.init();

	model.find_U();
	model.find_r();
	model.recalcDelta=0;
	model.init_uvw();
	model.TSolveWrapper();
	model.calc_maxFluxError();

	void * voidModel = &model;

	opt.set_min_objective(optiSolveObjFunc, voidModel);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.add_equality_constraint(F_constraint, voidModel, 1e-10);
	opt.add_equality_constraint(M_constraint, voidModel, 1e-10);
	opt.add_inequality_constraint(F_pos, voidModel, 1e-15);

	opt.set_xtol_rel(1e-4);

	std::vector<double> x(3);
	x[0] = model.U0; x[1] = model.r; x[2] = model.delta->average();
	/* x[0] = log(model.U0); x[1] = -1/model.r; x[2] = model.delta->average(); */

	double mfe;
	nlopt::result result = opt.optimize(x, mfe);

	model.printOutputs();
	/* Plot plot; */
	/* model.w->print(); */

	std::cout << "RESULT: " << x[0] << "\t" << x[1] <<"\t" << mfe <<"\n";
	std::cout << "#eval = " << count << "\n";

}

//////////////////////////////////////////


double deltaObjFunc(const std::vector<double> &x, std::vector<double> &grad, void * func_data)
{
	Model * model = (Model *) func_data;

	auto qsba = *(*(model->qNorth->differentiate(CX2, X)) + *(model->qNorth->differentiate(CX2, Y))) / model->qNorth->average();
	auto slope = *qsba * x[0];
	std::unique_ptr<Field> xy = std::make_unique<Field>(*model->delta);
	xy->set("x+y", model->variables);
	model->delta = *(*slope * *xy) + x[0];

	/* model->delta = *model->delta * (model->delta->average() / x[0]); */
	std::cout << x[0] << "\n";

	model->find_U();
	model->find_r();
	model->init_uvw();
	model->TSolveWrapper();
	model->calc_maxFluxError();

	return model->maxFluxError;

}

void optiSolveDelta()
{
	nlopt::opt opt(nlopt::LN_COBYLA, 1);

	std::vector<double> lb(1), ub(1);
	lb[0] = 1e-7; ub[0] = 1e-3;

	Model model("inputs.ini");
	model.init();

	model.find_U();
	model.find_r();
	model.recalcDelta=0;
	model.init_uvw();
	model.TSolveWrapper();
	model.calc_maxFluxError();

	void * voidModel = &model;

	opt.set_min_objective(deltaObjFunc, voidModel);
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	opt.set_xtol_rel(1e-4);

	std::vector<double> x(1);
	x[0] = model.delta->average();
	/* x[0] = log(model.U0); x[1] = -1/model.r; x[2] = model.delta->average(); */


	double mfe;
	nlopt::result result = opt.optimize(x, mfe);

	model.printOutputs();
	/* Plot plot; */
	/* model.w->print(); */

	std::cout << "RESULT: " << x[0] <<"\t" << mfe <<"\n";
	std::cout << "#eval = " << count << "\n";

}
