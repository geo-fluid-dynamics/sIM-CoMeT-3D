#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include<vector>
#include<string>

double obj_func_xy(const std::vector<double> &x, std::vector<double> &grad, void * func_data);
double c1_xy(const std::vector<double> &x);
void optimize_xy();

std::string circular(double r, double offsetTheta, int num);
double obj_func_r(const std::vector<double> &x, std::vector<double> &grad, void * func_data);
void optimize_r();

double optiSolveObjFunc(const std::vector<double> &x, std::vector<double> &grad, void * func_data);
void optiSolve();
double F_constraint(const std::vector<double> &x, std::vector<double> &grad, void * func_data);
double M_constraint(const std::vector<double> &x, std::vector<double> &grad, void * func_data);

double deltaObjFunc(const std::vector<double> &x, std::vector<double> &grad, void * func_data);
void optiSolveDelta();

#endif
