#include "packages/exprtk/exprtk.hpp"

void Field::set(std::string expression_string)
{

	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	/* std::string expression_string = "clamp(-1.0,sin(2 * pi * x) + cos(x / 2 * pi),+1.0)"; */

	double x;
	double y;
	double Lx = this->Lx;
	double Ly = this->Ly;

	symbol_table_t symbol_table;
	symbol_table.add_variable("x",x);
	symbol_table.add_variable("y",y);
	symbol_table.add_constants("Lx", Lx);
	symbol_table.add_constants("Ly", Ly);
	symbol_table.add_constants();

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	parser.compile(expression_string,expression);

	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
		{
			x = xVal(i);
			y = yVal(j);
			Tw->set(i,j,0, expression.value());
		}


}
