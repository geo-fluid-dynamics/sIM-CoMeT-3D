#ifndef SOLVER_H
#define SOLVER_H
#include "sparse.hpp"
#include "PDE.hpp"

#include <vector>

class Solver {

	public:
		Sparse * A;
		std::vector<double> * b;
		std::vector<double> * x;

		PDE * pde;
		Field * field;
		Field * BCFlag;

		int lda;

		Solver(PDE * ipde, Field * field, Field * bcflag);
		~Solver();
		void init();
		/* void solve(); */
		void solve(Field * x);
		void setA(Sparse iA);
		void buildb();

};
#endif
