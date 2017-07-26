#ifndef SPARSE_H
#define SPARSE_H
#include <vector>
#include "stencil.hpp"

class Sparse {

	public:
		std::vector<int> Ai;
		std::vector<int> Aj;
		std::vector<double> values;

		void push(int index1, int index2, double value);
		/* void pushStencil(Stencil& stencil, int row); */
		void pushStencil(Stencil& stencil, int row);
		bool isFinite();
		void print();
};
#endif
