#include "sparse.hpp"

void Sparse::push(int index1, int index2, double value)
{
	Ai.push_back(index1);
	Aj.push_back(index2);
	values.push_back(value);
}

void Sparse::pushStencil(Stencil& stencil, int row)
{
	int stride_x = 1;
	int stride_y = stencil.nx;
	int stride_z = stencil.nx*stencil.ny;

	for(int i=0; i<5; i++)
		for(int j=0; j<5; j++)
			for(int k=0; k<5; k++)
			{
				if(stencil._stencil[i][j][k] != 0)
					this->push(row, row + (i-2)*stride_x + (j-2)*stride_y + (k-2)*stride_z, stencil._stencil[i][j][k]);
			}
}

bool Sparse::isFinite()
{
	for(std::vector<int>::iterator i = Ai.begin(); i != Ai.end(); ++i)
	{
		if(!std::isfinite(*i))
			return false;
	}

	for(std::vector<int>::iterator i = Aj.begin(); i != Aj.end(); ++i)
	{
		if(!std::isfinite(*i))
			return false;
	}

	for(std::vector<double>::iterator i = values.begin(); i != values.end(); ++i)
	{
		if(!std::isfinite(*i))
			return false;
	}

	return true;

}

void Sparse::print()
{
	unsigned int size = Ai.size();

	for(unsigned int i = 0 ; i<size; i++)
	{
		printf("(%d,\t%d)\t%e\n", Ai[i], Aj[i], values[i]);
	}

}
