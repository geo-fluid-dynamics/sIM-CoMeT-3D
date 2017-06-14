#include <vector>
#include <cmath>

class Field {

	private:
		int nx;
		int ny;
		int nz;

		double Lx;
		double Ly;
		double Lz;

		double dx;
		double dy;
		double dz;

		std::vector<double> values;

	public:
		Field(Model&  model)
		{
			nx = model.nx;
			ny = model.ny;
			nz = model.nz;

			Lx = model.Lx;
			Ly = model.Ly;
			Lz = model.Lz;

			dx = -Lx + 2*Lx/(nx-1);
			dy = -Ly + 2*Ly/(ny-1);
			dz = -Lz + 2*Lz/(nz-1);

			values.reserve(nx*ny*nz);
		}

		Field(int inx, int iny, int inz, double iLx, double iLy, double iLz)
		{
			nx = inx;
			ny = iny;
			nz = inz;

			Lx = iLx;
			Ly = iLy;
			Lz = iLz;

			dx = -Lx + 2*Lx/(nx-1);
			dy = -Ly + 2*Ly/(ny-1);
			dz = -Lz + 2*Lz/(nz-1);

			values.reserve(nx*ny*nz);
		}

		void set(int i, int j, int k, double value)
		{
			values[i + nx*j + nx*ny*k] = value;
		}

		double get(int i, int j, int k)
		{
			return values[i + nx*j + nx*ny*k];
		}

		Field getSubfield(Field& field, int i1, int i2, int j1, int j2, int k1, int k2)
		{
			Field subfield(i2-i1+1, j2-j1+1, k2-k1+1, 0, 0, 0);

			/* for(std::vector<double>::iterator i = field.values.begin(); i != field.values.end(); ++i) */

			for(int i = i1; i <= i2; i++)
				for(int j = j1; j <= j2; j++)
					for(int k = k1; k <= k2; k++)
					{
						subfield.set(i-i1, j-j1, k-k1, field.get(i,j,k));
					}

			return subfield;
		}

		bool isFinite(Field& field)
		{
			for(std::vector<double>::iterator i = field.values.begin(); i != field.values.end(); ++i)
			{
				std::isfinite(*i);
			}
		}

		Field differentiate(Field& field, dMode mode, dDir dir)
		{

		}

		double integrate(Field& field)
		{

		}

		void add(Field& field)
		{

		}

		void multiply(double factor)
		{

		}

		void multiply(Field& field)
		{

		}

		Field copy(Field& field)
		{

		}

};
