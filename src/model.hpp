#ifndef MODEL_H
#define MODEL_H
#include <string>
#include "enum.hpp"
#include "field.hpp"

class Model {
	public:
		double mu   = 0.001;
		double Tm   = 0;
		double hm   = 3.337e5;
		double cpL  = 4222.22;
		double cpS  = 2049.41;
		double rhoS = 920;
		double rhoL = 1000;
		double kL   = 0.57;

		double Lx     = 0.075;
		double Ly     = 0.075;
		double Fscrew = 375;

		std::string TwExp = "30+10*x/Lx";
		std::string qwExp = "10";

		double Tinf = -20;

		int nx = 20;
		int ny = 1;
		int nz = 20;

		double theta = 0;
		double radianTheta = theta*3.1415927/180;

		int recalcDelta = 1;

		int maxLoopIter = 100;
		double MTol = 1e-10;
		double FTol = 1e-6;
		double allowedMaxFluxError = 100;
		double allowedRelativeMFE = 1e-6;
		double deltaCoeffMin = 0.7;
		double deltaCoeffMax = 1.3;
		double deltaRelax = 0.05;

		boundary southBC = DIRICHLET;
		boundary sidesBC = NEUMANN;

		Field * Tw;
		Field * qw;
		Field * delta;
		Field * p;
		Field * qNorth;
		Field * qStefan;
		Field * bcSouth;

		Field * T;
		Field * u;
		Field * v;
		Field * w;

		double U0;
		double r;
		double Mx;
		double My;
		double Mtheta;
		double F;

		double dx = 2.0*Lx/(nx-1);
		double dy = 2.0*Ly/(ny-1);
		double dz = 1.0/(nz - 1);

		double hmStar = hm + cpS*(Tm-Tinf);
		double alpha = kL/(rhoL*cpL);

		double maxFluxError;
		double relativeMFE;

		Model();
		Model(std::string iniPath);
		~Model();
		void solve();
		void print();
		void combinedUpdate();
		void combinedUpdate2();

	private:
		double xVal(int i);
		double yVal(int j);
		double zVal(int i, int j, int k);
		double UVal(int i, int j);

		void init_fields();
		void update_fields();

		void find_U();
		void find_r();
		void init_uvw();
		void calc_maxFluxError();
		void adjustDelta();

		void TSolveWrapper();
		void PSolveWrapper();

		void dumper();

};
#endif
