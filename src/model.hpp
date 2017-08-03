#ifndef MODEL_H
#define MODEL_H
#include <string>
#include "enum.hpp"
#include "field.hpp"

class Model {

	public:
		double mu;
		double Tm;
		double hm;
		double cpL;
		double cpS;
		double rhoS;
		double rhoL;
		double kL;
		double Lx;
		double Ly;
		double Fscrew;

		std::string TwExp;
		std::string qwExp;

		double Tinf;

		int nx;
		int ny;
		int nz;

		double theta;
		double radianTheta;

		int recalcDelta = 1;
		int maxMainIter;
		int maxFindIter;

		double MTol;
		double FTol;
		double allowedMaxFluxError;
		double allowedRelativeMFE;
		double deltaCoeffMin;
		double deltaCoeffMax;
		double deltaRelax;

		boundary southBC;
		boundary sidesBC;

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

		/* Model(); */
		Model(std::string iniPath);
		~Model();
		void solve();
		void print();
		void combinedUpdate();
		void combinedUpdate2();

	private:
		double xVal(int i) ;
		double yVal(int j)  ;
		double zVal(int i, int j, int k) ;
		double UVal(int i, int j) ;

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
