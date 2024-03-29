#ifndef MODEL_H
#define MODEL_H

#include <string>
#include "enum.hpp"
#include "field.hpp"
#include <memory>
#include <map>

class Model {

	public:
		std::map<std::string, double> variables;
		std::string modelID;

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

		std::unique_ptr<Field> Tw;
		std::unique_ptr<Field> qw;
		std::unique_ptr<Field> delta;
		std::unique_ptr<Field> p;
		std::unique_ptr<Field> qNorth;
		std::unique_ptr<Field> qStefan;
		std::unique_ptr<Field> bcSouth;

		std::unique_ptr<Field> T;
		std::unique_ptr<Field> u;
		std::unique_ptr<Field> v;
		std::unique_ptr<Field> w;
		std::unique_ptr<Field> vec;
		std::unique_ptr<Field> U;

		std::unique_ptr<Field> sensorFilter;

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
		double relativeU0;

		int exitVarFlag;
		std::map< std::string, Field * > fieldMap;

		std::unique_ptr<Field> radianThetaField;
		std::unique_ptr<Field> vecField;

		/* Model(); */
		Model(std::string iniPath);
		~Model();
		void solve();
		void printInputs();
		void printOutputs();
		void combinedUpdate();
		void combinedUpdate2();
		void Plot_FU();
		void Plot_Mr();

	private:
		double xVal(int i);
		double yVal(int j)  ;
		double zVal(int i, int j, int k) ;
		double UVal(int i, int j) ;

		void update_fields();

		void find_U();
		void find_r();
		void init_uvw();
		void calc_maxFluxError();
		void adjustDelta();

		void TSolveWrapper();
		void PSolveWrapper();

		void updateMap();

		void parseINI(std::string iniPath);

		std::unique_ptr<Field> get_rhsDdelta();


};
#endif
