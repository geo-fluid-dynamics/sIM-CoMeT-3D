class Model {
	public:
		double mu = 0.001;
		double Tm = 0;
		double hm = 3.337e5;
		double cpL = 4222.22;
		double cpS = 2049.41;
		double rhoS = 920;
		double kL = 0.57;

		double Lx = 0.075;
		double Ly = 0.075;
		double Fscrew = 800;
		string TwExp = "30+10*x/L_x";
		string qwExp = "10";

		/* vector<double> Tw(nx*ny); */
		/* vector<double> qw(nx*ny); */

		double Tinf = -20;

		int nx = 30;
		int ny = 30;
		int nz = 20;

		int maxLoopIter = 100;
		double MTol = 1e-10;
		double FTol = 1e-6;
		double maxFluxError = 100;
		double deltaCoeffMin = 0.7;
		double deltaCoeffMax = 1.3;
		double deltaRelax = 0.01;

		boundary southBC = DIRICHLET;
		boundary sidesBC = NEUMANN;

};
