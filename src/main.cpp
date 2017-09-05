#include "model.hpp"
#include <iostream>

int main(int argc, char * argv[])
{
	Model model("inputs.ini");
	model.solve();
	/* model.combinedUpdate(); */
	/* model.combinedUpdate2(); */
	/* model.Plot_FU(); */
	/* model.Plot_Mr(); */
}

