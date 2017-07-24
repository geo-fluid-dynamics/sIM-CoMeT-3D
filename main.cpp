#include "model.hpp"
#include <iostream>

int main(int argc, char * argv[])
{
	Model model("inputs.ini");
	model.print();
	/* model.legacySolve(); */
	model.solve();
}

