#include "ode.h"
#include "../modules/uint.h"
#include <cmath>

tensor<double> dotFuncs(const tensor<double>);


int main(){
	tensor<double> in(3,1,0);
	//in.set(1,1, 13);
	//in.set(2,1, 7);
	in.print();
	in.printLinear();
	in[1];
	dotFuncs(in);
	tensor<double> (*dotRulePtr)(tensor<double>);
	dotRulePtr = &dotFuncs;
	rungeKutta4<double> solver = rungeKutta4<double>(in, 0.1, 1, dotRulePtr);
	ode<double>* solverPtr = &solver;

	solverPtr -> solve();
	return 0;
}	

tensor<double> dotFuncs(const tensor<double> state)
{
	double phi 	= state[1].get();
	double y 	= state[2].get();
	double x 	= state[3].get();
	
	tensor <double> dotTensor(state, -999);
       dotTensor = {
	       2 * x * x + phi + y,
	       phi,
	       1
	};
       return dotTensor;
}
