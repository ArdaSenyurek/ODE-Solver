#include "ode.h"
#include "../modules/uint.h"
#include <cmath>

tensor<double> dotFuncs(const tensor<double>);


int main(){
	tensor<double> in(2,1,0);
	in.set(1,1, 5);
	in.print();
	in.printLinear();
	in[1];
	dotFuncs(in);
	tensor<double> (*dotRulePtr)(tensor<double>);
	dotRulePtr = &dotFuncs;
	euler<double> solver = euler<double>(in, 1.5, 3, dotRulePtr);
	ode<double>* solverPtr = &solver;

	solverPtr -> solve();
	return 0;
}	

tensor<double> dotFuncs(const tensor<double> state)
{
	double phi 	= state[1].get();
	double y 	= state[2].get();
	double t 	= state[3].get();
	
	tensor <double> dotTensor(state, -999);
       dotTensor = {
	
		2 * t,
		phi,
		1
	};
       return dotTensor;
/*
	double y 	= state[1].get();
	double t 	= state[2].get();
	
	tensor <double> dotTensor(state, -999);
       dotTensor = {
		3 * exp(-t) - 0.4 * y,
		1
	};
       return dotTensor;
	*/
}
