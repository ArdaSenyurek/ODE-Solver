#include "ode.h"
#include "../modules/uint.h"
#include <cmath>

tensor<float> dotFuncs(const tensor<float>);


int main(){
	tensor<float> in(3,1,0);
	in.print();
	in.printLinear();
	in[1];
	dotFuncs(in);
	tensor<float> (*dotRulePtr)(tensor<float>);
	dotRulePtr = &dotFuncs;
	euler<float> solver = euler<float>(in, 0.01, 0.05, dotRulePtr);
	ode<float>* solverPtr = &solver;

	solverPtr -> solve();
	return 0;
}	

tensor<float> dotFuncs(const tensor<float> state)
{
	float phi 	= state[1].get();
	float y 	= state[2].get();
	float t 	= state[3].get();
	
	tensor <float> dotTensor(state, -999);
       dotTensor = {
	
		2 * t,
		phi,
		1
	};
       return dotTensor;
}
