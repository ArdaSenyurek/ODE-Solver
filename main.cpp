#include "ode.h"
#include "../modules/uint.h"
#include <cmath>

tensor<float> dotFuncs(const tensor<float>);


int main(){
	tensor<float> in(3,1);
	in.setInOrder();
	in.set(2,0, 123);
	in.print();
	in.printLinear();
	in[1];
	tensor<float> (*dotRulePtr)(const tensor<float>);
	dotRulePtr = &dotFuncs;
	dotRulePtr(in).print();
	//euler<float>(1, 1, 0.01, diffEq);
	return 0;
}	

tensor<float> dotFuncs(const tensor<float> state)
{
	float phi 	= state[1].get();
	float y 	= state[2].get();
	float t 	= state[3].get();
	
	tensor <float> dotTensor(3,1);
       dotTensor = {
	
		2 * t,
		phi,
		1
	};
       return dotTensor;
}
