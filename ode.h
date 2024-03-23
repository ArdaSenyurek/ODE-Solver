#include "../modules/tensor.h"
#include "../modules/uint.h"

template <class qnt>
class ode
{
	protected: 
		uint rank_;
		float endTime_;
		float dT_;
		tensor<qnt>* inits_;
		tensor<qnt>* current_;
		qnt  (*rule_)(tensor<qnt>*);
		//bus* wire_;

	public:
		ode(int rank,
		    float end,
		    float dt,
		    tensor<qnt>* inVal,
		    qnt (*funcPtr)(tensor<qnt>*)
		    //bus* wire
		    ) 

		    :
		    rank_(rank),
		    endTime_(end),
		    dT_(dt),
		    inits_(inVal),
		    current_(inits_),
		    rule_(funcPtr)
		    //wire_(wire)
		{}

		virtual void solveStep() = 0;
		
		virtual tensor<qnt> dot(tensor<qnt>& ModifiedCurrent)
		{
			tensor<qnt> res(rank_ + 2, 1);
			res(0,0)  = rule_(&ModifiedCurrent);

			for(int x = 1; x < rank_ + 2; x++)
				// TODO: cant modify the tensor externally. Tweak the operator.
				res(x, 0) = ModifiedCurrent[x - 1];
			return res;
		}
		
		virtual void solve()
		{
			for(uint i = 0; i < dT_; i++)
			{
				solveStep();
				//wire_ -> update(current_);
			}
		}
};
template<class qnt>
class rungeKutta4 : ode<qnt>
{
	void solveStep() override
	{
		tensor<qnt> a = this -> dot(*(this -> current_));
		tensor<qnt> b = this -> dot(*(this -> current_) +  1 /(this -> dT_ / 2) / a);
		tensor<qnt> c = this -> dot(*(this -> current_) +  1 /(this -> dT_ / 2) / b);
		tensor<qnt> d = this -> dot(*(this -> current_) +  1 /this ->dT_ / c);
		                                                            
		this -> current_ += this -> dT_ /6 * (a + 2 * b + 2 * c + d);
	}
};

template<class qnt>
class euler : ode<qnt>
{
	public:
		euler(int rank,
		    float end,
		    float dt,
		    tensor<qnt>* inVal,
		    qnt(*funcPtr)(tensor<qnt>*)
		    //bus* wire
		    ) 

		    :

		    ode<qnt>(rank,
			     end,
			     dt,
			     inVal,
			     funcPtr
			     //bus* wire
		    ) 
		    {}

	void solveStep() override
	{
		tensor<qnt>& currentRef = *(this -> current_);
		currentRef =  currentRef + this -> dot(*(this -> current_)) / (1/(this -> dT_));
	}
};
