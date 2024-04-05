#include "../modules/tensor.h"
#include "../modules/uint.h"

template <class qnt>
class ode
{
	protected: 
		float endTime_;
		float dT_;
		tensor<qnt> inits_;
		tensor<qnt> (*rule_)(tensor<qnt>);
		tensor<qnt>* current_;
		uint rank_;
		//bus* wire_;

	public:
		ode(tensor<qnt> inVal,
				float dt,
				float end,
				tensor<qnt> (*funcPtr)(tensor<qnt>)
				) 
		    :
		    endTime_(end),
		    dT_(dt),
		    inits_(inVal),
		    rule_(funcPtr),
		    current_(&inits_),
		    rank_ (inits_.size())
		{
				
		}

		virtual void solveStep() = 0;
		
		tensor<qnt> dot(tensor<qnt>& state)
		{
			return rule_(*current_);
		}
		
		void solve()
		{
			for(float i = 0; i < endTime_; i += dT_)
			{
				solveStep();
				current_ -> print();
				//wire_ -> update(current_);
			}
		}
};
template<class qnt>
class rungeKutta4 : public ode<qnt>
{
	public:
	void solveStep() override
	{
		tensor<qnt> a = this -> dot(*(this -> current_));
		tensor<qnt> b = this -> dot(*(this -> current_) +  1 /(this -> dT_ / 2) / a);
		tensor<qnt> c = this -> dot(*(this -> current_) +  1 /(this -> dT_ / 2) / b);
		tensor<qnt> d = this -> dot(*(this -> current_) +  1 /this -> dT_ / c);
		                                                            
		this -> current_ += this -> dT_ /6 * (a + 2 * b + 2 * c + d);
	}
};

template<class qnt>
class euler : public ode<qnt>
{
	public:
		euler(tensor<qnt> inVal,
				float dt,
				float end,
				tensor<qnt> (*funcPtr)(tensor<qnt>)) 
			:
			ode<qnt>(inVal,dt,end,funcPtr) 
		    {
			    
		    }

	void solveStep() override
	{
		tensor<qnt>& currentRef = *(this -> current_);
		currentRef =  currentRef + this -> dot(currentRef) * (this -> dT_);
	}
};
