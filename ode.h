#include "tensor.h"
#include <math>

template <class qnt>
class ode
{
	protected: 
		uint rank_;
		float endTime_;
		float dT_;
		qnt* inits_;
		qnt* current_;
		qnt  (*rule_)(tensor<qnt>);
		bus* wire_;

	public:
		ode(uint rank,
		    float end,
		    float dt,
		    qnt* inVal,
		    tensor<qnt> (*funcPtr)(tensor<qnt>),
		    bus* wire) 

		    :
		    rank_(rank),
		    endTime_(end),
		    dT_(dt),
		    inits_(inVal),
		    current_(inits_),
		    wire_(wire)
		{}

		virtual void solveStep() = 0;
		
		virtual tensor<qnt> dot()
		{
			tensor<qnt> res(rank_ + 2, 1);
			res[0]  = rule_(current);

			for(int x = 1; x < rank_ + 2; x++)
				res[x] = current_[x - 1];
			return res;
		}
		
		virtual void solve()
		{
			for(uint i = 0; i < dT_; i++)
			{
				solveStep();
				wire_ -> update(current_);
			}
		}
};

template <class qnt>
class rungeKutta4 : public ode
{
	void solveStep() override
	{
		tensor<qnt> a = dot(current_);
		tensor<qnt> b = dot(current_ + dT_ / 2 * a);
		tensor<qnt> c = dot(current_ + dT_ / 2 * b);
		tensor<qnt> d = dot(current_ + dT_ * c);
		
		current_ += h/6 * (a + 2 * b + 2 * c + d);
	}
};

class euler : public ode
{

	void solveStep() override
	{
		current_ += dT_ * dot(current_);
	}
};
