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
		
		tensor<qnt> dot(const tensor<qnt> state)
		{
			return rule_(state);
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
		rungeKutta4(tensor<qnt> inVal,
				float dt,
				float end,
				tensor<qnt> (*funcPtr)(tensor<qnt>)) 
			:
			ode<qnt>(inVal,dt,end,funcPtr),
			a(this -> rank_, 1),
			b(this -> rank_, 1),
			c(this -> rank_, 1),
			d(this -> rank_, 1)
		    {
			    
		    }
		void solveStep() override
		{
			tensor<qnt>& s = *(this -> current_);
			const float dT = this -> dT_;
			a = this -> dot(s);
			b = this -> dot(s + a * dT/2);
			c = this -> dot(s + b * dT/2);
			d = this -> dot(s + c* dT);
			//a.print();
			//b.print();
			//c.print();
			//d.print();

			s +=  (a + b * 2 + c * 2 + d) * dT /6;
		}
	protected:



		tensor<qnt> a;
		tensor<qnt> b;
		tensor<qnt> c;
		tensor<qnt> d;

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
