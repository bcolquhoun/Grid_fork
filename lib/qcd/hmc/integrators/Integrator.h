//--------------------------------------------------------------------
/*! @file Integrator.h
 * @brief Classes for the Molecular Dynamics integrator
 *
 * @author Guido Cossu
 * Time-stamp: <2015-07-30 16:21:29 neo>
 */
//--------------------------------------------------------------------

#ifndef INTEGRATOR_INCLUDED
#define INTEGRATOR_INCLUDED

//class Observer;

#include <memory>

namespace Grid{
  namespace QCD{

    struct IntegratorParameters{

      int Nexp;
      int MDsteps;  // number of outer steps
      RealD trajL;  // trajectory length 
      RealD stepsize;

      IntegratorParameters(int Nexp_,
			   int MDsteps_, 
			   RealD trajL_):
        Nexp(Nexp_),
	MDsteps(MDsteps_),
	trajL(trajL_),
	stepsize(trajL/MDsteps)
        {
	  // empty body constructor
	};

    };

    // Should match any legal (SU(n)) gauge field
    // Need to use this template to match Ncol to pass to SU<N> class
    template<int Ncol,class vec> void generate_momenta(Lattice< iVector< iScalar< iMatrix<vec,Ncol> >, Nd> > & P,GridParallelRNG& pRNG){
      typedef Lattice< iScalar< iScalar< iMatrix<vec,Ncol> > > > GaugeLinkField;
      GaugeLinkField Pmu(P._grid);
      Pmu = zero;
      for(int mu=0;mu<Nd;mu++){
	SU<Ncol>::GaussianLieAlgebraMatrix(pRNG, Pmu);
	PokeIndex<LorentzIndex>(P, Pmu, mu);
      }
    }

    template<class GaugeField> struct ActionLevel{
      public:
      
	typedef Action<GaugeField>*  ActPtr; // now force the same colours as the rest of the code

	int multiplier;

	std::vector<ActPtr> actions;

        ActionLevel(int mul = 1) : multiplier(mul) {
	  assert (mul > 0);
	};

	void push_back(ActPtr ptr){
	  actions.push_back(ptr);
	}
    };

    template<class GaugeField> using ActionSet = std::vector<ActionLevel< GaugeField > >;

    /*! @brief Class for Molecular Dynamics management */   
    template<class GaugeField, class Algorithm >
    class Integrator : public Algorithm {

    private:

      IntegratorParameters Params;

      GridParallelRNG pRNG;        // Store this somewhere more sensible and pass as reference
      const ActionSet<GaugeField> as;

      int levels;              //
      double t_U;              // Track time passing on each level and for U and for P
      std::vector<double> t_P; //

      GaugeField P;  // is a pointer really necessary?

      //ObserverList observers; // not yet
      //      typedef std::vector<Observer*> ObserverList;
      //      void register_observers();
      //      void notify_observers();


      void update_P(GaugeField&U, int level,double ep){
	t_P[level]+=ep;
	for(int a=0; a<as[level].actions.size(); ++a){
	  GaugeField force(U._grid);
	  as[level].actions.at(a)->deriv(U,force);
	  P = P - force*ep;
	}

	std::cout<<GridLogMessage;
	for(int l=0; l<level;++l) std::cout<<"   ";	    
	std::cout<<"["<<level<<"] P " << " dt "<< ep <<" : t_P "<< t_P[level] <<std::endl;

      }

      void update_U(GaugeField&U, double ep){
	//rewrite exponential to deal automatically  with the lorentz index?
	//	GaugeLinkField Umu(U._grid);
	//	GaugeLinkField Pmu(U._grid);
	for (int mu = 0; mu < Nd; mu++){
	  auto Umu=PeekIndex<LorentzIndex>(U, mu);
	  auto Pmu=PeekIndex<LorentzIndex>(P, mu);
	  Umu = expMat(Pmu, ep, Params.Nexp)*Umu;
	  PokeIndex<LorentzIndex>(U, Umu, mu);
	}

	t_U+=ep;
	int fl = levels-1;
	std::cout<<GridLogMessage<<"   ";
	for(int l=0; l<fl;++l) std::cout<<"   ";	    
	std::cout<<"["<<fl<<"] U " << " dt "<< ep <<" : t_U "<< t_U <<std::endl;

      }
      
      friend void Algorithm::step (GaugeField& U, 
				   int level, 
				   std::vector<int>& clock,
				   Integrator<GaugeField,Algorithm>* Integ);

    public:

      Integrator(GridBase* grid, 
		 IntegratorParameters Par,
		 ActionSet<GaugeField> & Aset):
          Params(Par),
    	  as(Aset),
	  P(grid),
	  pRNG(grid),
	  levels(Aset.size())
      {
	std::vector<int> seeds({1,2,3,4,5}); // Fixme; Pass it the RNG as a ref
	pRNG.SeedFixedIntegers(seeds);

	t_P.resize(levels,0.0);
	t_U=0.0;
      };
      
      ~Integrator(){}

      //Initialization of momenta and actions
      void init(GaugeField& U){
	std::cout<<GridLogMessage<< "Integrator init\n";
	generate_momenta(P,pRNG);
	for(int level=0; level< as.size(); ++level){
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    as[level].actions.at(actionID)->init(U, pRNG);
	  }
	}
      }

      // Calculate action
      RealD S(GaugeField& U){

	LatticeComplex Hloc(U._grid);	Hloc = zero;
	// Momenta
	for (int mu=0; mu <Nd; mu++){
	  auto Pmu = PeekIndex<LorentzIndex>(P, mu);
	  Hloc -= trace(Pmu*Pmu);
	}
	Complex Hsum = sum(Hloc);
	
	RealD H = Hsum.real();
	RealD Hterm;
	std::cout<<GridLogMessage << "Momentum action H_p = "<< H << "\n";

	// Actions
	for(int level=0; level<as.size(); ++level){
	  for(int actionID=0; actionID<as[level].actions.size(); ++actionID){
	    Hterm = as[level].actions.at(actionID)->S(U);
	    std::cout<<GridLogMessage << "Level "<<level<<" term "<<actionID<<" H = "<<Hterm<<std::endl;
	    H += Hterm;
	  }
	}
	std::cout<<GridLogMessage << "Total action H = "<< H << "\n";
	
	return H;
      }

      void integrate(GaugeField& U){

	std::vector<int> clock;

	clock.resize(as.size(),0);

	// All the clock stuff is removed if we pass first, last to the step down the way
	for(int step=0; step< Params.MDsteps; ++step){   // MD step
	  Algorithm::step(U,0,clock, (this));
	}

	// Check the clocks all match
	for(int level=0; level<as.size(); ++level){
	  assert(fabs(t_U - t_P[level])<1.0e-4);
	  std::cout<<GridLogMessage<<" times["<<level<<"]= "<<t_P[level]<< " " << t_U <<std::endl;
	}	

      }
    };
    
  }
}
#endif//INTEGRATOR_INCLUDED
