#ifndef Hadrons_MSource_FermionSmear_hpp_
#define Hadrons_MSource_FermionSmear_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         FermionSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class FermionSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionSmearPar,
                                    std::string, source,
				    std::string, gauge,
				    unsigned int, N,
				    double, alpha);
};

template <typename FImpl>
class TFermionSmear: public Module<FermionSmearPar>
{
public:
  FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TFermionSmear(const std::string name);
    // destructor
    virtual ~TFermionSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionSmear, TFermionSmear<FIMPL>, MSource);

/******************************************************************************
 *                 TFermionSmear implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TFermionSmear<FImpl>::TFermionSmear(const std::string name)
: Module<FermionSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TFermionSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TFermionSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFermionSmear<FImpl>::setup(void)
{
  envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TFermionSmear<FImpl>::execute(void)
{
  envGetTmp(FermionField, source);
  envGetTmp(FermionField, tmp);
  auto &SmrSrc = envGet(PropagatorField,getName());
  auto &fullSrc = envGet(PropagatorField,par().source);
  auto &U = envGet(LatticeGaugeField, par().gauge);
  //LaplacianField<FImpl> LaplaceOperator(env().getGrid());
  //LaplaceOperator.ImportGauge(U);
  double prefactor = par().alpha/(double)(par().N);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_FermionSmear_hpp_
