/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSource/LaplacianSmear.hpp

Copyright (C) 2015-2018

Author: Guido Cossu <guido.cossu@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MSource_LaplacianSmear_hpp_
#define Hadrons_MSource_LaplacianSmear_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Laplacian smearing of an input source
 ------------
 * src_out = (1 + alpha Laplacian)^n * src_in 
 * where Laplacian is the laplacian operator 

 * options:
 - alpha: smearing strength
 - n: smearing levels
 
 */

/******************************************************************************
 *                                  TPoint                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class LaplacianSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LaplacianSmearPar,
                                    std::string, gauge,
                                    std::string, source,
                                    double, alpha,
                                    unsigned int, n,
                                    bool, spatial);
};

template <typename FImpl>
class TLaplacianSmear: public Module<LaplacianSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLaplacianSmear(const std::string name);
    // destructor
    virtual ~TLaplacianSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LaplacianSmear,       TLaplacianSmear<FIMPL>,        MSource);


/******************************************************************************
 *                       TPoint template implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLaplacianSmear<FImpl>::TLaplacianSmear(const std::string name)
: Module<LaplacianSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLaplacianSmear<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge, par().source};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLaplacianSmear<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLaplacianSmear<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envTmpLat(FermionField, "src_sc");
    envTmpLat(FermionField, "tmp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLaplacianSmear<FImpl>::execute(void)
{
    LOG(Message) << "Laplacian smearing source named "<< par().source << " with parameters alpha: " << par().alpha
                 << " and n:" << par().n << std::endl;

    auto        &U            = envGet(LatticeGaugeField, par().gauge);
    auto        &src          = envGet(PropagatorField, par().source);
    auto        &smeared_src  = envGet(PropagatorField, getName());

    envGetTmp(FermionField, tmp);
    envGetTmp(FermionField, src_sc);

    LaplacianFundamental<FImpl> lapl(U._grid, par().spatial);
    lapl.ImportGauge(U);

    smeared_src = zero;

    for (unsigned int s = 0; s < Ns; ++s)
      for (unsigned int c = 0; c < FImpl::Dimension; ++c){
        PropToFerm<FImpl>(src_sc, src, s, c);
        for (uint n = 0; n < par().n; n++){
            lapl.M(src_sc, tmp );
            src_sc += par().alpha * tmp; 
        }
        FermToProp<FImpl>(smeared_src, src_sc, s, c);
      }

    

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_LaplacianSmear_hpp_
