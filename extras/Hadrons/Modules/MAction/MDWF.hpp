/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MAction/MDWF.hpp

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>

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

#ifndef Hadrons_MAction_MDWF_hpp_
#define Hadrons_MAction_MDWF_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Domain wall quark action                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class MDWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MDWFPar,
                                    std::string, gauge,
                                    unsigned int, Ls,
                                    double      , mass,
                                    double      , M5,
				    double      , scale,
                                    std::string , boundary);
};

template <typename FImpl>
class TMDWF: public Module<MDWFPar>
{
public:
    FGS_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TMDWF(const std::string name);
    // destructor
    virtual ~TMDWF(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(MDWF, TMDWF<FIMPL>, MAction);

/******************************************************************************
 *                        MDWF template implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMDWF<FImpl>::TMDWF(const std::string name)
: Module<MDWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TMDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMDWF<FImpl>::setup(void)
{
    unsigned int size;
    
    size = 2*env().template lattice4dSize<typename FImpl::DoubledGaugeField>();
    env().registerObject(getName(), size, par().Ls);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMDWF<FImpl>::execute(void)
{
    LOG(Message) << "Setting up scaled Shamir fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << ", Ls= "
                 << par().Ls << " and scale=" << par().scale 
		 << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    LOG(Message) << "Fermion boundary conditions: " << par().boundary 
                 << std::endl;
    env().createGrid(par().Ls);
    auto &U      = *env().template getObject<LatticeGaugeField>(par().gauge);
    auto &g4     = *env().getGrid();
    auto &grb4   = *env().getRbGrid();
    auto &g5     = *env().getGrid(par().Ls);
    auto &grb5   = *env().getRbGrid(par().Ls);
    std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
    typename ScaledShamirFermion<FImpl>::ImplParams implParams(boundary);
    FMat *fMatPt = new ScaledShamirFermion<FImpl>(U, g5, grb5, g4, grb4,
						  par().mass, par().M5,
						  par().scale,
						  implParams);
    env().setObject(getName(), fMatPt);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_MDWF_hpp_
