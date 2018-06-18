/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MAction/MobiusDWF.hpp

Copyright (C) 2015-2018

Author: Brian Colquhoun <brian.colquhoun@kek.jp>
based on ZMobiusDWF Module by
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
#ifndef Hadrons_MAction_MobiusDWF_hpp_
#define Hadrons_MAction_MobiusDWF_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         MobiusDWF                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class MobiusDWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MobiusDWFPar,
                                    std::string                      , gauge,
                                    unsigned int                     , Ls,
                                    double                           , mass,
                                    double                           , M5,
                                    double                           , scale,
                                    std::string                      , boundary);
};

template <typename FImpl>
class TMobiusDWF: public Module<MobiusDWFPar>
{
public:
    FGS_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TMobiusDWF(const std::string name);
    // destructor
    virtual ~TMobiusDWF(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(MobiusDWF, TMobiusDWF<FIMPL>, MAction);

/******************************************************************************
 *                 TMobiusDWF implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMobiusDWF<FImpl>::TMobiusDWF(const std::string name)
: Module<MobiusDWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMobiusDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TMobiusDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMobiusDWF<FImpl>::setup(void)
{
    LOG(Message) << "Setting up Mobius domain wall fermion with m= "
                 << par().mass << ", M5= " << par().M5 << ", Ls= " << par().Ls 
                 << ", and scale= " << par().scale
                 << " using gauge field '" << par().gauge << "'"
                 << std::endl;
    LOG(Message) << "Fermion boundary conditions: " << par().boundary
                 << std::endl;

    env().createGrid(par().Ls);
    auto &U    = envGet(LatticeGaugeField, par().gauge);
    auto &g4   = *env().getGrid();
    auto &grb4 = *env().getRbGrid();
    auto &g5   = *env().getGrid(par().Ls);
    auto &grb5 = *env().getRbGrid(par().Ls);
    std::vector<Complex> boundary = strToVec<Complex>(par().boundary);
    typename MobiusFermion<FImpl>::ImplParams implParams(boundary);
    envCreateDerived(FMat, ScaledShamirFermion<FImpl>, getName(), par().Ls, U, g5,
                     grb5, g4, grb4, par().mass, par().M5,
                     par().scale, implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMobiusDWF<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_MobiusDWF_hpp_
