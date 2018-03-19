/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSink/Point.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MSink_Point_hpp_
#define Hadrons_MSink_Point_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                   Point                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSink)

class PointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PointPar,
                                    std::string, mom);
};

template <typename FImpl>
class TPoint: public Module<PointPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SINK_TYPE_ALIASES();
public:
    // constructor
    TPoint(const std::string name);
    // destructor
    virtual ~TPoint(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasPhase_{false}; 
    std::string momphName_;
};

MODULE_REGISTER_NS(Point,       TPoint<FIMPL>,        MSink);
MODULE_REGISTER_NS(ScalarPoint, TPoint<ScalarImplCR>, MSink);

/******************************************************************************
 *                          TPoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPoint<FImpl>::TPoint(const std::string name)
: Module<PointPar>(name)
, momphName_ (name + "_momph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPoint<FImpl>::setup(void)
{
    envTmpLat(LatticeComplex, "coor");
    envCacheLat(LatticeComplex, momphName_);
    envCreate(SinkFn, getName(), 1, nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPoint<FImpl>::execute(void)
{   
    LOG(Message) << "Setting up point sink function for momentum ["
                 << par().mom << "]" << std::endl;

    auto &ph = envGet(LatticeComplex, momphName_);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
	std::vector<LatticeComplex> saved_coor;
        ph = zero;
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
	    std::cout << "mu: " << mu << "  p[mu]: " << p[mu] << std::endl ;
	    saved_coor.push_back(coor) ;
	    //std::cout << "coor: " << coor << std::endl ; 
	    std::cout << "fdimensions[mu]: " << env().getGrid()->_fdimensions[mu] << std::endl ; 
            ph = ph + (p[mu]/env().getGrid()->_fdimensions[mu])*coor;
        }
	std::cout << "ph (before exp): " << ph << std::endl ; 
        ph = exp((Real)(2*M_PI)*i*ph);
	std::cout << "saved_coor[0] " << saved_coor[0] << std::endl ; 
	std::cout << "saved_coor[1] " << saved_coor[1] << std::endl ; 
	std::cout << "saved_coor[2] " << saved_coor[2] << std::endl ; 
	std::cout << "saved_coor[3] " << saved_coor[3] << std::endl ; 

        hasPhase_ = true;
    }
    auto sink = [&ph](const PropagatorField &field)
    {
        SlicedPropagator res;
	std::cout << "ph: " << ph << std::endl ; 
        PropagatorField  tmp = ph*field;
        
        sliceSum(tmp, res, Tp);
        
        return res;
    };
    envGet(SinkFn, getName()) = sink;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSink_Point_hpp_
