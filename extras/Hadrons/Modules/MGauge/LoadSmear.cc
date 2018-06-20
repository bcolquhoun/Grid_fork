/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MGauge/LoadSmear.cc

Copyright (C) 2015-2018

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

#include <Grid/Hadrons/Modules/MGauge/LoadSmear.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGauge;

/******************************************************************************
*                            TLoadSmear implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
TLoadSmear::TLoadSmear(const std::string name)
: Module<LoadSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
std::vector<std::string> TLoadSmear::getInput(void)
{
     std::vector<std::string> in;
     return in;
}

std::vector<std::string> TLoadSmear::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
void TLoadSmear::setup(void)
{
    envCreateLat(LatticeGaugeField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
void TLoadSmear::execute(void)
{
  FieldMetaData header;
  std::string fileName = par().file + "."
    + std::to_string(vm().getTrajectory());

  GridCartesian *Ugrid = env().getGrid();
  LatticeGaugeField U(Ugrid);

  LOG(Message) << "Reading from a file...." << std::endl;
  
  //  auto &U = envGet(LatticeGaugeField, getName());
  //SU3::ColdConfiguration(*env().get4dRng(), U);
  
  NerscIO::readConfiguration(U, header, fileName);
  LOG(Message) << "NERSC header:" << std::endl;
  dump_meta_data(header, LOG(Message));

  RealD plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(U);
  LOG(Message) << "Average thin plaquette: " << plaq << std::endl ;

  // Smearing, assuming PeriodicGauge Implementation policy                                                                                                                                                      
  LOG(Message) << "Smearing configuration with Nsmear = " << par().Nsmear << " rho = " << par().rho
	       << std::endl;
  Smear_Stout<PeriodicGimplR> Stout(par().rho);
  SmearedConfiguration<PeriodicGimplR> SmearingPolicy(Ugrid, par().Nsmear, Stout);
  SmearingPolicy.set_Field(U);

  // Problem with the following line. 
  LatticeGaugeField &SmearedU = *env().getObject<LatticeGaugeField>(getName());
  // This next line was an attempt at an update but isn't correct
  //LatticeGaugeField & SmearedU = *envCreateLat(LatticeGaugeField, getName());


  SmearedU = SmearingPolicy.get_SmearedU();
  RealD sm_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(SmearedU);
  LOG(Message) << "Average Stout smeared plaquette: " << sm_plaq << std::endl;

}