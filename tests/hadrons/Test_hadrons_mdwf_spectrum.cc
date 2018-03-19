/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_spectrum.cc
 
 Copyright (C) 2015
 
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
 
 See the full license in the file "LICENSE" in the top level distribution
 directory.
 *******************************************************************************/

#include <Grid/Hadrons/Application.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> flavour = {"b"};
    std::vector<double>      mass    = {.019};
    std::vector<double>      tol     = {1e-11};
    std::string srcName;
    std::string lapName;
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1010;
    globalPar.trajCounter.end   = 1020;
    globalPar.trajCounter.step  = 10;
    globalPar.seed              = "1 2 3 4";
    application.setPar(globalPar);

    // gauge field
    application.createModule<MGauge::Random>("gauge");
    // MGauge::Load::Par LPar;
    // LPar.file = "myfile";
    // application.createModule<MGauge::Load>("gauge",LPar);

    // MGauge::LoadSmear::Par LPar;
    // LPar.file = "myfile";
    // LPar.Nsmear = 3;
    // LPar.rho = 0.1;
    // application.createModule<MGauge::LoadSmear>("gauge",LPar);

    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);

    // MSource::Z2:: Par z2Par;
    // z2Par.tA = 0;
    // z2Par.tB = 0;
    // srcName = "z2_";
    // application.createModule<MSource::Z2>(srcName, z2Par);

    // MSource::LaplaceSmearing::Par LapPar;
    // LapPar.N = 20;
    // LapPar.alpha = 0.1;
    // LapPar.source = srcName;
    // LapPar.gauge = "gauge";
    // lapName= "z2smr_";
    // application.createModule<MSource::LaplaceSmearing>(lapName,LapPar);
    
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
 
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 1";
    // set source - smear the b
    // std::vector<std::string> smearing = {"z2smr"}; 

    for (unsigned int i = 0; i < flavour.size() ; ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12;
        actionPar.M5    = 1.0;
        actionPar.mass  = mass[i];
	//actionPar.scale = 2.0;
        actionPar.boundary = boundary;
        application.createModule<MAction::DWF>("MDWF_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action   = "MDWF_" + flavour[i];
	solverPar.residual = tol[i];
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
	if (flavour[i]=="s") {
	    quarkPar.source = "pt";
	  }
	else
	  {
	    quarkPar.source = "pt";
	  }
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);
    }

    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::Meson::Par mesPar;
        
        mesPar.output  = "mesons/pt_" + flavour[i] + flavour[j];
        mesPar.q1      = "Qpt_" + flavour[i];
        mesPar.q2      = "Qpt_" + flavour[j];
        mesPar.gammas  = "all";
        mesPar.sink    = "sink";
        application.createModule<MContraction::Meson>("meson_pt_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);
	
    }

    // execution
    application.saveParameterFile("spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
