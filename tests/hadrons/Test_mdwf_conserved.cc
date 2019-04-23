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
#include <Grid/Hadrons/Modules.hpp>

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
    std::vector<std::string> flavour = {"l"};
    std::vector<double>      mass    = {.015  };
    std::vector<double>      tol     = {1e-11 };
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.seed              = "1 2 3 4";
    application.setPar(globalPar);

    // gauge field
    // MGauge::LoadSmear::Par LPar;
    // LPar.file="myfile";
    // LPar.Nsmear=3;
    // LPar.rho=0.1;
    
    //application.createModule<MGauge::LoadSmear>("gauge",LPar);
    application.createModule<MGauge::Unit>("gauge");
    
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);

    // sources
    // MSource::Z2::Par z2Par;
    // z2Par.tA = 0;
    // z2Par.tB = 0;
    // application.createModule<MSource::Z2>("z2", z2Par);
    
    // MSource::LaplacianSmear::Par lapPar;
    // lapPar.gauge = "gauge";
    // lapPar.source = "pt";
    // lapPar.alpha = 20.;
    // lapPar.n = 200;
    // lapPar.spatial = true; // 3d smearing 
    // application.createModule<MSource::LaplacianSmear>("lap", lapPar);
    
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 1";

    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::MobiusDWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12; 
        actionPar.M5    = 1.0;
        actionPar.mass  = mass[i];
	actionPar.scale = 2.0;
        actionPar.boundary = boundary;
        application.createModule<MAction::MobiusDWF>("MDWF_" + flavour[i], actionPar);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        //solverPar.action       = "DWF_" + flavour[i];
        solverPar.residual     = tol[i];
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i], solverPar);
        
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.action = "MDWF_" + flavour[i];
        quarkPar.solver = "CG_" + flavour[i];

        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i], quarkPar);

        // quarkPar.source = "lap";
        // application.createModule<MFermion::GaugeProp>("Qlap_" + flavour[i], quarkPar);

	// // sink smearing: simply apply smearing to solved propagator
	// lapPar.source = "Qlap_" + flavour[i] ;
	// application.createModule<MSource::LaplacianSmear>("Qsm_"+flavour[i], lapPar);
    }


    for (unsigned int i = 0; i < flavour.size(); ++i)
    for (unsigned int j = i; j < flavour.size(); ++j)
    {
        MContraction::WardIdentity::Par mesPar;
	MContraction::Meson::Par mesParContract;
        
        //mesPar.output  = "mesons/pt_" + flavour[i] + flavour[j];
	mesPar.q      = "Qpt_" + flavour[i]+"_5d";
	mesPar.action = "MDWF_" + flavour[i];
	mesPar.mass = mass[i];
	mesPar.test_axial = true;
        //mesPar.q1      = "Qpt_" + flavour[j];
        //mesPar.gammas  = "all";
        //mesPar.sink    = "sink";
        application.createModule<MContraction::WardIdentity>("meson_wi_"
                                                      + flavour[i] + flavour[j],
                                                      mesPar);


        mesParContract.output  = "mesons/pt_" + flavour[i] + flavour[j];
	mesParContract.q1      = "Qpt_" + flavour[i];
	mesParContract.q2      = "Qpt_" + flavour[j];
	mesParContract.gammas = "all";
	mesParContract.sink = "sink";
        application.createModule<MContraction::Meson>("meson_ll_"
                                                      + flavour[i] + flavour[j],
                                                      mesParContract);

        // mesPar.output  = "mesons/ss_" + flavour[i] + flavour[j];
        // mesPar.q1      = "Qpt_" + flavour[i];
        // mesPar.q2      = "Qsm_" + flavour[j];
        // mesPar.gammas  = "all";
        // mesPar.sink    = "sink";
        // application.createModule<MContraction::Meson>("meson_ss_"
        //                                               + flavour[i] + flavour[j],
        //                                               mesPar);

    }
    
    // execution
    application.saveParameterFile("spectrum.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}