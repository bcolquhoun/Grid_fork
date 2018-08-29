/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_meson_3pt.cc
 
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
    std::vector<std::string> flavour = {"b"};
    std::vector<double>      mass    = {.68808  };
    std::vector<double>      lmass   = {.5  };
    unsigned int             nt      = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.seed                 = "1 2 3 4";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);
    
    // gauge field
    MGauge::LoadSmear::Par LPar;
    LPar.file="myfile";
    LPar.Nsmear=3;
    LPar.rho=0.1;
    
    application.createModule<MGauge::LoadSmear>("gauge", LPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 1";

    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "1. -1. 0.";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
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
        application.createModule<MAction::MobiusDWF>("MobiusDWF_" + flavour[i], actionPar);

        // actions
        MAction::MobiusDWF::Par actionPar_light;
        actionPar_light.gauge = "gauge";
        actionPar_light.Ls    = 12;
        actionPar_light.M5    = 1.0;
        actionPar_light.mass  = lmass[i];
	actionPar_light.scale = 2.0;
        actionPar_light.boundary = boundary;
        application.createModule<MAction::MobiusDWF>("MobiusDWF_light_" + flavour[i], actionPar_light);
        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.residual     = 1.0e-16;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
    }
    for (unsigned int t = 0; t < 1; t += 1)
    {
        std::string                           srcName;
	std::string                           smrSrcName;
	std::string                           smrSnkName;
        std::vector<std::string>              qName;
	std::string                           ptName;
        std::vector<std::vector<std::string>> seqName;

	// Point Source
	MSource::Point::Par pointPar;
	pointPar.position = "0 0 0 0";
	srcName = "pt_" + std::to_string(t);
	application.createModule<MSource::Point>(srcName, pointPar);


	// Set up quark smearing parameters
	MSource::LaplacianSmear::Par lapPar;
	lapPar.gauge = "gauge";
	lapPar.source = srcName;
	lapPar.alpha = 20.;
	lapPar.n = 200;
	lapPar.spatial = true;
	
	// smear source
	smrSrcName = "smrSrc_" + std::to_string(t) ;
	application.createModule<MSource::LaplacianSmear>(smrSrcName,lapPar);
        
        // Z2 source
        // MSource::Z2::Par z2Par;
        // z2Par.tA = t;
        // z2Par.tB = t;
        // srcName  = "z2_" + std::to_string(t);
        // application.createModule<MSource::Z2>(srcName, z2Par);


        for (unsigned int i = 0; i < flavour.size(); ++i)
        {
	  
	  // propagators
	  MFermion::GaugeProp::Par quarkPar;
	  qName.push_back("Qpt_" + flavour[i] + "_" + std::to_string(t));
	  // using light propagator action
	  quarkPar.action = "MobiusDWF_light_" + flavour[i];
	  quarkPar.solver = "CG_" + flavour[i];
	  quarkPar.source = smrSrcName ;

	  // propagator form smeared source
	  application.createModule<MFermion::GaugeProp>(qName[i], quarkPar);
	  ptName = "ptName";

	  // smear propagator at sink
	  lapPar.source = qName[i];
	  smrSnkName = "smrSnk" ;
	  application.createModule<MSource::LaplacianSmear>(smrSnkName,lapPar);
	  
	  // sequential sources
	  MSource::SeqGamma::Par seqPar;
	  // create from smeared-smeared propagator
	  seqPar.q   = smrSnkName;
	  seqPar.tA  = 28 % nt;
	  seqPar.tB  = 28 % nt;
	  seqPar.mom = "0. 0. 0. 0.";
	  seqName.push_back(std::vector<std::string>(Nd));
	  
	  for (unsigned int mu = 0; mu < 1; ++mu)
            {
	      // hex for gamma (0x0 is -gamma5 to agree with IroIro +gamma5)
	      seqPar.gamma = 0x0; 
	      seqName[i][mu] = "G" + std::to_string(seqPar.gamma)
	    	+ "_" + std::to_string(seqPar.tA) + "-"
	    	+ qName[i];
	      application.createModule<MSource::SeqGamma>(seqName[i][mu], seqPar);
            }
	  
	  // create propagator from sequential source with heavier quark
	  for (unsigned int mu = 0; mu < 1; ++mu)
            {
	      quarkPar.source = seqName[i][mu];
	      quarkPar.action = "MobiusDWF_" + flavour[i] ;
	      seqName[i][mu]  = "Q_" + flavour[i] + "-" + seqName[i][mu];
	      application.createModule<MFermion::GaugeProp>(seqName[i][mu], quarkPar);
            }

	  // point propagator 
	  quarkPar.action = "MobiusDWF_light_" + flavour[i] ;
	  quarkPar.source =  srcName;
	  application.createModule<MFermion::GaugeProp>(ptName, quarkPar);
        }
	        
        // contractions
        MContraction::Meson::Par mesPar;
        for (unsigned int i = 0; i < flavour.size(); ++i)
	  {
	    for (unsigned int j = i; j < flavour.size(); ++j)
	      {
		// smeared-local meson
		mesPar.output = "3pt/sl_" + flavour[i] + flavour[j];
		mesPar.q2     = ptName;
		mesPar.q1     = qName[j];
		mesPar.gammas = "all";
		mesPar.sink   = "sink";
		application.createModule<MContraction::Meson>("meson_sl_"
							      + std::to_string(t)
							      + "_"
							      + flavour[i]
							      + flavour[j],
							      mesPar);

		// smeared-smeared meson
		mesPar.output = "3pt/ss_" + flavour[i] + flavour[j];
		mesPar.q1     = smrSnkName;
		mesPar.q2     = ptName;
		mesPar.gammas = "all";
		mesPar.sink   = "sink";
		application.createModule<MContraction::Meson>("meson_ss_"
							      + std::to_string(t)
							      + "_"
							      + flavour[i]
							      + flavour[j],
							      mesPar);


	      }
	  }

	for (unsigned int i = 0; i < flavour.size(); ++i)
	  {
	    for (unsigned int j = 0; j < flavour.size(); ++j)
	      {
		for (unsigned int mu = 0; mu < 1; ++mu)
		  {
		    MContraction::Meson::Par mesPar;
		    
		    mesPar.output = "3pt/3pt_ss_" + flavour[i] + flavour[j] + "_"
		      + std::to_string(mu);
		    mesPar.q1     = seqName[j][mu];
		    mesPar.q2     = ptName;
		    mesPar.gammas = "all";
		    mesPar.sink   = "sink";
		    application.createModule<MContraction::Meson>("3pt_pt_"
								  + std::to_string(t)
								+ "_"
								  + flavour[i]
								  + flavour[j]
								  + "_"
								  + std::to_string(mu),
								  mesPar);
		}
	    }
	  }
    }    
    // execution
    application.saveParameterFile("meson3pt.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
