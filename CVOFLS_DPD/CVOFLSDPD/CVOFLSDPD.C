/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    CVOFLSDPD

Description
    Solver for 2 compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.
    The momentum and other fluid properties are of the "mixture" and single
    momentum,species i and energy equation are solved. 
    - IsoAdvector is used in order to advect the interface.
    - OpenSMOKE++ library for thermal and transport properties.
    - Evaporation rate calculated from diffusion flux
    - VOF solver coupled with LS 
    - CFD solver coupled with DPD
	- Runs in parallel

\*---------------------------------------------------------------------------*/


// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"
#include "math.h"
#include "math/PhysicalConstants.h"
#include "math/OpenSMOKEUtilities.h"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
//#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// ODE solvers
//#include "math/native-ode-solvers/MultiValueSolver"
//#include "math/external-ode-solvers/ODE_Parameters.h"

// Customized radiation model
//#include "OpenSMOKEradiationModel.H"

// OpenFOAM headers
//#include "meshSearch.H"
#include "isoAdvection.H"
#include "fvCFD.H"
#include "CMULES.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "myimmiscibleSemiCompressibleTwoPhaseMixtureSSF.H"
#include "fvIOoptionList.H"
#include "fvcSmooth.H"
#include "pimpleControl.H"
#include "IOMRFZoneList.H"
#include "liquidProperties.H"
#include "fixedFluxPressureFvPatchScalarField.H"
//#include "OFstream.H"

#include "cfdemCloud.H"
#include "implicitCouple.H"
#include "forceModel.H"
#include "smoothingModel.H"

//BzzLibraries
//#include "BzzMath.hpp"

// Homogeneous reactors for Operator Splitting Approach
//#include "DRG.h"
//#include "BatchReactorHomogeneousConstantPressure.H"
//#include "BatchReactorHomogeneousConstantPressure_ODE_Interface.H"
//#include "BatchReactorHomogeneousConstantVolume.H"
//#include "BatchReactorHomogeneousConstantVolume_ODE_Interface.H"

// Soot
//#include "sootUtilities.H"

// Utilities
#include "utilities.H"

//Spark
//#include "sparkModel.H"



// Parallel
#include "IPstream.H"
#include "OPstream.H"
#include <mpi.h>
#include "parallelClass.H"

// InterfaceClass
#include "interfaceSerialClass.H"
#include "interfaceParallelClass.H"



//template<typename Solver, typename OdeBatch>
//void SolveOpenSourceSolvers(OdeBatch& ode, const double t0, const double tf, const OpenSMOKE::OpenSMOKEVectorDouble& y0, OpenSMOKE::OpenSMOKEVectorDouble& yf, const OpenSMOKE::ODE_Parameters& parameters)
//{
//	Solver o(ode);
//	o.SetDimensions(y0.Size());
//	o.SetAbsoluteTolerance(parameters.absolute_tolerance());
//	o.SetRelativeTolerance(parameters.relative_tolerance());
//	o.SetAnalyticalJacobian(false);
//	o.SetInitialValues(t0, y0.GetHandle());
//	o.Solve(tf);
//	o.Solution(yf.GetHandle());
//} 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
//  #include "createControl.H"

    #include "createBasicFields.H"
    #include "readOptions.H"
    #include "createLiquidFields.H"
    #include "createGasFields.H"
    #include "createLiquidChemicalFields.H"
    #include "createGasChemicalFields.H" 
    #include "createChemicalFluxes.H"
//    #include "createDRGFields.H"
 // #include "createSurfaceTensionFields.H"

    #include "createFvOptions.H"
//    #include "createRadiationModel.H"
    #include "memoryAllocation.H"

    #include "GasProperties.H"  
    #include "LiquidProperties.H"
    pimpleControl pimple(mesh);
    #include "createMixtureFields.H"
    #include "createIsoAdvection.H"
    #include "createAdditionalFields.H"

    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        #include "compressibleCourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

       
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
   
	{
		// do particle stuff
		Info << "- evolve()" << endl;
		particleCloud.evolve(voidfraction,Us,U);
		Info << "Done_evolve()" << endl;

		//Ksl.internalField() = particleCloud.momCoupleM(0).impMomSource();
		//particleCloud.smoothingM().smoothen(Ksl);
		//Ksl.correctBoundaryConditions();

		//#include "solverDebugInfo.H"

		// get scalar source from DEM        
		particleCloud.forceM(1).manipulateScalarField(Tsource);
		Tsource.correctBoundaryConditions();

		// solve scalar transport equation
		//phi = fvc::interpolate(U*voidfraction) & mesh.Sf();

	}

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
             #include "updateEvaporationFlux.H" 
             #include "alphaEvaporation.H" 
             #include "alphaGeometricAdvection.H" 
             #include "GasProperties.H"                              
             #include "LiquidProperties.H"         
             #include "updateAlphaTransportProperties.H"
             #include "DensityEqn.H"     
             #include "UEqn.H"

             #include "TEqn.H"
             #include "YiEqn.H"


            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

           #include "solveAdvectionVariables.H"
           #include "conserveMass.H"
	   	   #include "solveLSFunction.H"
		   #include "calcDdtU.H"

        }


        runTime.write();
     

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
