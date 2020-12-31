/*

// CALCULATION OF NORMAL VECTORS TO THE INTERFACE
dimensionedScalar  epsilon
(
  dimensionedScalar("epsilon",dimless/dimLength, 1e-12)
);

gradalpha = fvc::grad(alphasmooth_aux);
normal = gradalpha/( mag(gradalpha) + epsilon);



// CALCULATION OF TANGENT VECTORS TO THE INTERFACE
forAll (CurvatureLSM.internalField(), celli)
{

  tangent.internalField()[celli].x() = normal.internalField()[celli].y();
  tangent.internalField()[celli].y() = -normal.internalField()[celli].x();

}

tangent.correctBoundaryConditions();



// CALCULATION OF CURVILINEAR COORDINATE ALONG THE INTERFACE

CurvCoord *= 0.0;
CurvCoord.correctBoundaryConditions();

   fvScalarMatrix CurvCoordEqn
   (
     fvm::laplacian(CurvCoord) == fvc::div(tangent)
   );


    CurvCoordEqn.relax();
    CurvCoordEqn.solve();


// FILTERING AND CLEANING THE CURVILINEAR COORDINATE 

std::vector<double> CurvCoordVector;
int NcellsCurvCoord = 0.0;

forAll (CurvatureLSM.internalField(), celli)
{
  if (alpha1smooth.internalField()[celli] > 0.4 && alpha1smooth.internalField()[celli] < 0.6)
  {
      NcellsCurvCoord++;
  }

  else if (alpha1smooth.internalField()[celli] < 0.2 || alpha1smooth.internalField()[celli] > 0.8)
  {
      CurvCoord.internalField()[celli] = 0.0;
  }

}

CurvCoord.correctBoundaryConditions();

CurvCoordVector.resize(NcellsCurvCoord);

NcellsCurvCoord = 0.0;

 forAll (CurvatureLSM.internalField(), celli)
 {

  if (alpha1smooth.internalField()[celli] > 0.4 && alpha1smooth.internalField()[celli] < 0.6)
  {
    NcellsCurvCoord++;
    CurvCoordVector[NcellsCurvCoord-1] = CurvCoord.internalField()[celli];
  }

 }

double CurvCoordMaximo = *max_element(CurvCoordVector.begin(), CurvCoordVector.end());
double CurvCoordMinimo = *min_element(CurvCoordVector.begin(), CurvCoordVector.end());


//#include "RegressionInterfaceReconstruction.H";
//#include "SplineInterfaceReconstruction.H";
//#include "SplineInterfaceReconstructionImplicit.H";
#include "BezierInterfaceReconstruction.H";


//CURVATURE
//CurvatureLSM =  -fvc::div(normal);


//#include "filterCurvature.H";


alpha1pc = (1./(1.-0.99))*(min(max(alpha1, 0.99/2.), 1.-0.99/2.)-0.99/2.);

//fcf = fvc::interpolate(surfaceTension* 4000./dimChange)*fvc::snGrad(alpha1pc);
fcf = fvc::interpolate(surfaceTension*CurvatureLSM)*fvc::snGrad(alpha1pc);



/*
TensionForce = fvc::div( surfaceTension*( Ei-normal*normal )*mag(fvc::grad(alpha1)) );

TensionForce = mag(TensionForce)*normal;

fcf = fvc::interpolate(TensionForce) & ( mesh.Sf()/mesh.magSf() );
*/


/*
   // Solve for capillary pressure pc
  surfaceScalarField   tension  =   fvc::interpolate(surfaceTension*CurvatureLSM*DiracDelta)
        * fvc::snGrad(lsf)
//* fvc::snGrad(alpha1)
        * mesh.magSf();
*/


 //TensionForce = fvc::reconstruct(fcf*mesh.magSf() );
/*
 TensionForce = mag(TensionForce)*normal;

fcf = fvc::interpolate(TensionForce) & ( mesh.Sf()/mesh.magSf() );





*/


            fvScalarMatrix pcEqn
            (
                fvm::laplacian(pc) == fvc::div(mixture.surfaceTensionForce()*mesh.magSf())
            );


    pcEqn.relax();
    pcEqn.solve();








       
