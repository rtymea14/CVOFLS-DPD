
/*---------------------------------------------------------------------------------------------------------------------------*\
Description

  The interface is here advected. IsoAdvector library by Roenby e Jasak (2016) is used, advecting the interface in a geometric 
  approach. The velocity and the fluxes needed for the advection are not U and phi, but UAdv and phiAdv, which don't take into
  account the vaporization velocity, which is too big and tends to destroy the interface.

  The user can also choose to use the MULES algorithm developed by Weller.

  The interface is advected only if gravity and the magnetic field are activated. In microgravity condition the interface cannot
  be displaced by the velocity field by definition.

\*----------------------------------------------------------------------------------------------------------------------------*/


// If gravity is not present the advection of the interface is suppressed

if ( (mag(g)).value() != 0.)
{

    #include "alphaControls.H"

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(mixture.cAlpha()*phic, max(phic));



    volScalarField divU(fvc::div(fvc::absolute(phi, U)));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", rhoPhi.dimensions(), 0)
        );

        for
        (
          subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaAdvEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
         
    }
    else
    {
          #include "alphaAdvEqn.H"
    }


}




