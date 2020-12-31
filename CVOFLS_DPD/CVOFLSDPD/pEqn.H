{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); 
    //volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    volVectorField HbyA("HbyA", U); 
    HbyA = rAU*UEqn.H();

    // First guess Flux calculated ignoring gravity and surface tension
    surfaceScalarField phiHbyA       
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    MRF.makeRelative(phiHbyA);
  // adjustPhi(phiHbyA, U, p_rgh);    


// Capillary fluxes
surfaceScalarField phi_c
    (
        IOobject
        (
            "phi_c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("phi_c", dimLength*dimLength*dimLength/dimTime, 0.)
    );


//Raeini et.al. 2012
/*

phi_c = ( fcf - 0.0*fvc::snGrad(pc) )*rAUf*mesh.magSf();  

phic = fvc::average(phi_c);
surfaceScalarField  phi_classic (  fcf*rAUf*mesh.magSf() );  
//surfaceScalarField  phi_classic (  phi_c );  

const dimensionedScalar  dummyFlux("dummyFlux", dimensionSet(0,3,-1,0,0,0,0), 1.0);
scalar Cflux = readScalar (LiquidPropertiesDictionary.lookup("Cflux") );   	    

    dimensionedScalar phi_c_thresh
    ( 
      Cflux*dummyFlux*gSum( mag(phi_classic.field()) ) / gSum( SMALL + pos( mag(phi_classic.field()) - SMALL ) ) 
    );

      phi_c = phi_c - max( min(phi_c, phi_c_thresh), -phi_c_thresh );


phicFiltered = fvc::average(phi_c);





//Brackbill et.al. 1992
//mixture.correct();
//phi_c = mixture.surfaceTensionForce()*rAUf*mesh.magSf(); 
//phi_c = mixture.phi_c(rAUf);



// Sussman et.al. 1994 (Level Of Set Method)

phi_c =   fvc::interpolate(surfaceTension*CurvatureLSM*DiracDelta)
        * fvc::snGrad(lsf)
       * rAUf*mesh.magSf();



// FixedCurvature

volScalarField alpha1pc = (1./(1.-0.999))*(min(max(alpha1, 0.999/2.), 1.-0.999/2.)-0.999/2.);

Curvature = 4.0/(Diameter);
Info << "Curvature =  " << Curvature.value() << " 1/m " <<  endl;

phi_c =   surfaceTension*Curvature*fvc::snGrad(alpha1pc)*rAUf*mesh.magSf();
*/



// Gravitational fluxes
surfaceScalarField phi_rgh
    (
        IOobject
        (
            "phi_rgh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("phi_rgh", dimLength*dimLength*dimLength/dimTime, 0.)
    );

     phi_rgh = -ghf*fvc::snGrad(rho)*rAUf*mesh.magSf();



     dimensionedScalar  Mass = fvc::domainIntegrate(rho1*alpha1);


    surfaceScalarField phig
    (
        phi_rgh
      + phi_c  
    );


    phiHbyA += phig;


// Update the pressure BCs to ensure flux consistency
   setSnGrad<fixedFluxPressureFvPatchScalarField>    
       (    
           p_rgh.boundaryField(),    
           (    
               phiHbyA.boundaryField()    
             - MRF.relative(mesh.Sf().boundaryField() & U.boundaryField())    
           )/(mesh.magSf().boundaryField()*rAUf.boundaryField())    
       );   
//   constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;


    if (pimple.transonic())
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

        p_rghEqnComp1 =
            fvc::ddt(rho1) + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1)
          + correction
            (
                psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp1().faceFluxCorrectionPtr());
        p_rghEqnComp1().relax();

        p_rghEqnComp2 =
            fvc::ddt(rho2) + fvc::div(phi, rho2) - fvc::Sp(fvc::div(phi), rho2)
          + correction
            (
                psi2*fvm::ddt(p_rgh)
              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp2().faceFluxCorrectionPtr());
        p_rghEqnComp2().relax();
    }


    else
    {
       // Total derivative of rho1
        p_rghEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rgh))
          + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1);

       // Total derivative of rho2
        p_rghEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rgh))
          + fvc::div(phi, rho2) - fvc::Sp(fvc::div(phi), rho2);
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0 (p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

  //      p_rghEqnIncomp.setReference(pRefCell, getRefCellValue(p, pRefCell));


    solve
      (
            (
                (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1()
              + (max(alpha2, scalar(0))/rho2)*p_rghEqnComp2()
            )

           + p_rghEqnIncomp
           - EvapCorrection*interfacialAreaCoefficient*mEvapTot*(1./rho1 - 1./rho2),

            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
       );



        if (pimple.finalNonOrthogonalIter())
        {
            p = max( p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh  ;



            dgdt =
            (
                pos(alpha2)*(p_rghEqnComp2 & p_rgh)/rho2
              - pos(alpha1)*(p_rghEqnComp1 & p_rgh)/rho1
            );


            phi = phiHbyA + p_rghEqnIncomp.flux();
            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }


     Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;

    // Update rho1,2 density from change in p_rgh (Neglected the difference rho_1*gh - rho_0*gh)


    rho1 += psi1*(p_rgh - p_rgh_0); 
   
    rho2 += psi2*(p_rgh - p_rgh_0);

    rho = alpha1*rho1 + alpha2*rho2;

    rhoPhi = fvc::interpolate(rho)*phi;



    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh ;
    p_rgh.correctBoundaryConditions();



/*
// Check of continuity errors for liquid phase
{
    scalar cumulativeContErrLiquid = 0;
    dimensionedScalar totalLiquidPhaseMass = fvc::domainIntegrate(alpha1*rho1);

    scalar sumLocalContErrLiquid =
        (fvc::domainIntegrate(mag(alpha1*(rho1 - p*psi1)))/totalLiquidPhaseMass).value();

    scalar globalContErrLiquid =
        (fvc::domainIntegrate(alpha1*(rho1 - p*psi1))/totalLiquidPhaseMass).value();

    cumulativeContErrLiquid += globalContErrLiquid;

    Info<< "time step continuity errors for liquid phase : sum local = " << sumLocalContErrLiquid
        << ", global = " << globalContErrLiquid
        << ", cumulative = " << cumulativeContErrLiquid
        << endl;
}

// Check of continuity errors for gas phase
{
    scalar cumulativeContErrGas = 0;
    dimensionedScalar totalGasPhaseMass = fvc::domainIntegrate(alpha2*rho2);

    scalar sumLocalContErrGas =
        (fvc::domainIntegrate(mag(alpha2*(rho2 - p*psi2)))/totalGasPhaseMass).value();

    scalar globalContErrGas =
        (fvc::domainIntegrate(alpha2*(rho2 - p*psi2))/totalGasPhaseMass).value();

    cumulativeContErrGas += globalContErrGas;

    Info<< "time step continuity errors for gas phase : sum local = " << sumLocalContErrGas
        << ", global = " << globalContErrGas
        << ", cumulative = " << cumulativeContErrGas
        << endl;
}
*/
// Lagrangian derivative of pressure (for energy equation)
    DpDt = fvc::DDt(phi, p);

}
