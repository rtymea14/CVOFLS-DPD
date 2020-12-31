{
    volScalarField rAU("rAU", 1.0/UAdvEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); 
    //volVectorField HbyA(constrainHbyA(rAU*UAdvEqn.H(), UAdv, p_rghAdv));
    volVectorField HbyA("HbyA", UAdv); 
    HbyA = rAU*UAdvEqn.H();

    // First guess Flux calculated ignoring gravity and surface tension
    surfaceScalarField phiHbyA       
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(UAdv, phiAdv)
    );
    MRF.makeRelative(phiHbyA);  


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
    );


    phiHbyA += phig;


// Update the pressure BCs to ensure flux consistency
   setSnGrad<fixedFluxPressureFvPatchScalarField>    
       (    
           p_rghAdv.boundaryField(),    
           (    
               phiHbyA.boundaryField()    
             - MRF.relative(mesh.Sf().boundaryField() & UAdv.boundaryField())    
           )/(mesh.magSf().boundaryField()*rAUf.boundaryField())    
       );   
//   constrainPressure(p_rghAdv, UAdv, phiHbyA, rAUf, MRF);

    tmp<fvScalarMatrix> p_rghAdvEqnComp1;
    tmp<fvScalarMatrix> p_rghAdvEqnComp2;


    if (pimple.transonic())
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phiAdv);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phiAdv);

        p_rghAdvEqnComp1 =
            fvc::ddt(rho1) + fvc::div(phiAdv, rho1) - fvc::Sp(fvc::div(phiAdv), rho1)
          + correction
            (
                psi1*fvm::ddt(p_rghAdv)
              + fvm::div(phid1, p_rghAdv) - fvm::Sp(fvc::div(phid1), p_rghAdv)
            );
        deleteDemandDrivenData(p_rghAdvEqnComp1().faceFluxCorrectionPtr());
        p_rghAdvEqnComp1().relax();

        p_rghAdvEqnComp2 =
            fvc::ddt(rho2) + fvc::div(phiAdv, rho2) - fvc::Sp(fvc::div(phiAdv), rho2)
          + correction
            (
                psi2*fvm::ddt(p_rghAdv)
              + fvm::div(phid2, p_rghAdv) - fvm::Sp(fvc::div(phid2), p_rghAdv)
            );
        deleteDemandDrivenData(p_rghAdvEqnComp2().faceFluxCorrectionPtr());
        p_rghAdvEqnComp2().relax();
    }


    else
    {
       // Total derivative of rho1
        p_rghAdvEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rghAdv))
          + fvc::div(phiAdv, rho1) - fvc::Sp(fvc::div(phiAdv), rho1);

       // Total derivative of rho2
        p_rghAdvEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rghAdv))
          + fvc::div(phiAdv, rho2) - fvc::Sp(fvc::div(phiAdv), rho2);
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rghAdv_0 (p_rghAdv);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghAdvEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rghAdv)
        );



    solve
      (
            (
                (max(alpha1, scalar(0))/rho1)*p_rghAdvEqnComp1()
              + (max(alpha2, scalar(0))/rho2)*p_rghAdvEqnComp2()
            )

           + p_rghAdvEqnIncomp,

            mesh.solver(p_rghAdv.select(pimple.finalInnerIter()))
       );



        if (pimple.finalNonOrthogonalIter())
        {
            pAdv = max( p_rghAdv + (alpha1*rho1 + alpha2*rho2)*gh, pMin);
            p_rghAdv = pAdv - (alpha1*rho1 + alpha2*rho2)*gh  ;

            phiAdv = phiHbyA + p_rghAdvEqnIncomp.flux();
            p_rghAdv.relax();

            UAdv = HbyA + rAU*fvc::reconstruct((phig + p_rghAdvEqnIncomp.flux())/rAUf);
            UAdv.correctBoundaryConditions();
            fvOptions.correct(UAdv);
        }
    }


     Info << " UAdvmin: " << min(UAdv).value() << " UAdvmax: " << max(UAdv).value() << endl;

    // Update rho1,2 density from change in p_rgh (Neglected the difference rho_1*gh - rho_0*gh)


    volScalarField rho1Adv = rho1 + psi1*(p_rghAdv - p_rghAdv_0); 
   
    volScalarField rho2Adv = rho2 + psi2*(p_rghAdv - p_rghAdv_0); 

    volScalarField rhoAdv = alpha1*rho1Adv + alpha2*rho2Adv;

    rhoPhiAdv = fvc::interpolate(rhoAdv)*phiAdv;



    // Correct p_rgh for consistency with p and the updated densities
    p_rghAdv = pAdv - rhoAdv*gh ;
    p_rghAdv.correctBoundaryConditions();




}
