/*---------------------------------------------------------------------------------------------------------------------------*\
Description

 The user can choose whether to use MULES or isoAdvector. The latter is highly recommended.

\*----------------------------------------------------------------------------------------------------------------------------*/


// Looking up which alpha1 solver to use
const dictionary& alpha1Dict(mesh.solverDict(alpha1.name()));

word interfaceMethod
(
    alpha1Dict.lookupOrDefault<word>("interfaceMethod", "MULES")
);

// Advancing alpah1 field
if (interfaceMethod == "MULES")
{
    #include "alphaEqnMULES.H"
}



else if (interfaceMethod == "isoAdvector")
{
    // If there are more than one outer corrector, we use a mixture of old and new U and phi for propagating alpha1 in all but the first outer iteration
    if (!pimple.firstIter())
    {
        // We are recalculating alpha1 from the its old time value
        alpha1 = alpha1.oldTime();
        // Temporarily storing new U and phi values in prevIter storage
        UAdv.storePrevIter();
        phiAdv.storePrevIter();

        // Overwriting new U and phi values with mixture of old and new values
        phiAdv = 0.5*(phiAdv + phiAdv.oldTime());
        UAdv = 0.5*(UAdv + UAdv.oldTime());
    }

    // Update alpha1
    advector.advect();



    // Update rhoPhi
    surfaceScalarField rho1f (fvc::interpolate(rho1));
    surfaceScalarField rho2f (fvc::interpolate(rho2)); 

    rhoPhi = advector.getRhoPhi(rho1f, rho2f);
    rhoPhiAdv = advector.getRhoPhi(rho1f, rho2f);



    alpha2 = 1.0 - alpha1;

    if (!pimple.firstIter())
    {
        // Restoring new U and phi values temporarily saved in prevIter() above
        UAdv = UAdv.prevIter();
        phiAdv = phiAdv.prevIter();
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
        << endl;
}

//U not solved but prescribed and optionally periodic
if (prescribedUAdv)
{

    //Setting velocity field and face fluxes for next time step
    const scalar t(runTime.time().value());
    const scalar dt(runTime.deltaT().value());
    if (reverseTime > 0.0 && t >= reverseTime)
    {
        Info<< "Reversing flow" << endl;
        phiAdv *= -1.0;
        UAdv *= -1.0;
        reverseTime = -1.0;
    }
    if (period > 0.0)
    {
        const scalar uFactor
        (
            .5*(Foam::cos(2.0*M_PI*t/period)
                + Foam::cos(2.0*M_PI*(t + dt)/period))
        );
        phiAdv = uFactor*phi0Adv;
        UAdv = uFactor*U0Adv;
    }
}




