{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phir(phic*mixture.nHatf());
   //   surfaceScalarField phir ( phic* (fvc::interpolate(normal) & mesh.Sf()) );

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dgdt.dimensions(), 0.0)
        );

        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            //  EvapCorrection*mEvapTot*alpha1/rho1
                EvapCorrection*mEvapTot*Foam::pow(Foam::mag(alpha1),0.5)/rho1
        );

        forAll(dgdt, celli)
        {
            if (dgdt[celli] > 0.0 && alpha1[celli] > 0.0)
            {
                Sp[celli] -= dgdt[celli]*alpha1[celli];
                Su[celli] += dgdt[celli]*alpha1[celli];
            }
            else if (dgdt[celli] < 0.0 && alpha1[celli] < 1.0)
            {
                Sp[celli] += dgdt[celli]*(1.0 - alpha1[celli]);
            }
        }


/*
        forAll(mevap, celli)
        {

            if (mEvapTot[celli] > 0.0)
            {
              //  Sp[celli] -= mevap[celli]/rho1[celli];
                Su[celli] -= mEvapTot[celli]*( alpha1[celli]*(1./rho1[celli] - 1./rho2[celli]) );
            }
            else if (mEvapTot[celli] < 0.0)
            {

               Sp[celli] -= mEvapTot[celli]*(1./rho1[celli] - 1./rho2[celli]); //alpha1 sottinteso
              Su[celli] += mEvapTot[celli]/rho1[celli];
             //   Sp[celli] += mEvapTot[celli]*(1./rho2[celli] - 1./rho1[celli]); //alpha1 sottinteso
             //  Su[celli] += mEvapTot[celli]/rho2[celli];
            }
        }
*/

        surfaceScalarField alphaPhi1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        MULES::explicitSolve
        (
            geometricOneField(),
             alpha1,
            phi,
            alphaPhi1,
            Sp,
            Su,
            1,
            0
        );


// Normalization of alpha

forAll(alpha1.internalField(), celli)
{

   if (alpha1.internalField()[celli] > 1.0)
   {
       alpha1.internalField()[celli] = 1.0;
   }

   if (alpha1.internalField()[celli] < 0.0)
   {
       alpha1.internalField()[celli] = 0.0;
   }

}



forAll(alpha1.boundaryField(), patchi)
{
    forAll(alpha1.boundaryField()[patchi], facei)
    {
        if (alpha1.boundaryField()[patchi][facei] > 1.0)
        {
           alpha1.boundaryField()[patchi][facei] = 1.0;
        }

   if (alpha1.boundaryField()[patchi][facei] < 0.0)
   {
       alpha1.boundaryField()[patchi][facei] = 0.0;
   }

    }
}


        surfaceScalarField rho1f(fvc::interpolate(rho1));
        surfaceScalarField rho2f(fvc::interpolate(rho2));
        rhoPhi = alphaPhi1*(rho1f - rho2f) + phi*rho2f;

        alpha2 = scalar(1) - alpha1;
    }

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << endl;
}
