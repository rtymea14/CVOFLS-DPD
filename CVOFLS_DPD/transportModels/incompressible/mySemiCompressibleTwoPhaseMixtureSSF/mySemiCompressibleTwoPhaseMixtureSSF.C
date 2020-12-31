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

\*---------------------------------------------------------------------------*/

#include "mySemiCompressibleTwoPhaseMixtureSSF.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::incompressibleTwoPhaseMixture::calcmu()
{
      const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );


         mu_=    limitedAlpha1*rho1_*nu1_ + (scalar(1) - limitedAlpha1)*rho2_*nu2_;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseMixture::incompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const volScalarField& rho1,
    const volScalarField& rho2,
    const volScalarField& nu1,
    const volScalarField& nu2,
        const volScalarField& cp1,
        const volScalarField& cp2,
        const volScalarField& lambda1,
        const volScalarField& lambda2
  //  const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),
/*
    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),


    rho1_("rho", dimDensity, nuModel1_->viscosityProperties().lookup("rho")),



    cp1_("cp", dimSpecificHeatCapacity, nuModel1_->viscosityProperties().lookup("cp")),
    cp2_("cp", dimSpecificHeatCapacity, nuModel2_->viscosityProperties().lookup("cp")),

    lambda1_("lambda", dimPower/(dimTemperature*dimLength), nuModel1_->viscosityProperties().lookup("lambda")),
    lambda2_("lambda", dimPower/(dimTemperature*dimLength), nuModel2_->viscosityProperties().lookup("lambda")),

*/

    U_(U),
    rho1_(rho1),
    rho2_(rho2),
    nu1_(nu1),
    nu2_(nu2),
    cp1_(cp1),
    cp2_(cp2),
    lambda1_(lambda1),
    lambda2_(lambda2),

 //   phi_(phi),

    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("mu", dimViscosity*dimDensity, 0.),
        calculatedFvPatchScalarField::typeName
    )

{
    calcmu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

/*
Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (

           "mu",
            limitedAlpha1*rho1_*nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()

            "mu",
            limitedAlpha1*rho1_*nu1_
          + (scalar(1) - limitedAlpha1)*rho2_*nu2_

        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (

            "muf",
            alpha1f*fvc::interpolate(rho1_)*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*fvc::interpolate(rho2_)*fvc::interpolate(nuModel2_->nu())

            "muf",
            alpha1f*fvc::interpolate(rho1_)*fvc::interpolate(nu1_)
          + (scalar(1) - alpha1f)*fvc::interpolate(rho2_)*fvc::interpolate(nu2_)

        )
    );
}

*
Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (
                alpha1f*fvc::interpolate(rho1_)*fvc::interpolate(nuModel1_->nu())
              + (scalar(1) - alpha1f)*fvc::interpolate(rho2_)*fvc::interpolate(nuModel2_->nu())
            )/(alpha1f*fvc::interpolate(rho1_) + (scalar(1) - alpha1f)*fvc::interpolate(rho2_))
        )
    );
}


*/

bool Foam::incompressibleTwoPhaseMixture::read()
{
/*
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
           // nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
           // nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

            nuModel1_->viscosityProperties().lookup("cp") >> cp1_;
            nuModel2_->viscosityProperties().lookup("cp") >> cp2_;

            nuModel1_->viscosityProperties().lookup("lambda") >> lambda1_;
            nuModel2_->viscosityProperties().lookup("lambda") >> lambda2_;

*/

            return true;
/*
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
*/

}



// ************************************************************************* //
