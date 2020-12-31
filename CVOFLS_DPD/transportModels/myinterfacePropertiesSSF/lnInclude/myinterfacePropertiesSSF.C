/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "myinterfacePropertiesSSF.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::GeometricBoundaryField& nHatb,
    surfaceVectorField::GeometricBoundaryField& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::GeometricBoundaryField& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Alpha filtering

if (NAlphaFilters_ == 0)

       {
         alpha1smooth_ = alpha1_; 
       }
else

       {
           alpha1smooth_ = Csk_*fvc::average(fvc::interpolate(alpha1_)) + (1.-Csk_)*alpha1_;
                 
           for(unsigned int i=0;i<NAlphaFilters_-1;i++)

                 {
                   alpha1smooth_ = Csk_*fvc::average(fvc::interpolate(alpha1smooth_)) + (1.-Csk_)*alpha1smooth_;

                 }

       }


    // Cell gradient of alpha
    //volVectorField gradAlpha(fvc::grad(alpha1smooth_, "nHat"));    // Provare a togliere nHat
      volVectorField gradAlpha( fvc::grad(alpha1smooth_,"nHat") );

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

    //gradAlphaf -=
    //    (mesh.Sf()/mesh.magSf())
    //   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlphaf/(Foam::mag(gradAlphaf) + deltaN_));

    // Correct Boundary Angle
//    correctContactAngle(nHatfv.boundaryField(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Expression for curvature
    Kvol_ = -fvc::div(nHatf_);


    // Curvature filtering

       if (NCurvatureFilters_ == 0)
       {
        Kf_ = (fvc::interpolate(Kvol_)); 
       }

       else

      {

          volScalarField  walpha_  (alpha1_*(1.-alpha1_) );

          volScalarField w ( sqrt( mag(walpha_) ) );
        
          volScalarField wsmooth (fvc::average(fvc::interpolate(w)));
        
          K_filtered_ = 2.*w*Kvol_ + (1.-2.*w)*( (fvc::average(fvc::interpolate(Kvol_*w)))/(wsmooth+1e-6) );
  
             for(unsigned int i=0;i<NCurvatureFilters_-1;i++)
             {
               K_filtered_ =  2.*w*Kvol_ + (1.-2.*w)*( (fvc::average(fvc::interpolate(K_filtered_*w)))/(wsmooth+1e-6) );                     
             }

         
         // Interpolated weighted value of curvature   
         Kf_ = (fvc::interpolate(K_filtered_*w))/(fvc::interpolate(w)+1e-6);

      }
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() 
{
      // Curtailing alpha field

       alpha1pc_ = (1./(1.-Cpc_))*(min(max(alpha1_, Cpc_/2.), 1.-Cpc_/2.)-Cpc_/2.);


  if (TypeOfGradient_ == 1)    // Classic method
  {
       return sigma_*Kf_*fvc::snGrad(alpha1pc_); 
    // return fvc::interpolate(sigma_*Kvol_)*fvc::snGrad(alpha1_); 
  }


}


Foam::tmp<Foam::surfaceScalarField> 
Foam::interfaceProperties::phi_c(const surfaceScalarField& rAUf_) const
{

          const fvMesh& mesh = alpha1_.mesh();

          const surfaceScalarField phi_c_i
         ( 
          fvc::snGrad(pc_)*rAUf_* mesh.magSf() 
         );

    const dimensionedScalar  dummyFlux("dummyFlux", dimensionSet(0,3,-1,0,0,0,0), 1.0);

    dimensionedScalar phi_c_thresh
    ( 
          Cflux_*dummyFlux*gSum( mag(phi_c_i.field()) ) / gSum( SMALL + pos( mag(phi_c_i.field()) - SMALL ) ) 
    );

    // Return filtered capillary flux
    return     phi_c_i - max( min(phi_c_i, phi_c_thresh), -phi_c_thresh );
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos(alpha1_ - 0.01)*pos(0.99 - alpha1_);
}

void Foam::interfaceProperties::calculatePsi0()
{
    //psi0_ == (double(5.0)*mixture_.alpha1() - double(4.0))*gamma_;
    psi0_ == (double(2.0)*alpha1_ - double(1.0))*gamma_;
}

void Foam::interfaceProperties::calculateDelta()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(mag(psi_[celli]) > epsilon_.value())
          delta_[celli] = double(0.0);
       else
          delta_[celli] = double(1.0)/(double(2.0)*epsilon_.value())*(double(1.0)+cos(M_PI*psi_[celli]/epsilon_.value()));
    }
}

void Foam::interfaceProperties::calculateH()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          H_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          H_[celli] = double(1.0);
       else
          H_[celli] = double(1.0)/double(2.0)*(double(1.0)+psi_[celli]/epsilon_.value()+sin(M_PI*psi_[celli]/epsilon_.value())/M_PI);
    }
}

void Foam::interfaceProperties::calculateHscale()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          Hscale_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          Hscale_[celli] = double(1.0);
       else
          Hscale_[celli] = double(1.0)/double(2.0)*(double(1.0)/double(2.0)+psi_[celli]/epsilon_.value()+psi_[celli]*psi_[celli]/(double(2.0)*epsilon_.value()*epsilon_.value())-(cos(double(2.0)*M_PI*psi_[celli]/epsilon_.value())-double(1.0))/(double(4.0)*M_PI*M_PI)+sin(M_PI*psi_[celli]/epsilon_.value())*(epsilon_.value()+psi_[celli])/(M_PI*epsilon_.value()));
    }
}

void Foam::interfaceProperties::calculateDeltaScale()
{
    deltaScale_ == double(2.0)*H_*delta_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const volScalarField& pc,
    const volScalarField& Psii,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    psi_(Psii),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),

    Cpc_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("Cpc")
        )
    ),

    Csk_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("Csk")
        )
    ),

    Cflux_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("Cflux")
        )
    ),

    NAlphaFilters_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("NAlphaFilters")
        )
    ),

    NCurvatureFilters_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("NCurvatureFilters")
        )
    ),


    TypeOfGradient_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("TypeOfGradient")
        )
    ),
    sigma_(dict.lookup("sigma")),

    deltaN_
    (
        "deltaN",
        1e-6/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    alpha1_(alpha1),

    U_(U),
   
    pc_(pc),

    alpha1smooth_
    (
        IOobject
        (
            "interfaceProperties:alpha1smooth_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
       ),
        alpha1_.mesh(),
        dimensionedScalar("alpha1smooth_", dimless, 0.0)
    ),
/*
    pcsmooth_
    (
        IOobject
        (
            "interfaceProperties:pcsmooth_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
       ),
        alpha1_.mesh(),
        dimensionedScalar("pcsmooth_", dimPressure, 0.0)
    ),
*/

    alpha1pc_
    (
        IOobject
        (
            "interfaceProperties:alpha1pc_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
       alpha1_.mesh(),
        dimensionedScalar("alpha1pc_", dimless, 0.0)
    ),


    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    Kf_
    (
        IOobject
        (
            "interfaceProperties:Kf_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf_", dimless/dimLength, 0.0)
    ),

    Kvol_
    (
        IOobject
        (
            "interfaceProperties:Kvol_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kvol_", dimless/dimLength, 0.0)
    ),

    K_filtered_
    (
        IOobject
        (
            "interfaceProperties:K_filtered_",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K_filtered_", dimless/dimLength, 0.0)
    ),

    deltaX_(dict.lookup("deltaX")),

    gamma_
    (
        "gamma",
        deltaX_*double(0.75)
    ),

    epsilon_
    (
        "epsilon",
        deltaX_*double(1.5)
    ),

    //deltaTau_(dict.lookup("deltaTau")),
    deltaTau_
    (
        "deltaTau",
        deltaX_*double(0.1)
    ),

    psi0_
    (
        IOobject
        (
            "psi0",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("psi0", dimless, 0.0),
        alpha1_.boundaryField().types()
    ),

    delta_
    (
        IOobject
        (
            "delta",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("delta", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    H_
    (
        IOobject
        (
            "H",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("H", dimless, 0.0),
        alpha1_.boundaryField().types()
    ),

    Hscale_
    (
        IOobject
        (
            "Hscale",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Hscale", dimless, 0.0),
        alpha1_.boundaryField().types()
    ),

    deltaScale_
    (
        IOobject
        (
            "deltaScale",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("deltaScale", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    )
{
    calculateK();
}


// ************************************************************************* //
