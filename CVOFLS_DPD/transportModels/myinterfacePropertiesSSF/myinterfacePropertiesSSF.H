/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Alpha curtailing coefficient
        scalar Cpc_;

        //- Alpha filtering coefficient
        scalar Csk_;

        //- Capillary flux error criterion
        scalar Cflux_;

        //- Number of alpha filters
         const unsigned int NAlphaFilters_;

        //- Number of curvature filters
         const unsigned int NCurvatureFilters_;

        //- Type of gradient calculation (classic vs Rho_Corrected)
         const int TypeOfGradient_;
        
	//- Surface tension
        dimensionedScalar sigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;
        const volScalarField& alpha1_;
        const volVectorField& U_;
        const volScalarField& pc_;
              volScalarField alpha1smooth_;
         //   volScalarField pcsmooth_;
              volScalarField alpha1pc_;
              surfaceScalarField nHatf_;
              surfaceScalarField Kf_;
              volScalarField Kvol_;
              volScalarField K_filtered_;
        //- Parameter for S-CLSVOF method (deltaX,gamma,epsilon,deltaTau)
        dimensionedScalar deltaX_;
        dimensionedScalar gamma_;
        dimensionedScalar epsilon_;
        dimensionedScalar deltaTau_;

        const volScalarField& psi_;
        volScalarField psi0_;
        volScalarField delta_;
        volScalarField H_;
        volScalarField Hscale_;
        volScalarField deltaScale_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfaceProperties(const interfaceProperties&);
        void operator=(const interfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::GeometricBoundaryField& nHat,
            surfaceVectorField::GeometricBoundaryField& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();

        //- Re-calculate the delta function and Heaviside function
        void calculateDelta();
        void calculateH();
        void calculateHscale();
        void calculateDeltaScale();

        //- Re-calculate the psi0
        void calculatePsi0();

public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const volScalarField& pc,
	    const volScalarField& Psii, 
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        scalar Cpc() const
        {
            return Cpc_;
        }

        scalar Csk() const
        {
            return Csk_;
        }

        scalar Cflux() const
        {
            return Cflux_;
        }

        scalar NAlphaFilters() const
        {
            return NAlphaFilters_;
        }

        scalar NCurvatureFilters() const
        {
            return NCurvatureFilters_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const surfaceScalarField& Kf() const
        {
            return Kf_;
        }

         volScalarField& alpha1smooth() 
        {
            return alpha1smooth_;
        }
/*
         volScalarField& pcsmooth() 
        {
            return pcsmooth_;
        }
*/
         volScalarField& alpha1pc() 
        {
            return alpha1pc_;
        }

         volScalarField& Kvol() 
        {
            return Kvol_;
        }

         volScalarField& K_filtered() 
        {
            return K_filtered_;
        }


        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        tmp<surfaceScalarField> surfaceTensionForce() ;

        tmp<surfaceScalarField> phi_c(const surfaceScalarField& rAUf_) const ;

        //- Indicator of the proximity of the interface : Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        void correct()
        {
            calculateK();
        }
        
	const dimensionedScalar& epsilon() const
        {
            return epsilon_;
        }

        const dimensionedScalar& deltaTau() const
        {
            return deltaTau_;
        }

        const volScalarField& psi0() const
        {
            return psi0_;
        }

        const volScalarField& delta() const
        {
            return delta_;
        }

        const volScalarField& H() const
        {
            return H_;
        }

        const volScalarField& Hscale() const
        {
            return Hscale_;
        }

        const volScalarField& deltaScale() const
        {
            return deltaScale_;
        }

        void correctPsi0()
        {
            calculatePsi0();
            psi0_.correctBoundaryConditions();
        }

        void correctDeltaH()
        {
            calculateDelta();
            calculateH();
            delta_.correctBoundaryConditions();
            H_.correctBoundaryConditions();
        }

        void correctDeltaHscale()
        {
            calculateDelta();
            calculateH();
            calculateHscale();
            calculateDeltaScale();
            delta_.correctBoundaryConditions();
            H_.correctBoundaryConditions();
            deltaScale_.correctBoundaryConditions();
            Hscale_.correctBoundaryConditions();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
