/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::incompressibleTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    incompressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseMixture_H
#define incompressibleTwoPhaseMixture_H

//#include "incompressible/transportModel/transportModel.H"
//#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
//#include "transportModel.H"
//#include "viscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseMixture
:
    public IOdictionary,
    public twoPhaseMixture
{
protected:

// Constructor data
        const volVectorField& U_;
        const volScalarField& rho1_;
        const volScalarField& rho2_;
        const volScalarField& nu1_;
        const volScalarField& nu2_;
        const volScalarField& cp1_;
        const volScalarField& cp2_;
        const volScalarField& lambda1_;
        const volScalarField& lambda2_;

        // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcmu();


public:

    TypeName("incompressibleTwoPhaseMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseMixture
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

        );

            volScalarField mu_;


    //- Destructor
    virtual ~incompressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 density
        const volScalarField& rho1() const
        {
            return rho1_;
        }


        //- Return const-access to phase2 density
        const volScalarField& rho2() const
        {
            return rho2_;
        }

        //- Return access to rho
        const volScalarField rho() 
        {
            return alpha1_*rho1_ + alpha2_*rho2_;
        }

        //- Return const-access to phase1 viscosity
        const volScalarField& nu1() const
        {
            return nu1_;
        }


        //- Return const-access to phase2 viscosity
        const volScalarField& nu2() const
        {
            return nu2_;
        }

        //- Return access to nu
       const  volScalarField nu() 
        {
            return alpha1_*nu1_ + alpha2_*nu2_;
        }


        //- Return const-access to phase1 cp
        const volScalarField& cp1() const
        {
            return cp1_;
        }

        //- Return const-access to phase2 cp
        const volScalarField& cp2() const
        {
            return cp2_;
        }

        //- Return access to cp
        const volScalarField cp() 
        {
            return (alpha1_*rho1_*cp1_ + alpha2_*rho2_*cp2_)/(alpha1_*rho1_ + alpha2_*rho2_);
        }

        //- Return const-access to phase1 lambda
        const volScalarField& lambda1() const
        {
            return lambda1_;
        }

        //- Return const-access to phase2 lambda
        const volScalarField& lambda2() const
        {
            return lambda2_;
        }

        //- Return access to lambda
         const volScalarField lambda() 
        {
            return alpha1_*lambda1_ + alpha2_*lambda2_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }


        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> mu() const
        {
            return mu_;
        }


        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcmu();
        }


        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
