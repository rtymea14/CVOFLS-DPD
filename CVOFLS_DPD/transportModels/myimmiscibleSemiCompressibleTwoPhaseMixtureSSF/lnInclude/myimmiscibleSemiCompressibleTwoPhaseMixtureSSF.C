/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "myimmiscibleSemiCompressibleTwoPhaseMixtureSSF.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immiscibleIncompressibleTwoPhaseMixture::
immiscibleIncompressibleTwoPhaseMixture
(
        const volVectorField& U,
        const volScalarField& rho1,
        const volScalarField& rho2,
        const volScalarField& nu1,
        const volScalarField& nu2,
        const volScalarField& cp1,
        const volScalarField& cp2,
        const volScalarField& lambda1,
        const volScalarField& lambda2,
        const volScalarField& pc,
        const volScalarField& Psii
)
:
    incompressibleTwoPhaseMixture(U, rho1, rho2, nu1, nu2, cp1, cp2, lambda1, lambda2),
    interfaceProperties(alpha1(), U, pc, Psii, *this)
{}


// ************************************************************************* //
