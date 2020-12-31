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

#ifndef parallelClass_H
#define parallelClass_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class parallelClass
{
       const fvMesh& mesh_;


public:

    // Constructors

        //- Construct from mesh
        parallelClass(const fvMesh& mesh);


    // Member Functions

       void exchange(List<scalarField>& fieldList);
       void exchange(List<tensorField>& fieldList);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parallelClass::parallelClass(const fvMesh& mesh):

     mesh_(mesh)
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::parallelClass::exchange(List<scalarField>& fieldList)
{  
		
	// Send all the fields from the slaves to the master
	if (Pstream::master())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		IPstream fromSlave(Pstream::commsTypes::scheduled, jSlave);
		fromSlave >> fieldList[jSlave]; 
	    }
	}

	else
	{
	    OPstream toMaster (Pstream::commsTypes::scheduled, Pstream::masterNo());			
	    toMaster << fieldList[Pstream::myProcNo()];	
	}


	// Distribute all the fields from the master to the slaves
	if (Pstream::master())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		OPstream toSlave (Pstream::commsTypes::scheduled, jSlave);

		forAll(fieldList, iProcess)
		{
		     toSlave << fieldList[iProcess];
		}
	    }
	}

	else
	{
		IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
		forAll(fieldList, iProcess)
	    	{
		   fromMaster >> fieldList[iProcess];
		}
	}
}




void Foam::parallelClass::exchange(List<tensorField>& fieldList)
{  
		
	// Send all the fields from the slaves to the master
	if (Pstream::master())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		IPstream fromSlave(Pstream::commsTypes::scheduled, jSlave);
		fromSlave >> fieldList[jSlave]; 
	    }
	}

	else
	{
	    OPstream toMaster (Pstream::commsTypes::scheduled, Pstream::masterNo());			
	    toMaster << fieldList[Pstream::myProcNo()];	
	}


	// Distribute all the fields from the master to the slaves
	if (Pstream::master())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		OPstream toSlave (Pstream::commsTypes::scheduled, jSlave);

		forAll(fieldList, iProcess)
		{
		     toSlave << fieldList[iProcess];
		}
	    }
	}

	else
	{
		IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
		forAll(fieldList, iProcess)
	    	{
		   fromMaster >> fieldList[iProcess];
		}
	}
}







#endif

}

// ************************************************************************* //
