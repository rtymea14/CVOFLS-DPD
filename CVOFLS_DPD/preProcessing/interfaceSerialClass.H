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


\*---------------------------------------------------------------------------*/

#ifndef interfaceSerialClass_H
#define interfaceSerialClass_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceSerialClass Declaration
\*---------------------------------------------------------------------------*/

class interfaceSerialClass
{
         volScalarField& alpha1_;


public:

    // Constructors

        //- Construct from volume fraction 
        interfaceSerialClass(volScalarField& alpha1);


    // Member Functions

       double GetGasValue(volScalarField& field, int celli);
      // std::vector<double> GetGasValue(volVectorField& field, int celli);
       double GetLiquidValue(volScalarField& field, int celli);
       int GetGasCell(int celli);
       int GetLiquidCell(int celli);
       bool checkInterface(volScalarField& gradAlpha);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceSerialClass::interfaceSerialClass(volScalarField& alpha1):

    alpha1_(alpha1)
{

}


// ********************************************************************* MEMBER FUNCTIONS ************************************************************************************* // 



// ******************************************************************************************************** // 
// ******************************************  CHECK INTERFACE ******************************************** // 
// ******************************************************************************************************** // 

bool Foam::interfaceSerialClass::checkInterface(volScalarField& gradAlpha)
{

     int DetectInterface = 0;

     forAll(alpha1_.internalField(), celli)
     {
       if (gradAlpha.internalField()[celli] > 1e-3 )
        {            
           DetectInterface++; 
           break;                 
        }
     }


       if (DetectInterface != 0)
       {
          return true;
       }

       else
       {
          return false;
       }

}


// ******************************************************************************************************** // 
// ******************************************  GET GAS VALUE ********************************************** // 
// ******************************************************************************************************** // 

double Foam::interfaceSerialClass::GetGasValue(volScalarField& field, int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];
              label marker = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue > 1e-3 )
              {
                 NumberOfIterations++;
                 std::vector<double> alphaValues(mesh.cellCells()[marker].size());

                 forAll (mesh.cellCells()[marker], ncelli)
                 {
                    label cellNeighbour = mesh.cellCells()[marker][ncelli];   
                    alphaValues[ncelli] = alpha1_.internalField()[cellNeighbour];   
                 }

                  double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                  double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 0.0;
                   marker=-2000;
                 }
                 
                 else
                 {

                   forAll (mesh.cellCells()[marker], ncelli)
                   {
                     label cellNeighbour = mesh.cellCells()[marker][ncelli];

                     if (alpha1_.internalField()[cellNeighbour] == minAlpha)
                     {
                       thresholdAlphaValue = minAlpha;
                       marker = cellNeighbour; 
                       break;                  
                     } 

                   }
                 }

               }

                 if (marker != -2000)
                 {
                   double GasValue = field.internalField()[marker];
                   return GasValue;
                 }

                 else
                 {
                   return -1.0;
                 }

}


// ******************************************************************************************************** // 
// ******************************************  GET GAS VALUE ********************************************** // 
// ******************************************************************************************************** // 

/*
std::vector<double> Foam::interfaceSerialClass::GetGasValue(volVectorField& field, int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];
              label marker = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue > 1e-3 )
              {
                 NumberOfIterations++;
                 std::vector<double> alphaValues(mesh.cellCells()[marker].size());

                 forAll (mesh.cellCells()[marker], ncelli)
                 {
                    label cellNeighbour = mesh.cellCells()[marker][ncelli];   
                    alphaValues[ncelli] = alpha1_.internalField()[cellNeighbour];   
                 }

                  double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                  double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 0.0;
                   marker=-2000;
                 }
                 
                 else
                 {

                   forAll (mesh.cellCells()[marker], ncelli)
                   {
                     label cellNeighbour = mesh.cellCells()[marker][ncelli];

                     if (alpha1_.internalField()[cellNeighbour] == minAlpha)
                     {
                       thresholdAlphaValue = minAlpha;
                       marker = cellNeighbour; 
                       break;                  
                     } 

                   }
                 }

               }

                 if (marker != -2000)
                 {
                   std::vector<double> GasValue = field.internalField()[marker];
                   return GasValue;
                 }

                 else
                 {
                   return -1.0;
                 }

}
*/



// ******************************************************************************************************** // 
// ******************************************  GET GAS CELL *********************************************** // 
// ******************************************************************************************************** // 



int Foam::interfaceSerialClass::GetGasCell(int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];
              label marker = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue > 1e-3 )
              {
                 NumberOfIterations++;
                 std::vector<double> alphaValues(mesh.cellCells()[marker].size());

                 forAll (mesh.cellCells()[marker], ncelli)
                 {
                    label cellNeighbour = mesh.cellCells()[marker][ncelli];   
                    alphaValues[ncelli] = alpha1_.internalField()[cellNeighbour];   
                 }

                  double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                  double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 0.0;
                   marker=-2000;
                 }
                 
                 else
                 {

                   forAll (mesh.cellCells()[marker], ncelli)
                   {
                     label cellNeighbour = mesh.cellCells()[marker][ncelli];

                     if (alpha1_.internalField()[cellNeighbour] == minAlpha)
                     {
                       thresholdAlphaValue = minAlpha;
                       marker = cellNeighbour; 
                       break;                  
                     } 

                   }
                 }

               }

                 if (marker != -2000)
                 {
                   return marker;
                 }

                 else
                 {
                   return -1;
                 }

}


// ******************************************************************************************************** // 
// ******************************************  GET LIQUID VALUE ******************************************* // 
// ******************************************************************************************************** // 



double Foam::interfaceSerialClass::GetLiquidValue(volScalarField& field, int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];
              label marker = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue < 0.999 )
              {
                 NumberOfIterations++;
                 std::vector<double> alphaValues(mesh.cellCells()[marker].size());

                 forAll (mesh.cellCells()[marker], ncelli)
                 {
                    label cellNeighbour = mesh.cellCells()[marker][ncelli];   
                    alphaValues[ncelli] = alpha1_.internalField()[cellNeighbour];   
                 }

                  double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());
                  double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 1.0;
                   marker=-2000;
                 }
                 
                 else
                 {

                   forAll (mesh.cellCells()[marker], ncelli)
                   {
                     label cellNeighbour = mesh.cellCells()[marker][ncelli];

                     if (alpha1_.internalField()[cellNeighbour] == maxAlpha)
                     {
                       thresholdAlphaValue = maxAlpha;
                       marker = cellNeighbour; 
                       break;                  
                     } 

                   }
                 }

               }

                 if (marker != -2000)
                 {
                   double LiquidValue = field.internalField()[marker];
                   return LiquidValue;
                 }

                 else
                 {
                   return -1.0;
                 }
}





// ******************************************************************************************************** // 
// ******************************************  GET LIQUID CELL ******************************************** // 
// ******************************************************************************************************** // 


int Foam::interfaceSerialClass::GetLiquidCell(int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];
              label marker = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue < 0.999 )
              {
                 NumberOfIterations++;
                std::vector<double> alphaValues(mesh.cellCells()[marker].size());

                 forAll (mesh.cellCells()[marker], ncelli)
                 {
                    label cellNeighbour = mesh.cellCells()[marker][ncelli];   
                    alphaValues[ncelli] = alpha1_.internalField()[cellNeighbour];   
                 }

                  double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                  double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 1.0;
                   marker=-2000;
                 }
                 
                 else
                 {

                   forAll (mesh.cellCells()[marker], ncelli)
                   {
                     label cellNeighbour = mesh.cellCells()[marker][ncelli];

                     if (alpha1_.internalField()[cellNeighbour] == maxAlpha)
                     {
                       thresholdAlphaValue = maxAlpha;
                       marker = cellNeighbour; 
                       break;                  
                     } 

                   }
                 }

               }

                 if (marker != -2000)
                 {
                   return marker;
                 }

                 else
                 {
                   return -1;
                 }

}

// ********************************************************************* END OF MEMBER FUNCTIONS ************************************************************************************* // 


#endif

}

// ************************************************************************* //
