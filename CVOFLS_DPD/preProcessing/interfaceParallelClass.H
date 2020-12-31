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


SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceParallelClass_H
#define interfaceParallelClass_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceParallelClass
{
         volScalarField& alpha1_;
         List<scalarField>& IDsList_;
         List<scalarField>& alphaList_;
         List<tensorField>& NeighboursList_;


public:

    // Constructors

        //- Construct from volume fraction 
        interfaceParallelClass(volScalarField& alpha1, List<scalarField>& IDsList, List<scalarField>& alphaList, List<tensorField>& NeighboursList);


    // Member Functions

       double GetGasValue(List<scalarField>& fieldList, int celli);
       double GetLiquidValue(List<scalarField>& fieldList, int celli);
       int GetGasCell(int celli);
       int GetLiquidCell(int celli);
       std::vector<double> markers(int celli);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceParallelClass::interfaceParallelClass(volScalarField& alpha1, List<scalarField>& IDsList, List<scalarField>& alphaList, List<tensorField>& NeighboursList):

    alpha1_(alpha1),
    IDsList_(IDsList),
    alphaList_(alphaList),
    NeighboursList_(NeighboursList)
{

}


// ********************************************************************* MEMBER FUNCTIONS ************************************************************************************* // 





// ******************************************************************************************************** // 
// ******************************************  GET GAS VALUE ********************************************** // 
// ******************************************************************************************************** // 

double Foam::interfaceParallelClass::GetGasValue(List<scalarField>& fieldList, int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];

              int markerProcessor = Pstream::myProcNo();
              int markerCell = celli;
              int NumberOfIterations = 0;

              int MaxNumberOfIterations = 10;

              while ( thresholdAlphaValue > 1e-3 )
              {
                     NumberOfIterations++;

                     int NeighboursSize = 0;
                     for (unsigned int i=0; i<9; i++)
                     { 
                       if (NeighboursList_[markerProcessor][markerCell].component(i) != -1)
                       NeighboursSize++;
                     }


                 std::vector<double> alphaValues(NeighboursSize);
                 std::vector<int> nProcessor(NeighboursSize);            // Processor at which I have that neighbour cell
                 std::vector<int> ncelli(NeighboursSize);                // Cell at which I have that neighbour cell

                 for (int i=0;i<NeighboursSize;i++)
                 {
                   int pointedNeighbour = NeighboursList_[markerProcessor][markerCell].component(i);

                   int markerCellInSameProcessor = 0;

	             forAll(IDsList_[markerProcessor], cellj)
                     {
                        if (IDsList_[markerProcessor][cellj] == pointedNeighbour)
                        {   
                           markerCellInSameProcessor++;                       
                           alphaValues[i] = alphaList_[markerProcessor][cellj];  
                           nProcessor[i] = markerProcessor; 
                           ncelli[i] = cellj; 
                           break; 
                        }
                     } 


                   if (markerCellInSameProcessor == 0)
                   {   
	              forAll(IDsList_, jProc)
	              {
                        if(jProc != markerProcessor)
                        {
	                   forAll(IDsList_[jProc], cellj)
                           {
                              if (IDsList_[jProc][cellj] == pointedNeighbour)
                              {
                                alphaValues[i] = alphaList_[jProc][cellj];  
                                nProcessor[i] = jProc; 
                                ncelli[i] = cellj; 
                                break; 
                              }
                           }
                        }

                       } 
 
                    }
                  }



                 double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                 double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>MaxNumberOfIterations)
                 {
                   thresholdAlphaValue = 0.0;
                   markerProcessor=-2000;
                 }


                 else
                 {
                   thresholdAlphaValue = minAlpha;

                   int minimumAlphaPosition = std::distance(alphaValues.begin(), min_element(alphaValues.begin(), alphaValues.end()));
                   markerProcessor = nProcessor[minimumAlphaPosition];
                   markerCell = ncelli[minimumAlphaPosition];
                 }


               }


                 if (markerProcessor != -2000)
                 {
                   double GasValue = fieldList[markerProcessor][markerCell];
                   return GasValue;
                 }


                 else
                 {
                   return -1.0;
                 }

}



// ******************************************************************************************************** // 
// ******************************************  GET GAS CELL *********************************************** // 
// ******************************************************************************************************** // 



int Foam::interfaceParallelClass::GetGasCell(int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];

              int markerProcessor = Pstream::myProcNo();
              int markerCell = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue > 1e-3 )
              {
                     NumberOfIterations++;

                     int NeighboursSize = 0;
                     for (unsigned int i=0; i<9; i++)
                     { 
                       if (NeighboursList_[markerProcessor][markerCell].component(i) != -1)
                       NeighboursSize++;
                     }


                 std::vector<double> alphaValues(NeighboursSize);
                 std::vector<int> nProcessor(NeighboursSize);            // Processor at which I have that neighbour cell
                 std::vector<int> ncelli(NeighboursSize);                // Cell at which I have that neighbour cell

                 for (int i=0;i<NeighboursSize;i++)
                 {
                   int pointedNeighbour = NeighboursList_[markerProcessor][markerCell].component(i);

                   int markerCellInSameProcessor = 0;

	             forAll(IDsList_[markerProcessor], cellj)
                     {
                        if (IDsList_[markerProcessor][cellj] == pointedNeighbour)
                        {   
                           markerCellInSameProcessor++;                       
                           alphaValues[i] = alphaList_[markerProcessor][cellj];  
                           nProcessor[i] = markerProcessor; 
                           ncelli[i] = cellj; 
                           break; 
                        }
                     } 


                   if (markerCellInSameProcessor == 0)
                   {   
	              forAll(IDsList_, jProc)
	              {
                        if(jProc != markerProcessor)
                        {
	                   forAll(IDsList_[jProc], cellj)
                           {
                              if (IDsList_[jProc][cellj] == pointedNeighbour)
                              {
                                alphaValues[i] = alphaList_[jProc][cellj];  
                                nProcessor[i] = jProc; 
                                ncelli[i] = cellj; 
                                break; 
                              }
                           }
                        }

                       } 
 
                    }
                  }

                 double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                 double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());


                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 0.0;
                   markerProcessor = -2000;
                 }

                 else
                 {
                   thresholdAlphaValue = minAlpha;

                   int minimumAlphaPosition = std::distance(alphaValues.begin(), min_element(alphaValues.begin(), alphaValues.end()));
                   markerProcessor = nProcessor[minimumAlphaPosition];
                   markerCell = ncelli[minimumAlphaPosition];
                 }

               }


                 if (markerProcessor != -2000)
                 {
                   double GasCell = IDsList_[markerProcessor][markerCell];
                   return GasCell;
                 }

                 else
                 {
                   return -1;
                 }
  
}




// ******************************************************************************************************** // 
// ******************************************  GET GAS MARKERS ********************************************* // 
// ******************************************************************************************************** // 


std::vector<double> Foam::interfaceParallelClass::markers(int celli)
{
              const fvMesh& mesh = alpha1_.mesh();


              double thresholdAlphaValue = alpha1_.internalField()[celli];

              int markerProcessor = Pstream::myProcNo();
              int markerCell = celli;
              int NumberOfIterations = 0;

              std::vector<double> markers(10);

                 for (int i=0;i<10;i++)
                   markers[i] = -2000.0;



              while ( thresholdAlphaValue > 1e-3 )
              {
                     NumberOfIterations++;

                     int NeighboursSize = 0;
                     for (unsigned int i=0; i<9; i++)
                     { 
                       if (NeighboursList_[markerProcessor][markerCell].component(i) != -1)
                       NeighboursSize++;
                     }



                 std::vector<double> alphaValues(NeighboursSize);
                 std::vector<int> nProcessor(NeighboursSize);            // Processor at which I have that neighbour cell
                 std::vector<int> ncelli(NeighboursSize);                // Cell at which I have that neighbour cell

                 for (int i=0;i<NeighboursSize;i++)
                 {
                   int pointedNeighbour = NeighboursList_[markerProcessor][markerCell].component(i);

                   int markerCellInSameProcessor = 0;

	             forAll(IDsList_[markerProcessor], cellj)
                     {
                        if (IDsList_[markerProcessor][cellj] == pointedNeighbour)
                        {   
                           markerCellInSameProcessor++;                       
                           alphaValues[i] = alphaList_[markerProcessor][cellj];  
                           nProcessor[i] = markerProcessor; 
                           ncelli[i] = cellj; 
                           break; 
                        }
                     } 


                   if (markerCellInSameProcessor == 0)
                   {   
	              forAll(IDsList_, jProc)
	              {
                        if(jProc != markerProcessor)
                        {
	                   forAll(IDsList_[jProc], cellj)
                           {
                              if (IDsList_[jProc][cellj] == pointedNeighbour)
                              {
                                alphaValues[i] = alphaList_[jProc][cellj];  
                                nProcessor[i] = jProc; 
                                ncelli[i] = cellj; 
                                break; 
                              }
                           }
                        }

                       } 
 
                    }
                  }

                 double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                 double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());


                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 0.0;
                   markerProcessor = -2000;
                 }

                 else
                 {
                   thresholdAlphaValue = minAlpha;

                   int minimumAlphaPosition = std::distance(alphaValues.begin(), min_element(alphaValues.begin(), alphaValues.end()));
                   markerProcessor = nProcessor[minimumAlphaPosition];
                   markerCell = ncelli[minimumAlphaPosition];
                 }

                  markers[NumberOfIterations-1] = markerCell;

               }

                   return markers;

}


// ******************************************************************************************************** // 
// ******************************************  GET LIQUID VALUE ******************************************* // 
// ******************************************************************************************************** // 



double Foam::interfaceParallelClass::GetLiquidValue(List<scalarField>& fieldList, int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];

              int markerProcessor = Pstream::myProcNo();
              int markerCell = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue < 0.999 )
              {
                     NumberOfIterations++;

                     int NeighboursSize = 0;
                     for (unsigned int i=0; i<9; i++)
                     { 
                       if (NeighboursList_[markerProcessor][markerCell].component(i) != -1)
                       NeighboursSize++;
                     }


                 std::vector<double> alphaValues(NeighboursSize);
                 std::vector<int> nProcessor(NeighboursSize);            // Processor at which I have that neighbour cell
                 std::vector<int> ncelli(NeighboursSize);                // Cell at which I have that neighbour cell

                 for (int i=0;i<NeighboursSize;i++)
                 {
                   int pointedNeighbour = NeighboursList_[markerProcessor][markerCell].component(i);

                   int markerCellInSameProcessor = 0;

	             forAll(IDsList_[markerProcessor], cellj)
                     {
                        if (IDsList_[markerProcessor][cellj] == pointedNeighbour)
                        {   
                           markerCellInSameProcessor++;                       
                           alphaValues[i] = alphaList_[markerProcessor][cellj];  
                           nProcessor[i] = markerProcessor; 
                           ncelli[i] = cellj; 
                           break; 
                        }
                     } 


                   if (markerCellInSameProcessor == 0)
                   {   
	              forAll(IDsList_, jProc)
	              {
                        if(jProc != markerProcessor)
                        {
	                   forAll(IDsList_[jProc], cellj)
                           {
                              if (IDsList_[jProc][cellj] == pointedNeighbour)
                              {
                                alphaValues[i] = alphaList_[jProc][cellj];  
                                nProcessor[i] = jProc; 
                                ncelli[i] = cellj; 
                                break; 
                              }
                           }
                        }

                       } 
 
                    }
                  }

                 double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                 double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 1.0;
                   markerProcessor=-2000;
                 }

                 else
                 {
                   thresholdAlphaValue = maxAlpha;

                   int minimumAlphaPosition = std::distance(alphaValues.begin(), max_element(alphaValues.begin(), alphaValues.end()));
                   markerProcessor = nProcessor[minimumAlphaPosition];
                   markerCell = ncelli[minimumAlphaPosition];
                 }

               }


                 if (markerProcessor != -2000)
                 {
                   double LiquidValue = fieldList[markerProcessor][markerCell];
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


int Foam::interfaceParallelClass::GetLiquidCell(int celli)
{
              const fvMesh& mesh = alpha1_.mesh();

              double thresholdAlphaValue = alpha1_.internalField()[celli];

              int markerProcessor = Pstream::myProcNo();
              int markerCell = celli;
              int NumberOfIterations = 0;

              while ( thresholdAlphaValue < 0.999 )
              {
                     NumberOfIterations++;

                     int NeighboursSize = 0;
                     for (unsigned int i=0; i<9; i++)
                     { 
                       if (NeighboursList_[markerProcessor][markerCell].component(i) != -1)
                       NeighboursSize++;
                     }


                 std::vector<double> alphaValues(NeighboursSize);
                 std::vector<int> nProcessor(NeighboursSize);            // Processor at which I have that neighbour cell
                 std::vector<int> ncelli(NeighboursSize);                // Cell at which I have that neighbour cell

                 for (int i=0;i<NeighboursSize;i++)
                 {
                   int pointedNeighbour = NeighboursList_[markerProcessor][markerCell].component(i);

                   int markerCellInSameProcessor = 0;

	             forAll(IDsList_[markerProcessor], cellj)
                     {
                        if (IDsList_[markerProcessor][cellj] == pointedNeighbour)
                        {   
                           markerCellInSameProcessor++;                       
                           alphaValues[i] = alphaList_[markerProcessor][cellj];  
                           nProcessor[i] = markerProcessor; 
                           ncelli[i] = cellj; 
                           break; 
                        }
                     } 


                   if (markerCellInSameProcessor == 0)
                   {   
	              forAll(IDsList_, jProc)
	              {
                        if(jProc != markerProcessor)
                        {
	                   forAll(IDsList_[jProc], cellj)
                           {
                              if (IDsList_[jProc][cellj] == pointedNeighbour)
                              {
                                alphaValues[i] = alphaList_[jProc][cellj];  
                                nProcessor[i] = jProc; 
                                ncelli[i] = cellj; 
                                break; 
                              }
                           }
                        }

                       } 
 
                    }
                  }

                 double minAlpha = *min_element(alphaValues.begin(), alphaValues.end());
                 double maxAlpha = *max_element(alphaValues.begin(), alphaValues.end());

                 if (minAlpha == maxAlpha || NumberOfIterations>10)
                 {
                   thresholdAlphaValue = 1.0;
                   markerProcessor=-2000;
                 }

                 else
                 {
                   thresholdAlphaValue = maxAlpha;

                   int minimumAlphaPosition = std::distance(alphaValues.begin(), max_element(alphaValues.begin(), alphaValues.end()));
                   markerProcessor = nProcessor[minimumAlphaPosition];
                   markerCell = ncelli[minimumAlphaPosition];
                 }

               }


                 if (markerProcessor != -2000)
                 {
                   double LiquidCell = IDsList_[markerProcessor][markerCell];
                   return LiquidCell;
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
