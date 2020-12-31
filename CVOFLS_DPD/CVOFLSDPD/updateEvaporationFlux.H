
/*--------------------------------------------------------------------------------------------------------------------------*\
Description

The evaporation flux is here evaluated. It is calculated based on the fuel molar fraction at the interface:

 m = - rho*D*(Mi/MW)*gradX / (1- Wi)                   

with the correction of the molecular way (see MicrogravityDroplet solver) and the presence of the Stephan Flow

\*--------------------------------------------------------------------------------------------------------------------------*/



if(evaporation == true)
{

   Info<< "Evaporating flux evaluation... "<< endl;

   mEvapTot *= 0.0;  
   Marker *= 0.0;    

              

   for (unsigned int i=0;i<NLS;i++)              
     mevap[i] *= 0.0;                                


  volScalarField gradAlpha = Foam::mag(fvc::grad(alpha1));

  double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();


if(evaporationMechanism == "FicksLaw")
{
     interfaceSerialClass interfaceSerial(alpha1);

     volScalarField heatFlux =  fvc::laplacian(lambda,T) - cp*fvc::div(rhoPhi, T) ;

     for (unsigned int i=0;i<NLS;i++)
     {              
       volScalarField gradX = Foam::mag(fvc::grad(Xsat[i]));   
       volScalarField Ysati = Ysat[i]; 
       volScalarField Xsati = Xsat[i]; 
       double MWi = thermodynamicsMapXML->MW(LiquidSpeciesIndices(i));

        double rhoDmixiFixed = 0.0;
        double gradXFixed = 0.0;
        double Wfixed = 0.0;
        double Xfixed = 0.0;
        double MWfixed = 0.0;
        FugacityCoefficientMixFixed = 0.0;
        FugacityCoefficientPureFixed = 0.0;
        Zfixed = 0.0;
        Poyntingfixed = 0.0;
        gammafixed = 0.0;
   
        int nCellsFOUND = 0;



// Calculate the number of found gas cells 

         forAll(alpha1.internalField(), celli)
         {
            if (gradAlpha.internalField()[celli] > 0.)
            {  
                           
               double rhoDmixi = interfaceSerial.GetGasValue(rhoDmix[LiquidSpeciesIndices(i)], celli);
               double gradXi = interfaceSerial.GetGasValue(gradX, celli);
               double Wfixedi = interfaceSerial.GetLiquidValue(Ysati, celli);
               double Xfixedi = interfaceSerial.GetLiquidValue(Xsati, celli);
               double MWfixedi = interfaceSerial.GetLiquidValue(MWsat, celli);
               double FugacityCoefficientMixFixedi = 0.0;
               double FugacityCoefficientPureFixedi = 0.0;
               double Zfixedi = 0.0;
               double Poyntingfixedi = 0.0;
               double gammafixedi = 0.0;

               if (fugacity == "PengRobinson")
               {
                 FugacityCoefficientMixFixedi = interfaceSerial.GetGasValue(FugacityCoefficientMix[LiquidSpeciesIndices(i)], celli);
                 FugacityCoefficientPureFixedi = interfaceSerial.GetLiquidValue(FugacityCoefficientPure[i], celli);
               }

               if ( EoSforGasDensity == "PengRobinson")
                 Zfixedi = interfaceSerial.GetGasValue(Zgas, celli);

               if (poyntingCorrection == true)
                 Poyntingfixedi = interfaceSerial.GetGasValue(Poynting[i], celli);

               if (activity == "Unifac")
                 gammafixedi = interfaceSerial.GetLiquidValue(gamma[i], celli);

               if (rhoDmixi != -1. && Wfixedi != -1.)
               {
                 nCellsFOUND++;
                 rhoDmixiFixed += rhoDmixi;
                 gradXFixed += gradXi;
                 Xfixed += Xfixedi;
                 Wfixed += Wfixedi;
                 MWfixed += MWfixedi;

                 if (fugacity == "PengRobinson")
                 {
                   FugacityCoefficientMixFixed += FugacityCoefficientMixFixedi;
                   FugacityCoefficientPureFixed += FugacityCoefficientPureFixedi;
                 }

                 if ( EoSforGasDensity == "PengRobinson")
                   Zfixed += Zfixedi;

                 if (poyntingCorrection == true)
                   Poyntingfixed += Poyntingfixedi;

                 if (activity == "Unifac")
                   gammafixed += gammafixedi;
               }
            }
         }

// Calculate the average properties at the interface cells

           reduce(nCellsFOUND, sumOp<int>());
           reduce(rhoDmixiFixed, sumOp<double>());
           reduce(gradXFixed, sumOp<double>());
           reduce(Wfixed, sumOp<double>());
           reduce(Xfixed, sumOp<double>());
           reduce(MWfixed, sumOp<double>());

           if (fugacity == "PengRobinson")
           {
             reduce(FugacityCoefficientMixFixed, sumOp<double>());
             reduce(FugacityCoefficientPureFixed, sumOp<double>());
           }

           if ( EoSforGasDensity == "PengRobinson")
             reduce(Zfixed, sumOp<double>());

           if (poyntingCorrection == true)
             reduce(Poyntingfixed, sumOp<double>());

           if (activity == "Unifac")
             reduce(gammafixed, sumOp<double>());

           rhoDmixiFixed /= (double)nCellsFOUND;
           gradXFixed /= (double)nCellsFOUND;
           Wfixed /= (double)nCellsFOUND;
           Xfixed /= (double)nCellsFOUND;
           MWfixed /= (double)nCellsFOUND;

           if (fugacity == "PengRobinson")
           {
             FugacityCoefficientMixFixed /= (double)nCellsFOUND;
             FugacityCoefficientPureFixed /= (double)nCellsFOUND;
           }

           if ( EoSforGasDensity == "PengRobinson")
             Zfixed /=(double)nCellsFOUND;

           if (poyntingCorrection == true)
             Poyntingfixed /= (double)nCellsFOUND;

           if (activity == "Unifac")
             gammafixed /= (double)nCellsFOUND;

            Info<< "rhoDmixiFixed = " << rhoDmixiFixed << endl;
            Info<< "gradXFixed = " << gradXFixed << endl;
            Info<< "Xfixed = " << Xfixed << endl;
            Info<< "Wfixed = " << Wfixed << endl;
            Info<< "MWfixed = " << MWfixed << endl;

               if (fugacity == "PengRobinson")
               {
                   Info<< "FugacityCoefficientMixFixed = " << FugacityCoefficientMixFixed << endl;
                   Info<< "FugacityCoefficientPureFixed = " << FugacityCoefficientPureFixed << endl;
               }

               if ( EoSforGasDensity == "PengRobinson")
                  Info<< "Zfixed = " << Zfixed << endl;

               if (poyntingCorrection == true)
                  Info<< "Poyntingfixed = " << Poyntingfixed << endl;

               if (activity == "Unifac")
                  Info<< "gammafixed = " << gammafixed << endl;


// Calculate evaporation flux  

        int nCellsInterface = 0;     
        int nEvaporatingCells = 0;   
        int nBoilingCells = 0; 


         forAll(alpha1.internalField(), celli)
         {
            if (gradAlpha.internalField()[celli] > 0.)
            {   
               nCellsInterface++;    
               double gradAlphai = gradAlpha.internalField()[celli]; 

              //   if (T.internalField()[celli] < 0.999*Tboiling)
              //   {  
                     nEvaporatingCells++;
                     mevap[i].internalField()[celli] = -rhoDmixiFixed*(MWi/MWfixed)*gradXFixed*gradAlphai; 
                     mevap[i].internalField()[celli] /= ( 1. - Wfixed ) ;    
                     mEvapTot.internalField()[celli] += mevap[i].internalField()[celli];  
             //    }

/*
                 else
                 {
                    if (dhev[i].internalField()[celli] != 0.)
                    {
                      nBoilingCells++;
                      mevap[i].internalField()[celli] = -Foam::mag( heatFlux.internalField()[celli]/dhev[i].internalField()[celli] );
                      mEvapTot.internalField()[celli] += mevap[i].internalField()[celli];
                    }
                 }
*/

            }

         }

/*
           reduce(nCellsInterface, sumOp<int>());
           reduce(nEvaporatingCells, sumOp<int>());
           reduce(nBoilingCells, sumOp<int>());


          if (nCellsInterface != 0)
          {
            Info<< "Evaporation = " << 100.*((double)nEvaporatingCells/(double)nCellsInterface) << " %" ;
            Info<< "  Boiling  = " << 100.*((double)nBoilingCells/(double)nCellsInterface) << " %" << endl;
          }
*/

      }


}





    else if(evaporationMechanism  == "BoilingFlux")
    {
       for (unsigned int i=0;i<NLS;i++)
       {   
         volScalarField gradT = Foam::mag(fvc::grad(T));

         mevap[i] = -lambda*gradT*gradAlpha / dhev[i];
         mEvapTot += mevap[i]; 
       }

    }


    else if(evaporationMechanism  == "FixedFlux")
    {

      dimensionedScalar LiquidSurface = fvc::domainIntegrate(gradAlpha);
      mEvapTot = -(FixedFlux/LiquidSurface)*gradAlpha; 

    }



    else
    {
       Info << "Unknown evaporation mechanism!" << endl;
       abort();
    }



     mEvapTot.correctBoundaryConditions();


     dimensionedScalar EvaporationFlux = fvc::domainIntegrate(mEvapTot);
     Info << " m = " << EvaporationFlux.value() <<" kg/s";


     //double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
     //Info << "      done in:  " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl; 
}





/*--------------------------------------------------------------------------------------------------------------------------*\
Description

The interfacialAreaCoefficient is here evaluated. This is needed to distribute the evaporation flux among the interface cells.
Three choices are available:

- Linear: The source is multiplied by alpha1, in order to keep it zero when alpha1 goes to zero (and avoid alpha1<0 regions).
- Parabola: The source term is multiplied by alpha1*(1-alpha1), getting to zero when alpha1= 0,1.
- Circle : The same as parabola, but goes to zero with infinite derivative (better)

\*--------------------------------------------------------------------------------------------------------------------------*/


 volScalarField interfacialAreaCoefficient = 0.0*alpha1;
 volScalarField compressibleAreaCoefficient = 0.0*alpha1;

   if(interfacialAreaCoefficientMethod == "Linear")
   {
      interfacialAreaCoefficient = alpha1;
      compressibleAreaCoefficient = alpha2;
   }


   else if(interfacialAreaCoefficientMethod == "Parabola")
   {
      interfacialAreaCoefficient = alpha1*(1.0-alpha1);
      compressibleAreaCoefficient = interfacialAreaCoefficient;
      
      forAll(alpha1.internalField(), celli)
      {
          if(interfacialAreaCoefficient.internalField()[celli] == 0.0)
           interfacialAreaCoefficient.internalField()[celli] = alpha1.internalField()[celli];

          if(compressibleAreaCoefficient.internalField()[celli] == 0.0)
           compressibleAreaCoefficient.internalField()[celli] = alpha2.internalField()[celli];
      }
      
   }


   else if(interfacialAreaCoefficientMethod == "Circle")
   {
      interfacialAreaCoefficient = Foam::pow( alpha1*(1.0-alpha1), 0.5);
      compressibleAreaCoefficient = interfacialAreaCoefficient;

      forAll(alpha1.internalField(), celli)
      {
          if(interfacialAreaCoefficient.internalField()[celli] == 0.0)
           interfacialAreaCoefficient.internalField()[celli] = alpha1.internalField()[celli];

          if(compressibleAreaCoefficient.internalField()[celli] == 0.0)
           compressibleAreaCoefficient.internalField()[celli] = alpha2.internalField()[celli];
      }

   }


   else 
   {
      Info << "Unknown interfacialAreaCoefficientMethod!" << endl;
      abort();
   }





/*--------------------------------------------------------------------------------------------------------------------------*\
Description

The EvapCorrection factor is here evaluated. This is needed because multiplying the source by the interfacial coefficient,
the evaporation flux is modified. This coefficient, based on mass conservation, ensures the correct global evaporation flux,
though differently distributed among interfacial cells

\*--------------------------------------------------------------------------------------------------------------------------*/


  // Evaluate correction for the evaporation flux based on mass conservation

   totIntegral = fvc::domainIntegrate(mEvapTot + massLiquidCorrection);

   totIntegralAlpha = fvc::domainIntegrate(
                                              (mEvapTot + massLiquidCorrection)*interfacialAreaCoefficient
                                          );


   if(totIntegralAlpha.value() != 0.0)
   EvapCorrection = totIntegral/totIntegralAlpha;

   else
   EvapCorrection = 2.0;


   Info << "Evaporation Correction factor: " <<  EvapCorrection.value() << endl;


  // Evaluate correction for the compressibility source based on mass conservation

dimensionedScalar totIntegralCompressible = 0.0*totIntegral;
dimensionedScalar totIntegralAlphaCompressible = 0.0*totIntegralAlpha;

   totIntegralCompressible = fvc::domainIntegrate(compressibilitySourceTerm);

   totIntegralAlphaCompressible = fvc::domainIntegrate(
                                                         compressibilitySourceTerm*compressibleAreaCoefficient
                                                      );


   if(totIntegralAlphaCompressible.value() != 0.0)
   CompressibleCorrection = totIntegralCompressible/totIntegralAlphaCompressible;

   else
   CompressibleCorrection = 2.0;

   Info << "Compressible Correction factor: " << CompressibleCorrection.value() << endl;


