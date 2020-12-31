/*---------------------------------------------------------------------------------------------------------------------------*\
Description

The mixture transport properties are updated base on new liquid and gas values and the new alpha1 distribution.

\*----------------------------------------------------------------------------------------------------------------------------*/


mu = alpha1*mu1 + alpha2*mu2;

nu = alpha1*nu1 + alpha2*nu2;

rho = alpha1*rho1 + alpha2*rho2;

cTot = cTot1*alpha1 + cTot2*alpha2;

rhoPhi = fvc::interpolate(rho)*phi; 

rhoPhiAdv = fvc::interpolate(rho)*phiAdv; 
 
cp = alpha1*cp1 + cp2*alpha2;

psi =  alpha1*psi1 + alpha2*psi2;

lambda = lambda1*alpha1 + lambda2*alpha2;



for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{

    Dmix[i] = Dmix2[i]*alpha2 +  Dmix1[i]*alpha1 ;
    rhoDmix[i] = rho2*Dmix2[i]*alpha2 +  rho1*Dmix1[i]*alpha1 ;

}


// Update saturation molecular weight

sumXiMWi *= 0.0;
sumXi *= 0.0;
sumYi *= 0.0;

for (unsigned int i=0;i<NLS;i++)
{  
   dimensionedScalar MWi( "Mi", dimensionSet(1,0,0,0,-1,0,0),thermodynamicsMapXML->MW(LiquidSpeciesIndices(i)) );

   sumXiMWi += MWi*Xsat[i];
   sumXi += Xsat[i];
   sumYi += Ysat[i];
}

MWsat = sumXiMWi + MWinert*(1. - sumXi);


/*
if (homogeneousReactions == true)
{ 
    std::vector<double> GaseousInterfaceMoleFractions(NGS);
   // Redistribution of inert mass fraction among the various inert species

   // All gas species apart from nitrogen
    for (label i=0; i<NGS; i++)
    {
       volScalarField& Yi = Y[GasSpeciesIndices(i)];
       volScalarField& Xi = X[GasSpeciesIndices(i)];

                interfaceSerialClass interfaceSerial(alpha1);
                volScalarField gradAlpha = Foam::mag(fvc::grad(alpha1));

                double YGasSurfaceFixed = 0.0;
                double XGasSurfaceFixed = 0.0;
                int nCellsFOUND = 0;

                forAll(alpha1.internalField(), celli)
                {
                    if (gradAlpha.internalField()[celli] > 0.)
                    {                             
                       double YGasSurfacei = interfaceSerial.GetGasValue(Yi, celli);  
                       double XGasSurfacei = interfaceSerial.GetGasValue(Xi, celli);  

                       if (YGasSurfacei != -1.)
                       {                     
                         nCellsFOUND++;
                         YGasSurfaceFixed += YGasSurfacei;
                         XGasSurfaceFixed += XGasSurfacei;
                       }

                    }
                }

                reduce(nCellsFOUND, sumOp<int>());
                reduce(YGasSurfaceFixed, sumOp<double>());
                reduce(XGasSurfaceFixed, sumOp<double>());

                YGasSurfaceFixed /= (double)nCellsFOUND;
                XGasSurfaceFixed /= (double)nCellsFOUND;

                GaseousInterfaceComposition[i] =   YGasSurfaceFixed;
                GaseousInterfaceMoleFractions[i] =   XGasSurfaceFixed;
    }

   // Nitrogen
    {
                volScalarField& Yi = Y[inertIndex];
                volScalarField& Xi = X[inertIndex];

                interfaceSerialClass interfaceSerial(alpha1);
                volScalarField gradAlpha = Foam::mag(fvc::grad(alpha1));

                double YInertSurfaceFixed = 0.0;
                double XInertSurfaceFixed = 0.0;
                int nCellsFOUND = 0;

                forAll(alpha1.internalField(), celli)
                {
                    if (gradAlpha.internalField()[celli] > 0.)
                    {                             
                       double YInertSurfacei = interfaceSerial.GetGasValue(Yi, celli);  
                       double XInertSurfacei = interfaceSerial.GetGasValue(Xi, celli);  

                       if (YInertSurfacei != -1.)
                       {                     
                         nCellsFOUND++;
                         YInertSurfaceFixed += YInertSurfacei;
                         XInertSurfaceFixed += XInertSurfacei;
                       }

                    }
                }

                reduce(nCellsFOUND, sumOp<int>());
                reduce(YInertSurfaceFixed, sumOp<double>());
                reduce(XInertSurfaceFixed, sumOp<double>());

                YInertSurfaceFixed /= (double)nCellsFOUND;
                XInertSurfaceFixed /= (double)nCellsFOUND;

                MWinert *= 0.0;

                for (label i=0; i<NGS; i++)
                    MWinert.value() += GaseousInterfaceMoleFractions[i]*thermodynamicsMapXML->MW(GasSpeciesIndices(i));

                    MWinert.value() += XInertSurfaceFixed*thermodynamicsMapXML->MW(inertIndex);

                    Info<< "MWinert = " << MWinert.value() << endl;
                
   }

}

*/








