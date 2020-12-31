/*
   // CREATE THE TWO-PHASE MIXTURE OBJECT
    immiscibleIncompressibleTwoPhaseMixture mixture(U, rho1, rho2, nu1, nu2, cp1, cp2, lambda1, lambda2, pc);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());
*/


//- Creating field rho [kg/m3]
Info<< "Creating field rho\n" << endl;
volScalarField rho
(
   IOobject
   (
	"rho",
         runTime.timeName(),
         mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
   ),
      alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Setting pressure field
p = p_rgh + rho*gh;


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}


//mesh.setFluxRequired(p_rgh.name());
//mesh.setFluxRequired(alpha1.name());


// Setting pressureAdv field
pAdv = p_rghAdv + rho*gh;


label pAdvRefCell = 0;
scalar pAdvRefValue = 0.0;
setRefCell
(
    pAdv,
    p_rghAdv,
    pimple.dict(),
    pAdvRefCell,
    pAdvRefValue
);

if (p_rghAdv.needReference())
{
    pAdv += dimensionedScalar
    (
        "pAdv",
        pAdv.dimensions(),
        pAdvRefValue - getRefCellValue(pAdv, pAdvRefCell)
    );
    p_rghAdv = pAdv - rho*gh;
}


//mesh.setFluxRequired(p_rghAdv.name());


dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));


//- Allocating memory for total concentration [kmol/m3]
Info<< "Creating field cTot\n" << endl;
volScalarField cTot
(
    IOobject
    (
        "mixture_cTot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
      cTot1*alpha1 + cTot2*alpha2
);


//- Allocating memory for mixture dynamic viscosity [kg/m/s]
Info<< "Creating field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "mixture_mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputMixtureProperties
    ),
      alpha1*mu1 + alpha2*mu2
);

//- Allocating memory for mixture cinematic viscosity [kg/m/s]
Info<< "Creating field nu\n" << endl;
volScalarField nu
(
    IOobject
    (
        "mixture_nu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputMixtureProperties
    ),
      alpha1*nu1 + alpha2*nu2
);

//- Allocating memory for mixture thermal conductivity [W/m/K]
Info<< "Creating field lambda mixture\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "mixture_lambda",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputMixtureProperties
    ),
      alpha1*lambda1 + alpha2*lambda2
);



//- Allocating memory for mixture mean specific heat [J/kg/K]
Info<< "Creating field Cp mixture\n" << endl;
volScalarField cp
(
    IOobject
    (
        "mixture_Cp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputMixtureProperties
    ),
      alpha1*cp1 + alpha2*cp2
);

//- Allocating memory for mixture mean specific heat [J/kg/K]
Info<< "Creating field psi mixture\n" << endl;
volScalarField psi
(
    IOobject
    (
        "mixture_psi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
      alpha1*psi1 + alpha2*psi2
);


PtrList<volScalarField> Dmix(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> rhoDmix(thermodynamicsMapXML->NumberOfSpecies());

Info << "   Creating species diffusivities: " << endl;
for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT
	);

	{
		Dmix.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"mixture_Dmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("Dmix", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);


		rhoDmix.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"mixture_rhoDmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("rhoDmix", dimensionSet(1, -1, -1, 0, 0), 1e-10)
			)
		);

         }

}
           
         for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
         {
           Dmix[i] = Dmix2[i]*alpha2 +  Dmix1[i]*alpha1 ;
           rhoDmix[i] = rho2*Dmix2[i]*alpha2 +  rho1*Dmix1[i]*alpha1 ;
         }


         for (unsigned int j=0;j<NLS;j++)
         { 
           Ysat[j] = Keq[j]*alpha1 + Ysat[j]*alpha2;
           Xsat[j] = Keq[j]*alpha1 + Xsat[j]*alpha2;

           double MWj = thermodynamicsMapXML->MW(LiquidSpeciesIndices(j));
           Xsat[j] = Ysat[j]/( MWj * (Ysat[j]/MWj + (1.-Ysat[j])/MWinert.value()) );
         }


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



// This vector contains the gas phase mass fractions at the interface of the non-soluble species (O2, CO, CO2 ...)

    std::vector<double> GaseousInterfaceComposition(NGS);
    std::vector<double> GaseousInterfaceRatio(NGS);

    for (label i=0; i<NGS; i++)
    {
       volScalarField& Yi = Y[GasSpeciesIndices(i)];

                interfaceSerialClass interfaceSerial(alpha1);
                volScalarField gradAlpha = Foam::mag(fvc::grad(alpha1));


                double YGasSurfaceFixed = 0.0;
                int nCellsFOUND = 0;

                forAll(alpha1.internalField(), celli)
                {
                    if (gradAlpha.internalField()[celli] > 0.)
                    {                             
                       double YGasSurfacei = interfaceSerial.GetGasValue(Yi, celli);  

                       if (YGasSurfacei != -1.)
                       {                     
                         nCellsFOUND++;
                         YGasSurfaceFixed += YGasSurfacei;
                       }

                    }
                }


                reduce(nCellsFOUND, sumOp<int>());
                reduce(YGasSurfaceFixed, sumOp<double>());



                YGasSurfaceFixed /= (double)nCellsFOUND;


                GaseousInterfaceComposition[i] =   YGasSurfaceFixed;
    }

// Find the inert mole fraction at the interface
    {
                volScalarField& Yi = Y[inertIndex];

                interfaceSerialClass interfaceSerial(alpha1);
                volScalarField gradAlpha = Foam::mag(fvc::grad(alpha1));

                double YInertSurfaceFixed = 0.0;
                int nCellsFOUND = 0;

                forAll(alpha1.internalField(), celli)
                {
                    if (gradAlpha.internalField()[celli] > 0.)
                    {                             
                       double YInertSurfacei = interfaceSerial.GetGasValue(Yi, celli);  

                       if (YInertSurfacei != -1.)
                       {                     
                         nCellsFOUND++;
                         YInertSurfaceFixed += YInertSurfacei;
                       }

                    }
                }

                reduce(nCellsFOUND, sumOp<int>());
                reduce(YInertSurfaceFixed, sumOp<double>());

                YInertSurfaceFixed /= (double)nCellsFOUND;


                if (YInertSurfaceFixed == 0.)
                    Info << " Warning: the nitrogen is zero in the gas phase!! " << endl;

                for (label i=0; i<NGS; i++)
                    GaseousInterfaceRatio[i] =  GaseousInterfaceComposition[i]/YInertSurfaceFixed;

                
   }




    for (label i=0; i<NGS; i++)
    {
       volScalarField& Ygasi = Ygas[i];

       Ygasi = alpha1*(1.0-sumYi)*GaseousInterfaceRatio[i]/(1.0+GaseousInterfaceRatio[i]) + alpha2*Ygasi;  
    }





// Create surface mass flux [kg/m2 s]
Info<< "Create surface mass flux...\n" << endl;
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);
rhoPhi.oldTime();

// Create surface mass flux [kg/m2 s]
Info<< "Create surface mass flux...\n" << endl;
surfaceScalarField rhoPhiAdv
(
    IOobject
    (
        "rhoPhiAdv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phiAdv
);
rhoPhiAdv.oldTime();


//Create alpha flux (MULES flux from previous time-step)
Info<< "Create alpha flux...\n" << endl;
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);
alphaPhi.oldTime();


// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;


#include "createMRF.H"


