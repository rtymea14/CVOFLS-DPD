
/*---------------------------------------------------------------------------------------------------------------------------*\
Description

 Based on Banerjee et.al. (2013) the saturation mass fraction equation is solved. It is assumed that all the liquid phase has the
 saturation interfacial concentration and it diffuses from all the volume. Afterwards the liquid concentration is posed = 1 
 (for monocomponent cases). 

\*----------------------------------------------------------------------------------------------------------------------------*/



// Mass Fluxes
#include "correctDiffusionFluxes.H"


if(speciesEquations == true)
{

 double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();


 // Convection discretization schemes
 tmp<fv::convectionScheme<scalar> > mvConvection
 (
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        rhoPhi,
        mesh.divScheme("div(rhoPhi,Yi)")
    )
 );



    
  volScalarField Yt = 0.0*Y[0];


// Solving transport equations for the saturation mass fraction of liquid species Ysat

    for (label i=0; i<NLS; i++)
    {
           	volScalarField& Yi = Y[LiquidSpeciesIndices(i)];
	     	volScalarField& rhoDmixi = rhoDmix[LiquidSpeciesIndices(i)];
	    	volScalarField& Keqi = Keq[i];
           	volScalarField& Ysati = Ysat[i];
           	volScalarField& Xsati = Xsat[i];
                dimensionedScalar MWi( "Mi", dimensionSet(1,0,0,0,-1,0,0),thermodynamicsMapXML->MW(LiquidSpeciesIndices(i)) );

                Ysati = alpha1*Keqi + alpha2*Ysati;


/*---------------------------------------------------------------------------------------------------------------------------*\
Description

 The diffusion fluxes must be corrected with the molecular weight, because the diffusivity coefficients are molar-based. Having 
 mass-fraction equation instead of a molar-fraction one, a MW correction is needed ( Transport Phenomena, p.534 )

\*----------------------------------------------------------------------------------------------------------------------------*/


              volScalarField correctionMolarFluxes = rhoDmixi;

              if(mwCorrectionInDiffusionFluxes == true)
              {
                  correctionMolarFluxes = rhoDmixi*
                                          1./(MWsat*MWinert)*
                                          1./ Foam::pow( Ysati/MWi+(1.-Ysati)/MWinert , 2.0);
              }

     


// Solve liquid species equation

			fvScalarMatrix YsatiEqn
			(
				fvm::ddt(rho, Ysati)
			      + mvConvection->fvmDiv(rhoPhi, Ysati)
			      - fvm::laplacian(correctionMolarFluxes, Ysati)
				== 
                             // - fvm::div(Jc,Ysati, "div(Jc,Ysati)")
			       fvOptions(rho, Ysati)
			);


			// Solve
			YsatiEqn.relax();
			fvOptions.constrain(YsatiEqn);
			YsatiEqn.solve(mesh.solver("T"));
			fvOptions.correct(Ysati);	


                Ysati = alpha1*Keqi + alpha2*Ysati;

                Xsati = Ysati/( MWi * (Ysati/MWi + (1.-Ysati)/MWinert) );


                Yi = alpha1 + Ysati*alpha2;
                Yi.correctBoundaryConditions();

		Yi.max(0.0);
		Yt += Yi;
                     
    }





// Solving transport equations for the gas species (not entering the liquid phase, like O2, CO2, reaction products...)

    for (label i=0; i<NGS; i++)
    {
       volScalarField& Yi = Y[GasSpeciesIndices(i)];
       volScalarField& rhoDmixi = rhoDmix[GasSpeciesIndices(i)];
       volScalarField& Ygasi = Ygas[i];


        Ygasi = alpha1*(1.0-sumYi)*GaseousInterfaceRatio[i]/(1.0+GaseousInterfaceRatio[i]) + alpha2*Ygasi;                        


			fvScalarMatrix YgasiEqn
			(
				fvm::ddt(rho, Ygasi)
			      + mvConvection->fvmDiv(rhoPhi, Ygasi)
			      - fvm::laplacian(rhoDmixi, Ygasi)
				== 
			       fvOptions(rho, Ygasi)
			);


			// Solve
			YgasiEqn.relax();
			fvOptions.constrain(YgasiEqn);
			YgasiEqn.solve(mesh.solver("Yi"));
			fvOptions.correct(Ygasi);	


           Ygasi = alpha1*(1.0-sumYi)*GaseousInterfaceRatio[i]/(1.0+GaseousInterfaceRatio[i]) + alpha2*Ygasi;   

                 Yi =  alpha2*Ygasi;
                 Yi.correctBoundaryConditions();                        


			
			// Sum of mass fractions
		    	Yi.max(0.0);
		   	Yt += Yi;


    }




    Info << "Inert species is " << Y[inertIndex].name() << " with local index equal to " << inertIndex << endl;
    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
	
    Info << "Transport equations of species solved in " << tEnd - tStart << " s " << endl;


}




