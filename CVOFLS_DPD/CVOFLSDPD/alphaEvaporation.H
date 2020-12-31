/*---------------------------------------------------------------------------------------------------------------------------*\
Description

The source term for alpha equation are here included in explicit way. The sources are:

  - mEvapTot for the liquid evaporation
  - massLiquidCorrection, which accounts for the eventual variation of thickness of the interface
  - compressibilitySourceTerm, which closes the mass conservation accounting for liquid rho variation and eventual errors.

\*----------------------------------------------------------------------------------------------------------------------------*/

		fvScalarMatrix alphaSourceEqn
		(
			fvm::ddt(alpha1)	    
			==
	    		EvapCorrection*interfacialAreaCoefficient*( 
                                                                    ( mEvapTot + massLiquidCorrection ) / rho1
                                                                  )
                      - CompressibleCorrection*compressibleAreaCoefficient*compressibilitySourceTerm / rho1
                       
		);




		alphaSourceEqn.relax();
		fvOptions.constrain(alphaSourceEqn);
		alphaSourceEqn.solve(mesh.solver("T"));
		fvOptions.correct(alpha1);


/*---------------------------------------------------------------------------------------------------------------------------*\
Description

  - If gravity is zero: it simply cleans eventual regions with alpha>1 or alpha<0
  - If gravity is not zero: the alpha1 field is further cleaned with max and min thresholds 0.99 and 0.01.

  Of course the mass is modified. It is recovered and stored in massLiquidCorrection.

\*----------------------------------------------------------------------------------------------------------------------------*/



dimensionedScalar InitialLiquidMass = fvc::domainIntegrate(rho1*alpha1);

if ( (mag(g)).value() == 0.)
{

  forAll(alpha1.internalField(), celli)
  {

   if (alpha1.internalField()[celli] > 1.0)
       alpha1.internalField()[celli] = 1.0;

   if (alpha1.internalField()[celli] < 0.0)
       alpha1.internalField()[celli] = 0.0;

  }

}



else
{

  forAll(alpha1.internalField(), celli)
  {

   if (alpha1.internalField()[celli] > 0.99)
       alpha1.internalField()[celli] = 1.0;

   if (alpha1.internalField()[celli] < 0.01)
       alpha1.internalField()[celli] = 0.0;

  }

}


alpha1.correctBoundaryConditions();


// Evaluation of massLiquidCorrection

dimensionedScalar FinalLiquidMass = fvc::domainIntegrate(rho1*alpha1);
dimensionedScalar LiquidMassDifference = InitialLiquidMass-FinalLiquidMass;
dimensionedScalar LiquidSurface = fvc::domainIntegrate( Foam::mag( fvc::grad(alpha1) ) );

  massLiquidCorrection = (
                           LiquidMassDifference/(LiquidSurface*runTime.deltaT())
                         )*Foam::mag(fvc::grad(alpha1));

alpha1.correctBoundaryConditions();

alpha2 = scalar(1) - alpha1;





// Update Diameter

  Info<< "Diameter evaluation... " ;
       
  LiquidVolume = fvc::domainIntegrate(alpha1)*72.0;
  Diameter = 2.0*Foam::pow(1.5*LiquidVolume/PhysicalConstants::pi,1.0/3.0);
  //Diameter = Foam::pow(6.*LiquidVolume/PhysicalConstants::pi + pow(FiberDiameter, 3.0) , 1./3.);

  Info << "D =  " << Diameter.value()*1000. << " mm " <<  endl;


// Save the liquid mass value before advection & fluid properties evaluation

dimensionedScalar MassAfterEvaporation = fvc::domainIntegrate(alpha1*rho1);
