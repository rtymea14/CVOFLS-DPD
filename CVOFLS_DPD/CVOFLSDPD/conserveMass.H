
/*---------------------------------------------------------------------------------------------------------------------------*\
Description

 At the end of alphaEvaporation the liquid mass has been stored. After the alpha advection and the variation of transport
 properties (rhoL in particular) the mass is here stored again and the difference included in alpha sources to enforce mass liquid
 conservation. 

 An immediate effect is the varation of droplet size when liquid density varies.

\*----------------------------------------------------------------------------------------------------------------------------*/


  dimensionedScalar MassAfterConvection = fvc::domainIntegrate(alpha1*rho1);

  dimensionedScalar dmdt = (MassAfterConvection-MassAfterEvaporation)/runTime.deltaT();
  dimensionedScalar LiquidS = fvc::domainIntegrate( Foam::mag( fvc::grad(alpha1) ) );

  compressibilitySourceTerm = (
                                  dmdt/LiquidS
                              )*Foam::mag(fvc::grad(alpha1));
