
/*---------------------------------------------------------------------------------------------------------------------------*\
Description

 When the alphaAdvection is activated, the droplet surface tends to be destroyed by the high evaporation velocity. So a secondary
 velocity field is included and solved, where the vaporization velocity is neglected. This is used ONLY to advect the interface.
 Of course a pressure coupling is also added to enforce continuity satisfaction.

 For ALL the other fields ( temperature, species ...) the ordinary velocity (with vaporization velocity contribution) is used.

 This is used only in evaporating cases in gravity conditions (in Microgravity alpha1 is not advected)

\*----------------------------------------------------------------------------------------------------------------------------*/

  if ( (mag(g)).value() != 0. && evaporation == true )
  {
       #include "UAdvEqn.H"

       while (pimple.correct())
       {
          #include "pAdvEqn.H"
       }
  }



  else
  {
    UAdv = U;
    p_rghAdv = p_rgh;
    pAdv = p;
    phiAdv = phi;
    rhoPhiAdv = rhoPhi;
  }
