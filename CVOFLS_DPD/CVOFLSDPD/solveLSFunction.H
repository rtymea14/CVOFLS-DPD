// solve Level-Set function as the re-initialization equation
   Info<< "solve the reinitialization equation"     
       << nl << endl;
   mixture.interfaceProperties::correctPsi0();
   Psii == mixture.interfaceProperties::psi0();

   for (int corr=0; corr<int(mixture.interfaceProperties::epsilon().value()/mixture.interfaceProperties::deltaTau().value()); corr++)
   {
      Psii = Psii + mixture.interfaceProperties::psi0()/(mag(mixture.interfaceProperties::psi0())+SMALL)*(double(1)-mag(fvc::grad(Psii))*dimMeter)*mixture.interfaceProperties::deltaTau();
      Psii.correctBoundaryConditions();
   }

   mixture.interfaceProperties::correctDeltaHscale();
