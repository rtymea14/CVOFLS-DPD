
{

    if(energyEquation == true)
    {
		//radiation->correct();

		massDiffusionInEnergyEquation *= 0.;
		if (iMassDiffusionInEnergyEquation == true) 
		{
			for (label i=0; i<NLS; i++)
			       massDiffusionInEnergyEquation -= CpSpecies2[LiquidSpeciesIndices(i)]*( fvc::reconstruct(J[i]) & fvc::grad(T));

			       massDiffusionInEnergyEquation -= CpSpecies2[inertIndex]*( fvc::reconstruct(Jinert) & fvc::grad(T));
		}

                vapHeat *= 0.0;

                for (unsigned int i=0;i<NLS;i++)
                {                           
                    vapHeat += mevap[i]*dhev[i];                                        
                }



		
		fvScalarMatrix TEqn
		(
			  cp*fvm::ddt(rho,T)
			+ cp*fvm::div(rhoPhi, T)		    
			==
	    		fvm::laplacian(lambda,T) + 
	   	 	massDiffusionInEnergyEquation +
			DpDt +
                        vapHeat +
		//	radiation->divq(T) +
			//RT +
	    		cp*fvOptions(rho, T)
		);

		TEqn.relax();
		fvOptions.constrain(TEqn);
		TEqn.solve(mesh.solver("T"));
		fvOptions.correct(T);
		



// Adjust temperature field: the liquid temperature cannot overcome the boiling temperature



         forAll(alpha1.internalField(), celli)
         {
            if (alpha1.internalField()[celli] > 1e-1 && T.internalField()[celli] > 0.999*Tboiling)
                T.internalField()[celli] = 0.999*Tboiling;
         }


         T.correctBoundaryConditions();


   }


    else
    {
       if(setTemperatureField == true)
       {
          T= alpha1*Tdroplet + alpha2*Tambient;

          T.correctBoundaryConditions();
       }
    }





    Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;



    // Spark
//    spark.SetSpark(runTime.value(), T);
}


