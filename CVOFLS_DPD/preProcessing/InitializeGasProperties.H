{

  if(speciesEquations == true)
  {


	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info<< "Gas properties evaluation... " ;

        OpenSMOKE::OpenSMOKEVectorDouble CpVector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble Dmixvector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble tetamixvector(thermodynamicsMapXML->NumberOfSpecies());
	Eigen::VectorXd massFractionsEigen(thermodynamicsMapXML->NumberOfSpecies());
        volScalarField Xt = 0.0*X[0];

	// Internal fields
	
		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();

		scalarField& cTotCells = cTot2.internalField();
		scalarField& psiCells =  psi2.internalField();
		scalarField& muCells  =  mu2.internalField();
		scalarField& nuCells  =  nu2.internalField();
		scalarField& lambdaCells = lambda2.internalField();
		scalarField& cpCells = cp2.internalField();
		scalarField& MWmixCells = MWmix.internalField();


		forAll(TCells, celli)
		{

			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			
			transportMapXML->SetPressure(pCells[celli]);
			transportMapXML->SetTemperature(TCells[celli]);
				
	
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].internalField()[celli];

				
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),MWmixCells[celli],massFractions.GetHandle());

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                        {
			   X[i].internalField()[celli] = moleFractions[i+1];
                        }


			cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			psiCells[celli]  = cTotCells[celli]*MWmixCells[celli]/pCells[celli];


			//if (viscositySpecies == true)
			{
				 muCells[celli] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());
                                 nuCells[celli] = muCells[celli]/(pCells[celli]*psiCells[celli]);
			}




            		if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            		{
				lambdaCells[celli] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
               			cpCells[celli] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/kmol/K]
				cpCells[celli] = cpCells[celli]/MWmixCells[celli];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector.GetHandle());

					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                               CpSpecies2[i].internalField()[celli] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
	
				}
			}
			

			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(), moleFractions.GetHandle());

				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix2[i].internalField()[celli] = Dmixvector[i+1];

			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient = lambdaCells[celli]/(pCells[celli]*psiCells[celli])/cpCells[celli];
				
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix2[i].internalField()[celli] = coefficient/LewisNumbers(i);
			}
/*

			if (physicalSootDiffusivity == true)
			{
				const double DmixReference = Dmix2[physicalSootDiffusivityReferenceIndex].internalField()[celli];

				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix2[physicalSootDiffusivityCorrectionIndex[i]].internalField()[celli] = DmixReference*physicalSootDiffusivityCorrection[i];
			}






			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());		
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];

					Dsoret[i].internalField()[celli] =Dmix2[index-1].internalField()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/MWmixCells[celli];
				} 
			}
*/
		
    }

	// Boundaries
	forAll(T.boundaryField(), patchi)
	{

		fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];

		fvPatchScalarField& pcTot = cTot2.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsi  = psi2.boundaryFieldRef()[patchi];
		fvPatchScalarField& pmu   = mu2.boundaryFieldRef()[patchi];
		fvPatchScalarField& pnu   = nu2.boundaryFieldRef()[patchi];
		fvPatchScalarField& plambda = lambda2.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcp = cp2.boundaryFieldRef()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryFieldRef()[patchi];
        
		forAll(pT, facei)
		{

			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);

			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),pMWmix[facei],massFractions.GetHandle());

			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				X[i].boundaryFieldRef()[patchi][facei] = moleFractions[i+1];


			pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
			ppsi[facei]  = pcTot[facei]*pMWmix[facei]/pp[facei];


			

			//if (viscositySpecies == true)
			{
			        pmu[facei] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());
                                pnu[facei] = pmu[facei]/(pp[facei]*ppsi[facei]);
			}


						
           		if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            		{
				plambda[facei] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
				
				pcp[facei] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
				pcp[facei] = pcp[facei]/pMWmix[facei];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->cpMolar_Species(CpVector.GetHandle());

					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                                                CpSpecies2[i].boundaryFieldRef()[patchi][facei] = CpVector[i+1] / thermodynamicsMapXML->MW(i);
				}
			} 
			



	if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(),moleFractions.GetHandle());

				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix2[i].boundaryFieldRef()[patchi][facei] = Dmixvector[i+1];

			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient =  plambda[facei]/(pp[facei]*ppsi[facei])/pcp[facei];

				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix2[i].boundaryFieldRef()[patchi][facei] = coefficient/LewisNumbers(i);

			}
/*
			if (physicalSootDiffusivity == true)
			{
				const double DmixReference = Dmix2[physicalSootDiffusivityReferenceIndex].boundaryField()[patchi][facei];

				for(int i=0;i<physicalSootDiffusivityCorrectionIndex.size();i++)
					Dmix2[physicalSootDiffusivityCorrectionIndex[i]].boundaryFieldRef()[patchi][facei] = DmixReference*physicalSootDiffusivityCorrection[i];

			}



			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];

					Dsoret[i].boundaryFieldRef()[patchi][facei] = Dmix2[index-1].boundaryField()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/pMWmix[facei];

				}
			}
*/
                }

	}

	 double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	 Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

    }

    else
    {
       Info << "Species equations must be activated for the evaluation of gas and liquid properties!" << endl;
       abort();
    }


}




// Update gas density
rho2 = Rho(p,psi2);




