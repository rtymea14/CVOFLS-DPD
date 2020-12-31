

Info<< " * Reading chemical species mole fractions or mass fractions..." << endl;

PtrList<volScalarField> Y(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> X(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> Dmix2(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> FugacityCoefficientMix(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> CpSpecies2(thermodynamicsMapXML->NumberOfSpecies());
//PtrList<volScalarField> FormationRates(outputFormationRatesIndices.size());


#if STEADYSTATE == 1
PtrList<volScalarField> sourceImplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
PtrList<volScalarField> sourceExplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
#endif


PtrList<volScalarField> Dsoret(transportMapXML->iThermalDiffusionRatios().size());
if (soretEffect == true)
{
	for(unsigned int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
	{ 
		unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
		indexSoret[index-1] = i;
		soretEffectList[index-1] = true;

		std::cout << "Soret effect: species " 	<< i+1 << "/" << transportMapXML->iThermalDiffusionRatios().size() 
							<< " " << thermodynamicsMapXML->NamesOfSpecies()[index-1] << " (" << index << ")" << std::endl;

		Dsoret.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"gas_Dsoret_"+ thermodynamicsMapXML->NamesOfSpecies()[index-1],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("Dsoret", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);
	}
}

// Loop over allthe species in the kinetic mechanism
for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{
	Info << "   Reading species: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " > ";

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// Sets of diffusivities and specific heats of species
	{
		Dmix2.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"gas_Dmix2_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("Dmix2", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);

		FugacityCoefficientMix.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"gas_FugacityCoefficientMix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("FugacityCoefficientMix", dimensionSet(0, 0, 0, 0, 0), 1.)
			)
		);


		CpSpecies2.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo_Cp2_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::NO_WRITE
                                ),
                                mesh,
		             	dimensionedScalar("Cp2", dimensionSet(0, 2, -2, -1, 0), 1e-10)
                        )
                );

		X.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "X_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::AUTO_WRITE
                                ),
                                mesh,
                                dimensionSet(0, 0, 0, 0, 0)
                        )
                );

	}

	// check if field exists and can be read


	//if (header.headerOk())
	if (header.headerOk())
	{
		Info << thermodynamicsMapXML->NamesOfSpecies()[i] << endl;
		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE //AUTO
				),
				mesh
			)
		);
	}

	else
	{
		Info << "Ydefault" << endl;

		volScalarField Ydefault
		(
			IOobject
			(
				"Ydefault",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh
		);

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE //auto
				),
				Ydefault
			)
		);
	}
}




// Calculate inert molecular mass

{
  MWinert *= 0.0;
  int counter = 0;

          forAll(T.internalField(), celli)
          {
             if (alpha1.internalField()[celli] == 0.)
             {
               counter++;

               double sumYiDividedMWi = 0.;
               std::vector<double> moleFractions(thermodynamicsMapXML->NumberOfSpecies());

               for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                 sumYiDividedMWi += Y[i].internalField()[celli]/thermodynamicsMapXML->MW(i);                

               for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                 moleFractions[i] = (Y[i].internalField()[celli]/thermodynamicsMapXML->MW(i))/sumYiDividedMWi;

               for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
                 MWinert.value() += moleFractions[i]*thermodynamicsMapXML->MW(i);  

               break;              
             }
             
          } 

       reduce(counter, sumOp<int>());
       reduce(MWinert.value(), sumOp<double>());

       MWinert.value() /= counter;

  Info<< "MWinert = " << MWinert.value() << endl;

}












