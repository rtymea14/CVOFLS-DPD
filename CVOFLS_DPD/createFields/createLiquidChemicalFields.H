
// Create gradient fields of liquid species for Stephan Maxwell diffusivities calculation and evaporating mass flux update
PtrList<volScalarField> Ysat(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
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


		Ysat.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)] + "sat",
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


// Create gradient fields of liquid species for Stephan Maxwell diffusivities calculation and evaporating mass flux update
PtrList<volScalarField> Xsat(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
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


		Xsat.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)] + "Xsat",
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

// Evaporating fluxes for all liquid species
PtrList<volScalarField> mevap(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
		mevap.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"mevap_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh,
                        dimensionSet(1, -3, -1, 0, 0)
			)
		);
         }

      mevap[i] *= 0.0;
  }


// Poynting
PtrList<volScalarField> Poynting(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
		Poynting.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"Poynting_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh,
                        dimensionSet(0, 0, 0, 0, 0)
			)
		);
         }

      Poynting[i] *= 0.0;
  }


// Fugacity coefficients for pure species
PtrList<volScalarField> FugacityCoefficientPure(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
		FugacityCoefficientPure.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"FugacityCoefficientPure_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh,
                        dimensionSet(0, 0, 0, 0, 0)
			)
		);
         }

      FugacityCoefficientPure[i] *= 0.0;
  }



// Create objects for evaluation of liquid properties
//speciesMap* species_map;
//mixtureL*   liquid_mixtureCells;
//mixtureG*   gas_mixtureCells;

std::vector<std::string> listLiquidSpecies_;
std::vector<std::string> listGasSpecies_;

listLiquidSpecies_.resize(listLiquidSpecies.size());
listGasSpecies_.resize(thermodynamicsMapXML->NumberOfSpecies());

                          for (unsigned int j=0;j<listLiquidSpecies.size();j++)
                          {                          
                            listLiquidSpecies_[j] = listLiquidSpecies[j];
                          }

                          for (unsigned int j=0;j<listGasSpecies_.size();j++)
                          {                          
                            listGasSpecies_[j] = thermodynamicsMapXML->NamesOfSpecies()[j];
                          }


//species_map = new speciesMap(liquid_folder);
//
//liquid_mixtureCells = new mixtureL(listLiquidSpecies_, *species_map); 
//gas_mixtureCells = new mixtureG(listGasSpecies_, *species_map); 



   // CREATE THE TWO-PHASE MIXTURE OBJECT
    immiscibleIncompressibleTwoPhaseMixture mixture(U, rho1, rho2, nu1, nu2, cp1, cp2, lambda1, lambda2, pc, Psii);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());
