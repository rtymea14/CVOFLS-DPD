Info<< "Reading/creating Liquid transportProperties...\n" << endl;


//- Allocating memory for total liquid concentration [kmol/m3]
Info<< "Creating field cTot1\n" << endl;
volScalarField cTot1
(
    IOobject
    (
        "liquid_cTot1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cTot1", dimensionSet(0, -3, 0, 0, 1), 1e-10)    
);

// Create psi1 = rho1/p field [s2/m2]
Info<< "Creating field psi1...\n" << endl;
volScalarField psi1
(
    IOobject
    (
        "psi1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
        mesh,       
        dimensionedScalar("psi1", dimensionSet(0, -2, 2, 0, 0), 0.0)
);
psi1.oldTime();

//- Creating field rho1 [kg/m3]
Info<< "Creating field rho1\n" << endl;
volScalarField rho1
(
   IOobject
   (
	"rho1",
         runTime.timeName(),
         mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
   ),
    mesh,
   dimensionedScalar("rho1", dimensionSet(1, -3, 0, 0, 0), 1000.)
);

rho1.oldTime();


//- Allocating memory for liquid dynamic viscosity [kg/m/s]
Info<< "Creating field mu liquid\n" << endl;
volScalarField mu1
(
    IOobject
    (
        "liquid_mu1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputLiquidProperties
    ),
    mesh,  
    dimensionedScalar("mu1", dimensionSet(1, -1, -1, 0, 0), 1e-4)
);

//- Allocating memory for liquid cinematic viscosity [kg/m/s]
Info<< "Creating field nu liquid\n" << endl;
volScalarField nu1
(
    IOobject
    (
        "liquid_nu1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputLiquidProperties
    ),
    mesh,
    dimensionedScalar("nu1", dimensionSet(0, 2, -1, 0, 0), 1.0e-6)
);

//- Allocating memory for liquid thermal conductivity [W/m/K]
Info<< "Creating field lambda liquid\n" << endl;
volScalarField lambda1
(
    IOobject
    (
        "liquid_lambda1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputLiquidProperties
    ),
    mesh,  
    dimensionedScalar("lambda1", dimensionSet(1, 1, -3, -1, 0), 1.0)
);



//- Allocating memory for liquid mean specific heat [J/kg/K]
Info<< "Creating field cp1\n" << endl;
volScalarField cp1
(
    IOobject
    (
        "liquid_cp1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputLiquidProperties
    ),
    mesh,
    dimensionedScalar("cp1", dimensionSet(0, 2, -2, -1, 0), 4219.)
);



//- Allocating memory for total evaporation power [J/m3]
Info<< "Creating field vapHeat\n" << endl;
volScalarField vapHeat
(
    IOobject
    (
        "liquid_vapHeat",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("vapHeat", dimensionSet(1, -1, -3, 0, 0), 1e-10)     
);



PtrList<volScalarField> Cp1Species(thermodynamicsMapXML->NumberOfSpecies());

for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	{

		Cp1Species.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo_Cp1Species_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                                ),
                                mesh,
			        dimensionedScalar("Cp1Species", dimensionSet(0, 2, -2, -1, 0), 4219.0)
                        )
                );

	}
}


PtrList<volScalarField> Dmix1(thermodynamicsMapXML->NumberOfSpecies());

for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	{
		Dmix1.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"liquid_Dmix1_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("Dmix1", dimensionSet(0, 2, -1, 0, 0), 1.0e-8)
			)
		);

	}
}



//- Create equilibrium constants for evaporation
PtrList<volScalarField> Keq(NLS);

Info << "   Creating species L-V equilibrium: " <<  endl;
for(int i=0;i<NLS;i++)
{


	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::READ_IF_PRESENT
	);

	{
		Keq.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"TD_Keq_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("Keq", dimensionSet(0, 0, 0, 0, 0), 0.0)
			)
		);
         }

}

// Heat of vaporatization for all liquid species
PtrList<volScalarField> dhev(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
		dhev.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"dhev_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
                        dimensionedScalar("dhev_", dimensionSet(0, 2, -2, 0, 0), 2260.0e3)
			)
		);
         }
  }


//- Create activity coefficients
PtrList<volScalarField> gamma(NLS);

  for(int i=0;i<NLS;i++)
  {
	{
		gamma.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"gamma_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
                        dimensionSet(0, 0, 0, 0, 0)
			)
		);
         }

    gamma[i] *= 0.0;
  }


