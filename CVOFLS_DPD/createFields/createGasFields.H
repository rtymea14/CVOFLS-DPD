Info<< "Reading/creating Gas transportProperties...\n" << endl;


//- Allocating memory for total gas concentration [kmol/m3]
Info<< "Creating field cTot2\n" << endl;
volScalarField cTot2
(
    IOobject
    (
        "gas_cTot2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cTot2", dimensionSet(0, -3, 0, 0, 1), 1e-10)
    
);


// Create psi2 = rho2/p field [s2/m2]
Info<< "Creating field psi2...\n" << endl;
volScalarField psi2
(
    IOobject
    (
        "psi2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
        mesh,       
        dimensionedScalar("psi2", dimensionSet(0, -2, 2, 0, 0), 1e-5)
);
psi2.oldTime();

//- Creating field rho2 [kg/m3]
Info<< "Creating field rho2\n" << endl;
volScalarField rho2
(
   IOobject
   (
	"rho2",
         runTime.timeName(),
         mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
   ),
  //   Rho(p,psi2)
        mesh,       
        dimensionedScalar("rho2", dimDensity, 1.)
);

rho2.oldTime();


//- Allocating memory for gas dynamic viscosity [kg/m/s]
Info<< "Creating field mu gas\n" << endl;
volScalarField mu2
(
    IOobject
    (
        "gas_mu2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputGasProperties
    ),
    mesh,
    dimensionedScalar("mu2", dimensionSet(1, -1, -1, 0, 0), 1e-5)
);


//- Allocating memory for gas cinematic viscosity [kg/m/s]
Info<< "Creating field nu gas\n" << endl;
volScalarField nu2
(
    IOobject
    (
        "gas_nu2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputGasProperties
    ),
    mesh, 
    dimensionedScalar("nu2", dimensionSet(0, 2, -1, 0, 0), 1.48e-5)
);

//- Allocating memory for gas thermal conductivity [W/m/K]
Info<< "Creating field lambda gas\n" << endl;
volScalarField lambda2
(
    IOobject
    (
        "gas_lambda2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputGasProperties
    ),
    mesh,   
    dimensionedScalar("lambda2", dimensionSet(1, 1, -3, -1, 0), 0.025)
);


//- Allocating memory for gas mean specific heat [J/kg/K]
Info<< "Creating field Cp gas\n" << endl;
volScalarField cp2
(
    IOobject
    (
        "gas_Cp2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputGasProperties
    ),
    mesh,   
    dimensionedScalar("cp2", dimensionSet(0, 2, -2, -1, 0), 1000.)
);

//- Allocating memory for Z1
Info<< "Creating field Zgas\n" << endl;
volScalarField Zgas
(
    IOobject
    (
        "gas_Zgas",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        outputGasProperties
    ),
    mesh,   
    dimensionedScalar("Zgas", dimless, 1.)
);



// Creation of insoluble gas species (O2, CO2, CO ecc...)

// Create gradient fields of liquid species for Stephan Maxwell diffusivities calculation and evaporating mass flux update

PtrList<volScalarField> Ygas(NGS);

  for(int i=0;i<NGS;i++)
  {

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[GasSpeciesIndices(i)],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);


	if (header.headerOk())
	{

		volScalarField Ypresent
		(
			IOobject
			(
				thermodynamicsMapXML->NamesOfSpecies()[GasSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			mesh
		);


		Ygas.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[GasSpeciesIndices(i)] + "gas",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE //AUTO
				),
				Ypresent
			)
		);
	}

        else
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


		Ygas.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[GasSpeciesIndices(i)] + "gas",
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

