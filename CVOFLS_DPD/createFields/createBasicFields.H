Info<< "Reading field psi\n" << endl;
volScalarField Psii
(
    IOobject
    (
        "Psii",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read pressure [Pa]
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
p_rgh.oldTime();



// Read pressure [Pa]
Info<< "Reading field p_rghAdv\n" << endl;
volScalarField p_rghAdv
(
    IOobject
    (
        "p_rghAdv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
p_rghAdv.oldTime();



// Read velocity field [m/s]
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
U.oldTime();


// Read velocity field Uadv [m/s]
Info<< "Reading field UAdv\n" << endl;
volVectorField UAdv
(
    IOobject
    (
        "UAdv",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
UAdv.oldTime();

Info<< "\nReading momentum exchange field Ksl\n" << endl;
    volScalarField Ksl
    (
        IOobject
        (
            "Ksl",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    );

    Info<< "\nReading voidfraction field voidfraction = (Vgas/Vparticle)\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading particle velocity field Us\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    volVectorField gradPsi
    (   
        IOobject
        (
            "gradPsi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,-1,0,0,0), vector::zero)  // 1/s
    );
    
    //volVectorField gradAlpha
    //(   
    //    IOobject
    //    (
    //        "gradAlpha",
    //        mesh.time().timeName(),
    //        mesh,
    //        IOobject::NO_READ,
    //        IOobject::AUTO_WRITE
    //    ),
    //    mesh,
    //    dimensionedVector("zero", dimensionSet(0,-1,0,0,0), vector::zero)  // 1/s
    //);
    
    volVectorField ddtU_
    (   
        IOobject
        (
            "ddtU_",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-2,0,0), vector::zero)  // 1/s
    );
    
    volTensorField gradU_
    (   
        IOobject
        (
            "gradU_",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimensionSet(0,0,-1,0,0), tensor::zero)  // 1/s
    );

// Read temperature field [K]
Info<< "Reading field T...\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
T.oldTime();

//- Create Dp/Dt term [Pa]
volScalarField DpDt
(
    IOobject
    (
        "DpDt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DpDt", dimensionSet(1, -1, -3, 0, 0), 0.0)
);
DpDt.oldTime();

/*
// Create initial condition pressure field [Pa] (in order to initialize gas density field)
Info<< "Create initial pressure p_0...\n" << endl;
  volScalarField p_0
    (
        IOobject
        (
            "p_0",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
       ),
       mesh,
       dimensionedScalar("p_0", dimPressure, 101325.)
    );
    
p_0.oldTime();
*/


#include "createPhi.H"
#include "createPhiAdv.H"



//mesh.setFluxRequired(p_rgh.name());

//- Allocating memory for radiative heat transfer term [W/m3]
volScalarField Qrad
(
    IOobject
    (
        "rad_Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rad_Qrad", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

//- Allocating memory for dummy field useful only to tune the radiative heat transfer model
volScalarField dQ
(
	IOobject
	(
	    "dQ",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

//- Allocating memory for reaction heat [W/m3]
volScalarField Q
(
    IOobject
    (
        "gas_Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gas_Q", dimensionSet(1, -1, -3, 0, 0), 0.0)
); 

//- Allocating memory for mean molecular weight [kg/mol]
Info<< "Creating field MWmix\n" << endl;
volScalarField MWmix
(
    IOobject
    (
        "gas_MW",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 0, 0, 0, -1)
);



//- Allocating memory for mean molecular weight [kg/mol]
Info<< "Creating field MWsat\n" << endl;
volScalarField MWsat
(
    IOobject
    (
        "MWsat",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 0, 0, 0, -1)
);


// Create sumXiMWi
Info<< "Creating sumXiMWi...\n" << endl;
volScalarField sumXiMWi
(
    IOobject
    (
        "sumXiMWi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
        dimensionedScalar("sumXiMWi",  dimensionSet(1, 0, 0, 0, -1), 0.0)
);

// Create sumXi
Info<< "Creating sumXi...\n" << endl;
volScalarField sumXi
(
    IOobject
    (
        "sumXi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
        dimensionedScalar("sumXi",  dimless, 0.0)
);

// Create sumYi
Info<< "Creating sumYi...\n" << endl;
volScalarField sumYi
(
    IOobject
    (
        "sumYi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
        dimensionedScalar("sumYi",  dimless, 0.0)
);



dimensionedScalar MWinert("MWinert", dimensionSet(1,0,0,0,-1,0,0), 28.84000);
double FugacityCoefficientMixFixed = 1.0;
double FugacityCoefficientPureFixed = 1.0;
double Zfixed = 1.0;
double Poyntingfixed = 1.0;
double gammafixed = 1.0;

Info<< "Creating massDiffusionInEnergyEquation\n" << endl;
volScalarField massDiffusionInEnergyEquation
(
    IOobject
    (
		"gas_massDiffusionInEnergyEquation",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.)
);



#include "readGravitationalAcceleration.H"
//#include "readhRef.H"
//#include "gh.H"
Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());
surfaceScalarField ghf("ghf", g & mesh.Cf());

/*
// Create pressure [Pa]
Info<< "Create pressure...\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_0    
);
p.oldTime();


// Create pressure [Pa]
Info<< "Create pressure...\n" << endl;
volScalarField pAdv
(
    IOobject
    (
        "pAdv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_0    
);
pAdv.oldTime();
*/

// Read pressure [Pa]
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p_rgh
);
p.oldTime();

// Read pressure [Pa]
Info<< "Reading field pAdv\n" << endl;
volScalarField pAdv
(
    IOobject
    (
        "pAdv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    p
);
pAdv.oldTime();



// Read capillary pressure  [Pa]
Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
pc.oldTime();



// Create compressibilitySourceTerm
Info<< "Creating compressibilitySourceTerm...\n" << endl;
volScalarField compressibilitySourceTerm
(
    IOobject
    (
        "compressibilitySourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
        dimensionedScalar("compressibilitySourceTerm",  dimDensity/dimTime, 0.0)
);
compressibilitySourceTerm.oldTime();


// Constants for Diameter evaluation
dimensionedScalar  Diameter("Diameter", dimLength, 0.001);
dimensionedScalar  LiquidVolume("LiquidVolume", dimVolume, 0.0);
dimensionedScalar  FiberDiameter("FiberDiameter", dimLength, 1e-4);

dimensionedScalar dimMeter
(
    "dimMeter",
    dimLength,
    1.0
);

// create CFDEM Cloud
	cfdemCloud particleCloud(mesh);
	#include "checkModelType.H"
