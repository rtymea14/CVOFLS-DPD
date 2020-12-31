
/*
// Read capillary pressure  [Pa]
Info<< "Reading field pc\n" << endl;
volScalarField pc
(
    IOobject
    (
        "pc",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
pc.oldTime();


label pcRefCell = 0;
scalar pcRefValue = 0.0;
setRefCell
(
    pc,
    pimple.dict(),
    pcRefCell,
    pcRefValue
);



// Create smoothed alpha field [-]
Info<< "Create alpha1smooth...\n" << endl;
  volScalarField alpha1smooth
    (
        IOobject
        (
            "alpha1smooth",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedScalar("alpha1smooth", dimless, 0.0)
    );
    
alpha1smooth.oldTime();


// Create curtailed alpha field [-]
Info<< "Create alpha1pc...\n" << endl;
  volScalarField alpha1pc
    (
        IOobject
        (
            "alpha1pc",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedScalar("alpha1pc", dimless, 0.0)
    );
    
alpha1pc.oldTime();


// Create curvature field [1/m]
Info<< "Create curvature...\n" << endl;
  volScalarField Kvol
    (
        IOobject
        (
            "Kvol",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedScalar("Kvol", dimless/dimLength, 0.0)
    );
    
Kvol.oldTime();


// Create filtered curvature field [1/m]
Info<< "Create filtered curvature...\n" << endl;
  volScalarField K_filtered
    (
        IOobject
        (
            "K_filtered",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedScalar("K_filtered", dimless/dimLength, 0.0)
    );
    
K_filtered.oldTime();


// Create surface tension force field [N/m^3]
Info<< "Create surface tension force...\n" << endl;
  volVectorField TensionForce
    (
        IOobject
        (
            "TensionForce",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedVector("TensionForce", dimForce/dimVol, Foam::vector(0,0,0))
    );
    
TensionForce.oldTime();

// Create evaporation Cells Marker
Info<< "filteredTensionForce..\n" << endl;
volVectorField FilteredTensionForce
(
    IOobject
    (
        "FilteredTensionForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
        dimensionedVector("FilteredTensionForce", dimForce/dimVol, Foam::vector(0,0,0))
);
FilteredTensionForce.oldTime();

// Create capillary pressure gradient field 
Info<< "Create gradalpha force...\n" << endl;
  volVectorField gradalpha
    (
        IOobject
        (
            "gradalpha",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
         dimensionedVector("gradalpha", dimless/dimLength, Foam::vector(0,0,0))
    );
    
gradalpha.oldTime();


// Create filtered surface tension force
Info<< "Create filtered surface tension force...\n" << endl;
  surfaceScalarField fcf
    (
        IOobject
        (
            "fcf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
         dimensionedScalar("fcf", dimForce/dimVol, 0.0)
    );
    
fcf.oldTime();


// CurvatureLSM
Info<< "Reading field CurvatureLSM\n" << endl;
volScalarField CurvatureLSM
(
    IOobject
    (
        "CurvatureLSM",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh, 
    dimensionedScalar("CurvatureLSM", dimless/dimLength  , 0.0)//,
    //zeroGradientFvPatchScalarField::typeName
);


// Create surfaceTension
Info<< "Create surfaceTension...\n" << endl;
  volScalarField surfaceTension
    (
        IOobject
        (
            "surfaceTension",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedScalar("surfaceTension", dimForce/dimLength, 0.073)
    );
surfaceTension.oldTime();


// Create surfaceTension
Info<< "Create alphasmooth_aux...\n" << endl;
  volScalarField alphasmooth_aux
    (
        IOobject
        (
            "alphasmooth_aux",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
       ),
       mesh,
        dimensionedScalar("alphasmooth_aux", dimless, 0.0)
    );
alphasmooth_aux.oldTime();




// CurvatureLSM
Info<< "Reading field CurvatureLSM\n" << endl;
volVectorField normal
(
    IOobject
    (
        "normal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("normal", dimless, Foam::vector(0,0,0))//,
    //zeroGradientFvPatchVectorField::typeName
);


// CurvatureLSM
Info<< "Reading field tangent\n" << endl;
volVectorField tangent
(
    IOobject
    (
        "tangent",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("tangent", dimless, Foam::vector(0,0,0))

);


// Create surfaceTension
Info<< "CurvCoord...\n" << endl;
  volScalarField CurvCoord
    (
        IOobject
        (
            "CurvCoord",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
       ),
    mesh
    );
CurvCoord.oldTime();
*/


dimensionedScalar Curvature("Curvature", dimless/dimLength, 0.0);
dimensionedScalar surfaceTension("surfaceTension", dimForce/dimLength, 0.0);





