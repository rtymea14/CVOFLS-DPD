// Interpolation scheme
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
	fields.add(Y[i]);
}

// Create total mEvapTot
Info<< "Creating field mEvapTot...\n" << endl;
volScalarField mEvapTot
(
    IOobject
    (
        "mEvapTot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
	mesh,
        dimensionedScalar("mEvapTot",  dimDensity/dimTime, 0.0)
);
mEvapTot.oldTime();


// Create massLiquidCorrection
Info<< "Creating massLiquidCorrection...\n" << endl;
volScalarField massLiquidCorrection
(
    IOobject
    (
        "massLiquidCorrection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
        dimensionedScalar("massLiquidCorrection",  dimDensity/dimTime, 0.0)
);
massLiquidCorrection.oldTime();

// Create dgdt for comprimibility effects on alpha
volScalarField dgdt
(
   
    IOobject
    (
        "dgdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
     pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001))
);
dgdt.oldTime();


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture.lookup("pMin")
);



// Spark
//spark.Setup(mesh);


dimensionedScalar  dimChange
(
  dimensionedScalar("dimChange",dimLength , 1.0)
);


//Info<< "Creating DeltaT field" << endl;
//volScalarField DeltaT
//(
//	IOobject
//	(
//		"DeltaT",
//		runTime.timeName(),
//		mesh,
//		IOobject::NO_READ,
//		IOobject::AUTO_WRITE
//	),
//	mesh,
//	dimensionedScalar("one", dimTime, 0),
//
//	extrapolatedCalculatedFvPatchScalarField::typeName
//
//);


// Magnetic fields

// Create magnetic potential  
Info<< "Reading field magnetic potential  \n" << endl;
volScalarField MagneticPotential  
(
    IOobject
    (
        "MagneticPotential",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("MagneticPotential", dimLength*dimLength/dimTime/dimTime, 0.0)
);
MagneticPotential.oldTime();


// Create Marker 
Info<< "Reading field Marker  \n" << endl;
volScalarField Marker
(
    IOobject
    (
        "Marker",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Marker", dimless, 0.0)
);
Marker.oldTime();


dimensionedScalar  Mass0(  dimensionedScalar("Mass0", dimMass, 0.) );

  Mass0 = fvc::domainIntegrate(rho1*alpha1);


  dimensionedScalar averageLiquidRho = Mass0/fvc::domainIntegrate(alpha1);

// Create YdefaultBackup
volScalarField YdefaultBackup
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





/*
// Create IDs
Info<< "IDs  \n" << endl;
volScalarField IDs
(
    IOobject
    (
        "IDs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //mesh,
   // dimensionedScalar("IDs", dimless, 0.0)
);
IDs.oldTime();


forAll(IDs.internalField(),celli)
{
   IDs[celli] = celli;
}



// Create Neighbours
Info<< "Neighbours  \n" << endl;
volTensorField Neighbours
(
    IOobject
    (
        "Neighbours",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //mesh,
    //dimensionSet(0, 0, 0, 0, 0) 
);
Neighbours.oldTime();



    List<scalarField> IDsList;
    List<scalarField> alphaList;
    List<List<scalarField>> rhoDmixiList;
    List<tensorField> NeighboursList;

    parallelClass initializeParallel(mesh);

      if (Pstream::parRun())
      {
	IDsList.resize(Pstream::nProcs());
	IDsList[Pstream::myProcNo()] = IDs.internalField();    
 
	alphaList.resize(Pstream::nProcs());
	alphaList[Pstream::myProcNo()] = alpha1.internalField();

	rhoDmixiList.resize(NLS);

        for (unsigned int i=0;i<NLS;i++)
        rhoDmixiList[i].resize(Pstream::nProcs());

           for (unsigned int i=0;i<NLS;i++)
           {
             rhoDmixiList[i][Pstream::myProcNo()] = rhoDmix[LiquidSpeciesIndices(i)].internalField();
             initializeParallel.exchange(rhoDmixiList[i]);
           }


	NeighboursList.resize(Pstream::nProcs());
	NeighboursList[Pstream::myProcNo()] = Neighbours.internalField(); 

        initializeParallel.exchange(IDsList);
        initializeParallel.exchange(alphaList);
        initializeParallel.exchange(NeighboursList);
      }
*/

// Constants for evaporating flux evaluation
dimensionedScalar totIntegral("totIntegral", dimMass/dimTime, 1e-10);
dimensionedScalar totIntegralAlpha("totIntegralAlpha", dimMass/dimTime, 1e-10);
dimensionedScalar EvapCorrection("EvapCorrection", dimless, 1e-10);
dimensionedScalar CompressibleCorrection("CompressibleCorrection", dimless, 1e-10);
dimensionedScalar FixedFlux("FixedFlux", dimMass/dimTime, 1e-8);



// Compute the boiling temperature at the given pressure

    int NTemperaturePoints = 50000;
    std::vector<double> TemperatureResiduals (NTemperaturePoints);

    double Tmin = 300.;
    double Tmax = 700.;

    double intervalT = (Tmax-Tmin)/((double)NTemperaturePoints-1.);
    double pVap = 0.0;
    double Tazz = Tmin;
    double Tboiling = Tmin;


     for (unsigned int j=0;j< NTemperaturePoints;j++)
     { 
        Tazz = Tmin+(double)j*intervalT;

        pVap = Liquid->pv(thermodynamicPressure, Tazz); //species_map->pVap(listLiquidSpecies_[0], Tazz, thermodynamicPressure);

        if( Foam::mag(pVap/thermodynamicPressure - 1.0) < 1e-3 )
        {
           Tboiling =  Tazz;
           break;
        }

     }


    Info<< "Boiling temperature at " << thermodynamicPressure << " Pa = " << Tboiling << " K "<< endl;



  dimensionedScalar Tdroplet("Tdroplet", dimTemperature, 0.0);
  dimensionedScalar Tambient("Tambient", dimTemperature, 0.0);

  Tdroplet.value() = readScalar(physicalModelDictionary.lookup("Tdroplet"));
  Tambient.value() = readScalar(physicalModelDictionary.lookup("Tambient"));

/*
// Set temperature field equal to the initial one, but accounting for alpha field 

     if(energyEquation == false && setTemperatureField == true)
     {
         double onlyGas = 0.;
         double onlyLiquid = 0.;

         forAll(alpha1.internalField(), celli)
         {
            if (alpha1.internalField()[celli] == 1.0 )
              onlyLiquid++;

            else
              onlyGas++;             
         }


         if (onlyLiquid == 0.0)
            Tambient.value() = T.internalField()[0];

         if (onlyGas == 0.0)
            Tdroplet.value() = T.internalField()[0];

         else
         {

           forAll(alpha1.internalField(), celli)
           {
             if (alpha1.internalField()[celli] == 1.0 )
             {
               Tdroplet.value() = T.internalField()[celli];
               break;
             }
             
           }


           forAll(alpha1.internalField(), celli)
           {
             if (alpha1.internalField()[celli] == 0.0 )
             {
               Tambient.value() = T.internalField()[celli];
               break;
             }
             
           }

         }
     }
*/



