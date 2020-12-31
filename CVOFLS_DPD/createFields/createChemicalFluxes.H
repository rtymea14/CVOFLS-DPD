
PtrList<surfaceScalarField> J(NLS);

for (label i=0; i<NLS; i++)
{
	J.set
	(
		i,
		new surfaceScalarField
		(
			IOobject
			(
				"J_"+ thermodynamicsMapXML->NamesOfSpecies()[LiquidSpeciesIndices(i)],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionSet(1, 0, -1, 0, 0)
		)
	);
}

surfaceScalarField Jinert
(
    IOobject
    (
	"Jinert",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Jinert", dimensionSet(1, 0, -1, 0, 0), 0.0)
);


surfaceScalarField Jc
(
    IOobject
    (
	"fluxes_Jc",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("J_Jc", dimensionSet(1, 0, -1, 0, 0), 0.0)
);
