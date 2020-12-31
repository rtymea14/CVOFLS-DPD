
sumXi *= 0.0;

    for (label i=0; i<NLS; i++)
    {
	volScalarField& Xsati = Xsat[i];
	volScalarField& Ysati = Ysat[i];
	volScalarField& rhoDmixi = rhoDmix[LiquidSpeciesIndices(i)];
	dimensionedScalar MWi("MWi", dimensionSet(1,0,0,0,-1,0,0),scalar(thermodynamicsMapXML->MW(LiquidSpeciesIndices(i))) ); 

        sumXi += Ysati;

	J[i] = ( fvc::interpolate ( -rhoDmixi*MWi/MWsat ) ) * ( fvc::interpolate ( fvc::grad(Xsati) ) ) & mesh.Sf() ;

    }

     Jinert = ( fvc::interpolate ( -rhoDmix[inertIndex]*MWinert/MWsat ) ) * ( fvc::interpolate ( fvc::grad(1.-sumXi) ) ) & mesh.Sf() ;
       


     // Correction flux

    Jc *= 0.;
    for (label i=0; i<NLS; i++)
    {
	Jc -= J[i];
    }

    Jc -=  Jinert;


     // Corrected fluxes

    for (label i=0; i<NLS; i++)
    {
	J[i] += Jc*fvc::interpolate(Ysat[i]);
    }

  Jinert += Jc*fvc::interpolate(1.-sumXi);
