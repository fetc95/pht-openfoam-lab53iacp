Info<< "Reading reaction properties\n" << endl;

IOdictionary reactionProperties
(
    IOobject
    (
        "reactionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//------------------------------------- fuel combustion -------------------------------------------------------

dictionary combDict(reactionProperties.subDict("heterogeneousCombustion"));
scalarList combustionCoeffs(combDict.get<scalarList>("coeffs"));

//Info<< combustionCoeffs << endl;
//Info<< combustionCoeffs.size() << endl;

if (combustionCoeffs.size() != 5)
{
	FatalIOErrorIn(args.executable().c_str(), reactionProperties)
        << "Entry combustionCoeffs must contain 5 values\n"
        << combDict << exit(FatalIOError);
}

dimensionedScalar Q_comb(dimEnergy/dimMass,combustionCoeffs[0]);
dimensionedScalar Acomb(dimless/dimTime,combustionCoeffs[1]);
dimensionedScalar EdivR(dimTemperature,combustionCoeffs[2]);
dimensionedScalar xi_g(dimless,combustionCoeffs[3]);
dimensionedScalar xi_p(dimless,combustionCoeffs[4]);
dimensionedScalar xi_gp(scalar(1)+xi_g-xi_p);

const word oxSpecie(combDict.get<word>("oxidiser"));
if (!composition.species().found(oxSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), combDict)
        << "Oxidiser " << oxSpecie << " not found in available gas species "
        << composition.species() << exit(FatalIOError);
}

const word gprodSpecie(combDict.get<word>("gasProduct"));
if (!composition.species().found(gprodSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), combDict)
        << "Gaseous reaction product " << gprodSpecie << " not found in available gas species "
        << composition.species() << exit(FatalIOError);
}

dimensionedScalar TCombInit(dimTemperature, combDict.getOrDefault<scalar>("Tcrit", 500));

volScalarField Rcomb
(
    IOobject
    (
        "Rcomb",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermoSFuel.rho()*Acomb*scalar(0)
);

//------------------------------------- metal evaporation/condensation -------------------------------------------------------

dictionary evpDict(reactionProperties.subDict("evaporationCondensation"));
scalarList evpcndCoeffs(evpDict.get<scalarList>("coeffs"));

if (evpcndCoeffs.size() != 2)
{
	FatalIOErrorIn(args.executable().c_str(), reactionProperties)
        << "Entry evpcndCoeffs must contain 2 values: Heat and temperature of evaporation/condensation\n"
        << evpcndCoeffs << exit(FatalIOError);
}

dimensionedScalar Q_evap(dimEnergy/dimMass,evpcndCoeffs[0]);
dimensionedScalar T_evap(dimTemperature,evpcndCoeffs[1]);
dimensionedScalar Aevp(dimless/dimTemperature/dimTime,  evpDict.getOrDefault<scalar>("Aevp",1));
dimensionedScalar deltaEvap(dimTemperature, evpDict.getOrDefault<scalar>("range",1));
dimensionedScalar T_evap1(T_evap-deltaEvap/2.);
dimensionedScalar T_evap2(T_evap+deltaEvap/2.);

const word gmetalSpecie(evpDict.get<word>("gasMetal"));
if (!composition.species().found(gmetalSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), evpDict)
        << "Metal specie " << gmetalSpecie << " not found in available gas species "
        << composition.species() << exit(FatalIOError);
}

volScalarField Revp
(
    IOobject
    (
        "Revp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermoSMetal.rho()*Aevp*T_evap1*scalar(0)
);

volScalarField Rcnd
(
    IOobject
    (
        "Rcnd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*Aevp*T_evap2*scalar(0)
);

// Total mass source
volScalarField RR
(
    IOobject
    (
        "RR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (1-xi_p) * Rcomb + Revp - Rcnd
);

//Info<< reactionProperties << endl;
// Out reaction Info

Info<< "Oxidasion reaction equation:\n"
    << "\tSolid fuel + " << xi_g.value() << " " << oxSpecie << " = "
	<< xi_gp.value() << " " << gprodSpecie << " + " << xi_p.value() << " Solid products"
	<< endl;

Info<< "\tHeat of reaction = " << Q_comb.value() << " J/kg\n"
	<< "\tReaction rate factor = " << Acomb.value() << " 1/s\n"
	<< "\tTa = " << EdivR.value() << " K\n"
	<< endl;
	
Info<< "Evaporation/condensation:\n"
	<< "\tHeat of evaporation = " << Q_evap.value() << " J/kg\n"
	<< "\tTemperature of evaporation = " << T_evap.value() << " K\n"
	<< endl;
	