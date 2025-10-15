// Reading thermophyscial properties of all solid components and phase change properties of the metal

Info<< "Reading thermophysical properties of the solid\n" << endl;

IOdictionary thermoSolidProperties
(
    IOobject
    (
        IOobject::groupName("thermophysicalProperties",solidPhaseName),
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//wordList solidSpecies(thermoSolidProperties.get<wordList>("species"));

const word inertSolidSpecieName("inert");
const word fuelSolidSpecieName("fuel");
const word productSolidSpecieName("product");
const word metalSolidSpecieName("metal");

/*const word inertSolidThermoDictName(IOobject::groupName(inertSolidSpecieName+"ThermophysicalProperties",solidPhaseName));
const word fuelSolidThermoDictName(IOobject::groupName(fuelSolidSpecieName+"ThermophysicalProperties",solidPhaseName));
const word productSolidThermoDictName(IOobject::groupName(productSolidSpecieName+"ThermophysicalProperties",solidPhaseName));
const word metalSolidThermoDictName(IOobject::groupName(metalSolidSpecieName+"ThermophysicalProperties",solidPhaseName));*/

dictionary inertSolidThermoDict(thermoSolidProperties.subDict(inertSolidSpecieName));
dictionary fuelSolidThermoDict(thermoSolidProperties.subDict(fuelSolidSpecieName));
dictionary productSolidThermoDict(thermoSolidProperties.subDict(productSolidSpecieName));
dictionary metalSolidThermoDict(thermoSolidProperties.subDict(metalSolidSpecieName));

volScalarField Tc
(
    IOobject
    (
        IOobject::groupName("T",solidPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*
volScalarField TcFictive
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	Tc
);*/

Info<< "Solid fuel\n" << endl;
//Create thermo instance of the component
autoPtr<solidThermo> sfThermo
(
//    solidThermo::New(mesh,solidPhaseName,fuelSolidThermoDictName)
	solidThermo::New(mesh,fuelSolidThermoDict,solidPhaseName)
);
solidThermo& thermoSFuel = sfThermo();

// Read effective rho of the specie if set and calculate its vol fraction.
// If not set, look up the vol fraction of the specie and calculate its effective rho.

volScalarField fuelRhoEff
(
    IOobject
    (
        IOobject::groupName(fuelSolidSpecieName+"RhoBulk", solidPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

if (!fuelRhoEff.headerOk())
{
    volScalarField fuelVol
	(
		IOobject
		(
			IOobject::groupName(fuelSolidSpecieName+"Vol", solidPhaseName),
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);
    fuelRhoEff = fuelVol * thermoSFuel.rho();
}

volScalarField fuelVol(IOobject::groupName(fuelSolidSpecieName+"Vol", solidPhaseName), fuelRhoEff/thermoSFuel.rho());

Info<< "Solid reaction products\n" << endl;
autoPtr<solidThermo> spThermo
(
//    solidThermo::New(mesh,solidPhaseName,productSolidThermoDictName)
    solidThermo::New(mesh,productSolidThermoDict,solidPhaseName)
);
solidThermo& thermoSProduct = spThermo();

volScalarField productRhoEff
(
    IOobject
    (
        IOobject::groupName(productSolidSpecieName+"RhoBulk", solidPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

if (!productRhoEff.headerOk())
{
    volScalarField productVol
	(
		IOobject
		(
			IOobject::groupName(productSolidSpecieName+"Vol", solidPhaseName),
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);
    productRhoEff = productVol * thermoSProduct.rho();
}

volScalarField productVol(IOobject::groupName(productSolidSpecieName+"Vol", solidPhaseName), productRhoEff/thermoSProduct.rho());

Info<< "Solid metal\n" << endl;
autoPtr<solidThermo> smThermo
(
//    solidThermo::New(mesh,solidPhaseName,metalSolidThermoDictName)
    solidThermo::New(mesh,metalSolidThermoDict,solidPhaseName)
);
solidThermo& thermoSMetal = smThermo();

volScalarField metalRhoEff
(
    IOobject
    (
        IOobject::groupName(metalSolidSpecieName+"RhoBulk", solidPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume, Zero),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

if (!metalRhoEff.headerOk())
{
    volScalarField metalVol
	(
		IOobject
		(
			IOobject::groupName(metalSolidSpecieName+"Vol", solidPhaseName),
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		),
		mesh
	);
    metalRhoEff = metalVol * thermoSMetal.rho();
}

volScalarField metalVol(IOobject::groupName(metalSolidSpecieName+"Vol", solidPhaseName), metalRhoEff/thermoSMetal.rho());
//metalRhoEff.write();

dimensionedScalar Q_melt("Hf",dimEnergy/dimMass,thermoSMetal.subDict("mixture").subDict("thermodynamics"));
dimensionedScalar T_melt("Tf",dimTemperature,thermoSMetal.subDict("mixture").subDict("thermodynamics"));

dimensionedScalar deltaMelt
(
    dimTemperature,
    thermoSMetal.subDict("mixture").subDict("thermodynamics").getOrDefault<scalar>("range",1)
);

volScalarField lFraction
(
    IOobject
    (
        "liquidMetalFraction",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    scalar(0)*a_g
);

volScalarField ddTlFraction("ddTLiquidFraction", scalar(0)*a_g/deltaMelt);

Info<< "Solid inert specie\n" << endl;
autoPtr<solidThermo> siThermo
(
//    solidThermo::New(mesh,solidPhaseName,inertSolidThermoDictName)
    solidThermo::New(mesh,inertSolidThermoDict,solidPhaseName)
);
solidThermo& thermoSInert = siThermo();
//thermoSInert.validate(args.executable(), "h", "e");

volScalarField inertVol
(
	IOobject::groupName(inertSolidSpecieName+"Vol", solidPhaseName), 
	scalar(1)-a_g-fuelVol-productVol-metalVol
);

volScalarField inertRhoEff
(
    IOobject
    (
        IOobject::groupName(inertSolidSpecieName+"RhoBulk", solidPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    inertVol * thermoSInert.rho(),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

volScalarField solidKappaEff
(
    IOobject
    (
        IOobject::groupName("kappaBulk", solidPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fuelVol*thermoSFuel.kappa() + productVol*thermoSProduct.kappa() + metalVol*thermoSMetal.kappa() + inertVol*thermoSInert.kappa(),
	extrapolatedCalculatedFvPatchScalarField::typeName
);

