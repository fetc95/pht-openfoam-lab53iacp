Info<< "Reading thermophysical properties of the fluid\n" << endl;
//fluidThermo psiReactionThermo
autoPtr<psiReactionThermo> pThermo
(
    psiReactionThermo::New(mesh,fluidPhaseName)
);
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

if (composition.species().size() != 4)
{
	FatalIOErrorIn(args.executable().c_str(), thermo)
        << "4 gas species must set "
        << composition.species() << exit(FatalIOError);
}

const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert gas specie " << inertSpecie << " not found in available gas species "
        << composition.species() << exit(FatalIOError);
}

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, i)
{
    fields.add(Y[i]);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);
