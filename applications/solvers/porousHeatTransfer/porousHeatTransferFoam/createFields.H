#include "createRDeltaT.H"

word fluidPhaseName = "fluid";
word solidPhaseName = "solid";

Info<< "Reading fluid thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh,fluidPhaseName)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

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

Info<< "Reading solid thermophysical properties\n" << endl;

autoPtr<solidThermo> sThermo
(
    solidThermo::New(mesh,solidPhaseName)
);
solidThermo& thermoS = sThermo();
thermoS.validate(args.executable(), "h", "e");
//Info<<thermoS.updateT()<<endl;
//exit(1);
/*volScalarField Tc
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermoS.T()//,mesh
);*/

Info<< "Reading phase properties\n" << endl;
#include "readTransportProperties.H"

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

volScalarField kappaSolid
(
    IOobject
    (
        "kappa.solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermoS.kappa()
);

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

//#include "createDpdt.H"

//#include "createK.H"

//#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"
