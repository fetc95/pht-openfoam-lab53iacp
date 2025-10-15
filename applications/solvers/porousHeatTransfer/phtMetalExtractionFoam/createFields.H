#include "createRDeltaT.H"

word fluidPhaseName("fluid");
word solidPhaseName("solid");

#include "createPhases/createFluidPhase.H"
#include "createPhases/createPorousMedium.H"
#include "createPhases/createSolidPhase.H"
#include "createPhases/readReactions.H"

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

#include "compressibleCreatePhi.H"

#include "readGravitationalAcceleration.H"
dimensionedVector gu(g.dimensions(),g.value());

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

//#include "createDpdt.H"

//#include "createK.H"

//#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"
