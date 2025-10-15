const volScalarField& psi = thermo.psi();
const label inertIndex(composition.species()[inertSpecie]);
const label oxidiserIndex(composition.species()[oxSpecie]);
const label gprodIndex(composition.species()[gprodSpecie]);
const label gmetalIndex(composition.species()[gmetalSpecie]);
