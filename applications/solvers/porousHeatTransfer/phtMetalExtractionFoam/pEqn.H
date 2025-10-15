if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
//surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(a_g*rho*rAU));
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(a_g*a_g*rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(a_g*rho)*fvc::flux(HbyA)
//	fvc::interpolate(rho)*fvc::flux(HbyA)
//  + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
//    + rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf)
);

//fvc::makeRelative(phiHbyA, rho, U);
//MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
//constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);
constrainPressure(p, rho, a_g*U, phiHbyA, rhorAUf);
/*
if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA) + fvm::div(phid, p)
     ==
        fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{*/
    fvScalarMatrix pDDtEqn
    (
//        ( a_g*fvc::ddt(rho) + rho*fvc::ddt(a_g)+ a_g * psi*correction(fvm::ddt(p)) )
//        fvc::ddt(a_g, rho) + a_g*psi*correction(fvm::ddt(p))
		fvm::ddt(psi,p)
      + fvc::div(phiHbyA) 
     ==
	    RR
      + fvOptions(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p));

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
		//	U = fvc::reconstruct(phi)/a_g/rho;
        }
    }
//}


// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - a_g*rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
//K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}

thermo.correctRho(psi*p - psip0, rhoMin, rhoMax) ;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

rho = thermo.rho();

// Correct rhoUf if the mesh is moving
//fvc::correctRhoUf(rhoUf, rho, U, phi);
/*
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}*/
