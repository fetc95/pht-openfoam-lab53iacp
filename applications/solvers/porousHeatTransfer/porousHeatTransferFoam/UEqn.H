// Solve the Momentum equation

//MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    addedMassDivA2 * a * fvm::ddt(rho, U) 
    + addedMassDivA2 * fvm::div(phi, U)
    + fvm::Sp(KD * thermo.mu() + rho * F * mag(U), U)
//  + MRF.DDt(rho, U)
//  - fvm::laplacian(mu, U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
//    K = 0.5*magSqr(U);
}
