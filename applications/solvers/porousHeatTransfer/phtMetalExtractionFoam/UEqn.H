// Solve the Momentum equation

//MRF.correctBoundaryVelocity(U);
/*
dimensionedVector g(dimLength/dimTime/dimTime,vector(0,-9.81,0));
KdField = 1. / (1./KD + (1./KdEnd-1./KD)*eta);*/
tmp<fvVectorMatrix> tUEqn
(
      
    fvm::ddt(a_g, rho, U) 
    + fvm::div(phi, U)
    + a_g * rho * chi * (1 - a_g) * (
    fvc::ddt(U) + (U & fvc::grad(U))
	)
//    + fvm::Sp(a_g*a_g*KdField * thermo.mu() + rho * F * mag(U), U)
	 + fvm::Sp
	 (
	     /*a_g*a_g*/(
		     (150*(1-a_g)*(1-a_g)/(prtclDmtr*prtclDmtr*a_g/*a_g*a_g*/)) * thermo.mu() 
		     + rho * (1.75*(1-a_g)/(prtclDmtr/*a_g*a_g*a_g*/)) * mag(/*a_g*/U)
			 ), 
		U
	 )
//  + MRF.DDt(rho, U)
//  - fvm::laplacian(mu, U)
 ==
    a_g*rho*gu
    + fvOptions(a_g, rho, U)
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -a_g*fvc::grad(p));
//    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
//    K = 0.5*magSqr(U);
}
