{
    HTcoef = HTcoef0 * pow((1 - a_g) / (1-a_g0),HTorder);
	
	volScalarField& he = thermo.he();
	volScalarField Qgdot(Q_evap*Rcnd);

    fvScalarMatrix EEqn
    (
//        a_g*rho*fvm::ddt(he)
//		+a_g*rho*(U & fvc::grad(he))
        fvm::ddt(a_g, rho, he) 
        + fvm::div(phi, he)
		-fvm::Sp(RR, he)
//		-(c_cf - mu_p*c_cp)*rho_cf0*W_rate*Tc
        - fvm::laplacian(a_g * thermo.alpha(), he)
//        + fvm::Sp(HTcoef / thermo.Cp(), he)
//        - HTcoef / thermo.Cp() * he
/*      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )*/
     ==
        HTcoef * (Tc - thermo.T())
		+ Qgdot
		+ fvOptions(a_g, rho, he)
    );

    EEqn.relax();
    fvOptions.constrain(EEqn);
    EEqn.solve();
    fvOptions.correct(he);
    thermo.correct();
	
	// Solid phase
	
	volScalarField solidRhoEffCp("solidRhoEffCp",fuelRhoEff*thermoSFuel.Cp() + productRhoEff*thermoSProduct.Cp() + metalRhoEff*thermoSMetal.Cp() + inertRhoEff*thermoSInert.Cp());
	solidKappaEff = fuelVol*thermoSFuel.kappa() + productVol*thermoSProduct.kappa() + metalVol*thermoSMetal.kappa() + inertVol*thermoSInert.kappa();
	volScalarField Qsdot(Q_comb * Rcomb - Q_evap * Revp);
	ddTlFraction = 4*exp(-pow(4*(Tc-T_melt)/deltaMelt,2)) / (Foam::sqrt(Foam::constant::mathematical::pi)*deltaMelt);

    fvScalarMatrix TcEqn
    (
        solidRhoEffCp*fvm::ddt(Tc)
		- fvm::laplacian(solidKappaEff, Tc)
        + fvm::Sp(HTcoef, Tc)
        - HTcoef * thermo.T()
		+ metalRhoEff * Q_melt * ddTlFraction * fvm::ddt(Tc)
        ==
		Qsdot
        + fvOptions(solidRhoEffCp, Tc)
    );

    TcEqn.relax();
    fvOptions.constrain(TcEqn);
    TcEqn.solve();
    fvOptions.correct(Tc);
	
	thermoSFuel.correct();
	thermoSProduct.correct();
	thermoSMetal.correct();
	thermoSInert.correct();
	
	lFraction = 0.5*(Foam::erf(4*(Tc-T_melt)/deltaMelt) + 1);

//    if(runTime.writeTime()){ solidRhoEffCp.write(); solidKappaEff.write(); 
//	volScalarField("inertCp",thermoSInert.Cp()).write(); volScalarField("inertKappa",thermoSInert.kappa()).write();}
}
