{
    volScalarField& heF = thermo.he();
    volScalarField& heS = thermoS.he();

    fvScalarMatrix EFEqn
    (
        a * fvm::ddt(rho, heF) 
        + fvm::div(phi, heF)
        - a * fvm::laplacian(thermo.alpha(), heF)
        - HTcoef * (thermoS.T() - thermo.T())

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
        fvOptions(rho, heF)
    );

    EFEqn.relax();
	fvOptions.constrain(EFEqn);
    EFEqn.solve();
    fvOptions.correct(heF);
	thermo.correct();
	
	volScalarField& Tc = thermoS.T();
	volScalarField rhoCp(thermoS.rho()*thermoS.Cp());

    fvScalarMatrix ESEqn
    (
      /*  (1-a) * fvm::ddt(rhoCp, Tc) 
		- (1-a) * fvm::Sp(fvc::ddt(rhoCp), Tc)*/
		(1-a) * rhoCp * fvm::ddt(Tc)
        - (1-a) * fvm::laplacian(thermoS.kappa(), Tc)
        + fvm::Sp(HTcoef, Tc) 
		- HTcoef * thermo.T()
     ==
        fvOptions(rhoCp, Tc)
    );

    ESEqn.relax();
    fvOptions.constrain(ESEqn);
    ESEqn.solve();
    fvOptions.correct(Tc);
    heS = thermoS.he(p, Tc);
//    Info << thermoS.updateT() << endl;
    thermoS.correct();
//    thermo.init(p, Tc, heF);
}
