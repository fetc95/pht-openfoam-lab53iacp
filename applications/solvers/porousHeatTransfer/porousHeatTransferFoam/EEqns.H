{
    volScalarField& heF = thermo.he();
    volScalarField& heS = thermoS.he();

    fvScalarMatrix EFEqn
    (
        a * fvm::ddt(rho, heF) 
        + fvm::div(phi, heF)
        - fvm::laplacian(a * thermo.alpha(), heF)
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

    fvScalarMatrix ESEqn
    (
        (1-a) * fvm::ddt(thermoS.rho(), heS) 
        - fvm::laplacian((1-a) * thermoS.alpha(), heS)
    //    - fvc::laplacian((1-a) * thermoS.kappa(), thermoS.T())
        + HTcoef * (thermoS.T() - thermo.T())
     ==
        fvOptions(thermoS.rho(), heS)
    );

    ESEqn.relax();

    fvOptions.constrain(EFEqn);
    EFEqn.solve();
    fvOptions.correct(heF);

    fvOptions.constrain(ESEqn);
    ESEqn.solve();
    fvOptions.correct(heS);

    thermo.correct();
    thermoS.correct();
}
