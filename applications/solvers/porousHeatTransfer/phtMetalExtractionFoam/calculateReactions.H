{
	// calculating fuel and product solid species
	
	volScalarField factorReact = Foam::pos0(Tc-TCombInit) * Y[oxidiserIndex]*Acomb*exp(-EdivR/Tc);//Y[oxidiserIndex]//Oxdsr

	Info<< "EulerImplicit: Solving for solid species effective densities" << endl;
/*    fvScalarMatrix fuelEqn
    (
        fvm::ddt(fuelRhoEff)
      + fvm::Sp(factorReact, fuelRhoEff)
    );

    fuelEqn.solve(mesh.solver("rhoEffi.solid"));*/
	fuelRhoEff = fuelRhoEff.oldTime() / (1 + runTime.deltaT()*factorReact);
	Rcomb = fuelRhoEff * factorReact;
	
/*	fvScalarMatrix productEqn
    (
        fvm::ddt(productRhoEff)
      - xi_p * Rcomb
    );
	
	productEqn.solve(mesh.solver("rhoEffi.solid"));*/
	productRhoEff = productRhoEff.oldTime() + runTime.deltaT()*xi_p*Rcomb;
	fuelVol = fuelRhoEff / thermoSFuel.rho();
	productVol = productRhoEff / thermoSProduct.rho();
	
	//calculating metal solid specie
	
	Rcnd = Aevp * a_g*rho*Y[gmetalIndex] * (T_evap2-thermo.T())*Foam::pos0(T_evap2-thermo.T());//Y[gmetalIndex]//MetalFluid
	volScalarField factorEvp = Aevp*(Tc-T_evap1)*Foam::pos0(Tc-T_evap1);

/*    fvScalarMatrix metalEqn
    (
        fvm::ddt(metalRhoEff)
      + fvm::Sp(factorEvp, metalRhoEff)
	  - Rcnd
    );

    metalEqn.solve(mesh.solver("rhoEffi.solid"));*/
	metalRhoEff = (metalRhoEff.oldTime() + runTime.deltaT()*Rcnd) / (1 + runTime.deltaT()*factorEvp);
	Revp = metalRhoEff*factorEvp;
	metalVol = metalRhoEff/thermoSMetal.rho();
	a_g = 1 - (inertVol + fuelVol + productVol + metalVol);
	
	//update overal mass source to fluid phase
	RR = (1 - xi_p) * Rcomb + Revp - Rcnd;
}
