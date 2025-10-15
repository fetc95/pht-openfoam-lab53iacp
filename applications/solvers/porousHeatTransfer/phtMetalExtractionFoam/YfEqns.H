{

//    D_g = D_g0 * pow(thermo.T()/T_norm, b);
    
	volScalarField factorReact(Foam::pos0(Tc-TCombInit) * fuelRhoEff * Acomb * exp(-EdivR/Tc));
	
	volScalarField& Oxdsr = Y[oxidiserIndex];
	
    fvScalarMatrix OxdsrEqn
    (
        fvm::ddt(a_g, rho, Oxdsr) 
        + fvm::div(phi, Oxdsr)
//        - fvm::laplacian(a_g * rho * D_g, Oxdsr)
        + fvm::Sp(xi_g * factorReact, Oxdsr)
     ==
        fvOptions(a_g, rho, Oxdsr)
    );

    OxdsrEqn.relax();
    fvOptions.constrain(OxdsrEqn);
    OxdsrEqn.solve(mesh.solver("Yi.fluid"));
	Oxdsr.max(0.0);
	
	volScalarField& gProd = Y[gprodIndex];
	
    fvScalarMatrix gProdEqn
    (
        fvm::ddt(a_g, rho, gProd) 
        + fvm::div(phi, gProd)
     ==
        xi_gp*Oxdsr*factorReact
        + fvOptions(a_g, rho, gProd)
    );

    gProdEqn.relax();
    fvOptions.constrain(gProdEqn);
    gProdEqn.solve(mesh.solver("Yi.fluid"));
	gProd.max(0.0);
	
	volScalarField factorCnd(Aevp*a_g*rho*(T_evap2-thermo.T())*Foam::pos0(T_evap2-thermo.T()));
	
	volScalarField& MetalFluid = Y[gmetalIndex];

	fvScalarMatrix MetalFluidEqn
	(
	    fvm::ddt(a_g, rho, MetalFluid) 
        + fvm::div(phi, MetalFluid)
		+ fvm::Sp(factorCnd, MetalFluid)
		==
		Revp
		+ fvOptions(a_g, rho, MetalFluid)
    );
	
	MetalFluidEqn.relax();
	fvOptions.constrain(MetalFluidEqn);
	MetalFluidEqn.solve(mesh.solver("Yi.fluid"));
	MetalFluid.max(0.0);
	
	Y[inertIndex] = scalar(1) - Y[oxidiserIndex] - Y[gprodIndex] - Y[gmetalIndex];
    Y[inertIndex].max(0.0);
}
