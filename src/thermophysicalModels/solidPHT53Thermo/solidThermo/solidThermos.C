/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "makeSolidThermo.H"
#include "solidThermo.H"
#include "heSolidThermo.H"
#include "hePhaseChangeSolidThermo.H"

#include "specie.H"
#include "rhoConst.H"
#include "icoPolynomial.H"
//#include "hConstThermo.H"
//#include "hPowerThermo.H"
//#include "hPolynomialThermo.H"
//#include "hNonsymErfThermo.H"
#include "hMeltingSolidificationCubicThermo.H"
#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"
//#include "exponentialSolidTransport.H"
//#include "polynomialSolidTransport.H"
#include "pureMixture.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
//#include "liquidFractionEnthalpy.H"
#include "thermo.H"

//#include "solidThermoPhysicsTypes.H"
#include "solidFetc95ThermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

/*makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hNonsymErfThermo,
    rhoConst,
    specie
);*/

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hMeltingSolidificationCubicThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    hePhaseChangeSolidThermo,
    pureMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hMeltingSolidificationCubicThermo,
    rhoConst,
    specie
);
/*
makeSolidThermo
(
    solidThermo,
    hePhaseChangeSolidThermo,
    pureMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hNonsymErfThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    hePhaseChangeSolidThermo,
    pureMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hMeltingSolidificationCubicThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    constAnIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    exponentialSolidTransport,
    sensibleEnthalpy,
    hPowerThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    polynomialSolidTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);


makeSolidThermoPhysicsType
(
    solidThermo,
    heSolidThermo,
    pureMixture,
    hTransportThermoPoly8SolidThermoPhysics
);*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
