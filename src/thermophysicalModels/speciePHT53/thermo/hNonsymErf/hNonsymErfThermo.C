/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "hNonsymErfThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hNonsymErfThermo<EquationOfState>::hNonsymErfThermo(const dictionary& dict)
:
    EquationOfState(dict),
    CpS_(dict.subDict("thermodynamics").get<scalar>("CpS")),
    CpL_(dict.subDict("thermodynamics").get<scalar>("CpL")),
    Hf_(dict.subDict("thermodynamics").get<scalar>("Hf")),
    Tph_(dict.subDict("thermodynamics").get<scalar>("Tph")),
    T2minusT1_(dict.subDict("thermodynamics").get<scalar>("T2minusT1")),
    addCoef_(dict.subDict("thermodynamics").get<scalar>("addCoef"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::hNonsymErfThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("thermodynamics");
        os.writeEntry("CpS", CpS_);
        os.writeEntry("CpL", CpL_);
        os.writeEntry("Hf", Hf_);
        os.writeEntry("Tph", Tph_);
        os.writeEntry("T2minusT1", T2minusT1_);
        os.writeEntry("addCoef", addCoef_);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hNonsymErfThermo<EquationOfState>& ct
)
{
    ct.write(os);
    return os;
}


// ************************************************************************* //
