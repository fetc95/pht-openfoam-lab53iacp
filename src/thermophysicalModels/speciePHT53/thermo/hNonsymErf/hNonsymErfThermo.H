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

Class
    Foam::hNonsymErfThermo

Group
    grpSpecieThermo

Description
    Package of the thermodynamic properties
    templated into the EquationOfState.
    The enthalpy changes as error function of the 
    piece-wise constant temperature function
    to model non-symmetric solid-liquid phase changes.

SourceFiles
    hNonsymErfThermoI.H
    hNonsymErfThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hNonsymErfThermo_H
#define hNonsymErfThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hNonsymErfThermo;

template<class EquationOfState>
inline hNonsymErfThermo<EquationOfState> operator+
(
    const hNonsymErfThermo<EquationOfState>&,
    const hNonsymErfThermo<EquationOfState>&
);

template<class EquationOfState>
inline hNonsymErfThermo<EquationOfState> operator*
(
    const scalar,
    const hNonsymErfThermo<EquationOfState>&
);

template<class EquationOfState>
inline hNonsymErfThermo<EquationOfState> operator==
(
    const hNonsymErfThermo<EquationOfState>&,
    const hNonsymErfThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hNonsymErfThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hNonsymErfThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hNonsymErfThermo
:
    public EquationOfState
{
    // Private data

        scalar CpS_;
        scalar CpL_;
        scalar Hf_;
        scalar Tph_;
        scalar T2minusT1_;
        scalar addCoef_;


    // Private Member Functions

        //- Construct from components
        inline hNonsymErfThermo
        (
            const EquationOfState& st,
            const scalar cps,
            const scalar cpl,
            const scalar hf,
            const scalar tph,
            const scalar t2minust1,
            const scalar addcoef
        );


public:

    // Constructors

        //- Construct from dictionary
        hNonsymErfThermo(const dictionary& dict);

        //- Construct as named copy
        inline hNonsymErfThermo(const word&, const hNonsymErfThermo&);

        //- Construct and return a clone
        inline autoPtr<hNonsymErfThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hNonsymErfThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hNonsymErf<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- liquid fraction
            inline scalar fL(const scalar T) const;

            //- liquid fraction temperature partial derivative
            inline scalar dfLdT(const scalar T) const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            #include "HtoEthermo.H"


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;



        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hNonsymErfThermo&);


    // Friend operators

        friend hNonsymErfThermo operator+ <EquationOfState>
        (
            const hNonsymErfThermo&,
            const hNonsymErfThermo&
        );

        friend hNonsymErfThermo operator* <EquationOfState>
        (
            const scalar,
            const hNonsymErfThermo&
        );

        friend hNonsymErfThermo operator== <EquationOfState>
        (
            const hNonsymErfThermo&,
            const hNonsymErfThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hNonsymErfThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hNonsymErfThermoI.H"

#ifdef NoRepository
    #include "hNonsymErfThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
