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
    Foam::hMeltingSolidificationCubicThermo

Group
    grpSpecieThermo

Description
    Package of the thermodynamic properties
    templated into the EquationOfState.
    The enthalpy changes as error function of the 
    piece-wise constant temperature function
    to model non-symmetric solid-liquid phase changes.

SourceFiles
    hMeltingSolidificationCubicThermoI.H
    hMeltingSolidificationCubicThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hMeltingSolidificationCubicThermo_H
#define hMeltingSolidificationCubicThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hMeltingSolidificationCubicThermo;

template<class EquationOfState>
inline hMeltingSolidificationCubicThermo<EquationOfState> operator+
(
    const hMeltingSolidificationCubicThermo<EquationOfState>&,
    const hMeltingSolidificationCubicThermo<EquationOfState>&
);

template<class EquationOfState>
inline hMeltingSolidificationCubicThermo<EquationOfState> operator*
(
    const scalar,
    const hMeltingSolidificationCubicThermo<EquationOfState>&
);

template<class EquationOfState>
inline hMeltingSolidificationCubicThermo<EquationOfState> operator==
(
    const hMeltingSolidificationCubicThermo<EquationOfState>&,
    const hMeltingSolidificationCubicThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hMeltingSolidificationCubicThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hMeltingSolidificationCubicThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hMeltingSolidificationCubicThermo
:
    public EquationOfState
{
    // Private data

        scalar CpS_;
        scalar CpL_;
        scalar Hf_;
        scalar Tph_;
        scalar T1_;
        scalar T2_;


    // Private Member Functions

        //- Construct from components
        inline hMeltingSolidificationCubicThermo
        (
            const EquationOfState& st,
            const scalar cps,
            const scalar cpl,
            const scalar hf,
            const scalar tph,
            const scalar t1,
            const scalar t2
        );


public:

    // Constructors

        //- Construct from dictionary
        hMeltingSolidificationCubicThermo(const dictionary& dict);

        //- Construct as named copy
        inline hMeltingSolidificationCubicThermo(const word&, const hMeltingSolidificationCubicThermo&);

        //- Construct and return a clone
        inline autoPtr<hMeltingSolidificationCubicThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hMeltingSolidificationCubicThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hMeltingSolidificationCubic<" + EquationOfState::typeName() + '>';
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
            inline scalar Ff(const scalar T) const;

            //- liquid fraction temperature partial derivative
            inline scalar dFf(const scalar T) const;

            //- liquid fraction integral
            inline scalar IFf(const scalar T) const;
            
            inline scalar ffa(const scalar T) const;

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

        inline void operator+=(const hMeltingSolidificationCubicThermo&);


    // Friend operators

        friend hMeltingSolidificationCubicThermo operator+ <EquationOfState>
        (
            const hMeltingSolidificationCubicThermo&,
            const hMeltingSolidificationCubicThermo&
        );

        friend hMeltingSolidificationCubicThermo operator* <EquationOfState>
        (
            const scalar,
            const hMeltingSolidificationCubicThermo&
        );

        friend hMeltingSolidificationCubicThermo operator== <EquationOfState>
        (
            const hMeltingSolidificationCubicThermo&,
            const hMeltingSolidificationCubicThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hMeltingSolidificationCubicThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hMeltingSolidificationCubicThermoI.H"

#ifdef NoRepository
    #include "hMeltingSolidificationCubicThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
