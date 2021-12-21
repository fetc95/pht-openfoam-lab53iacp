/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd
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
    Foam::hePhaseChangeSolidThermo

Description
    Energy for a solid mixture with liquid fraction

SourceFiles
    hePhaseChangeSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePhaseChangeSolidThermo_H
#define hePhaseChangeSolidThermo_H

#include "heThermo.H"
//#include "fvMesh.H"
#include "volFields.H"
#include "autoPtr.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class hePhaseChangeSolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicSolidThermo, class MixtureType>
class hePhaseChangeSolidThermo
:
    public heThermo<BasicSolidThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        hePhaseChangeSolidThermo(const hePhaseChangeSolidThermo<BasicSolidThermo, MixtureType>&);


protected:

    //- liquid fraction []
    volScalarField fL_;

//    bool fLOwner_;


public:

    //- Runtime type information
    TypeName("hePhaseChangeSolidThermo");


    // Constructors

        //- Construct from mesh and phase name
        hePhaseChangeSolidThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        hePhaseChangeSolidThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

         //- Construct from mesh, dictionary name and phase name
        hePhaseChangeSolidThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~hePhaseChangeSolidThermo();


    // Member functions

        //- Update properties
        virtual void correct();


        // Derived thermal properties

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const;

            //- Return true if thermal conductivity is isotropic
            virtual bool isotropic() const
            {
                return MixtureType::thermoType::isotropic;
            }


        // Per patch calculation

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<vectorField> Kappa(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hePhaseChangeSolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
