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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hNonsymErfThermo<EquationOfState>::hNonsymErfThermo
(
    const EquationOfState& st,
    const scalar cps,
    const scalar cpl,
    const scalar hf,
    const scalar tph,
    const scalar t2minust1,
    const scalar addcoef
)
:
    EquationOfState(st),
    CpS_(cps),
    CpL_(cpl),
    Hf_(hf),
    Tph_(tph),
    T2minusT1_(t2minust1),
    addCoef_(addcoef)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hNonsymErfThermo<EquationOfState>::hNonsymErfThermo
(
    const word& name,
    const hNonsymErfThermo& ct
)
:
    EquationOfState(name, ct),
    CpS_(ct.CpS_),
    CpL_(ct.CpL_),
    Hf_(ct.Hf_),
    Tph_(ct.Tph_),
    T2minusT1_(ct.T2minusT1_),
    addCoef_(ct.addCoef_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hNonsymErfThermo<EquationOfState>>
Foam::hNonsymErfThermo<EquationOfState>::clone() const
{
    return autoPtr<hNonsymErfThermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hNonsymErfThermo<EquationOfState>>
Foam::hNonsymErfThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hNonsymErfThermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return CpS_ + (CpL_-CpS_)*fL(T) + Hf_*dfLdT(T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return CpS_*T + (CpL_-CpS_)*fL(T)*(T-Tph_) + Hf_*fL(T) + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::Hc() const
{
    return Hf_;
}

template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::fL
(
    const scalar T
) const
{
    const scalar bb1 = addCoef_ / (addCoef_ + 1);
    const scalar nonsymX = 1 + (addCoef_-1)*Foam::pos0(T-Tph_);
    return bb1*Foam::erf(4*nonsymX*(T-Tph_)/T2minusT1_)/nonsymX + bb1;
}

template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::dfLdT
(
    const scalar T
) const
{
    const scalar bb1 = addCoef_ / (addCoef_ + 1);
    const scalar nonsymX = 1 + (addCoef_-1)*Foam::pos0(T-Tph_);
    return 8*bb1*exp(-pow(4*nonsymX*(T-Tph_)/T2minusT1_,2))/(Foam::sqrt(Foam::constant::mathematical::pi)*T2minusT1_);
}

template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return 0;//Cp_*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hNonsymErfThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hNonsymErfThermo<EquationOfState>::operator+=
(
    const hNonsymErfThermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();

        CpS_ = Y1*CpS_ + Y2*ct.CpS_;
        CpL_ = Y1*CpL_ + Y2*ct.CpL_;
        Hf_ = Y1*Hf_ + Y2*ct.Hf_;
        Tph_ = Y1*Tph_ + Y2*ct.Tph_;
        T2minusT1_ = Y1*T2minusT1_ + Y2*ct.T2minusT1_;
        addCoef_ = Y1*addCoef_ + Y2*ct.addCoef_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hNonsymErfThermo<EquationOfState> Foam::operator+
(
    const hNonsymErfThermo<EquationOfState>& ct1,
    const hNonsymErfThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return hNonsymErfThermo<EquationOfState>
        (
            eofs,
            ct1.CpS_,
            ct1.CpL_,
            ct1.Hf_,
            ct1.Tph_,
            ct1.T2minusT1_,
            ct1.addCoef_
        );
    }
    else
    {
        return hNonsymErfThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.CpS_ + ct2.Y()/eofs.Y()*ct2.CpS_,
            ct1.Y()/eofs.Y()*ct1.CpL_ + ct2.Y()/eofs.Y()*ct2.CpL_,
            ct1.Y()/eofs.Y()*ct1.Hf_ + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Y()/eofs.Y()*ct1.Tph_ + ct2.Y()/eofs.Y()*ct2.Tph_,
            ct1.Y()/eofs.Y()*ct1.T2minusT1_ + ct2.Y()/eofs.Y()*ct2.T2minusT1_,
            ct1.Y()/eofs.Y()*ct1.addCoef_ + ct2.Y()/eofs.Y()*ct2.addCoef_
        );
    }
}


template<class EquationOfState>
inline Foam::hNonsymErfThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hNonsymErfThermo<EquationOfState>& ct
)
{
    return hNonsymErfThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.CpS_,
        ct.CpL_,
        ct.Hf_,
        ct.Tph_,
        ct.T2minusT1_,
        ct.addCoef_
    );
}


template<class EquationOfState>
inline Foam::hNonsymErfThermo<EquationOfState> Foam::operator==
(
    const hNonsymErfThermo<EquationOfState>& ct1,
    const hNonsymErfThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hNonsymErfThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.CpS_ - ct1.Y()/eofs.Y()*ct1.CpS_,
        ct2.Y()/eofs.Y()*ct2.CpL_ - ct1.Y()/eofs.Y()*ct1.CpL_,
        ct2.Y()/eofs.Y()*ct2.Hf_ - ct1.Y()/eofs.Y()*ct1.Hf_,
        ct2.Y()/eofs.Y()*ct2.Tph_ - ct1.Y()/eofs.Y()*ct1.Tph_,
        ct2.Y()/eofs.Y()*ct2.T2minusT1_ - ct1.Y()/eofs.Y()*ct1.T2minusT1_,
        ct2.Y()/eofs.Y()*ct2.addCoef_ - ct1.Y()/eofs.Y()*ct1.addCoef_
    );
}


// ************************************************************************* //
