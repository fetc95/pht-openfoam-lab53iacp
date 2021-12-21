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
inline Foam::hMeltingSolidificationCubicThermo<EquationOfState>::hMeltingSolidificationCubicThermo
(
    const EquationOfState& st,
    const scalar cps,
    const scalar cpl,
    const scalar hf,
    const scalar tph,
    const scalar t1,
    const scalar t2
)
:
    EquationOfState(st),
    CpS_(cps),
    CpL_(cpl),
    Hf_(hf),
    Tph_(tph),
    T1_(t1),
    T2_(t2)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hMeltingSolidificationCubicThermo<EquationOfState>::hMeltingSolidificationCubicThermo
(
    const word& name,
    const hMeltingSolidificationCubicThermo& ct
)
:
    EquationOfState(name, ct),
    CpS_(ct.CpS_),
    CpL_(ct.CpL_),
    Hf_(ct.Hf_),
    Tph_(ct.Tph_),
    T1_(ct.T1_),
    T2_(ct.T2_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hMeltingSolidificationCubicThermo<EquationOfState>>
Foam::hMeltingSolidificationCubicThermo<EquationOfState>::clone() const
{
    return autoPtr<hMeltingSolidificationCubicThermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::hMeltingSolidificationCubicThermo<EquationOfState>>
Foam::hMeltingSolidificationCubicThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<hMeltingSolidificationCubicThermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
	return CpS_ + (CpL_-CpS_)*Ff(T-Tph_) + Hf_*dFf(T-Tph_)/**/ + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return CpS_*T + (CpL_-CpS_)*IFf(T-Tph_) + Hf_*Ff(T-Tph_)/**/ + EquationOfState::H(p, T);
//	Info << T << "\t" << dFf(T-Tph_) << endl;
}


template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::Hc() const
{
    return Hf_;
}

template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::Ff
(
    const scalar T
) const
{
//    const scalar Tr=T-Tph_;
    const scalar Tr1 = T1_-Tph_;
    const scalar Tr2 = T2_-Tph_;
    scalar rslt = 0;
    if (T<Tr1)
    {
        rslt = 0;
    } else if (Tr1<=T && T<0)
    {
        rslt = T*(pow(T,3)-2*pow(T,2)*Tr1+2*pow(Tr1,3))/((Tr2-Tr1)*pow(Tr1,3))-Tr1/(Tr2-Tr1);
    } else if (T>=0 && T<=Tr2)
    {
        rslt = T*(pow(T,3)-2*pow(T,2)*Tr2+2*pow(Tr2,3))/(pow(Tr2,3)*(Tr2-Tr1))-Tr1/(Tr2-Tr1);
    } else rslt = 1;
    
    return rslt;
}

template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::dFf
(
    const scalar T
) const
{
    const scalar Tr1 = T1_-Tph_;
    const scalar Tr2 = T2_-Tph_;
    scalar rslt = 0;
    if (Tr1<=T && T<0)
    {
        rslt = 4*pow(T,3)/((Tr2-Tr1)*pow(Tr1,3))-6*pow(T,2)/((Tr2-Tr1)*pow(Tr1,2))+2.0/(Tr2-Tr1);
    } else if (T>=0 && T<=Tr2)
    {
        rslt = 4*pow(T,3)/(pow(Tr2,3)*(Tr2-Tr1))-6*pow(T,2)/(pow(Tr2,2)*(Tr2-Tr1))+2.0/(Tr2-Tr1);
    } else rslt = 0.0;

    return rslt;
}

template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::IFf
(
    const scalar T
) const
{
    const scalar Tr1 = T1_-Tph_;
    const scalar Tr2 = T2_-Tph_;
    scalar rslt = 0;
    if (T<Tr1)
    {
        rslt = 0;
    } else if (Tr1<=T && T<0)
    {
        rslt = (3./10)*pow(Tr1,2)/(Tr2-Tr1)-(1./10)*T*(10*pow(Tr1,4)-2*pow(T,4)+5*pow(T,3)*Tr1-10*T*pow(Tr1,3))/((Tr2-Tr1)*pow(Tr1,3));
    } else if (T>=0 && T<=Tr2)
    {
        rslt = (3./10)*pow(Tr1,2)/(Tr2-Tr1)-(1./10)*T*(10*Tr1*pow(Tr2,3)-2*pow(T,4)+5*pow(T,3)*Tr2-10*T*pow(Tr2,3))/(pow(Tr2,3)*(Tr2-Tr1));
    } else rslt = T-3./10*(Tr2+Tr1);

    return rslt;
}

template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return 0;//Cp_*log(T/Tstd) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::dGdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return 0;
}

template<class EquationOfState>
inline Foam::scalar Foam::hMeltingSolidificationCubicThermo<EquationOfState>::ffa
(
    const scalar T
) const
{
    return Ff(T-Tph_);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::hMeltingSolidificationCubicThermo<EquationOfState>::operator+=
(
    const hMeltingSolidificationCubicThermo<EquationOfState>& ct
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
        T1_ = Y1*T1_ + Y2*ct.T1_;
        T2_ = Y1*T2_ + Y2*ct.T2_;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::hMeltingSolidificationCubicThermo<EquationOfState> Foam::operator+
(
    const hMeltingSolidificationCubicThermo<EquationOfState>& ct1,
    const hMeltingSolidificationCubicThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return hMeltingSolidificationCubicThermo<EquationOfState>
        (
            eofs,
            ct1.CpS_,
            ct1.CpL_,
            ct1.Hf_,
            ct1.Tph_,
            ct1.T1_,
            ct1.T2_
        );
    }
    else
    {
        return hMeltingSolidificationCubicThermo<EquationOfState>
        (
            eofs,
            ct1.Y()/eofs.Y()*ct1.CpS_ + ct2.Y()/eofs.Y()*ct2.CpS_,
            ct1.Y()/eofs.Y()*ct1.CpL_ + ct2.Y()/eofs.Y()*ct2.CpL_,
            ct1.Y()/eofs.Y()*ct1.Hf_ + ct2.Y()/eofs.Y()*ct2.Hf_,
            ct1.Y()/eofs.Y()*ct1.Tph_ + ct2.Y()/eofs.Y()*ct2.Tph_,
            ct1.Y()/eofs.Y()*ct1.T1_ + ct2.Y()/eofs.Y()*ct2.T1_,
            ct1.Y()/eofs.Y()*ct1.T2_ + ct2.Y()/eofs.Y()*ct2.T2_
        );
    }
}


template<class EquationOfState>
inline Foam::hMeltingSolidificationCubicThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const hMeltingSolidificationCubicThermo<EquationOfState>& ct
)
{
    return hMeltingSolidificationCubicThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct),
        ct.CpS_,
        ct.CpL_,
        ct.Hf_,
        ct.Tph_,
        ct.T1_,
        ct.T2_
    );
}


template<class EquationOfState>
inline Foam::hMeltingSolidificationCubicThermo<EquationOfState> Foam::operator==
(
    const hMeltingSolidificationCubicThermo<EquationOfState>& ct1,
    const hMeltingSolidificationCubicThermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return hMeltingSolidificationCubicThermo<EquationOfState>
    (
        eofs,
        ct2.Y()/eofs.Y()*ct2.CpS_ - ct1.Y()/eofs.Y()*ct1.CpS_,
        ct2.Y()/eofs.Y()*ct2.CpL_ - ct1.Y()/eofs.Y()*ct1.CpL_,
        ct2.Y()/eofs.Y()*ct2.Hf_ - ct1.Y()/eofs.Y()*ct1.Hf_,
        ct2.Y()/eofs.Y()*ct2.Tph_ - ct1.Y()/eofs.Y()*ct1.Tph_,
        ct2.Y()/eofs.Y()*ct2.T1_ - ct1.Y()/eofs.Y()*ct1.T1_,
        ct2.Y()/eofs.Y()*ct2.T2_ - ct1.Y()/eofs.Y()*ct1.T2_
    );
}


// ************************************************************************* //
