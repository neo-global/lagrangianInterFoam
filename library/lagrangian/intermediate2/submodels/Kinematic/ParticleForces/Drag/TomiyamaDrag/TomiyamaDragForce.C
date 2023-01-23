/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "TomiyamaDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::TomiyamaDragForce<CloudType>::CdRe(const scalar Re, const scalar Eo, const word cont) const
{
    if(cont=="none")
    {
       return Re*max(min(16.0/Re*(1+0.15*pow(Re,0.687)),48.0/Re), 8.0*Eo/(3.0*(Eo+4.0)));
    }
    else if(cont=="mild")
    {
       return Re*max(min(24.0/Re*(1+0.15*pow(Re,0.687)),72.0/Re), 8.0*Eo/(3.0*(Eo+4.0)));
    }
    else if(cont=="high")
    {
       return Re*max(24.0/Re*(1+0.15*pow(Re,0.687)), 8.0*Eo/(3.0*(Eo+4.0)));
    }
    else
    {
	
       return Re*max(min(16.0/Re*(1+0.15*pow(Re,0.687)),48.0/Re), 8.0*Eo/(3.0*(Eo+4.0)));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TomiyamaDragForce<CloudType>::TomiyamaDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, true),
    sigma_(readScalar(this->coeffs().lookup("sigma")))
{}


template<class CloudType>
Foam::TomiyamaDragForce<CloudType>::TomiyamaDragForce
(
    const TomiyamaDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df),
    sigma_(df.sigma_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TomiyamaDragForce<CloudType>::~TomiyamaDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::TomiyamaDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(Zero, 0.0);
    word contamination(word(this->coeffs().lookup("contamination")));

    value.Sp() = mass*0.75*muc*CdRe(Re,p.Eo(td, sigma_), contamination)/(p.rho()*sqr(p.d()));

    return value;
}


// ************************************************************************* //
