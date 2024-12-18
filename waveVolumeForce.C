/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "waveVolumeForce.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(waveVolumeForce, 0);
    addToRunTimeSelectionTable(option, waveVolumeForce, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::waveVolumeForce::waveVolumeForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    Um_(coeffs_.get<vector>("Um")),
    T_(coeffs_.get<scalar>("Period"))
{
    coeffs_.readEntry("fields", fieldNames_);

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::waveVolumeForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    dimensioned<vector> waveForce_
    (
        dimensionSet(0,1,-2,0,0),
        calVf(Um_,T_)
    );
    eqn += waveForce_;
    
}


void Foam::fv::waveVolumeForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    dimensioned<vector> waveForce_
    (
        dimensionSet(0,1,-2,0,0),
        calVf(Um_,T_)
    );
    eqn += rho*waveForce_;
}


Foam::vector Foam::fv::waveVolumeForce::calVf
(
    const vector& um,
    const scalar& pd
)
{   
    #define PI 3.1415926
    scalar SMALL = 1e-6;
    scalar t = mesh().time().value();
    vector umwave = 2*PI*um*cos(2*PI*t/(pd+SMALL))/(pd+SMALL);              //CHANGE HERE
    //Info << "Wave force at time " << t << ": " << umwave << endl;         //USE FOR DEBUG
    return umwave;
}

bool Foam::fv::waveVolumeForce::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
