/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "butterflyValveFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::butterflyValveFvPatchVectorField::
butterflyValveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    cyclicPatchName_(),
    cyclicPatchLabel_(-1),
    initWallSf_(0),
    initCyclicSf_(0),
    nbrCyclicSf_(0),
    openFraction_(0),
    openingTime_(0),
    openingRate_(0.1),
    curTimeIndex_(-1)
{}


Foam::butterflyValveFvPatchVectorField::
butterflyValveFvPatchVectorField
(
    const butterflyValveFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    openingRate_(ptf.openingRate_),
    curTimeIndex_(-1)
{}


Foam::butterflyValveFvPatchVectorField::
butterflyValveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    cyclicPatchName_(dict.lookup("cyclicPatch")),
    cyclicPatchLabel_(p.patch().boundaryMesh().findPatchID(cyclicPatchName_)),
    initWallSf_(p.Sf()),
    initCyclicSf_(p.boundaryMesh()[cyclicPatchLabel_].Sf()),
    nbrCyclicSf_
    (
        refCast<const cyclicFvPatch>
        (
            p.boundaryMesh()[cyclicPatchLabel_]
        ).neighbFvPatch().Sf()
    ),
    openFraction_(readScalar(dict.lookup("openFraction"))),
    openingTime_(readScalar(dict.lookup("openingTime"))),
    openingRate_(readScalar(dict.lookup("openingRate"))),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(Zero);
}


Foam::butterflyValveFvPatchVectorField::
butterflyValveFvPatchVectorField
(
    const butterflyValveFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    openingRate_(ptf.openingRate_),
    curTimeIndex_(-1)
{}


Foam::butterflyValveFvPatchVectorField::
butterflyValveFvPatchVectorField
(
    const butterflyValveFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    openFraction_(ptf.openFraction_),
    openingTime_(ptf.openingTime_),
    openingRate_(ptf.openingRate_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::butterflyValveFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    //- Note: cannot map field from cyclic patch anyway so just recalculate
    //  Areas should be consistent when doing autoMap except in case of
    //  topo changes.
    //- Note: we don't want to use Sf here since triggers rebuilding of
    //  fvMesh::S() which will give problems when mapped (since already
    //  on new mesh)
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::butterflyValveFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // See autoMap.
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::butterflyValveFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const fvPatch& cyclicPatch = patch().boundaryMesh()[cyclicPatchLabel_];
        const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
        (
            cyclicPatch
        ).neighbFvPatch();

        scalar t  =  this->db().time().timeOutputValue();
        if (t > openingTime_)
        {
            scalar dt =  this->db().time().deltaTValue();
            openFraction_ =
                max
                (
                    min
                    (
                        openFraction_+dt*openingRate_,
                        1 - 1e-6
                    ),
                    1e-6
                );
    
            Info<< "Valve time = " << t << ", openFraction = " << openFraction_ << endl;
        }

        vectorField::subField Sfw = this->patch().patch().faceAreas();
        const vectorField newSfw((1 - openFraction_)*initWallSf_);
        forAll(Sfw, facei)
        {
            Sfw[facei] = newSfw[facei];
        }
        const_cast<scalarField&>(patch().magSf()) = mag(patch().Sf());

        // Update owner side of cyclic
        const_cast<vectorField&>(cyclicPatch.Sf()) =
            openFraction_*initCyclicSf_;
        const_cast<scalarField&>(cyclicPatch.magSf()) =
            mag(cyclicPatch.Sf());
        // Update neighbour side of cyclic
        const_cast<vectorField&>(nbrPatch.Sf()) =
            openFraction_*nbrCyclicSf_;
        const_cast<scalarField&>(nbrPatch.magSf()) =
            mag(nbrPatch.Sf());

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::butterflyValveFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "cyclicPatch", cyclicPatchName_);
    writeEntry(os, "openFraction", openFraction_);
    writeEntry(os, "openingTime", openingTime_);
    writeEntry(os, "openingRate", openingRate_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        butterflyValveFvPatchVectorField
    );
}


// ************************************************************************* //
