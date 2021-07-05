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

Class
    Foam::butterflyValveFvPatchVectorField

Description
    This butterfly valve boundary condition simulates the opening of a baffle
    according to a prescribed opening profile, by merging the behaviours of 
    wall and cyclic conditions. The baffle joins two mesh regions, where the 
    open fraction determines the interpolation weights applied to each cyclic- 
    and neighbour-patch contribution.

    For a linear opening, the baffle has a constant opening or closing rate. 

        \f[
            x = x_{old} + R*frac{dt}
        \f]

    where

    \vartable
        x       | baffle open fraction [0-1]
        x_{old} | baffle open fraction on previous evaluation
        dt      | simulation time step
        R       | valve opening rate         
    \endvartable

    The open fraction is then applied to scale the patch areas.

Usage
    \table
        Property     | Description             | Required    | Default value
        cyclicPatch  | cylclic patch name      | yes         |
        openFraction | current opatch open fraction [0-1]| yes |
        openingTime  | time taken to open the baffle | yes |
        openingRate  | valve opening rate | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            butterflyValve;
        cyclicPatch     cyclic1;
        openFraction    0.01;
        openingTime     0.5;
        openingRate     0.1;     
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::cyclicFvPatchField

SourceFiles
    butterflyValveFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef butterflyValveFvPatchVectorField_H
#define butterflyValveFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class butterflyValveFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class butterflyValveFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the cyclic patch used when the active baffle is open
        word cyclicPatchName_;

        //- Index of the cyclic patch used when the active baffle is open
        label cyclicPatchLabel_;

        //- Initial wall patch areas
        vectorField initWallSf_;

        //- Initial this-side cyclic patch areas
        vectorField initCyclicSf_;

        //- Initial neighbour-side cyclic patch areas
        vectorField nbrCyclicSf_;

        //- Current fraction of the active baffle which is open
        scalar openFraction_;

        //- Time taken for the active baffle to open
        scalar openingTime_;

        //- Valve opening rate 
        scalar openingRate_;

        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("butterflyValve");


    // Constructors

        //- Construct from patch and internal field
        butterflyValveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        butterflyValveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given butterflyValveFvPatchVectorField
        //  onto a new patch
        butterflyValveFvPatchVectorField
        (
            const butterflyValveFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        butterflyValveFvPatchVectorField
        (
            const butterflyValveFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new butterflyValveFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        butterflyValveFvPatchVectorField
        (
            const butterflyValveFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new butterflyValveFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
