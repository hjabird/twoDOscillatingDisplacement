/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "twoDOscillatingDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    axis_(vector::zero),
    origin_(vector::zero),
    mean_pitch_angle_(0.0),
    pitch_amplitude_(0.0),
    omega_(0.0),
    p0_(p.localPoints()),
    heave_amplitude_(vector::zero),
    pitch_phase_offset_(0.0)     
{}


twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    axis_(dict.lookup("axis")),
    origin_(dict.lookup("origin")),
    mean_pitch_angle_(readScalar(dict.lookup("mean_pitch_angle"))),
    pitch_amplitude_(readScalar(dict.lookup("pitch_amplitude"))),
    omega_(readScalar(dict.lookup("omega"))),
    heave_amplitude_(dict.lookup("heave_amplitude")),
    pitch_phase_offset_(readScalar(dict.lookup("pitch_phase_offset_")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const twoDOscillatingDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    mean_pitch_angle_(ptf.mean_pitch_angle_),
    pitch_amplitude_(ptf.pitch_amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_, mapper),
    heave_amplitude_(ptf.heave_amplitude_),
    pitch_phase_offset_(ptf.pitch_phase_offset_)

{}


twoDOscillatingDisplacementPointPatchVectorField::
twoDOscillatingDisplacementPointPatchVectorField
(
    const twoDOscillatingDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    axis_(ptf.axis_),
    origin_(ptf.origin_),
    mean_pitch_angle_(ptf.mean_pitch_angle_),
    pitch_amplitude_(ptf.pitch_amplitude_),
    omega_(ptf.omega_),
    p0_(ptf.p0_),
    heave_amplitude_(ptf.heave_amplitude_),
    pitch_phase_offset_(ptf.pitch_phase_offset_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void twoDOscillatingDisplacementPointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);

    p0_.autoMap(m);
}


void twoDOscillatingDisplacementPointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const twoDOscillatingDisplacementPointPatchVectorField& aODptf =
        refCast<const twoDOscillatingDisplacementPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(aODptf, addr);

    p0_.rmap(aODptf.p0_, addr);
}


void twoDOscillatingDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();

    scalar angle = mean_pitch_angle_ 
		+ pitch_amplitude_ * sin(omega_ * t.value() + pitch_phase_offset_);
    vector axisHat = axis_/mag(axis_);
    vectorField p0Rel(p0_ - origin_);

    vectorField::operator=
    (
        p0Rel * (cos(angle) - 1)
      + (axisHat ^ p0Rel * sin(angle))
      + (heave_amplitude_ * sin(omega_ * t.value()))
      + (axisHat & p0Rel) * (1 - cos(angle))*axisHat
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void twoDOscillatingDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("axis")
        << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("origin")
        << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("mean_pitch_angle")
        << mean_pitch_angle_ << token::END_STATEMENT << nl;
    os.writeKeyword("pitch_amplitude")
        << pitch_amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega")
        << omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("heave_amplitude")
        << heave_amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("pitch_phase_offset")
        << pitch_phase_offset_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    twoDOscillatingDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
