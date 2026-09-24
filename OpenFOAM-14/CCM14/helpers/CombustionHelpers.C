/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2021 Shijie Xu, Shenghui Zhong
                             | Copyright (C) 2025 Yuchen Zhou
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM CCM extension.

    OpenFOAM CCM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM CCM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM CCM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Chemistry Coordinate Mapping (CCM) acceleration method for OpenFOAM-10

\*---------------------------------------------------------------------------*/

#include "CombustionHelpers.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::CombustionHelpers<ThermoType>::CombustionHelpers
(
    const fvMesh& mesh,
    const fluidMulticomponentThermo& thermo,
    const multicomponentMixture<ThermoType>& mixture,
    const PtrList<ThermoType>& specieThermos,
    const dictionary& CCMdict,
    const label& nSpecie,
    const PtrList<volScalarField>& Yvf,
    hashedWordList& updateVars,
    scalarList& JhCoeff,
    scalarList& JcCoeff,
    scalarList& JnCoeff,
    scalarList& JoCoeff,
    scalarList& Jh_h2oCoeff,
    scalarList& Jc_co2Coeff,
    scalarList& Jo_co2_h2oCoeff,
    volScalarField& Jh,
    volScalarField& Jc,
    volScalarField& Jn,
    volScalarField& Jo,
    volScalarField& Jh_h2o,
    volScalarField& Jc_co2,
    volScalarField& Jo_co2_h2o,
    volScalarField& J,
    volScalarField& phieq,
    volScalarField& chi,
    const Switch& oldStyleChi,
    const Switch& oldStylePhi,
    const scalar& fuelOtoC,
    hashedWordList& principalVars
)
:
    mesh_(mesh),
    thermo_(thermo),
    mixture_(mixture),
    specieThermos_(specieThermos),
    CCMdict_(CCMdict),
    nSpecie_(nSpecie),
    Yvf_(Yvf),
    updateVars_(updateVars),
    JhCoeff_(JhCoeff),
    JcCoeff_(JcCoeff),
    JnCoeff_(JnCoeff),
    JoCoeff_(JoCoeff),
    Jh_h2oCoeff_(Jh_h2oCoeff),
    Jc_co2Coeff_(Jc_co2Coeff),
    Jo_co2_h2oCoeff_(Jo_co2_h2oCoeff),
    Jh_(Jh),
    Jc_(Jc),
    Jn_(Jn),
    Jo_(Jo),
    Jh_h2o_(Jh_h2o),
    Jc_co2_(Jc_co2),
    Jo_co2_h2o_(Jo_co2_h2o),
    J_(J),
    phieq_(phieq),
    chi_(chi),
    oldStyleChi_(oldStyleChi),
    oldStylePhi_(oldStylePhi),
    fuelOtoC_(fuelOtoC),
    principalVars_(principalVars)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJc()
{
    Jc_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jc_ += JcCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJh()
{
    Jh_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jh_ += JhCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJo()
{
    Jo_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jo_ += JoCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJn()
{
    Jn_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jn_ += JnCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJh_h2o()
{
    Jh_h2o_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jh_h2o_ += Jh_h2oCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJc_co2()
{
    Jc_co2_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jc_co2_ += Jc_co2Coeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateJo_co2_h2o()
{
    Jo_co2_h2o_ = Zero;
    forAll(Yvf_, specieI)
    {
        Jo_co2_h2o_ += Jo_co2_h2oCoeff_[specieI] * Yvf_[specieI];
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updateChi()
{
    if (oldStyleChi_)
    {
        volVectorField gradMZ_0 = fvc::grad(phieq_);
        dimensionedScalar dimDiffutionCoeff(dimLength*dimLength/dimTime, 1.);
        chi_ = dimDiffutionCoeff*(gradMZ_0&gradMZ_0);
        forAll(chi_, i)
        {
            chi_[i] = Foam::exp(-chi_[i]);
        }
        forAll(chi_.boundaryField(), patchi)
        {
            forAll(chi_.boundaryField()[patchi], facei)
            {
                chi_.boundaryFieldRef()[patchi][facei] = Foam::exp(-chi_.boundaryField()[patchi][facei]);
            }
        }
    }
    else
    {
        volScalarField D
        (
            mesh_.lookupObjectRef<volScalarField>("thermo:kappa")
            /
            thermo_.Cp()/thermo_.rho()
        );
        const word& chiComponent = CCMdict_.lookup("chiComponent");
        volScalarField& Z(mesh_.objectRegistry::lookupObjectRef<volScalarField>(chiComponent));
        volVectorField gradZ(fvc::grad(Z));
        chi_ = 2*D*(gradZ&gradZ);
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updatePhi()
{
    if (oldStylePhi_)
    {
        phieq_ = 2*Jc_/atomicWeights["C"] + 0.5*Jh_/atomicWeights["H"]-fuelOtoC_*Jc_/atomicWeights["C"];
        phieq_ /= max(Jo_/atomicWeights["O"]-fuelOtoC_*Jc_/atomicWeights["C"], 1e-4);
    }
    else
    {
        phieq_ = 2*Jc_co2_/atomicWeights["C"] + 0.5*Jh_h2o_/atomicWeights["H"]
        -fuelOtoC_*Jc_co2_/atomicWeights["C"];
        phieq_ /= max(Jo_co2_h2o_/atomicWeights["O"]-fuelOtoC_*Jc_co2_/atomicWeights["C"], 1e-4);        
    }
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::updatePV()
{
    // Always update J* first because phieq and chi depend on them
    if (updateVars_.found("Jc"))
    {
        updateJc();
    }
    if (updateVars_.found("Jh"))
    {
        updateJh();
    }
    if (updateVars_.found("Jo"))
    {
        updateJo();
    }
    if (updateVars_.found("Jn"))
    {
        updateJn();
    }
    if (updateVars_.found("Jh_h2o"))
    {
        updateJh_h2o();
    }
    if (updateVars_.found("Jc_co2"))
    {
        updateJc_co2();
    }
    if (updateVars_.found("Jo_co2_h2o"))
    {
        updateJo_co2_h2o();
    }
    
    // Update J based on JElement setting
    if (updateVars_.found("J"))
    {
        const word& JElement = CCMdict_.lookup("JElement");
        if (JElement == "h")
        {
            J_ = Jh_;
        }
        else if (JElement == "c")
        {
            J_ = Jc_;
        }
        else if (JElement == "o")
        {
            J_ = Jo_;
        }
        else if (JElement == "n")
        {
            J_ = Jn_;
        }
        else
        {
            FatalErrorIn("CombustionHelpers<ThermoType>::updatePV()")
                << "Unknown JElement: " << JElement << exit(FatalError);
        }
    }
    
    // Update phieq (depends on J*)
    if (updateVars_.found("phieq"))
    {
        updatePhi();
    }
    
    // Update chi (depends on other variables)
    if (updateVars_.found("chi"))
    {
        updateChi();
    }
}



template<class ThermoType>
bool Foam::CombustionHelpers<ThermoType>::appendIfNotExisting
(
    const Foam::word& varName, 
    Foam::hashedWordList& listToBeAppended
)
{
    bool found = false;
    found = listToBeAppended.found(varName);
    if (!found)
    {
        listToBeAppended.append(varName);
    }
    return found;
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::outputHashedWordList(const hashedWordList& list) const
{
    forAll(list, i)
    {
        Info << list[i] << " ";
    }
    Info << endl;
}


template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::initPvLists()
{
    // Deal with variables to update
    if (principalVars_.found("phieq"))
    {
        if (oldStylePhi_)
        {
            appendIfNotExisting("Jc", updateVars_);
            appendIfNotExisting("Jh", updateVars_);
            appendIfNotExisting("Jo", updateVars_);
            appendIfNotExisting("phieq", updateVars_);
        }
        else
        {
            appendIfNotExisting("phieq", updateVars_);
            appendIfNotExisting("Jc_co2", updateVars_);
            appendIfNotExisting("Jh_h2o", updateVars_);
            appendIfNotExisting("Jo_co2_h2o", updateVars_);
        }
    }
    if (principalVars_.found("chi"))
    {
        appendIfNotExisting("chi", updateVars_);
        if (!CCMdict_.found("chiComponent"))
        {
            FatalErrorIn("CombustionHelpers<ThermoType>::initPvLists()")
                << "The chiComponent entry not found in CCM dictionary" 
                << exit(FatalError);
        }
    }
    if (principalVars_.found("J"))
    {
        appendIfNotExisting("J", updateVars_);
        if (!CCMdict_.found("JElement"))
        {
            FatalErrorIn("CombustionHelpers<ThermoType>::initPvLists()")
                << "The JElement entry not found in CCM dictionary"
                << exit(FatalError);
        }
        else
        {
            word JElement = CCMdict_.lookup("JElement");
            if (JElement == "C" || JElement == "c")
            {
                appendIfNotExisting("Jc", updateVars_);
            }
            else if (JElement == "H" || JElement == "h")
            {
                appendIfNotExisting("Jh", updateVars_);
            }
            else if (JElement == "N" || JElement == "n")
            {
                appendIfNotExisting("Jn", updateVars_);
            }
            else if (JElement == "O" || JElement == "o")
            {
                appendIfNotExisting("Jo", updateVars_);
            }
            else
            {
                FatalErrorIn("CombustionHelpers<ThermoType>::initPvLists()")
                    << "The JElement entry should be one of C, H, N, O"
                    << exit(FatalError);
            }
        }
    }


    Info << "Variables to be update: ";
    outputHashedWordList(updateVars_);
}





template<class ThermoType>
void Foam::CombustionHelpers<ThermoType>::initJCoeffs()
{
    label H2Oindex(-1);
    label CO2index(-1);

    if (mixture_.species().found("H2O"))
    {
        H2Oindex = mixture_.species()["H2O"];
    }
    if (mixture_.species().found("CO2"))
    {
        CO2index = mixture_.species()["CO2"];
    }
    
    // Initialize JcCoeff_, JnCoeff_, JoCoeff_, JhCoeff_
    const scalar MWh = atomicWeights["H"];
    const scalar MWc = atomicWeights["C"];
    const scalar MWn = atomicWeights["N"];
    const scalar MWo = atomicWeights["O"];

    for (label i=0; i<nSpecie_; i++)
    {
        const List<specieElement>& spComposition = mixture_.specieComposition(i);
        forAll(spComposition, elementI)
        {
            scalar MW = specieThermos_[i].W();
            if (spComposition[elementI].name() == "H")
            {
                JhCoeff_[i] = spComposition[elementI].nAtoms() * MWh / MW;
            }
            if (spComposition[elementI].name() == "C")
            {
                JcCoeff_[i] = spComposition[elementI].nAtoms() * MWc / MW;
            }
            if (spComposition[elementI].name() == "N")
            {
                JnCoeff_[i] = spComposition[elementI].nAtoms() * MWn / MW;
            }
            if (spComposition[elementI].name() == "O")
            {
                JoCoeff_[i] = spComposition[elementI].nAtoms() * MWo / MW;
            }
        }
    }

    // Initialize modified J coefficients
    Jh_h2oCoeff_ = JhCoeff_;
    Jc_co2Coeff_ = JcCoeff_;
    Jo_co2_h2oCoeff_ = JoCoeff_;

    if (H2Oindex != -1)
    {
        Jh_h2oCoeff_[H2Oindex] = 0.0;
        Jo_co2_h2oCoeff_[H2Oindex] = 0.0;
    }

    if (CO2index != -1)
    {
        Jc_co2Coeff_[CO2index] = 0.0;
        Jo_co2_h2oCoeff_[CO2index] = 0.0;
    }

}

// ************************************************************************* //