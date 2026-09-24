/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "CCM_chemistryModel.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "cpuLoad.H"
#include <iomanip>
#include <iostream>

// * * * * * * * * * * * * * * * * Helpers  * * * * * * * * * * * * * * //
// Smallest power of two >= n (n>=1). Used for MPHashTable bucket counts.
static inline Foam::label ccmNextPow2(Foam::label n)
{
    Foam::label p = 1;
    while (p < n) p <<= 1;
    return p;
}

// Loose elementwise comparison used by Phase 1+ FatalError validators.
static inline bool ccmNearlyEq(Foam::scalar a, Foam::scalar b, Foam::scalar tol = 1e-10)
{
    return Foam::mag(a - b) <= tol * (1.0 + Foam::mag(a) + Foam::mag(b));
}

static inline bool ccmNearlyEq
(
    const Foam::scalarField& a,
    const Foam::scalarField& b,
    Foam::scalar tol = 1e-10
)
{
    if (a.size() != b.size()) return false;
    forAll(a, i) if (!ccmNearlyEq(a[i], b[i], tol)) return false;
    return true;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModels::CCM<ThermoType>::CCM
(
    const fluidMulticomponentThermo& thermo
)
:
    standard(thermo),
    log_(this->lookupOrDefault("log", false)),
    cpuLoad_(this->lookupOrDefault("cpuLoad", false)),
    jacobianType_
    (
        this->found("jacobian")
      ? jacobianTypeNames.read(this->lookup("jacobian"))
      : jacobianType::fast
    ),
    mixture_
    (
        dynamicCast<const multicomponentMixture<ThermoType>>(this->thermo())
    ),
    specieThermos_(mixture_.specieThermos()),
    reactions_(thermo.species(), specieThermos_, this->mesh(), *this),
    RR_(nSpecie_),
    Y_(nSpecie_),
    c_(nSpecie_),
    YTpWork_(scalarField(nSpecie_ + 2)),
    YTpYTpWork_(scalarSquareMatrix(nSpecie_ + 2)),
    fuelOtoC_(this->subDict("CCM").lookupOrDefault<scalar>("ratioOxygenToCarbonElementInFuel",0.0)),
    JhCoeff_(scalarList(nSpecie_, 0.)),
    JcCoeff_(scalarList(nSpecie_, 0.)),
    JnCoeff_(scalarList(nSpecie_, 0.)),
    JoCoeff_(scalarList(nSpecie_, 0.)),
    Jh_h2oCoeff_(scalarList(nSpecie_, 0.)),
    Jc_co2Coeff_(scalarList(nSpecie_, 0.)),
    Jo_co2_h2oCoeff_(scalarList(nSpecie_, 0.)),
    CCMdict_(this->subDict("CCM")),
    ccmVars_(nullptr),
    ccmDebug_(nullptr),
    combustionHelpers_(nullptr),
    updateVars_(wordList()),
    nSlice_(CCMdict_.lookupOrDefault<label>("nSlice",50)),
    ignoreMin_(CCMdict_.lookupOrDefault<scalar>("ignoreMin",1e-8)),
    Jh_
    (
        volScalarField
        (
            IOobject
            (
                "Jh",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jc_
    (
        volScalarField
        (
            IOobject
            (
                "Jc",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jn_
    (
        volScalarField
        (
            IOobject
            (
                "Jn",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jo_
    (
        volScalarField
        (
            IOobject
            (
                "Jo",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jh_h2o_
    (
        volScalarField
        (
            IOobject
            (
                "Jh_h2o",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jc_co2_
    (
        volScalarField
        (
            IOobject
            (
                "Jc_co2",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    Jo_co2_h2o_
    (
        volScalarField
        (
            IOobject
            (
                "Jo_co2_h2o",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    J_
    (
        volScalarField
        (
            IOobject
            (
                "J",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    phieq_
    (
        volScalarField
        (
            IOobject
            (
                "phieq",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless, 0)
        )
    ),
    chi_
    (
        volScalarField
        (
            IOobject
            (
                "chi",
                this->mesh().time().name(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimless/dimTime, 0)
        )
    ),
    maxiRepresentation_(29728),
    redistributeEvery_(-1),
    redistributedThisStep_(true),
    loads_(Pstream::nProcs()),
    localCores_(-1),

    chemComm_(MPI_COMM_WORLD),
    oldStyleChi_(CCMdict_.lookupOrDefault<Switch>("oldStyleChi",false)),
    oldStylePhi_(CCMdict_.lookupOrDefault<Switch>("oldStylePhi",false)),
    // accelerationInfo_("CCMAccInfo"),
    Ymax_(nSpecie_),
    debugTime_(CCMdict_.lookupOrDefault<Switch>("debugTime", true)),
    keyArena_(),
    keyLen_(mesh().nCells()),
    keyStride_(0),
    zoneRemainder_(mesh().nCells()),
    zoneHash_(mesh().nCells()),
    stepTimes_(16, 0.0),
    stepNames_(),
    stepTimer_(),
    currentStepIndex_(0),
    numECVarsToAdd_(CCMdict_.subOrEmptyDict("ecMode").lookupOrDefault<label>("numECVarsToAdd", 3)),
    numECVarsToRemove_(CCMdict_.subOrEmptyDict("ecMode").lookupOrDefault<label>("numECVarsToRemove", 1)),
    ecUpdateFreq_(CCMdict_.subOrEmptyDict("ecMode").lookupOrDefault<label>("updateFreq", 10)),
    currentStep_(0),
    initMode_(CCMdict_.lookupOrDefault<word>("mode", "default")),
    highMach_(CCMdict_.lookupOrDefault<Switch>("highMach", true)),
    optimizedODE_(CCMdict_.lookupOrDefault<Switch>("optimizedODE", false)),
    fastChemistryPtr_(
        CCMdict_.lookupOrDefault<Switch>("disableFastChemAlloc", false)
      ? autoPtr<basicFastChemistryModel>()
      : basicFastChemistryModel::New(mesh())
    )
{

    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "chemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction() << endl;

    if (log_)
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }

  Info << nl << endl;
    Info << "================= CCM Initialization ===================" << endl;
    // Memory monitoring at constructor start
    if (Pstream::myProcNo() <= 1)
    {
        std::ifstream memStatus0("/proc/self/status");
        std::string memLine0;
        while (std::getline(memStatus0, memLine0))
        {
            if (memLine0.substr(0, 5) == "VmRSS")
            {
                Pout << "CCM constructor START: " << memLine0 << endl;
                break;
            }
        }
    }

    // Initialize CCM Variable management class
    ccmVars_.set
    (
        new CCMVars<ThermoType>
        (
            mixture_,
            specieThermos_,
            this->mesh(),
            CCMdict_
        )
    );
    principalVars_ = ccmVars_->principalVars();
    ecEnabled_ = ccmVars_->ecMode();



    if (oldStyleChi_)
    {
        Info << "Using old style scalar dissipation rate" << endl;
    }
    else
    {
        Info << "Using correct scalar dissipation rate" << endl;
    }

    if (oldStylePhi_)
    {
        Info << "Using the equivalence ratio that will NOT change during the combustion process" << endl;
    }
    else
    {
        Info << "Using the equivalence ratio that will change during the combustion process" << endl;
    }

    if (highMach_)
    {
        Info << "Treating pressure variations in high Mach cases" << endl;
    }

    if (optimizedODE_)
    {
        Info << "Using optimized ODE solver routines" << endl;
        GetRRFunc_ = &chemistryModels::CCM<ThermoType>::getRRGivenYTP_Optimized;
    }
    else
    {
        Info << "Using OpenFOAM-based ODE solver routines" << endl;
        GetRRFunc_ = &chemistryModels::CCM<ThermoType>::getRRGivenYTP_Basic;
    }



    // Initialize CCM Debug helper
    ccmDebug_.set
    (
        new CCMDebug<ThermoType>
        (
            nSpecie_,
            Ymax_,
            mixture_,
            stepTimes_,
            stepNames_,
            currentStepIndex_
        )
    );

    // Initialize CCM Combustion helper
    combustionHelpers_.set
    (
        new CombustionHelpers<ThermoType>
        (
            this->mesh(),
            this->thermo(),
            mixture_,
            specieThermos_,
            CCMdict_,
            nSpecie_,
            Yvf_,
            updateVars_,
            JhCoeff_,
            JcCoeff_,
            JnCoeff_,
            JoCoeff_,
            Jh_h2oCoeff_,
            Jc_co2Coeff_,
            Jo_co2_h2oCoeff_,
            Jh_,
            Jc_,
            Jn_,
            Jo_,
            Jh_h2o_,
            Jc_co2_,
            Jo_co2_h2o_,
            J_,
            phieq_,
            chi_,
            oldStyleChi_,
            oldStylePhi_,
            fuelOtoC_,
            principalVars_
        )
    );

    // Initialize Element mass fractions (Jc, Jh, Jn, Jo ...)
    combustionHelpers_->initJCoeffs();

    // Initialize variable lists, including the list of principal variables (PV) 
    // and the PVs to be updated
    combustionHelpers_->initPvLists();

    // update phi/J/chi
    combustionHelpers_->updatePV();


    // For all principal variables, set the write option to AUTO_WRITE
    forAll(principalVars_, pi)
    {
        volScalarField& field = this->mesh().objectRegistry::lookupObjectRef<volScalarField>(principalVars_[pi]);
        field.writeOpt() = IOobject::AUTO_WRITE;
    }





    // Create the fields for the chemistry sources
    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + Yvf_[fieldi].name(),
                    this->mesh().time().name(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.T().mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    Info<< "CCMchemistryModel: Number of species = " << nSpecie_
        << " and reactions = " << nReaction() << endl;

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialised (by default 'active' is true)



    // Initialize inter-core communication
    string mode = CCMdict_.subOrEmptyDict("communicator").lookupOrDefault<word>("mode", "global");
    if (mode == "global")
    {
        localCores_ = Pstream::nProcs();
    }
    else if (mode == "distributed")
    {
        localCores_ = CCMdict_.subDict("communicator").lookupOrDefault<label>("localCores", Pstream::nProcs());
        if (localCores_ < 2  || Pstream::nProcs() % localCores_ != 0)
        {
            FatalErrorIn("CCMchemistryModel constructor")
                << "For distributed mode, localCores must be at least 2 and divide the total number of processes."
                << " localCores: " << localCores_
                << " nProcs: " << Pstream::nProcs()
                << exit(FatalError);
        }
        
        // redistribution
        redistributeEvery_ = CCMdict_.subDict("communicator").lookupOrDefault<label>("redistributeEvery", 5);

        // mpi
        label worldRank = Pstream::myProcNo();
        label color = worldRank/localCores_;
        label key = worldRank%localCores_;
        MPI_Comm_split(MPI_COMM_WORLD, color, key, &chemComm_);
    }
    
    


    // MPHashTable allocation. Sized once here; per-step uses clear().
    // perBlock estimates (bytes):
    //   TCEntry  : header + Y + Ystd + 8 scalars + pad  ~= 18*nSp + 120
    //   rateEntry: header + RR + dt + blockSize         ~= 10*nSp + 48
    // Per-core pool   = (nCells/N) * perBlock           [mpToCore, mpFromCore,
    //                                                    mpReturnToCore, mpReceivedRR]
    // Whole-mesh pool =  nCells     * perBlock           [mpGathered, mpLocalRR]
    {
        const Foam::label nL  = this->mesh().nCells();
        const Foam::label N   = localCores_;
        const Foam::label per = max(nL / N, Foam::label(20));
        const Foam::label perBlockTC = 18 * nSpecie_ + 120;
        const Foam::label perBlockRR = 10 * nSpecie_ + 48;
        const Foam::label nbPer = ccmNextPow2(per);
        const Foam::label nbAll = ccmNextPow2(max(nL, Foam::label(1)));

        mpToCore_.resize(N);
        mpFromCore_.resize(N);
        mpReturnToCore_.resize(N);
        mpReceivedRR_.resize(N);
        for (Foam::label i = 0; i < N; ++i)
        {
            mpToCore_.set     (i, new MPHashTable(nbPer, per * perBlockTC));
            mpFromCore_.set   (i, new MPHashTable(nbPer, per * perBlockTC));
            mpReturnToCore_.set(i, new MPHashTable(nbPer, per * perBlockRR));
            mpReceivedRR_.set  (i, new MPHashTable(nbPer, per * perBlockRR));
        }
        mpGathered_.set(new MPHashTable(nbAll, nL * perBlockTC));
        mpLocalRR_.set (new MPHashTable(nbAll, nL * perBlockRR));

        // Send / recv staging buffers. Strict upper bound on per-rank
        // outgoing bytes is nL * perBlockTC (degenerate case: every cell
        // routes to the same peer). Same bound covers rateEntry (smaller
        // perBlock) and incoming bytes (peer can't send more cells than
        // it owns; assumed symmetric load).
        const size_t bufBytes = size_t(nL) * size_t(perBlockTC);
        mpSendBuf_.assign(bufBytes, 0);
        mpRecvBuf_.assign(bufBytes, 0);

        Info<< "MPHashTable allocated: nCells=" << nL
            << " nCores=" << N << " nSpecie=" << nSpecie_
            << " perBlockTC=" << perBlockTC
            << " perBlockRR=" << perBlockRR
            << " bucketsPer=" << nbPer << " bucketsAll=" << nbAll
            << endl;
    }

    if (log_)
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }

    // Memory monitoring at constructor end
    if (Pstream::myProcNo() <= 1)
    {
        std::ifstream memStatusEnd("/proc/self/status");
        std::string memLineEnd;
        while (std::getline(memStatusEnd, memLineEnd))
        {
            if (memLineEnd.substr(0, 5) == "VmRSS")
            {
                Pout << "CCM constructor END: " << memLineEnd << endl;
                break;
            }
        }
    }
    Info << "=============== End CCM Initialization =================" << endl;
    Info << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::chemistryModels::CCM<ThermoType>::~CCM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::derivatives
(
    const scalar time,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt
) const
{

    forAll(Y_, i)
    {
        Y_[i] = max(YTp[i], 0);
    }


    const scalar T = YTp[nSpecie_];
    const scalar p = YTp[nSpecie_ + 1];

    // Evaluate the mixture density
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]/specieThermos_[i].rho(p, T);
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        c_[i] = rhoM/specieThermos_[i].W()*Y_[i];
    }

    // Evaluate contributions from reactions
    dYTpdt = Zero;
    forAll(reactions_, ri)
    {
        reactions_[ri].dNdtByV
        (
            p,
            T,
            c_,
            li,
            dYTpdt,
            false, //reduction_,
            cTos_,
            0
        );
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTpdt[i];
        dYidt *= WiByrhoM;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp
    scalar CpM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        CpM += Y_[i]*specieThermos_[i].Cp(p, T);
    }

    // dT/dt
    scalar& dTdt = dYTpdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        dTdt -= dYTpdt[i]*specieThermos_[sToc(i)].ha(p, T);
    }
    dTdt /= CpM;

    // dp/dt = 0 (pressure is assumed constant)
    scalar& dpdt = dYTpdt[nSpecie_ + 1];
    dpdt = 0;
}


template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::jacobian
(
    const scalar t,
    const scalarField& YTp,
    const label li,
    scalarField& dYTpdt,
    scalarSquareMatrix& J
) const
{
    forAll(c_, i)
    {
        Y_[i] = max(YTp[i], 0);
    }

    const scalar T = YTp[nSpecie_];
    const scalar p = YTp[nSpecie_ + 1];

    // Evaluate the specific volumes and mixture density
    scalarField& v = YTpWork_[0];
    for (label i=0; i<Y_.size(); i++)
    {
        v[i] = 1/specieThermos_[i].rho(p, T);
    }
    scalar rhoM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        rhoM += Y_[i]*v[i];
    }
    rhoM = 1/rhoM;

    // Evaluate the concentrations
    for (label i=0; i<Y_.size(); i ++)
    {
        c_[i] = rhoM/specieThermos_[i].W()*Y_[i];
    }

    // Evaluate the derivatives of concentration w.r.t. mass fraction
    scalarSquareMatrix& dcdY = YTpYTpWork_[0];
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar rhoMByWi = rhoM/specieThermos_[sToc(i)].W();
        switch (jacobianType_)
        {
            case jacobianType::fast:
                {
                    dcdY(i, i) = rhoMByWi;
                }
                break;
            case jacobianType::exact:
                for (label j=0; j<nSpecie_; j++)
                {
                    dcdY(i, j) =
                        rhoMByWi*((i == j) - rhoM*v[sToc(j)]*Y_[sToc(i)]);
                }
                break;
        }
    }

    // Evaluate the mixture thermal expansion coefficient
    scalar alphavM = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        alphavM += Y_[i]*rhoM*v[i]*specieThermos_[i].alphav(p, T);
    }

    // Evaluate contributions from reactions
    dYTpdt = Zero;
    scalarSquareMatrix& ddNdtByVdcTp = YTpYTpWork_[1];
    for (label i=0; i<nSpecie_ + 2; i++)
    {
        for (label j=0; j<nSpecie_ + 2; j++)
        {
            ddNdtByVdcTp[i][j] = 0;
        }
    }
    forAll(reactions_, ri)
    {
        reactions_[ri].ddNdtByVdcTp
        (
            p,
            T,
            c_,
            li,
            dYTpdt,
            ddNdtByVdcTp,
            false,//reduction_,
            cTos_,
            0,
            nSpecie_,
            YTpWork_[1],
            YTpWork_[2]
        );
    }

    // Reactions return dNdtByV, so we need to convert the result to dYdt
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar WiByrhoM = specieThermos_[sToc(i)].W()/rhoM;
        scalar& dYidt = dYTpdt[i];
        dYidt *= WiByrhoM;

        for (label j=0; j<nSpecie_; j++)
        {
            scalar ddNidtByVdYj = 0;
            switch (jacobianType_)
            {
                case jacobianType::fast:
                    {
                        const scalar ddNidtByVdcj = ddNdtByVdcTp(i, j);
                        ddNidtByVdYj = ddNidtByVdcj*dcdY(j, j);
                    }
                    break;
                case jacobianType::exact:
                    for (label k=0; k<nSpecie_; k++)
                    {
                        const scalar ddNidtByVdck = ddNdtByVdcTp(i, k);
                        ddNidtByVdYj += ddNidtByVdck*dcdY(k, j);
                    }
                    break;
            }

            scalar& ddYidtdYj = J(i, j);
            ddYidtdYj = WiByrhoM*ddNidtByVdYj + rhoM*v[sToc(j)]*dYidt;
        }

        scalar ddNidtByVdT = ddNdtByVdcTp(i, nSpecie_);
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddNidtByVdcj = ddNdtByVdcTp(i, j);
            ddNidtByVdT -= ddNidtByVdcj*c_[sToc(j)]*alphavM;
        }

        scalar& ddYidtdT = J(i, nSpecie_);
        ddYidtdT = WiByrhoM*ddNidtByVdT + alphavM*dYidt;

        scalar& ddYidtdp = J(i, nSpecie_ + 1);
        ddYidtdp = 0;
    }

    // Evaluate the effect on the thermodynamic system ...

    // Evaluate the mixture Cp and its derivative
    scalarField& Cp = YTpWork_[3];
    scalar CpM = 0, dCpMdT = 0;
    for (label i=0; i<Y_.size(); i++)
    {
        Cp[i] = specieThermos_[i].Cp(p, T);
        CpM += Y_[i]*Cp[i];
        dCpMdT += Y_[i]*specieThermos_[i].dCpdT(p, T);
    }

    // dT/dt
    scalarField& ha = YTpWork_[4];
    scalar& dTdt = dYTpdt[nSpecie_];
    for (label i=0; i<nSpecie_; i++)
    {
        ha[sToc(i)] = specieThermos_[sToc(i)].ha(p, T);
        dTdt -= dYTpdt[i]*ha[sToc(i)];
    }
    dTdt /= CpM;

    // dp/dt = 0 (pressure is assumed constant)
    scalar& dpdt = dYTpdt[nSpecie_ + 1];
    dpdt = 0;

    // d(dTdt)/dY
    for (label i=0; i<nSpecie_; i++)
    {
        scalar& ddTdtdYi = J(nSpecie_, i);
        ddTdtdYi = 0;
        for (label j=0; j<nSpecie_; j++)
        {
            const scalar ddYjdtdYi = J(j, i);
            ddTdtdYi -= ddYjdtdYi*ha[sToc(j)];
        }
        ddTdtdYi -= Cp[sToc(i)]*dTdt;
        ddTdtdYi /= CpM;
    }

    // d(dTdt)/dT
    scalar& ddTdtdT = J(nSpecie_, nSpecie_);
    ddTdtdT = 0;
    for (label i=0; i<nSpecie_; i++)
    {
        const scalar dYidt = dYTpdt[i];
        const scalar ddYidtdT = J(i, nSpecie_);
        ddTdtdT -= dYidt*Cp[sToc(i)] + ddYidtdT*ha[sToc(i)];
    }
    ddTdtdT -= dTdt*dCpMdT;
    ddTdtdT /= CpM;

    // d(dTdt)/dp = 0 (pressure is assumed constant)
    scalar& ddTdtdp = J(nSpecie_, nSpecie_ + 1);
    ddTdtdp = 0;

    // d(dpdt)/dYiTp = 0 (pressure is assumed constant)
    for (label i=0; i<nSpecie_ + 2; i++)
    {
        scalar& ddpdtdYiTp = J(nSpecie_ + 1, i);
        ddpdtdYiTp = 0;
    }
}


template<class ThermoType>
Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::chemistryModels::CCM<ThermoType>::reactionRR
(
    const label reactioni
) const
{
    tmp<volScalarField::Internal> tRR =
        volScalarField::Internal::New
        (
            "RR:" + reactions_[reactioni].name(),
            this->mesh(),
            dimensionedScalar(dimMoles/dimVolume/dimTime, 0)
        );
    volScalarField::Internal& RR = tRR.ref();

    if (!this->chemistry_)
    {
        return tRR;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    reactionEvaluationScope scope(*this);

    const Reaction<ThermoType>& R = reactions_[reactioni];

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rho*Yi/specieThermos_[i].W();
        }

        scalar omegaf, omegar;

        RR[celli] =
            R.omega
            (
                p,
                T,
                c_,
                celli,
                omegaf,
                omegar
            );
    }

    return tRR;
}


template<class ThermoType>
Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::chemistryModels::CCM<ThermoType>::specieReactionRR
(
    const label reactioni
) const
{
    PtrList<volScalarField::Internal> RR(nSpecie_);
    for (label i=0; i<nSpecie_; i++)
    {
        RR.set
        (
            i,
            volScalarField::Internal::New
            (
                "RR:" + reactions_[reactioni].name() + ":" + Yvf_[i].name(),
                this->mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            ).ptr()
        );
    }

    if (!this->chemistry_)
    {
        return RR;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    scalarField& dNdtByV = YTpWork_[0];

    reactionEvaluationScope scope(*this);

    const Reaction<ThermoType>& R = reactions_[reactioni];

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rho*Yi/specieThermos_[i].W();
        }

        dNdtByV = Zero;

        R.dNdtByV
        (
            p,
            T,
            c_,
            celli,
            dNdtByV,
            false,//reduction_,
            cTos_,
            0
        );

        for (label i=0; i<nSpecie_; i++)
        {
            RR[i][celli] = dNdtByV[i]*specieThermos_[i].W();
        }
    }

    return RR;
}


template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::calculate()
{
    if (!this->chemistry_)
    {
        return;
    }

    if (!zone_.all())
    {
        forAll(RR_, fieldi)
        {
            RR_[fieldi] = Zero;
        }
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    scalarField& dNdtByV = YTpWork_[0];

    reactionEvaluationScope scope(*this);

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            const scalar Yi = Yvf_[i][celli];
            c_[i] = rho*Yi/specieThermos_[i].W();
        }

        dNdtByV = Zero;

        forAll(reactions_, ri)
        {
            reactions_[ri].dNdtByV
            (
                p,
                T,
                c_,
                celli,
                dNdtByV,
                false,//reduction_,
                cTos_,
                0
            );
        }

        for (label i=0; i<nSpecie_; i++)
        {
            RR_[sToc(i)][celli] = dNdtByV[i]*specieThermos_[sToc(i)].W();
        }
    }
}


template<class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::chemistryModels::CCM<ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
   Info << nl << endl;

   zone_.regenerate();
    Info << "==================== CCM Solution ======================" << endl;
    scalar deltaTMin = great;
    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    // Reset redistribute flag for this timestep
    redistributedThisStep_ = false;

    // update zoneIndex_ and zoneRemainder_ if mesh has changed
    updateCCM4MeshChange();



    // Reset step timing array
    if (debugTime_)
    {
        stepTimes_ = 0.0;

        // Reset automatic step counter
        currentStepIndex_ = 0;
        stepNames_.clear();

        // Reset timer
        stepTimer_.cpuTimeIncrement(); // Reset timer
    }


    
    // Derived fields (phieq/chi/J) depend only on oldTime primitives,
    // so updating once per solve() is enough -- distribute may run
    // multiple times inside an examine step.
    combustionHelpers_->updatePV();
    ccmVars_->updateMinMaxSpan();

    const bool examine = ecEnabled_ && (currentStep_ % ecUpdateFreq_ == 0);

    if (examine)
    {
        Info << "EC update at step " << currentStep_ << " with ecVars: ";
        ccmVars_->outputHashedWordList(ccmVars_->ecVars());

        
        this->distributeReactionEntry(ccmVars_->ecVars(), true);

        scalarField stdnorm = maxYstdNorm();
        ccmVars_->sortSpeciesByStdNorm(stdnorm);

        if (ccmVars_->ecOkay(stdnorm))
        {
            // First trial in control -> try trimming. removeN returns the
            // names it demoted; if the trim reshapes the clustering enough
            // to push another unused species above tolerance, just addVars
            // them back.
            wordList removed = ccmVars_->removeN(numECVarsToRemove_, stdnorm);

            if (removed.size() > 0)
            {
                this->distributeReactionEntry(ccmVars_->ecVars(), true);
                stdnorm = maxYstdNorm();
                ccmVars_->sortSpeciesByStdNorm(stdnorm);

                if (!ccmVars_->ecOkay(stdnorm))
                {
                    Info << "removeN broke ecOkay; rolling back: "
                         << removed << endl;
                    ccmVars_->addVars(removed);
                    this->distributeReactionEntry(ccmVars_->ecVars(), true);
                    stdnorm = maxYstdNorm();
                    ccmVars_->sortSpeciesByStdNorm(stdnorm);
                }
            }
        }
        else
        {
            // First trial not in control -> add until it is.
            while (!ccmVars_->ecOkay(stdnorm))
            {
                ccmVars_->addN(numECVarsToAdd_, stdnorm);
                this->distributeReactionEntry(ccmVars_->ecVars(), true);
                stdnorm = maxYstdNorm();
                ccmVars_->sortSpeciesByStdNorm(stdnorm);
            }
        }
    }
    else if (ecEnabled_)
    {
        this->distributeReactionEntry(ccmVars_->ecVars(), false);
    }
    else
    {
        this->distributeReactionEntry();
    }

    ccmVars_->reportMinMax();
    Info << "Executed ecVars:" << endl;
    ccmVars_->outputHashedWordList(ccmVars_->ecVars());


    updateReactionRate();

    distributeReactionRate();

    if (zone_.all())
    {
        deltaTMin = min(deltaTChem_).value();
    }
    else
    {
        deltaTMin = GREAT;
        for (label zci = 0; zci<zone_.nCells(); zci++)
        {
            const label celli = zone_.celli(zci);
            deltaTMin = min(deltaTMin, deltaTChem_[celli]);
        }
    }

    // Output timing analysis if enabled
    if (debugTime_)
    {
        ccmDebug_->outputTimingAnalysis(stepTimes_);
    }

    // Memory monitoring - output VmRSS for master and one slave
    if (Pstream::myProcNo() <= 1)
    {
        std::ifstream memStatus("/proc/self/status");
        std::string memLine;
        while (std::getline(memStatus, memLine))
        {
            if (memLine.substr(0, 5) == "VmRSS")
            {
                Pout << "Step " << currentStep_ << " " << memLine << endl;
                break;
            }
        }
    }

    // Increment step counter
    currentStep_++;

    // Restore default stream format before leaving solve()
    // (CCM uses fixed/scientific + setprecision which affects global std::cout)
    std::cout.unsetf(std::ios_base::fixed | std::ios_base::scientific);
    std::cout << std::setprecision(6);

    Info << "================== End CCM Solution ====================" << endl;
    Info << nl << endl;
    return deltaTMin;
}


template<class ThermoType>
Foam::scalar Foam::chemistryModels::CCM<ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ThermoType>
Foam::scalar Foam::chemistryModels::CCM<ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModels::CCM<ThermoType>::tc() const
{
    tmp<volScalarField> ttc
    (
        volScalarField::New
        (
            "tc",
            this->mesh(),
            dimensionedScalar(dimTime, small),
            extrapolatedCalculatedFvPatchScalarField::typeName
        )
    );
    scalarField& tc = ttc.ref();

    if (!this->chemistry_)
    {
        ttc.ref().correctBoundaryConditions();
        return ttc;
    }

    tmp<volScalarField> trhovf(this->thermo().rho());
    const volScalarField& rhovf = trhovf();

    const volScalarField& Tvf = this->thermo().T();
    const volScalarField& pvf = this->thermo().p();

    reactionEvaluationScope scope(*this);

    const label nZoneCells = zone_.nCells();
    for(label zci = 0; zci<nZoneCells; zci++)
    {
        const label celli = zone_.celli(zci);

        const scalar rho = rhovf[celli];
        const scalar T = Tvf[celli];
        const scalar p = pvf[celli];

        for (label i=0; i<nSpecie_; i++)
        {
            c_[i] = rho*Yvf_[i][celli]/specieThermos_[i].W();
        }

        // A reaction's rate scale is calculated as its molar
        // production rate divided by the total number of moles in the
        // system.
        //
        // The system rate scale is the average of the reactions' rate
        // scales weighted by the reactions' molar production rates. This
        // weighting ensures that dominant reactions provide the largest
        // contribution to the system rate scale.
        //
        // The system time scale is then the reciprocal of the system rate
        // scale.
        //
        // Contributions from forward and reverse reaction rates are
        // handled independently and identically so that reversible
        // reactions produce the same result as the equivalent pair of
        // irreversible reactions.

        scalar sumW = 0, sumWRateByCTot = 0;
        forAll(reactions_, i)
        {
            const Reaction<ThermoType>& R = reactions_[i];
            scalar omegaf, omegar;
            R.omega(p, T, c_, celli, omegaf, omegar);

            scalar wf = 0;
            forAll(R.rhs(), s)
            {
                wf += R.rhs()[s].stoichCoeff*omegaf;
            }
            sumW += wf;
            sumWRateByCTot += sqr(wf);

            scalar wr = 0;
            forAll(R.lhs(), s)
            {
                wr += R.lhs()[s].stoichCoeff*omegar;
            }
            sumW += wr;
            sumWRateByCTot += sqr(wr);
        }

        tc[celli] =
            sumWRateByCTot == 0 ? vGreat : sumW/sumWRateByCTot*sum(c_);
    }

    ttc.ref().correctBoundaryConditions();
    return ttc;
}


template<class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::chemistryModels::CCM<ThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (!this->chemistry_)
    {
        return tQdot;
    }

    reactionEvaluationScope scope(*this);

    scalarField& Qdot = tQdot.ref();

    forAll(Yvf_, i)
    {
        const label nZoneCells = zone_.nCells();
        for(label zci = 0; zci<nZoneCells; zci++)
        {
            const label celli = zone_.celli(zci);

            const scalar hi = specieThermos_[i].hf();
            Qdot[celli] -= hi*RR_[i][celli];
        }
    }

    return tQdot;
}


template<class ThermoType>
Foam::scalarField 
Foam::chemistryModels::CCM<ThermoType>::getRRGivenYTP
(
    scalarField Y,
    scalar& T,
    scalar& p,
    const scalar& deltaT,
    scalar& deltaTChem,
    const scalar& rho,
    const scalar& rho0
)
{

        // Initialise time progress
    scalar timeLeft = deltaT;

    scalarField Yupdate(Y);
    scalar Tupdate(T);
    scalar pupdate(p);

    // Calculate the chemical source terms
    while (timeLeft > small)
    {
        scalar dt = timeLeft;

        solve(pupdate, Tupdate, Yupdate, 0, dt, deltaTChem);

        timeLeft -= dt;
    }

    // update RR
    scalarField RR(nSpecie_);
    for (label i=0; i<nSpecie_; i++)
    {
        RR[i] = rho*(Yupdate[i] - Y[i])/deltaT;
    }

    return RR;
}

template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::solve
(
    scalar& p,
    scalar& T,
    scalarField& c,
    const label li,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    // Reset the size of the ODE system to the simplified size when mechanism
    // reduction is active
    if (odeSolver_->resize())
    {
        odeSolver_->resizeField(cTp_);
    }

    const label nSpecie = this->nSpecie();

    // Copy the concentration, T and P to the total solve-vector
    for (int i=0; i<nSpecie; i++)
    {
        cTp_[i] = c[i];
    }
    cTp_[nSpecie] = T;
    cTp_[nSpecie+1] = p;

    if (debug)
    {
        scalarField dcTp(this->nEqns(), rootSmall);
        dcTp[nSpecie] = T*rootSmall;
        dcTp[nSpecie+1] = p*rootSmall;
        this->check(0, cTp_, dcTp, li);
    }

    odeSolver_->solve(0, deltaT, cTp_, li, subDeltaT);

    for (int i=0; i<nSpecie; i++)
    {
        c[i] = max(0.0, cTp_[i]);
    }
    T = cTp_[nSpecie];
    p = cTp_[nSpecie+1];
}


template<class ThermoType>
Foam::scalarField Foam::chemistryModels::CCM<ThermoType>::maxYstdNorm() const
{
    scalarField Ymaxstd(nSpecie_, 0.0);
    {
        const char* src = mpGathered_->poolStart();
        const label nE  = mpGathered_->nEntries();
        for (label i = 0; i < nE; ++i)
        {
            const auto* hdr =
                reinterpret_cast<const MPHashTable::Entry*>(src);
            const label bs = hdr->blockSize;
            TCEntry v
            (
                const_cast<void*>(static_cast<const void*>(src)),
                nSpecie_
            );
            if (v.getN() > 1)
            {
                const scalarField yStd = v.getYstd();
                forAll(yStd, yi)
                {
                    Ymaxstd[yi] = max(Ymaxstd[yi], yStd[yi]);
                }
            }
            src += bs;
        }
    }
    Pstream::listCombineGather(Ymaxstd, maxEqOp<scalar>());
    Pstream::listCombineScatter(Ymaxstd);

    // Per-species global max mass fraction, computed locally so this
    // function is self-contained (no precondition on Ymax_).
    scalarField Ymaxvalue(nSpecie_, 0.0);
    forAll(Yvf_, i)
    {
        Ymaxvalue[i] = max(Yvf_[i].primitiveField());
    }
    Pstream::listCombineGather(Ymaxvalue, maxEqOp<scalar>());
    Pstream::listCombineScatter(Ymaxvalue);

    return Ymaxstd / (Ymaxvalue + 1e-8);
}

template<class ThermoType>
Foam::label Foam::chemistryModels::CCM<ThermoType>::getNextStepIndex() const
{
    label index = currentStepIndex_;
    currentStepIndex_++;
    
    // Ensure stepTimes_ array is large enough
    if (currentStepIndex_ >= stepTimes_.size())
    {
        // Double the size to accommodate more steps
        stepTimes_.setSize(stepTimes_.size()*2, 0.0);
    }
    
    return index;
}


template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::mpAlltoallvBytes
(
    int N,
    const void* const srcPtrs[],
    const int srcBytes[],
    const int srcEntries[],
    const char* outRecvPtrs[],
    int outRecvBytes[],
    int outRecvEntries[]
) const
{
    if (Pstream::parRun() && N > 1)
    {
        // Exchange per-peer byte counts and entry counts so each receiver
        // can size its segment before the payload Alltoallv.
        std::vector<int> sb(srcBytes,   srcBytes   + N);
        std::vector<int> se(srcEntries, srcEntries + N);
        std::vector<int> rb(N), re(N);
        MPI_Alltoall(sb.data(), 1, MPI_INT, rb.data(), 1, MPI_INT,
                     chemComm_);
        MPI_Alltoall(se.data(), 1, MPI_INT, re.data(), 1, MPI_INT,
                     chemComm_);

        // Prefix-sum to per-peer offsets into mpSendBuf_/mpRecvBuf_.
        std::vector<int> sd(N, 0), rd(N, 0);
        int sTot = 0, rTot = 0;
        for (int p = 0; p < N; ++p)
        {
            sd[p] = sTot;  sTot += sb[p];
            rd[p] = rTot;  rTot += rb[p];
        }
        if
        (
            size_t(sTot) > mpSendBuf_.size()
         || size_t(rTot) > mpRecvBuf_.size()
        )
        {
            FatalErrorInFunction
                << "mpAlltoallvBytes: staging buffer too small: sTot="
                << sTot << " rTot=" << rTot
                << " bufSize=" << mpSendBuf_.size()
                << exit(FatalError);
        }

        // Pack each peer's src into the staging buffer at its prefix slot.
        for (int p = 0; p < N; ++p)
        {
            if (sb[p] > 0)
            {
                std::memcpy
                (
                    mpSendBuf_.data() + sd[p],
                    srcPtrs[p],
                    size_t(sb[p])
                );
            }
        }

        MPI_Alltoallv
        (
            mpSendBuf_.data(), sb.data(), sd.data(), MPI_BYTE,
            mpRecvBuf_.data(), rb.data(), rd.data(), MPI_BYTE,
            chemComm_
        );

        for (int p = 0; p < N; ++p)
        {
            outRecvPtrs[p]    = mpRecvBuf_.data() + rd[p];
            outRecvBytes[p]   = rb[p];
            outRecvEntries[p] = re[p];
        }
    }
    else
    {
        // Serial / single-peer: route srcPtrs[0] through unchanged.
        outRecvPtrs[0]    = static_cast<const char*>(srcPtrs[0]);
        outRecvBytes[0]   = srcBytes[0];
        outRecvEntries[0] = srcEntries[0];
    }
}


template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::distributeReactionEntry
(
    const hashedWordList& encodingVars,
    bool debug
)
{

    // Running species/non-species sets live in CCMVars; both are kept in
    // principalVars_ order so the encoding ordering stays deterministic
    // across ranks. encodingVars is no longer consulted -- it was always
    // ccmVars_->ecVars() (or principalVars_ in non-ec mode, which equals
    // ecVars()), and CCMVars exposes the same partition directly.
    (void)encodingVars;

    const labelList&             runSp  = ccmVars_->runningSpIndices();
    const List<volScalarField*>& runNS  = ccmVars_->runningNonSpFields();
    const label nRunSp = runSp.size();
    const label nRunNS = runNS.size();
    const label nEnc   = nRunSp + nRunNS;

    // Pull per-var (min, max, span) from CCMVars (kept fresh by
    // ccmVars_->updateMinMaxSpan() in solve()), split by kind so the cell
    // loop indexes them directly without a varsToUse indirection.
    scalarList encMinSp (nRunSp);
    scalarList encSpanSp(nRunSp);
    scalarList encMinNS (nRunNS);
    scalarList encSpanNS(nRunNS);

    // Step 1: Initialization & Setup
    CCM_TIMING_START(Initialization, debugTime_, stepTimer_);
    // combustionHelpers_->updatePV() is hoisted to solve() (once per step).

    {
        const auto& vIdx = ccmVars_->varIndex();
        const scalarList& mn = ccmVars_->minvals();
        const scalarList& sp = ccmVars_->spans();
        const hashedWordList& species = mixture_.species();
        forAll(runSp, i)
        {
            const std::size_t pi = vIdx.at(species[runSp[i]]);
            encMinSp [i] = mn[pi];
            encSpanSp[i] = sp[pi];
        }
        forAll(runNS, i)
        {
            const std::size_t pi = vIdx.at(runNS[i]->name());
            encMinNS [i] = mn[pi];
            encSpanNS[i] = sp[pi];
        }
    }


    // CPU time analysis
    cpuTime solveCpuTime_;

    // basicChemistryModel::correct();

    // tmp<volScalarField> trhovf(this->thermo().rho());
    // const volScalarField& rhovf = trhovf();
    const volScalarField& rho0vf =
        this->mesh().template lookupObject<volScalarField>
        (
            this->thermo().phasePropertyName("rho")
        ).oldTime();
    const volScalarField& T0vf = this->thermo().T().oldTime();
    const volScalarField& p0vf = this->thermo().p().oldTime();
    reactionEvaluationScope scope(*this);
    scalarField Y0(nSpecie_);


    // chemistryCpuTime.reset();

    // Number of peers in the MPI_COMM_WORLD all-to-all. Pstream::nProcs()
    // returns 1 in serial mode and the world size in parallel — matches
    // the sizing of mpToCore_/mpFromCore_/mpReturnToCore_/mpReceivedRR_.

    // Clear mpToCore_ at step start. The other mp* tables are cleared
    // in the phase that fills them (mpFromCore_/mpGathered_ in the merge,
    // mpReturnToCore_/mpReceivedRR_ in distributeReactionRate).
    forAll(mpToCore_, i)
    {
        mpToCore_[i].clear();
    }

    CCM_TIMING_END(Initialization, debugTime_, stepTimer_, *this);

    // logP encoding setup: calculate log(P) min/max/span if enabled
    scalar pMin = -1, pMax = -1, pSpan = -1;
    bool useLogPOnThisRun = true;
    if (highMach_)
    {
        // Calculate log(P) range over all cells
        pMax = gMax(p0vf);
        pMin = max(gMin(p0vf),1);// 1 Pa minimum to avoid log(0)

        Info << nl << "High Mach Pressure: [" << pMin << ", " << pMax << "] Pa, ratio=" << pMax/pMin;

        if (pMax/pMin < 1.1)
        {
            useLogPOnThisRun = false;
            Info << "pMax = " << pMax << ", pMin = " << pMin << ", ratio < 1.1, log(P) encoding disabled" << nl;
        }
        else
        {
            Info << "log(P) encoding enabled" << nl;
            pMax = Foam::log(pMax);
            pMin = Foam::log(pMin);
            pSpan = (pMax - pMin) / nSlice_;
        }
    }


    // Step 2: Cell Grouping & Zone Calculation

    CCM_TIMING_START(Cell_Grouping, debugTime_, stepTimer_);

    // Phase 1.1: scratch buffer used by the merge branch of the MPHashTable
    // mirror to build a "single-cell" TCEntry view (count=1) without touching
    // the destination pool. Sized once at the worst-case keyLen for this call.
    const label scratchKeyLen   = nEnc * 2 + 2;  // listToBase256Word + up to 2 char logP
    const label scratchUserSize = TCEntry::userSize(scratchKeyLen, nSpecie_, debug);
    const label scratchTotal    = MPHashTable::headerBytes() + scratchUserSize;
    std::vector<char> mergeScratch(size_t(scratchTotal), 0);

    labelList posList(nEnc);

    // Loop tiling: process cells in blocks of GROUP_TILE. The per-field
    // values (Y for every specie, plus T/p/rho/rho0) are first streamed
    // into small tile buffers with one contiguous read per field, then
    // the per-cell encode + hash-table work reads back from those small,
    // cache-resident buffers. This turns the (4 + nSpecie_) simultaneous
    // strided field streams of the naive loop into one stream at a time.
    const label nGroupCells = rho0vf.size();
    const label GROUP_TILE   = 16;

    // Size the contiguous key arena for this step. keyStride_ is the
    // worst-case key length (nEnc*2 running-var bytes + up to 2 logP
    // bytes); every cell's key lives at celli*keyStride_ so no per-cell
    // offset table is needed. setSize only reallocates when nEnc grows.
    keyStride_ = nEnc * 2 + 2;
    const label arenaNeed = nGroupCells * keyStride_;
    if (keyArena_.size() < arenaNeed)
    {
        keyArena_.setSize(arenaNeed);
    }
    if (keyLen_.size() != nGroupCells)
    {
        keyLen_.setSize(nGroupCells);
    }

    if (zone_.all())
    {
        // species-major tile: Ytile[i*GROUP_TILE + t] = Y of specie i, cell t
        std::vector<scalar> Ytile(size_t(nSpecie_) * GROUP_TILE, 0.0);
        std::vector<scalar> Ttile (size_t(GROUP_TILE), 0.0);
        std::vector<scalar> Ptile (size_t(GROUP_TILE), 0.0);
        std::vector<scalar> rhoTile (size_t(GROUP_TILE), 0.0);
        std::vector<scalar> rho0Tile(size_t(GROUP_TILE), 0.0);

        for (label base = 0; base < nGroupCells; base += GROUP_TILE)
        {
            const label hi = min(base + GROUP_TILE, nGroupCells);
            const label nT = hi - base;

            // --- Tile load: one contiguous read per field ---
            for (label i = 0; i < nSpecie_; ++i)
            {
                const scalarField& Yi = Yvf_[i].oldTime();
                scalar* dst = &Ytile[size_t(i) * GROUP_TILE];
                for (label t = 0; t < nT; ++t)
                {
                    dst[t] = Yi[base + t];
                }
            }
            for (label t = 0; t < nT; ++t)
            {
                const label celli = base + t;
                Ttile  [t] = T0vf [celli];
                Ptile  [t] = p0vf [celli];
                rhoTile[t] = rho0vf[celli];
                rho0Tile[t] = rho0vf[celli];
            }

            // --- Per-cell compute from the cache-resident tile ---
            for (label t = 0; t < nT; ++t)
            {
                const label celli = base + t;

                const scalar rho  = rhoTile[t];
                const scalar rho0 = rho0Tile[t];
                scalar p = Ptile[t];
                scalar T = Ttile[t];

                for (label i = 0; i < nSpecie_; ++i)
                {
                    Y_[i] = Y0[i] = Ytile[size_t(i) * GROUP_TILE + t];
                }

                // Build this cell's zone key straight into its arena slot —
                // no std::string, no heap allocation. The running-var part is
                // a fixed 2 base-256 bytes per value (high, low), matching the
                // old listToBase256Word; the optional logP suffix matches the
                // old toBase256Word (1 byte, or 2 when posP > 222).
                char* keyPtr = &keyArena_[size_t(celli) * keyStride_];

                forAll(runSp, svi)
                {
                    const scalar value = Y_[runSp[svi]];
                    const scalar minvalue = encMinSp[svi];
                    const scalar span     = encSpanSp[svi];
                    const label pos = min(max(label(floor((value-minvalue)/span)),0),maxiRepresentation_);
                    posList[svi] = pos;
                }

                forAll(runNS, nvi)
                {
                    const scalar value = (*runNS[nvi])[celli];
                    const scalar minvalue = encMinNS[nvi];
                    const scalar span     = encSpanNS[nvi];
                    const label pos = min(max(label(floor((value-minvalue)/span)),0),maxiRepresentation_);
                    posList[nRunSp + nvi] = pos;
                }

                for (label i = 0; i < nEnc; ++i)
                {
                    const uint value = uint(posList[i]);
                    keyPtr[2*i+0] = char(value / 223 + 32);
                    keyPtr[2*i+1] = char(value % 223 + 32);
                }
                label keyLen = nEnc * 2;

                // Add logarithmic pressure suffix if enabled
                if (useLogPOnThisRun)
                {
                    scalar logP = Foam::log(max(p, 1.0));  // Avoid log(0) by clamping to 1 Pa minimum
                    label posP = min(max(label(floor((logP - pMin) / pSpan)), 0), maxiRepresentation_);
                    const uint pv = uint(posP);
                    if (pv > 222)
                    {
                        keyPtr[keyLen++] = char(pv / 223 + 32);
                        keyPtr[keyLen++] = char(pv % 223 + 32);
                    }
                    else
                    {
                        keyPtr[keyLen++] = char(pv % 223 + 32);
                    }
                }

                keyLen_[celli] = keyLen;

                // Hash the key once here and reuse it: for the round-robin
                // remainder, for the mpToCore_ find/insert below, and (stored
                // in zoneHash_) for the mpReceivedRR_ lookup in Rate_Assignment.
                uint encodeZoneIndex = encode(keyPtr, keyLen);
                zoneHash_[celli] = encodeZoneIndex;
                label remainder = encodeZoneIndex % localCores_;
                zoneRemainder_[celli] = remainder;

                // Insert / merge this cell's (Y, T, p, ...) into mpToCore_[remainder]:
                //   new key      -> insertManual + initialise
                //   existing key -> build a count=1 TCEntry on mergeScratch and
                //                   merge it into the existing block.
                {
                    MPHashTable& mpDst  = mpToCore_[remainder];
                    void* existing = mpDst.find(keyPtr, keyLen, encodeZoneIndex);
                    if (existing)
                    {
                        auto* hdr = reinterpret_cast<MPHashTable::Entry*>
                            (mergeScratch.data());
                        hdr->keyLength = keyLen;
                        hdr->blockSize = MPHashTable::headerBytes()
                                    + TCEntry::userSize(keyLen, nSpecie_, debug);
                        hdr->next      = -1;
                        hdr->pad_      = 0;
                        std::memcpy(hdr->data(), keyPtr, size_t(keyLen));
                        TCEntry incoming(mergeScratch.data(), nSpecie_);
                        incoming.initialise
                        (
                            Y0, T, p, deltaTChem_[celli], rho0, rho, 1, debug
                        );
                        TCEntry(existing, nSpecie_).mergeInPlace(incoming);
                    }
                    else
                    {
                        void* blk = mpDst.insertManual
                        (
                            keyPtr,
                            keyLen,
                            TCEntry::userSize(keyLen, nSpecie_, debug),
                            encodeZoneIndex
                        );
                        TCEntry(blk, nSpecie_).initialise
                        (
                            Y0, T, p, deltaTChem_[celli], rho0, rho, 1, debug
                        );
                    }
                }
            }
        }
    }
    else
    {
        label nZoneCells = zone_.nCells();
        scalarField noTileY(nSpecie_);
        scalar noTileT, noTileP, noTileRho;

        // initialize
        for (label zci = 0; zci < nZoneCells; zci++)
        {
            label celli = zone_.celli(zci);
            for (label i = 0; i < nSpecie_; ++i)
            {
                const scalarField& Yi = Yvf_[i].oldTime();
                noTileY[i] = Yi[celli];
            }
            noTileT = T0vf[celli];
            noTileP = p0vf[celli];
            noTileRho = rho0vf[celli];


            // encode
            forAll(runSp, svi)
            {
                const scalar value = noTileY[runSp[svi]];
                const scalar minvalue = encMinSp[svi];
                const scalar span     = encSpanSp[svi];
                const label pos = min(max(label(floor((value-minvalue)/span)),0),maxiRepresentation_);
                posList[svi] = pos;
            }

            forAll(runNS, nvi)
            {
                const scalar value = (*runNS[nvi])[celli];
                const scalar minvalue = encMinNS[nvi];
                const scalar span     = encSpanNS[nvi];
                const label pos = min(max(label(floor((value-minvalue)/span)),0),maxiRepresentation_);
                posList[nRunSp + nvi] = pos;
            }

            char* keyPtr = &keyArena_[size_t(celli) * keyStride_];
            for (label i = 0; i < nEnc; ++i)
            {
                const uint value = uint(posList[i]);
                keyPtr[2*i+0] = char(value / 223 + 32);
                keyPtr[2*i+1] = char(value % 223 + 32);
            }
            label keyLen = nEnc * 2;

            // Add logarithmic pressure suffix if enabled
            if (useLogPOnThisRun)
            {
                scalar logP = Foam::log(max(noTileP, 1.0));  // Avoid log(0) by clamping to 1 Pa minimum
                label posP = min(max(label(floor((logP - pMin) / pSpan)), 0), maxiRepresentation_);
                const uint pv = uint(posP);
                if (pv > 222)
                {
                    keyPtr[keyLen++] = char(pv / 223 + 32);
                    keyPtr[keyLen++] = char(pv % 223 + 32);
                }
                else
                {
                    keyPtr[keyLen++] = char(pv % 223 + 32);
                }
            }
            keyLen_[celli] = keyLen;

            // Hash the key once here and reuse it: for the round-robin
            // remainder, for the mpToCore_ find/insert below, and (stored
            // in zoneHash_) for the mpReceivedRR_ lookup in Rate_Assignment.
            uint encodeZoneIndex = encode(keyPtr, keyLen);
            zoneHash_[celli] = encodeZoneIndex;
            label remainder = encodeZoneIndex % localCores_;
            zoneRemainder_[celli] = remainder;

            // Insert / merge this cell's (Y, T, p, ...) into mpToCore_[remainder]:
            //   new key      -> insertManual + initialise
            //   existing key -> build a count=1 TCEntry on mergeScratch and
            //                   merge it into the existing block.
            {
                MPHashTable& mpDst  = mpToCore_[remainder];
                void* existing = mpDst.find(keyPtr, keyLen, encodeZoneIndex);
                if (existing)
                {
                    auto* hdr = reinterpret_cast<MPHashTable::Entry*>
                        (mergeScratch.data());
                    hdr->keyLength = keyLen;
                    hdr->blockSize = MPHashTable::headerBytes()
                                + TCEntry::userSize(keyLen, nSpecie_, debug);
                    hdr->next      = -1;
                    hdr->pad_      = 0;
                    std::memcpy(hdr->data(), keyPtr, size_t(keyLen));
                    TCEntry incoming(mergeScratch.data(), nSpecie_);
                    incoming.initialise
                    (
                        noTileY, noTileT, noTileP, deltaTChem_[celli], noTileRho, noTileRho, 1, debug
                    );
                    TCEntry(existing, nSpecie_).mergeInPlace(incoming);
                }
                else
                {
                    void* blk = mpDst.insertManual
                    (
                        keyPtr,
                        keyLen,
                        TCEntry::userSize(keyLen, nSpecie_, debug),
                        encodeZoneIndex
                    );
                    TCEntry(blk, nSpecie_).initialise
                    (
                        noTileY, noTileT, noTileP, deltaTChem_[celli], noTileRho, noTileRho, 1, debug
                    );
                }
            }
            
        }


       

            
        }
    



    CCM_TIMING_END(Cell_Grouping, debugTime_, stepTimer_, *this);

    CCM_TIMING_START(Load_Redistribution, debugTime_, stepTimer_);
    // Load balancing redistribution
    if (Pstream::parRun() && (redistributeEvery_) > 0 &&
        (currentStep_ % redistributeEvery_ == 0) && !redistributedThisStep_)
    {
        redistributeLoads();
        redistributedThisStep_ = true;
    }
    CCM_TIMING_END(Load_Redistribution, debugTime_, stepTimer_, *this);

    // Step 3: Round-Robin Distribution
    CCM_TIMING_START(Distribution, debugTime_, stepTimer_);

    // Manual MPI_Alltoallv: ship mpToCore_[p].pool to mpFromCore_[p].pool.
    // The 2x MPI_Alltoall (bytes + entries) + 1x MPI_Alltoallv pipeline is
    // packed into mpAlltoallvBytes; this block just gathers inputs from
    // mpToCore_ and feeds the recv pointers into mpFromCore_->decodeFrom.
    {
        const int numPeers = mpToCore_.size();
        forAll(mpFromCore_, peer) mpFromCore_[peer].clear();

        // Per-peer send descriptors: pool pointer, byte count, entry count.
        std::vector<const void*> sendPoolPtrs   (numPeers);
        std::vector<int>         sendByteCounts (numPeers);
        std::vector<int>         sendEntryCounts(numPeers);
        for (int peer = 0; peer < numPeers; ++peer)
        {
            sendPoolPtrs   [peer] = mpToCore_[peer].poolStart();
            sendByteCounts [peer] = int(mpToCore_[peer].poolUsed());
            sendEntryCounts[peer] = int(mpToCore_[peer].nEntries());
        }

        // Per-peer recv descriptors: filled by mpAlltoallvBytes.
        std::vector<const char*> recvPoolPtrs   (numPeers);
        std::vector<int>         recvByteCounts (numPeers);
        std::vector<int>         recvEntryCounts(numPeers);
        mpAlltoallvBytes
        (
            numPeers,
            sendPoolPtrs.data(),
            sendByteCounts.data(),
            sendEntryCounts.data(),
            recvPoolPtrs.data(),
            recvByteCounts.data(),
            recvEntryCounts.data()
        );

        // Decode each peer's received bytes into mpFromCore_[peer].
        for (int peer = 0; peer < numPeers; ++peer)
        {
            if (recvEntryCounts[peer] > 0)
            {
                mpFromCore_[peer].decodeFrom
                (
                    recvPoolPtrs[peer],
                    recvEntryCounts[peer]
                );
            }
        }
    }

    CCM_TIMING_END(Distribution, debugTime_, stepTimer_, *this);

    CCM_TIMING_START(Merging, debugTime_, stepTimer_);

    // Merge mpFromCore_[p] into mpGathered_ peer-by-peer.
    // Walk each peer's pool block-by-block (via the blockSize header in
    // every Entry, the same way decodeFrom does). On key collision the
    // existing block is updated by TCEntry::mergeInPlace; otherwise a
    // fresh block is allocated in mpGathered_'s pool and the source
    // userSize region is memcpy'd verbatim (layout-compatible).
    {
        const label N = mpFromCore_.size();
        mpGathered_->clear();
        for (label p = 0; p < N; ++p)
        {
            const char* src = mpFromCore_[p].poolStart();
            const label nE  = mpFromCore_[p].nEntries();
            for (label i = 0; i < nE; ++i)
            {
                const auto* hdr =
                    reinterpret_cast<const MPHashTable::Entry*>(src);
                const label bs = hdr->blockSize;

                TCEntry srcView
                (
                    const_cast<void*>(static_cast<const void*>(src)),
                    nSpecie_
                );
                const string key    = srcView.getKey();
                const label  keyLen = label(key.size());

                void* existing = mpGathered_->find(key);
                if (existing)
                {
                    TCEntry(existing, nSpecie_).mergeInPlace(srcView);
                }
                else
                {
                    void* blk = mpGathered_->insertManual
                    (
                        key,
                        TCEntry::userSize(keyLen, nSpecie_, debug)
                    );
                    auto* dstHdr =
                        reinterpret_cast<MPHashTable::Entry*>(blk);
                    std::memcpy
                    (
                        dstHdr->data(),
                        hdr->data(),
                        size_t(TCEntry::userSize(keyLen, nSpecie_, debug))
                    );
                }
                src += bs;
            }
        }
    }

    CCM_TIMING_END(Merging, debugTime_, stepTimer_, *this);

    // Step 3.3: Post-merge operations (statistics and debug)
    CCM_TIMING_START(Post_merge_Total, debugTime_, stepTimer_);

    label localPhaseSpaceSize = mpGathered_->nEntries();
    label phaseSpaceSize      = localPhaseSpaceSize;
    reduce(phaseSpaceSize, sumOp<label>());
    
    label meanPhaseSpaceSize = phaseSpaceSize/Pstream::nProcs();
    scalar unbalanceRatio = mag(scalar(localPhaseSpaceSize)/(max(1.0, meanPhaseSpaceSize))-1.0);
    reduce(unbalanceRatio, maxOp<scalar>());

    scalar totalCells = this->mesh().nCells();
    reduce(totalCells, sumOp<scalar>());


    scalar accRatio = totalCells/scalar(phaseSpaceSize);

    Info << nl << "CCM status" << endl;   
    Info << "====================================" << nl;
    Info << setw(20) << "acceleration ratio"
     << setw(20) << "phase space size"
     << setw(20) << "Unbalance ratio" << endl;

    Info << setw(20) << fixed << setprecision(2) << accRatio
    << setw(20) << fixed << setprecision(2) << phaseSpaceSize
    << setw(20) << fixed << setprecision(2) << unbalanceRatio << endl;
    Info << "====================================" << nl;
    CCM_TIMING_END(Post_merge_Total, debugTime_, stepTimer_, *this);

}

template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::updateReactionRate()
{

    // Step 4: Chemistry Computation
    CCM_TIMING_START(Chemistry_Computation, debugTime_, stepTimer_);
    // Calculate reaction rates for the reactionEntries assigned to this processor 
    // according to the remainders
    // Store the results in the reaction rate entries(rrEntries)


    // Walk mpGathered_'s pool block-by-block, compute the rate for each
    // zone via the configured GetRRFunc_, and store the result as a
    // rateEntry in mpLocalRR_.
    mpLocalRR_->clear();
    {
        const scalar& deltaT = mesh().time().deltaTValue();
        const char* src = mpGathered_->poolStart();
        const label nE  = mpGathered_->nEntries();
        for (label i = 0; i < nE; ++i)
        {
            const auto* hdr =
                reinterpret_cast<const MPHashTable::Entry*>(src);
            const label bs = hdr->blockSize;

            TCEntry v
            (
                const_cast<void*>(static_cast<const void*>(src)),
                nSpecie_
            );
            const string key  = v.getKey();
            scalarField Y     = v.getY();
            scalar      T     = v.getT();
            scalar      p     = v.getP();
            scalar      deltaTChem = v.getDtChem();
            const scalar rho  = v.getRho();
            const scalar rho0 = v.getRho0();

            scalarField cellRR = (this->*GetRRFunc_)
            (
                Y, T, p, deltaT, deltaTChem, rho, rho0
            );

            void* blk = mpLocalRR_->insertManual
            (
                key, rateEntry::userSize(label(key.size()), nSpecie_)
            );
            rateEntry(blk, nSpecie_).initialise(cellRR, deltaTChem);

            src += bs;
        }
    }

    CCM_TIMING_END(Chemistry_Computation, debugTime_, stepTimer_, *this);
}

template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::distributeReactionRate()
{
    // Step 6: pack per-peer return data. For each peer p, walk the keys
    // it sent us (carried in mpFromCore_[p]) and look up the rate we
    // computed for that zone in mpLocalRR_, then stuff it into
    // mpReturnToCore_[p] as a rateEntry.
    CCM_TIMING_START(Hash_Table_Construction, debugTime_, stepTimer_);
    {
        forAll(mpReturnToCore_, p) mpReturnToCore_[p].clear();
        const label N = mpFromCore_.size();
        for (label p = 0; p < N; ++p)
        {
            const char* src = mpFromCore_[p].poolStart();
            const label nE  = mpFromCore_[p].nEntries();
            for (label i = 0; i < nE; ++i)
            {
                const auto* hdr =
                    reinterpret_cast<const MPHashTable::Entry*>(src);
                const label bs     = hdr->blockSize;
                const label keyLen = hdr->keyLength;
                const string key
                (
                    static_cast<const char*>(hdr->data()),
                    size_t(keyLen)
                );

                void* rrBlk = mpLocalRR_->find(key);
                if (rrBlk)
                {
                    rateEntry rv(rrBlk, nSpecie_);
                    void* dst = mpReturnToCore_[p].insertManual
                    (
                        key, rateEntry::userSize(keyLen, nSpecie_)
                    );
                    rateEntry(dst, nSpecie_).initialise
                    (
                        rv.getRR(), rv.getDeltaTChem()
                    );
                }
                src += bs;
            }
        }
    }
    CCM_TIMING_END(Hash_Table_Construction, debugTime_, stepTimer_, *this);

    // Step 7: ship mpReturnToCore_[p].pool back to mpReceivedRR_[p].pool via
    // the same 2x MPI_Alltoall + 1x MPI_Alltoallv pipeline as Step 3,
    // mirroring the distribution but carrying rateEntry blocks instead of
    // TCEntry blocks.
    CCM_TIMING_START(Communication_Execution, debugTime_, stepTimer_);
    {
        const label numPeers = mpReturnToCore_.size();
        forAll(mpReceivedRR_, peer) mpReceivedRR_[peer].clear();

        // Per-peer send descriptors: pool pointer, byte count, entry count.
        std::vector<const void*> sendPoolPtrs   (numPeers);
        std::vector<int>         sendByteCounts (numPeers);
        std::vector<int>         sendEntryCounts(numPeers);
        for (label peer = 0; peer < numPeers; ++peer)
        {
            sendPoolPtrs   [peer] = mpReturnToCore_[peer].poolStart();
            sendByteCounts [peer] = int(mpReturnToCore_[peer].poolUsed());
            sendEntryCounts[peer] = int(mpReturnToCore_[peer].nEntries());
        }

        // Per-peer recv descriptors: filled by mpAlltoallvBytes.
        std::vector<const char*> recvPoolPtrs   (numPeers);
        std::vector<int>         recvByteCounts (numPeers);
        std::vector<int>         recvEntryCounts(numPeers);
        mpAlltoallvBytes
        (
            int(numPeers),
            sendPoolPtrs.data(),
            sendByteCounts.data(),
            sendEntryCounts.data(),
            recvPoolPtrs.data(),
            recvByteCounts.data(),
            recvEntryCounts.data()
        );

        // Decode each peer's received bytes into mpReceivedRR_[peer].
        for (label peer = 0; peer < numPeers; ++peer)
        {
            if (recvEntryCounts[peer] > 0)
            {
                mpReceivedRR_[peer].decodeFrom
                (
                    recvPoolPtrs[peer],
                    recvEntryCounts[peer]
                );
            }
        }
    }
    CCM_TIMING_END(Communication_Execution, debugTime_, stepTimer_, *this);

    // Step 8: write RR_/deltaTChem_ for every cell by looking up its zone
    // in the corresponding mpReceivedRR_[r] table.
    CCM_TIMING_START(Rate_Assignment, debugTime_, stepTimer_);
    {
        // Loop tiling, mirror of the grouping loop. Phase 1 finds each
        // cell's rate block (reusing the FNV hash stashed in zoneHash_ so
        // the lookup skips re-hashing the key) and copies its RR into a
        // small block-major tile buffer with one contiguous read per block.
        // Phase 2 scatters the tile out species-major so each RR_[i] field
        // is written as one contiguous stream instead of nSpecie_ strided
        // writes per cell.
        if (zone_.all())
        {
            const label nRateCells = zoneRemainder_.size();
            const label RATE_TILE   = 16;

            // block-major: rrTile[t*nSpecie_ + i]
            std::vector<scalar> rrTile(size_t(nSpecie_) * RATE_TILE, 0.0);
            std::vector<scalar> dtTile(size_t(RATE_TILE), 0.0);

            for (label base = 0; base < nRateCells; base += RATE_TILE)
            {
                const label hi = min(base + RATE_TILE, nRateCells);
                const label nT = hi - base;

                // --- Phase 1: find + gather into the tile ---
                for (label t = 0; t < nT; ++t)
                {
                    const label celli = base + t;
                    const label r = zoneRemainder_[celli];
                    const char* keyPtr = &keyArena_[size_t(celli) * keyStride_];
                    void* blk = mpReceivedRR_[r].find
                    (
                        keyPtr, keyLen_[celli], zoneHash_[celli]
                    );
                    rateEntry v(blk, nSpecie_);
                    const scalar* rr = v.rawRR();
                    scalar* dst = &rrTile[size_t(t) * nSpecie_];
                    for (label i = 0; i < nSpecie_; ++i)
                    {
                        dst[i] = rr[i];
                    }
                    dtTile[t] = min(v.getDeltaTChem(), deltaTChemMax_);
                }

                // --- Phase 2: scatter species-major (contiguous RR_ writes) ---
                for (label i = 0; i < nSpecie_; ++i)
                {
                    scalarField& RRi = RR_[i];
                    for (label t = 0; t < nT; ++t)
                    {
                        RRi[base + t] = rrTile[size_t(t) * nSpecie_ + i];
                    }
                }
                for (label t = 0; t < nT; ++t)
                {
                    deltaTChem_[base + t] = dtTile[t];
                }
            }
        }
        else
        {
            deltaTChem_ = dimensionedScalar(dimTime, GREAT);
            // set reaction rates to Zero 
            forAll(RR_, i)
            {
                RR_[i] = Zero;
            }
            
            // attribute zone reaction rates
            scalarField tempRR = scalarField(nSpecie_);
            const label nZoneCells = zone_.nCells();
            
            for (label zci = 0; zci < nZoneCells; zci++)
            {
                const label celli = zone_.celli(zci);
                const label r = zoneRemainder_[celli];
                char* keyPtr = &keyArena_[size_t(celli) * keyStride_];
                void* blk = mpReceivedRR_[r].find
                    (
                        keyPtr, keyLen_[celli], zoneHash_[celli]
                    );
                rateEntry v(blk, nSpecie_);

                // get RR and deltaTChem
                const scalar* rr = v.rawRR();
                std::memcpy(tempRR.data(), rr, size_t(nSpecie_) * sizeof(scalar));
                scalar dtChemTemp = min(v.getDeltaTChem(), deltaTChemMax_);
                for (label i = 0; i < nSpecie_; ++i)
                {
                    RR_[i][celli] = tempRR[i];
                }
                deltaTChem_[celli] = dtChemTemp;
            }
        }
    }
    CCM_TIMING_END(Rate_Assignment, debugTime_, stepTimer_, *this);
}

template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::updateCCM4MeshChange()
{
    if (mesh().dynamic())
    {
        label nCells = mesh().nCells();
        if (keyLen_.size() != nCells)
        {
            // keyArena_ is (re)sized to nCells*keyStride_ at the start of
            // the grouping loop, where keyStride_ for the step is known.
            keyLen_.setSize(nCells);
            zoneRemainder_.setSize(nCells);
            zoneHash_.setSize(nCells);
        }
    }
}


template<class ThermoType>
Foam::scalarField Foam::chemistryModels::CCM<ThermoType>::getRRGivenYTP_Basic
(
    scalarField& Y,
    scalar& T,
    scalar& p,
    const scalar& deltaT,
    scalar& deltaTChem,
    const scalar& rho,
    const scalar& rho0
)
{
    return getRRGivenYTP(Y, T, p, deltaT, deltaTChem, rho, rho0);
}

template<class ThermoType>
Foam::scalarField Foam::chemistryModels::CCM<ThermoType>::getRRGivenYTP_Optimized
(
    scalarField& Y,
    scalar& T,
    scalar& p,
    const scalar& deltaT,
    scalar& deltaTChem,
    const scalar& rho,
    const scalar& rho0
)
{
    return fastChemistryPtr_->getRRGivenYTP(Y, T, p, deltaT, deltaTChem, rho, rho0);
}

template<class ThermoType>
void Foam::chemistryModels::CCM<ThermoType>::redistributeLoads()
{

    // calculate curernt load
    label& localLoad = loads_[Pstream::myProcNo()];
    localLoad = 0;
    forAll(mpToCore_, i)
    {
        localLoad += mpToCore_[i].nEntries();
    }

    // Gather all loads
    Pstream::gatherList(loads_);
    Pstream::scatterList(loads_);

    Info << "loads_" << loads_ << endl;
    
    // sort loads
    labelList index(loads_.size());
    forAll(index, i){index[i] = i;}
    std::sort(index.begin(), index.end(), [this](label a, label b)
    {
        bool returnValue = (loads_[a] == loads_[b]) ? (a < b) : (loads_[a] < loads_[b]);
        return returnValue;
    });

    forAll(index,i)
    {
        Info << "rank " << index[i] << " load: " << loads_[index[i]] << endl;
    }

    label myLoc = 0;
    forAll(index, i)
    {
        if (index[i] == Pstream::myProcNo())
        {
            myLoc = i;
            break;
        }
    }

    label nGroups = Pstream::nProcs()/localCores_;
    label  color = myLoc%nGroups;
    label  key   = myLoc/nGroups;

    MPI_Comm_free(&chemComm_);
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &chemComm_);

}

// ************************************************************************* //
