/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "basicFastChemistryModel.H"
// #include "basicThermo.H"
#include "compileTemplate.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicFastChemistryModel> Foam::basicFastChemistryModel::New
(
    const fvMesh& mesh
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            "chemistryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );


    const Foam::dictionary& solverDict = chemistryDict.typeDict("ode");
    word solverName = solverDict.lookup("solver");
    if (solverName == "seulex" || solverName == "OptSeulex")
    {
        solverName = "OptSeulex";
    }
    else if (solverName == "rodas34" || solverName == "OptRodas34")
    {
        solverName = "OptRodas34";
    }
    else if (solverName == "Rosenbrock34" || solverName == "OptRosenbrock34")
    {
        solverName = "OptRosenbrock34";
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported Fast ODE solver name: " << solverName << endl << nl
            << "Available choices: seulx, Rodas34, Rosenbrock34" << exit(FatalError);
    }
    
    // initialize method name
    const word methodName("FastChemistryModel");



    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);

    Info<< "Selecting chemistry solver " << chemistryTypeDictNew << endl;

    // FastChemistryModel is no longer templated, so we only need solver<method>
    const word chemSolverNameName =
        solverName + '<' + methodName + '>';



    typename meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(chemSolverNameName);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        // Info<< "Looking for: " << chemSolverNameName << endl;
        // Info<< "Available solvers in table: " << meshConstructorTablePtr_->sortedToc() << endl;
        FatalErrorInFunction
            << "Unknown chemistry solver type "
            << chemSolverNameName << nl << nl
            << "Valid chemistry solver types are:" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<basicFastChemistryModel>(cstrIter()(mesh));
}


// ************************************************************************* //
