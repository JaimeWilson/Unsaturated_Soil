/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createConstants.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
         #include "CalPvandHw.H"	
         #include "CalHydraulicconductivity.H"
        //#include "lamEqn.H"
        //#include "KEqn.H"
	//#include "CtEqn.H"
	//#include "CEqn.H"

        volVectorField  A(fvc::grad(k));

        while (simple.correctNonOrthogonal())
        {

	// Diferential Equation Temperature

            fvScalarMatrix TEqn

            (
                fvm::ddt(Ct, T) - fvm::laplacian(lam, T) - fvc::laplacian(lam2, Theta) //- Foam::fvc::grad(k) // - fvm::laplacian(lam2, Theta)
             //==
                //fvOptions(T)*dimensionedScalar("scl1",dimensionSet(1,-1,-2,-1,0,0,0),1.)
            );

            //fvOptions.constrain(TEqn);

            TEqn.solve();

            //fvOptions.correct(T);

	// Diferential Equation Theta

            fvScalarMatrix ThetaEqn

            (
                fvm::ddt(C,Theta) - fvm::laplacian(K, Theta) - fvc::laplacian(K2, T) //- Foam::fvc::grad(k)  //fvm::laplacian(K2, T) 

             //==
                //fvOptions(Theta)*dimensionedScalar("scl2",dimensionSet(0,-1,0,0,0,0,0),1.)
            );

            //fvOptions.constrain(ThetaEqn);

            ThetaEqn.solve();

            //fvOptions.correct(Theta);

        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //
