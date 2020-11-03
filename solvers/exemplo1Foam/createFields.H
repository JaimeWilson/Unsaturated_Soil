Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Theta\n" << endl;

volScalarField Theta
(
    IOobject
    (
        "Theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field k\n" << endl;

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,dimensionedScalar("k", dimensionSet(0,1,1,0,0,0,0), 0.)
);

        //IOobject::READ_IF_PRESENT,
        //IOobject::AUTO_WRITE

//Info<< "Reading field lam\n" << endl;

//volScalarField lam
//(
    //IOobject
    //(
        //"lam",
        //runTime.timeName(),
        //mesh,
        //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
    //),
    //mesh, dimensionedScalar("lam", dimensionSet(1,1,-3,-1,0,0,0), 4e-01)
//);

//Info<< "Reading field K\n" << endl;

//volScalarField K
//(
    //IOobject
    //(
        //"K",
        //runTime.timeName(),
        //mesh,
        //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
    //),
    //mesh, dimensionedScalar("K", dimensionSet(0,1,-1,0,0,0,0), 4e-01)
//);

//Info<< "Reading heat capacity field Ct\n" << endl; // ML-1T−2Θ−1

//volScalarField Ct
//(
    //IOobject
    //(
        //"Ct",
        //runTime.timeName(),
        //mesh,
        //IOobject::MUST_READ,
        //IOobject::AUTO_WRITE
    //),
    //mesh, dimensionedScalar("Ct", dimensionSet(1,-1,-2,-1,0,0,0), 1.)
//);

Info<< "Reading hidraulic capacity field S\n" << endl; //M1 L-3

volScalarField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, dimensionedScalar("S", dimensionSet(0,0,0,0,0,0,0), 1.)
);

Info<< "Reading hidraulic capacity field Uw\n" << endl; //M1 L-3

volScalarField Uw
(
    IOobject
    (
        "Uw",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh, dimensionedScalar("Uw", dimensionSet(0,0,0,0,0,0,0), 1.)
);


#include "createFvOptions.H"
