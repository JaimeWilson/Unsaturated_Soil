volVectorField gradk(fvc::grad(k));

        volScalarField gradkx
        (
            IOobject
            (
                "gradkx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradk.component(vector::X)
        );

        volScalarField gradky
        (
            IOobject
            (
                "gradky",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradk.component(vector::Y)
        );

        volScalarField gradkz
        (
            IOobject
            (
                "gradkz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradk.component(vector::Z)
        );
