#include "../SETTINGS.h"

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field Umean\n" << endl;
    volVectorField Umean
    (
        IOobject
        (
            "Umean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading field nut\n" << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    dimensionedScalar N2
    (
        "N2",
        dimensionSet(0, 0, -2, 0, 0, 0, 0),
        N*N
    );
    dimensionedScalar VISC
    (
        "VISC",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
        nu
    );

        // DISSIPATION
        volTensorField      gradu( fvc::grad(U - Umean) );

        Info<< "Computing field DIS" << endl;
        volScalarField DIS
        (
            IOobject
            (
                "DIS",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (-(VISC+nut))*(gradu && gradu)  //Nomore sqrt(N2*N2*N2)
        );
