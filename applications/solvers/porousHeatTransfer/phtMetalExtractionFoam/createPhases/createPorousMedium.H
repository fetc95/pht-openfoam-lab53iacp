Info<< "Reading porosity\n" << endl;
volScalarField a_g
(
    IOobject
    (
        "Porosity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField a_g0
(
    IOobject
    (
        "Porosity0",
        runTime.timeName(0),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
//    a_g
	mesh
);

IOdictionary interphaseProperties
(
    IOobject
    (
        "interphaseProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Porous frame resistance

dimensionedScalar chi
(
    dimless,
    interphaseProperties.getOrDefault<scalar>("virtualMassCoeff", 0.5)
);

dimensionedScalar prtclDmtr
(
    "d",
    dimLength,
    interphaseProperties
);

// Heat transfer

dimensionedScalar HTcoef0
(
    "Kh0",
    dimensionSet(1, -1, -3, -1, 0, 0, 0),
    interphaseProperties
);

dimensionedScalar HTorder
(
    dimless,
    interphaseProperties.getOrDefault<scalar>("KhOrder", 1)
);

volScalarField HTcoef
(
    IOobject
    (
        "HTcoef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
	HTcoef0
);
