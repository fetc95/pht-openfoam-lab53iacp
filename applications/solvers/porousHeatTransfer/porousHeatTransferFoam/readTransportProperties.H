IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Porosity
dimensionedScalar a
(
    "a",
    dimless,
    transportProperties
);

// Porous resistance properties

dimensionedScalar chi
(
    "chi",
    dimless,
    transportProperties
);

dimensionedScalar KD
(
    "KD",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    transportProperties
);

dimensionedScalar F
(
    "F",
    dimensionSet(0, -1, 0, 0, 0, 0, 0),
    transportProperties
);

dimensionedScalar addedMassDivA2 = (1 + (1 - a) * chi) / pow(a, 2);

// Heat transfer coefficient

dimensionedScalar HTcoef
(
    "HTcoef",
    dimensionSet(1, -1, -3, -1, 0, 0, 0),
    transportProperties
);

