IOdictionary turbDict
(
    IOobject
    (
	"turbulenceProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

word turbType (turbDict.lookup("simulationType"));

volScalarField kSGS
(
    IOobject
    (
	"kSGS",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kSGS",sqr(dimVelocity),0.0)
);

volScalarField uSGS
(
    IOobject
    (
	"uSGS",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("uSGS",dimVelocity,0.0)
);

