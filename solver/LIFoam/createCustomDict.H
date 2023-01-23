IOdictionary customDict
(
    IOobject
    (
	"values",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);
