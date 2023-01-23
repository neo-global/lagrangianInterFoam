#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


volScalarField mu
(
    IOobject
    (
        "mu",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mixture.mu()
);


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
 surfaceScalarField alphaPhi1Un
 (
     IOobject
     (
         "alphaPhi1Un",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", phi.dimensions(), Zero)
 );

 surfaceScalarField alphaPhi2Un
 (
     IOobject
     (
         "alphaPhi2Un",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedScalar("zero", phi.dimensions(), Zero)
 );

 volVectorField Um
 (
     IOobject
     (
         "Um",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     U
 );

 volScalarField alpham
 (
     IOobject
     (
         "alpham",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     alpha1
 );


#include "createMRF.H"

//Initialize basicKinematicCollidingCloud

const word kinematicCloudName
(
    args.optionLookupOrDefault<word>("cloudName", "kinematicCloud")
);
Info<< "Constructing kinematic cloud "<< kinematicCloudName <<endl;

const volScalarField& gasPhase = mesh.lookupObject<volScalarField>("alpha.air");

volScalarField phaseSu
(
    IOobject
    (
        "phaseSu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0) 
);

volVectorField slipVel
(
    IOobject
    (
	"slipVel",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, Zero)
);

volScalarField Cd 
(
    IOobject
    (
        "Cd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) 
);

volScalarField nPInCell 
(
    IOobject
    (
        "nPInCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) 
);

volScalarField surfaceCell 
(
    IOobject
    (
        "surfaceCell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) 
);
volScalarField cellsWithConv 
(
    IOobject
    (
        "cellsWithConv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) 
);

volScalarField Db 
(
    IOobject
    (
        "Db",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0) 
);


basicKinematicCloud kinematicCloud
(
    kinematicCloudName,
    rho,
    U,
    mu,
    gasPhase,
    phaseSu,
    g
);

volScalarField voidage 
(
    IOobject
    (
        "voidage",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) 
);

#include "createFvOptions.H"
#include "createIsoAdvection.H"
#include "createCustomDict.H"
#include "createResIndex.H"
