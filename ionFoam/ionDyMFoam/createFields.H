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
Info<< "############################################################################" << endl;
Info<< "###############   SETTING UP RAD.-INDUCED VISCOSITY-MODEL   ################" << endl;
Info<< "################################# BY #######################################" << endl;
Info<< "##################### G. HOBLER & J. KIESENHOFER ###########################" << endl;
Info<< "############################################################################\n" << endl;

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const scalar& pi = Foam::constant::mathematical::pi;

volVectorField C = mesh.C();           // Cell center coordinates
//surfaceVectorField Cf = mesh.Cf();   // Face center coordinates
pointField points = mesh.points();     // Node coordinates

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


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

typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
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

volScalarField voidfraction  // volume fraction of fluid, the rest is solid phase within one cell
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
);

volScalarField xDir  
    (
        IOobject
        (
            "xDir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C.component(0)
);

volScalarField yDir  
    (
        IOobject
        (
            "yDir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C.component(1)
);

volScalarField zDir  
    (
        IOobject
        (
            "zDir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        C.component(2)
);

dimensionedScalar dim_meter
(
        "dim_meter",
        dimensionSet(0,1,0,0,0,0,0), // [kg m s K mol A cd]
        scalar(1.0)
);

dimensionedScalar dim_time
(
        "dim_time",
        dimensionSet(0,0,1,0,0,0,0), // [kg m s K mol A cd]
        scalar(1.0)
);

dimensionedScalar dim_density
(
        "dim_density",
        dimensionSet(1,-3,0,0,0,0,0), // [kg m s K mol A cd]
        scalar(1.0)
);

dimensionedScalar dim_mass
(
        "dim_mass",
        dimensionSet(1,0,0,0,0,0,0), // [kg m s K mol A cd]
        scalar(1.0)
);

dimensionedScalar dim_vol
(
        "dim_vol",
        dimensionSet(0,3,0,0,0,0,0), // [kg m s K mol A cd]
        scalar(1.0)
);

volScalarField phaseTest  
    (
        IOobject
        (
            "phaseTest",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
);

volScalarField phase 
    (
        IOobject
        (
            "phase",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
);

volScalarField solid_VOF
    (
        IOobject
        (
            "solid_VOF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
);

volScalarField contErr
    (
        IOobject
        (
            "contErr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
);

volVectorField nuEff
    (
        IOobject
        (
            "nuEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C()/mag(mesh.C())
);

volScalarField MC_source
    (
        IOobject
        (
            "MC_source",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
);

volScalarField MC
    (
        IOobject
        (
            "MC_source",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1
);

//label maxindex=0;

forAll(mesh.C(),i)
{
    //MC_source[i] = MC_source[i]+rand();
    MC_source[i] = MC_source[i]+rand();
    //if(MC_source[i]>-1.0)
    //{
    //if(mesh.C()[i].component(0) >  mesh.C()[maxindex].component(0)) 
    //      maxindex=i;
    //}
    Info<< "MC_source = " << MC_source[i] << endl;
}
MC_source = MC_source/max(MC_source)/10;
MC = MC_source;

volVectorField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField S
(
    IOobject
    (
        "S",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh.C()/mag(mesh.C())
);

dimensionedScalar xMax = max(points.component(0));
dimensionedScalar yMax = max(points.component(1));
dimensionedScalar zMax = max(points.component(2));

volScalarField xDir0 = xDir-xMax.value()/2*dim_meter;
volScalarField yDir0 = yDir-yMax.value()/2*dim_meter;
volScalarField zDir0 = zDir-zMax.value()/4*dim_meter;

volScalarField radiusC = sqrt(sqr(8*(xDir0))+sqr(zDir));
volScalarField radius0 = sqrt(sqr((xDir0))+sqr(zDir0));
phaseTest = cos(radiusC/dim_meter*300); //+cos(radiusC/dim_meter*300);

dimensionedScalar beamWaist = 2200*1e-9*dim_meter;
dimensionedScalar wavelength =  8*1e-9*dim_meter;

dimensionedScalar Rayleigh_length = pi*sqr(beamWaist)/wavelength;
dimensionedScalar wavevector = 2*pi/wavelength;
volScalarField Wz = beamWaist*sqrt(1+sqr(zDir0/Rayleigh_length));
volScalarField Rz = zDir0*(1+sqr(Rayleigh_length/zDir0));

volScalarField stokesDistr = exp(-2*(sqr(radius0/dim_meter)/sqr(beamWaist/dim_meter)));
stokesDistr = stokesDistr/max(stokesDistr);
MC_source = MC_source+stokesDistr;

phase = wavevector*(sqr(xDir0)/(2*Rz)+zDir0);
S = -fvc::grad(phase)*dim_meter;
S = S/mag(S);
volVectorField Z = -fvc::grad(MC_source)*dim_meter;
Z = S/mag(Z);

//dimensionedScalar speedOfIons  =      6000*dim_meter/dim_time;
//dimensionedScalar speedOfLight = 299792458*dim_meter/dim_time;
//dimensionedScalar kappa = speedOfIons/speedOfLight;
dimensionedScalar kappa = 1.0;

surfaceScalarField ionf ( "ionf" , fvc::interpolate(S) & mesh.Sf() * dim_meter / dim_time * kappa );
surfaceScalarField stokesf ( "stokesf" , fvc::interpolate(Z) & mesh.Sf() * dim_meter / dim_time * kappa );

// Mass flux pushed by beam
surfaceScalarField rhoIonf
(
    IOobject
    (
        "rhoIonf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*ionf
);

surfaceScalarField rhoStokesf
(
    IOobject
    (
        "rhoStokesf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*stokesf
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

#include "createMRF.H"
#include "createFvOptions.H"
