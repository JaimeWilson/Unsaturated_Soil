// This reads a dictionary file.
Info << "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
     IOobject
     (
	"transportProperties", // name of the dictionary
	runTime.constant(), // location in the case - this one is in constant
	mesh, // needs the mesh object reference to do some voodoo - unimportant now
	IOobject::MUST_READ_IF_MODIFIED, // the file will be re-read if it gets modified during time stepping
	IOobject::NO_WRITE // read-only
     )
);


Info<< "Reading Ct\n" << endl; 

dimensionedScalar Ct("Ct", dimensionSet(1,-1,-2,-1,0,0,0), transportProperties);

Info<< "Reading lam\n" << endl; 

dimensionedScalar lam("lam", dimensionSet(1,1,-3,-1,0,0,0), transportProperties);

Info<< "Reading lam2\n" << endl; 

dimensionedScalar lam2("lam2", dimensionSet(1,1,-3,0,0,0,0), transportProperties);

Info<< "Reading K2\n" << endl; 

dimensionedScalar K2("K2", dimensionSet(0,1,-1,-1,0,0,0), transportProperties);

Info<< "Reading K\n" << endl; 

dimensionedScalar K("K", dimensionSet(0,1,-1,0,0,0,0), transportProperties);

Info<< "Reading C\n" << endl; 

dimensionedScalar C("C", dimensionSet(0,-1,0,0,0,0,0), transportProperties);

//Info<< "Reading C2\n" << endl; 

//dimensionedScalar C2("C2", dimensionSet(1,-1,-2,0,0,0,0), transportProperties);

//Info<< "Reading Ct2\n" << endl; 

//dimensionedScalar Ct2("Ct2", dimensionSet(0,-1,0,-1,0,0,0), transportProperties);

//Teor volumetrico
const scalar Vws(0.49);
const scalar Vwr(0.08);
const scalar Vw1(0.24);
const scalar Nada1(1.10);
const scalar Sl1(700.0);

// Permeabilidade
const scalar Ks(0.000000012);
const scalar Kw1(0.000000000000012);
const scalar Nada2(1.25);
const scalar Sl2(40.0);


scalar temp1(0.0);
scalar temp2(0.0);
scalar isaturate(0.0);
scalar j(0.0);
scalar mark(0.0);








