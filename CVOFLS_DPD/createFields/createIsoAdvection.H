
/*


// Instantiating isoAdvection object for advancing alpha1
isoAdvection advector(alpha1, phi, U);

// Checking if U is to be treated as prescribed
const dictionary& UDict(mesh.solverDict(U.name()));

bool prescribedU
(
    UDict.lookupOrDefault<bool>("prescribedU", false)
);

// Checking if U is to be treated as periodic or reversing
const scalar period(UDict.lookupOrDefault<scalar>("period", 0));
scalar reverseTime(UDict.lookupOrDefault<scalar>("reverseTime", 0));

//Creating copies of initial U and phi for periodic flows
//Rewrite to  avoid these extra fields when predcribedU = false
const volVectorField U0(U);
const surfaceScalarField phi0(phi);

*/


// Instantiating isoAdvection object for advancing alpha1
isoAdvection advector(alpha1, phiAdv, UAdv);

// Checking if U is to be treated as prescribed
const dictionary& UAdvDict(mesh.solverDict(UAdv.name()));

bool prescribedUAdv
(
    UAdvDict.lookupOrDefault<bool>("prescribedUAdv", false)
);

// Checking if U is to be treated as periodic or reversing
const scalar period(UAdvDict.lookupOrDefault<scalar>("period", 0));
scalar reverseTime(UAdvDict.lookupOrDefault<scalar>("reverseTime", 0));

//Creating copies of initial U and phi for periodic flows
//Rewrite to  avoid these extra fields when predcribedU = false
const volVectorField U0Adv(UAdv);
const surfaceScalarField phi0Adv(phiAdv);
