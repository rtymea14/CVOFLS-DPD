const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

//bool MULESCorr(alphaControls.lookupOrDefault<Switch>("MULESCorr", false));

// Isotropic compression coefficient
scalar icAlpha
(
    alphaControls.lookupOrDefault<scalar>("icAlpha", 1)
);
