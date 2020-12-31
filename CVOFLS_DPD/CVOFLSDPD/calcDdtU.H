{
    gradPsi = fvc::grad(Psii);
    //gradAlpha = fvc::grad(voidfraction);
    ddtU_ = fvc::ddt(U);
    gradU_ = fvc::grad(U);
//fvScalarMatrix alphaSEqn
//(
//    fvm::Sp(scalar(1),alphaS) - fvm::laplacian(DAlpha,alphaS) == alpha1
//);
//alphaSEqn.solve();
}
