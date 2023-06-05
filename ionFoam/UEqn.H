    MRF.correctBoundaryVelocity(U);
    
    //scalar nu = 1e-3;
    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)  // or "fvm::div(rhoStokesf, U)"
      + MRF.DDt(rho, U)
      //+ turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
    
    
    //volTensorField gradU(fvc::grad(U));
    //volScalarField devUU(dev2(gradU));
    
    //volScalarField devUU = dev(fvc::grad(U));
    //fvVectorMatrix nuEffEqn
    //(
           //fvm::div(ionf, nuEff*devUU)
          
    //);
    
    //nuEff.relax();
    //solve(nuEffEqn);
    //fvOptions.correct(nuEff); 
