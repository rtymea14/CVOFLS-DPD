/*---------------------------------------------------------------------------------------------------------------------------*\
Description

  The liquid properties ( density, Keq, activity coefficients conductivity, diffusivities, specific heat and viscosity) 
  are evaluated based on the OpenSMOKE++ library

\*----------------------------------------------------------------------------------------------------------------------------*/




{

  if(speciesEquations == true)
  {
	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
	Info<< "Liquid properties evaluation... " ;


		const scalarField& TCells = T.internalField();
		const scalarField& pCells = p.internalField();
		scalarField& rhoLCells = rho1.internalField();
		scalarField& cTotLCells = cTot1.internalField();
		scalarField& psiLCells =  psi1.internalField();
		scalarField& muLCells  =  mu1.internalField();
		scalarField& nuLCells  =  nu1.internalField();
		scalarField& lambdaLCells = lambda1.internalField();
		scalarField& cpLCells = cp1.internalField();
		scalarField& MWmixCells = MWmix.internalField();

// Update gradient of liquid species
/*
                for (unsigned int i=0;i<NLS;i++)
                {                          
                  gradX[i] = fvc::grad(X[LiquidSpeciesIndices(i)]);  
                  gradY[i] = fvc::grad(Y[LiquidSpeciesIndices(i)]);                                 
                }
*/

               forAll(TCells, celli)
               {

                          // Setting vector of liquid species
                          std::vector<double> XL_aux(NLS);
                         
                         //  if(NLS == 1)  
                             XL_aux[0] = 1.0; 
/*
                           else
                           {
                            for (unsigned int j=0;j<NLS;j++)
                            {                          
                              XL_aux[j] = X[LiquidSpeciesIndices(j)].internalField()[celli];                                     
                            }
                           }
*/


                          // Setting liquid Cp and heat of vaporization of species
                          for (unsigned int j=0;j<NLS;j++)
                          {                           
                            Cp1Species[LiquidSpeciesIndices(j)].internalField()[celli] = Liquid->Cp(pCells[celli], TCells[celli]);  // for now just set a const. value for the liquid //species_map->cpL(listLiquidSpecies_[j], TCells[celli]); 
                            dhev[j].internalField()[celli] = Liquid->hl(pCells[celli], TCells[celli]);	   //species_map->deltaHv(listLiquidSpecies_[j], TCells[celli]);                                   
                          }
/*


                          // Correct Molar Fractions for the sum to be 1
                          double sumXelements = 0.0;                          
                          for (unsigned int j=0;j<NLS;j++)
                          {                          
                            sumXelements += XL_aux[j];                                         
                          }

                           if(sumXelements != 1.0)
                           {
                               unsigned int MaxIndex;
                               double MaxMolarFraction = *std::max_element(XL_aux.begin(), XL_aux.end()); 

                               for (unsigned int j=0;j<NLS;j++)
                               {
                                   if(XL_aux[j] == MaxMolarFraction)
                                   MaxIndex = j;
                               }

                               for (unsigned int j=0;j<NLS;j++)
                               {  
                                 if(j != MaxIndex)
                                 sumXelements += XL_aux[j];                                         
                               }

                               XL_aux[MaxIndex] = 1.0 - sumXelements;

                           }
*/



                          // Thermodynamic equilibrium
                          for (unsigned int j=0;j<NLS;j++)
                          {  

                            double pVap = Liquid->pv(pCells[celli], TCells[celli]); 
                            double pTot = pCells[celli];

                            if (fugacity == "Raoult")
                             Keq[j].internalField()[celli] = pVap/pTot;

                            else if (fugacity == "PengRobinson")
                            {
                                FugacityCoefficientPure[j].internalField()[celli] =  Liquid->fugacityCoefficientPure(TCells[celli], pVap);
                                double FugacityCoefficientPurej = FugacityCoefficientPure[j].internalField()[celli];

                                Keq[j].internalField()[celli] = (pVap*FugacityCoefficientPurej) /(pTot*FugacityCoefficientMixFixed);
                            }

                            else
                            {
                             Info << "Available fugacity models are -Raoult- and -PengRobinson-" << endl;
			     abort();
                            } 


                            if (activity == "Ideal")
                              Keq[j].internalField()[celli] *= 1.0 ;

                            else if (activity == "Unifac")
                            {
                              gamma[j].internalField()[celli] =  Liquid->gamma(TCells[celli], pCells[celli])[j];
                              double gammaj = gamma[j].internalField()[celli];

                              Keq[j].internalField()[celli] *= gammaj;
                            }

                            else
                            {
                             Info << "Available activity models are -Ideal- and -Unifac-" << endl;
			     abort();
                            } 


                            if (poyntingCorrection == true)
                            {
                               double molarVolumej = 1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(j)) / Liquid->rho(pCells[celli], TCells[celli]);  
                               Poynting[j].internalField()[celli] = std::exp(molarVolumej * (pCells[celli] - pVap) /(8.3144621 * TCells[celli]));

                               Keq[j].internalField()[celli] *= Poyntingfixed ;
                            } 

                            else  
                            Keq[j].internalField()[celli] *= 1.0 ;

                          }


			    // Mixing properties of liquid phase
                              cpLCells[celli] = Liquid->Cp(pCells[celli], TCells[celli]); 

                            if (NLS == 1)
                              rhoLCells[celli] = Liquid->rho(pCells[celli], TCells[celli]); 

                            else
                              rhoLCells[celli] = Liquid->rho(pCells[celli], TCells[celli]); 

                              muLCells[celli] = Liquid->mu(pCells[celli], TCells[celli]); 
                              nuLCells[celli] = muLCells[celli]/rhoLCells[celli];
                              psiLCells[celli] = 0.0;
		              		  cTotLCells[celli] = rhoLCells[celli]/MWmixCells[celli];
                              lambdaLCells[celli] = Liquid->K(pCells[celli], TCells[celli]); 

              }


	// BOUNDARIES
	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pT = T.boundaryField()[patchi];
		fvPatchScalarField& pp = p.boundaryField()[patchi];
		fvPatchScalarField& pcTotL = cTot1.boundaryField()[patchi];
		fvPatchScalarField& ppsiL  = psi1.boundaryField()[patchi];
		fvPatchScalarField& prhoL  = rho1.boundaryField()[patchi];
		fvPatchScalarField& pmuL   = mu1.boundaryField()[patchi];
		fvPatchScalarField& pnuL   = nu1.boundaryField()[patchi];
		fvPatchScalarField& plambdaL = lambda1.boundaryField()[patchi];
		fvPatchScalarField& pcpL = cp1.boundaryField()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryField()[patchi];
               

		forAll(pT, facei)
		{

                          std::vector<double> XL_aux(NLS);
                       //    if(NLS == 1)  
                             XL_aux[0] = 1.0; 
/*
                           else
                           {
                            for (unsigned int j=0;j<NLS;j++)
                            {                          
                              XL_aux[j] = X[LiquidSpeciesIndices(j)].boundaryField()[patchi][facei];                                     
                            }
                           }
*/


                        for (unsigned int j=0;j<NLS;j++)
                        {                                        
                          Cp1Species[LiquidSpeciesIndices(j)].boundaryField()[patchi][facei] = Liquid->Cp(pp[facei], pT[facei]); 
                          dhev[j].boundaryField()[patchi][facei] = Liquid->hl(pp[facei], pT[facei]);  
                        }
/*
                          // Correct Molar Fractions for the sum to be 1
                          double sumXelements = 0.0;                          
                          for (unsigned int j=0;j<NLS;j++)
                          {                          
                            sumXelements += XL_aux[j];                                         
                          }

                           if(sumXelements != 1.0)
                           {
                               unsigned int MaxIndex;
                               double MaxMolarFraction = *std::max_element(XL_aux.begin(), XL_aux.end()); 
                               for (unsigned int j=0;j<NLS;j++)
                               {
                                   if(XL_aux[j] == MaxMolarFraction)
                                   MaxIndex = j;
                               }

                               for (unsigned int j=0;j<NLS;j++)
                               {  
                                 if(j != MaxIndex)
                                 sumXelements += XL_aux[j];                                         
                               }

                               XL_aux[MaxIndex] = 1.0 - sumXelements;
                           }
*/


                        // Thermodynamic equilibrium
                        for (unsigned int j=0;j<NLS;j++)
                        { 
                            double pVap = Liquid->pv(pp[facei], pT[facei]); 
                            double pTot = pp[facei];


                            if (fugacity == "Raoult")
                             Keq[j].boundaryField()[patchi][facei]  = pVap/pTot;

                            else if (fugacity == "PengRobinson")
                            {
                              FugacityCoefficientPure[j].boundaryField()[patchi][facei] =  Liquid->fugacityCoefficientPure(pT[facei], pVap);
                              double FugacityCoefficientPurej = FugacityCoefficientPure[j].boundaryField()[patchi][facei];

                              Keq[j].boundaryField()[patchi][facei]  = (pVap*FugacityCoefficientPurej) /(pTot*FugacityCoefficientMixFixed);
                            }

                            else
                            {
                             Info << "Available fugacity models are -Raoult- and -PengRobinson-" << endl;
			     abort();
                            } 



                            if (activity == "Ideal")
                             Keq[j].boundaryField()[patchi][facei]  *= 1.0 ;

                            else if (activity == "Unifac")
                            {
                              gamma[j].boundaryField()[patchi][facei] = Liquid->gamma(pT[facei], pp[facei])[j];
                              double gammaj = gamma[j].boundaryField()[patchi][facei];

                              Keq[j].boundaryField()[patchi][facei]  *= gammaj;
                            }

                            else
                            {
                             Info << "Available activity models are -Ideal- and -Unifac-" << endl;
			     abort();
                            } 


                            if (poyntingCorrection == true)
                            {
                              double molarVolumej = 1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(j)) / Liquid->rho(pp[facei], pT[facei]); 
                              Poynting[j].boundaryField()[patchi][facei] = std::exp(molarVolumej * (pp[facei] - pVap) /(8.3144621 * pT[facei]));

                              Keq[j].boundaryField()[patchi][facei] *= Poyntingfixed ;
                            }

                            else  
                             Keq[j].boundaryField()[patchi][facei] *= 1.0 ;

                        }


			    // Mixing rules
                              pcpL[facei] = Liquid->Cp(pp[facei], pT[facei]); 

                            if (NLS == 1)
                              prhoL[facei] = Liquid->rho(pp[facei], pT[facei]); 

                            else
                              prhoL[facei] = Liquid->rho(pp[facei], pT[facei]); 


                              pmuL[facei] = Liquid->mu(pp[facei], pT[facei]);   
                              pnuL[facei] = pmuL[facei]/prhoL[facei];
                              ppsiL[facei] = 0.0;
		              pcTotL[facei] = prhoL[facei]/pMWmix[facei];
                              plambdaL[facei] = Liquid->K(pp[facei], pT[facei]); 

            }

         }






/*---------------------------------------------------------------------------------------------------------------------------*\
Description

  Based on Banerjee at.al (2013), the saturation concentration calculated in molar fractions is here converted in mass fractions

\*----------------------------------------------------------------------------------------------------------------------------*/



        volScalarField sumXiMWi = 0.0*Keq[0];
        volScalarField sumXi = 0.0*Keq[0];

        for (unsigned int j=0;j<NLS;j++)
        {  
          double MWj = thermodynamicsMapXML->MW(LiquidSpeciesIndices(j));

          sumXiMWi += MWj*Keq[j];
          sumXi += Keq[j];
        }



        for (unsigned int j=0;j<NLS;j++)
        {  
          double MWj = thermodynamicsMapXML->MW(LiquidSpeciesIndices(j));
          Keq[j] =  Keq[j]*MWj /( sumXiMWi + MWinert.value()*(1.0-sumXi) );
        }




	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

   }

    else
    {
       Info << "Species equations must be activated for the evaluation of gas and liquid properties!" << endl;
       abort();
    }



}
