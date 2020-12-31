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
                            Cp1Species[LiquidSpeciesIndices(j)].internalField()[celli] = species_map->cpL(listLiquidSpecies_[j], TCells[celli]); 
                            dhev[j].internalField()[celli] = species_map->deltaHv(listLiquidSpecies_[j], TCells[celli]);                                   
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
                            gamma[j].internalField()[celli] =  liquid_mixtureCells->gamma(TCells[celli], pCells[celli], XL_aux)[j];

                            double pVap = species_map->pVap(listLiquidSpecies_[j], TCells[celli], pCells[celli]);
                            double pTot = pCells[celli];
                            double gammaj = gamma[j].internalField()[celli];


                            if (fugacity == "Ideal")
                            Keq[j].internalField()[celli] = (pVap) /pTot;
                      
                            else if (fugacity == "NonIdeal")
                            Keq[j].internalField()[celli] = (pVap*gammaj) /pTot;

                            else
                            {
                             Info << "Fugacity model not known!" << endl;
			     abort();
                            } 



                          }


			    // Mixing properties of liquid phase
                              cpLCells[celli] = liquid_mixtureCells->cpL_mix(TCells[celli], XL_aux);
                              rhoLCells[celli] = liquid_mixtureCells->rhoL_mix(TCells[celli], pCells[celli], XL_aux);
                              muLCells[celli] = liquid_mixtureCells->etaL_mix(TCells[celli], XL_aux);
                              nuLCells[celli] = muLCells[celli]/rhoLCells[celli];
                              psiLCells[celli] = 0.0;
		              cTotLCells[celli] = rhoLCells[celli]/MWmixCells[celli];
                              lambdaLCells[celli] = liquid_mixtureCells->lambdaL_mix(TCells[celli], XL_aux);

/*
		
                          // Liquid phase diffusivities (Stephan-Maxwell->pseudo Fick approach): PROBLEMA: si tiene conto anche dei gradienti l-v all'interfaccia
                          if(NLS>1)
                          {
                            std::vector< std::vector<double> > gradX_aux(NLS-1);
                            for (unsigned int i=0;i<NLS-1;i++)
                            {
                              gradX_aux[i].resize(3);
                            }

                            Eigen::MatrixXd  gradY_aux(NLS, 3);
                            Eigen::MatrixXd MolarDiffusFluxes(NLS-1, 3);

                            for (unsigned int i=0;i<NLS-1;i++)
                            {
                              for (unsigned int j=0; j<3; j++)
                              {                          
                                gradX_aux[i][j] =  gradX[i].internalField()[celli][j];    
                                gradY_aux(i,j) =   gradY[i].internalField()[celli][j];                             
                              }
                            }

                              for (unsigned int j=0; j<3; j++)
                              {                            
                                gradY_aux(NLS-1,j) =   gradY[NLS-1].internalField()[celli][j];                             
                              }

                            MolarDiffusFluxes = liquid_mixtureCells[celli]->molar_diffusion_fluxes_vectors(TCells[celli], pCells[celli], XL_aux, cTotLCells[celli], gradX_aux);


                            for (unsigned int j=0;j<NLS-1;j++)
                            {                                                             
                              Dmix1[LiquidSpeciesIndices(j)].internalField()[celli] = ( MolarDiffusFluxes.row(j).norm() *1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(j)) ) 
                                                                         / ( rhoLCells[celli] * gradY_aux.row(j).norm() + 1e-6) + 1e-20;                                 
                            }
                          
                            Eigen::VectorXd sumLiquidFluxes(3);

                            for (unsigned int j=0;j<NLS-1;j++)
                            {                          
                              sumLiquidFluxes += MolarDiffusFluxes.row(j)*1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(j));                             
                            }
          
                             Dmix1[LiquidSpeciesIndices(NLS-1)].internalField()[celli] = sumLiquidFluxes.norm() / ( rhoLCells[celli] * gradY_aux.row(NLS-1).norm() + 1e-6) + 1e-20; 
                         }



                    if(NLS == 1)
                    LiquidDiffFlux[LiquidSpeciesIndices(0)].internalField()[celli] *= 0.0;

                    else
                    {
                       Eigen::MatrixXd MolarDiffusFluxes(NLS-1, 3);

                       std::vector< std::vector<double> > gradX_aux(NLS-1);
                       for (unsigned int i=0;i<NLS-1;i++)
                       {
                          gradX_aux[i].resize(3);
                       }

                       for (unsigned int i=0;i<NLS-1;i++)
                       {
                         for (unsigned int j=0; j<3; j++)
                         {                          
                           gradX_aux[i][j] =  gradX[i].internalField()[celli][j];                               
                         }
                       }

                       for (unsigned int i=0;i<NLS-1;i++)
                       {
                         MolarDiffusFluxes = liquid_mixtureCells[celli]->molar_diffusion_fluxes_vectors(TCells[celli], pCells[celli], XL_aux, cTotLCells[celli], gradX_aux);

                         for (unsigned int j=0; j<3; j++)
                         {                          
                           LiquidDiffFlux[LiquidSpeciesIndices(i)].internalField()[celli][j] = MolarDiffusFluxes(i,j);
                           LiquidDiffFlux[LiquidSpeciesIndices(i)].internalField()[celli][j] *= 1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(i));                               
                         }

                       }
                            
                       Eigen::VectorXd sumLiquidFluxes(3);

                       for (unsigned int i=0;i<NLS-1;i++)
                       {  
                         for (unsigned int j=0; j<3; j++)
                         {                          
                           sumLiquidFluxes[j] += LiquidDiffFlux[LiquidSpeciesIndices(i)].internalField()[celli][j];        
                                                 
                         }                        
                                                     
                       }

                       for (unsigned int j=0; j<3; j++)
                       {                          
                         LiquidDiffFlux[LiquidSpeciesIndices(NLS-1)].internalField()[celli][j] = -sumLiquidFluxes[j];                                                         
                       }   
  
                     }
*/
              }



	// BOUNDARIES
	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
		fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcTotL = cTot1.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsiL  = psi1.boundaryFieldRef()[patchi];
		fvPatchScalarField& prhoL  = rho1.boundaryFieldRef()[patchi];
		fvPatchScalarField& pmuL   = mu1.boundaryFieldRef()[patchi];
		fvPatchScalarField& pnuL   = nu1.boundaryFieldRef()[patchi];
		fvPatchScalarField& plambdaL = lambda1.boundaryFieldRef()[patchi];
		fvPatchScalarField& pcpL = cp1.boundaryFieldRef()[patchi];
		fvPatchScalarField& pMWmix = MWmix.boundaryFieldRef()[patchi];
               

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
                              XL_aux[j] = X[LiquidSpeciesIndices(j)].boundaryFieldRef()[patchi][facei];                                     
                            }
                           }
*/


                        for (unsigned int j=0;j<NLS;j++)
                        {                                        
                          Cp1Species[LiquidSpeciesIndices(j)].boundaryFieldRef()[patchi][facei] = species_map->cpL(listLiquidSpecies_[j], pT[facei]);
                          dhev[j].boundaryFieldRef()[patchi][facei] = species_map->deltaHv(listLiquidSpecies_[j], pT[facei]); 
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
                            gamma[j].boundaryFieldRef()[patchi][facei] = liquid_mixtureCells->gamma(pT[facei], pp[facei], XL_aux)[j];

                            double pVap = species_map->pVap(listLiquidSpecies_[j], pT[facei], pp[facei]);
                            double pTot = pp[facei];
                            double gammaj = gamma[j].boundaryFieldRef()[patchi][facei];


                          if (fugacity == "Ideal")
                          Keq[j].boundaryFieldRef()[patchi][facei] = (pVap) /pTot;
                      
                          else if (fugacity == "NonIdeal")
                          Keq[j].boundaryFieldRef()[patchi][facei] = (pVap*gammaj) /pTot;

                          else
                           {
                            Info << "Fugacity model not known!" << endl;
			    abort();
                           } 


                        }


			    // Mixing rules
                              pcpL[facei] = liquid_mixtureCells->cpL_mix(pT[facei], XL_aux);
                              prhoL[facei] = liquid_mixtureCells->rhoL_mix(pT[facei], pp[facei], XL_aux);
                              pmuL[facei] = liquid_mixtureCells->etaL_mix(pT[facei], XL_aux);
                              pnuL[facei] = pmuL[facei]/prhoL[facei];
                              ppsiL[facei] = 0.0;
		              pcTotL[facei] = prhoL[facei]/pMWmix[facei];
                              plambdaL[facei] = liquid_mixtureCells->lambdaL_mix(pT[facei], XL_aux);
/*

                          // Liquid phase diffusivities (Stephan-Maxwell->pseudo Fick approach): PROBLEMA: si tiene conto anche dei gradienti l-v all'interfaccia
                          if(NLS>1)
                          {
                            std::vector< std::vector<double> > gradX_aux(NLS-1);
                            for (unsigned int i=0;i<NLS-1;i++)
                            {
                              gradX_aux[i].resize(3);
                            }

                            Eigen::MatrixXd  gradY_aux(NLS, 3);
                            Eigen::MatrixXd MolarDiffusFluxes(NLS-1, 3);

                            for (unsigned int i=0;i<NLS-1;i++)
                            {
                              for (unsigned int j=0; j<3; j++)
                              {                          
                                gradX_aux[i][j] =  gradX[i].boundaryFieldRef()[patchi][facei][j];    
                                gradY_aux(i,j) =   gradY[i].boundaryFieldRef()[patchi][facei][j];                             
                              }
                            }

                              for (unsigned int j=0; j<3; j++)
                              {                            
                                gradY_aux(NLS-1,j) =   gradY[NLS-1].boundaryFieldRef()[patchi][facei][j];                             
                              }

                            MolarDiffusFluxes = liquid_mixturePatch[facei]->molar_diffusion_fluxes_vectors(pT[facei], pp[facei], XL_aux, pcTotL[facei], gradX_aux);

                            for (unsigned int j=0;j<NLS-1;j++)
                            {                                                             
                              Dmix1[LiquidSpeciesIndices(j)].boundaryFieldRef()[patchi][facei] = ( MolarDiffusFluxes.row(j).norm() *1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(j)) ) 
                                                                                               / ( prhoL[facei] * gradY_aux.row(j).norm() + 1e-6);                                 
                            }
                          
                            Eigen::VectorXd sumLiquidFluxes(3);

                            for (unsigned int j=0;j<NLS-1;j++)
                            {                          
                              sumLiquidFluxes += MolarDiffusFluxes.row(j)*1e-3*thermodynamicsMapXML->MW(LiquidSpeciesIndices(j));                             
                            }
                             Dmix1[LiquidSpeciesIndices(NLS-1)].boundaryFieldRef()[patchi][facei] = sumLiquidFluxes.norm() / ( prhoL[facei] * gradY_aux.row(NLS-1).norm() + 1e-6); 
                         }
*/

            }

         }




// Transform saturation molar fractions into saturation mass fractions (Banerjee et.al. 2013)

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
          //dimensionedScalar MWair("MWair", dimless, 100.204998);
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
