//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// Geant4 class G4ChemDissociationChannels_option1
//
// H. Tran 16.12.2022
//
#include <cmath> //Multiple powers
#include "G4ChemDissociationChannels_option1.hh"
#include "G4EmDNAChemistry_option3.hh"  // Enger lab
#include "G4DNAChemistryManager.hh" //Enger lab

#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4Electron_aq.hh"
#include "G4FakeMolecule.hh"
#include "G4H2.hh"
#include "G4H2O.hh"
#include "G4H2O2.hh"
#include "G4H3O.hh"
#include "G4HO2.hh"
#include "G4Hydrogen.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4O2.hh"
#include "G4O3.hh"
#include "G4OH.hh"
#include "G4Oxygen.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ChemDissociationChannels_option1::G4ChemDissociationChannels_option1(double initialTemperature) //Constructor that takes temperature as input
 : Temperature(initialTemperature) {}

G4ChemDissociationChannels_option1::G4ChemDissociationChannels_option1()
    : Temperature(300.0)
{
    G4cout << "Default constructor called. Temperature = " << Temperature << G4endl;
}                                                                                                //to keep the constructor without parameters.

void G4ChemDissociationChannels_option1::ConstructMolecule(double Temperature) //Instantiation of all molecules with their appropriate temperature-dependent diffusion coefficient value 
{

  G4cout << "[DEBUG] fDissociationChannel initialized with T = " << Temperature << G4endl;//for debugging

  //-----------------------------------
  //  G4Electron::Definition(); // safety

  //-----------------------------------
  // Create the definition

  G4H2O::Definition();
  G4Hydrogen::Definition();
  G4H3O::Definition();
  G4OH::Definition();
  G4Electron_aq::Definition();
  G4H2O2::Definition();
  G4H2::Definition();

  G4O2::Definition();
  G4HO2::Definition();
  G4Oxygen::Definition();
  G4O3::Definition();

  //-----------------------------------
  // Temperature dependent polynomials of diffusion coefficient of different species.

  G4double t = Temperature; //Kelvin
  G4double T = t - 273.15; //Celsius

  G4double p = 0.999 + 1.094e-4*T - 7.397e-6*pow(T,2) + 2.693e-8*pow(T,3) - 4.714e-11*pow(T,4); //density

  G4double D_oh = (0.00408*pow(t,1.38897)+p*(1.0/pow(t,2)-174613.7394/t+847.03131-0.71041*t)+pow(p,2)*log(p)*(1.0/pow(t,2)-153864.1999/t+782.41564-0.64852*t)+pow(p,2)*(1.0/pow(t,2)+163191.6423/t-795.71944+0.63849*t))*1e-9;

  G4double D_h2o2 = (0.0046471*pow(t,1.2939)+p*(1.05646e8/pow(t,2)-511648/t+691.903-0.1729*t)+pow(p,2)*log(p)*(8.88978e7/pow(t,2)-358381/t+334.773+0.0213276*t)+pow(p,2)*(-1.01281e8/pow(t,2)+475755/t -608.627+0.116203*t))*1e-9;


  G4double D_o2 = (1.82779*pow(t,0.422868)+p*(1.0/pow(t,2)-102443/t+334.021-0.119239*t)+pow(p,2)*log(p)*(1.0/pow(t,2)-102959/t+334.195-0.117517*t)+pow(p,2)*(1.0/pow(t,2)+100433/t-347.059+0.125558*t))*1e-9;


  G4double D_h2 = (116.211*pow(t,0.0538917)+p*(1.0/pow(t,2)+372312/t-1794.21+1.42193*t)+pow(p,2)*log(p)*(1.0/pow(t,2)+476427/t -1880.99+1.62233*t)+pow(p,2)*(1.0/pow(t,2)-377905/t+1654.15-1.39764*t))*1e-9;


  G4double D_om = (1.82779*pow(t,0.422868)+p*(1.0/pow(t,2)-102443/t+334.021-0.119239*t)+pow(p,2)*log(p)*(1.0/pow(t,2)-102959/t+334.195-0.117517*t)+pow(p,2)*(1.0/pow(t,2)+100433/t-347.059+0.125558*t)-0.55152409)*1e-9;

  G4double D_o3m = (1.82779*pow(t,0.422868)+p*(1.0/pow(t,2)-102443/t+334.021-0.119239*t)+pow(p,2)*log(p)*(1.0/pow(t,2)-102959/t+334.195-0.117517*t)+pow(p,2)*(1.0/pow(t,2)+100433/t-347.059+0.125558*t)-0.55152409)*1e-9;

  G4double D_o3 = D_o3m;

  G4double D_oxy = (1.82779*pow(t,0.422868)+p*(1.0/pow(t,2)-102443/t+334.021-0.119239*t)+pow(p,2)*log(p)*(1.0/pow(t,2)-102959/t+334.195-0.117517*t)+pow(p,2)*(1.0/pow(t,2)+100433/t-347.059+0.125558*t)-0.55152409)*1e-9;

  G4double D_o2m = (1.82779*pow(t,0.422868)+p*(1.0/pow(t,2)-102443/t+334.021-0.119239*t)+pow(p,2)*log(p)*(1.0/pow(t,2)-102959/t+334.195-0.117517*t)+pow(p,2)*(1.0/pow(t,2)+100433/t-347.059+0.125558*t)-0.55152409-0.25)*1e-9;


  G4double D_ho2 = D_h2o2;


  G4double D_ho2m = (0.0046471*pow(t,1.2939)+p*(1.05646e8/pow(t,2)-511648/t+691.903-0.1729*t)+pow(p,2)*log(p)*(8.88978e7/pow(t,2)-358381/t+334.773+0.0213276*t)+pow(p,2)*(-1.01281e8/pow(t,2)+475755/t -608.627+0.116203*t)-1.034362355)*1e-9;

  G4double log10D_eaq = (-0.7638-1058.18/t-pow((254.92/t),40));

  G4double D_eaq = pow(10,log10D_eaq)*1e-4;

//  G4double D_eaq = pow(10,(-0.7638-1058.18/T-pow((254.92/T),40)))*1e-4;


  G4double D_ohm = 2.666e-9+9.769e-11*T+3.303e-13*pow(T,2)-7.295e-16*pow(T,3);

  G4double D_hp = 5.361e-9 + 1.659e-10*T-7.48e-14*pow(T,2)-1.0018e-16*pow(T,3);

  G4double D_h = 5.361e-9 + 1.659e-10*T-7.48e-14*pow(T,2)-1.0018e-16*pow(T,3)-2.46e-09;


  auto G4OHm = new G4MoleculeDefinition("OH",/*mass*/ 17.00734 * g / Avogadro * c_squared,
                                        2.8e-9 * (m * m / s), -1,
                                          5, 0.958 * angstrom, // radius
                                          2 // number of atoms
                                        );

  auto G4HO2m = new G4MoleculeDefinition("HO_2", 33.0034 * g / Avogadro * c_squared,
                                         2.3e-9 * (m * m / s), -1, 0,
                                        2.1 * angstrom, 3);
  auto G4Om = new G4MoleculeDefinition("O", 15.99773 * g / Avogadro * c_squared,
                                       2.0e-9 * (m * m / s), 0, 0,
                                       2.0 * angstrom, 1);
  //____________________________________________________________________________
  auto molTable = G4MoleculeTable::Instance();
  molTable->CreateConfiguration("H3Op", G4H3O::Definition());
  molTable->GetConfiguration("H3Op")->SetDiffusionCoefficient(D_hp
                                                                                 * (m2 / s));
  molTable->GetConfiguration("H3Op")->SetVanDerVaalsRadius(0.25 * nm);

  molTable->CreateConfiguration("°OH", G4OH::Definition());
  molTable->GetConfiguration("°OH")->SetDiffusionCoefficient(D_oh * (m2 / s));
  molTable->GetConfiguration("°OH")->SetVanDerVaalsRadius(0.22 * nm);

  G4MolecularConfiguration* OHm =
    molTable->CreateConfiguration("OHm",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                  G4OHm,
                                                     -1,  // charge
                                                    D_ohm  * (m2 / s));
  OHm->SetMass(17.0079 * g / Avogadro * c_squared);
  OHm->SetVanDerVaalsRadius(0.33 * nm);

  molTable->CreateConfiguration("e_aq", G4Electron_aq::Definition());
  molTable->GetConfiguration("e_aq")->SetDiffusionCoefficient(D_eaq * (m2 / s));
  molTable->GetConfiguration("e_aq")->SetVanDerVaalsRadius(0.50 * nm);

  molTable->CreateConfiguration("H", G4Hydrogen::Definition());
  molTable->GetConfiguration("H")->SetDiffusionCoefficient(D_h * (m2 / s));
  molTable->GetConfiguration("H")->SetVanDerVaalsRadius(0.19 * nm);

  molTable->CreateConfiguration("H2", G4H2::Definition());
  molTable->GetConfiguration("H2")->SetDiffusionCoefficient(D_h2 * (m2 / s));
  molTable->GetConfiguration("H2")->SetVanDerVaalsRadius(0.14 * nm);

  molTable->CreateConfiguration("H2O2", G4H2O2::Definition());
  molTable->GetConfiguration("H2O2")->SetDiffusionCoefficient(D_h2o2 * (m2 / s));
  molTable->GetConfiguration("H2O2")->SetVanDerVaalsRadius(0.21 * nm);

  // molecules extension (RITRACKS)

  molTable->CreateConfiguration("HO2°", G4HO2::Definition());
  molTable->GetConfiguration("HO2°")->SetDiffusionCoefficient(D_ho2 * (m2 / s));
  molTable->GetConfiguration("HO2°")->SetVanDerVaalsRadius(0.21 * nm);

  G4MolecularConfiguration* HO2m =
    molTable->CreateConfiguration("HO2m",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                        G4HO2m,
                                                     -1,  // charge
                                                     D_ho2m * (m2 / s));
  HO2m->SetMass(33.00396 * g / Avogadro * c_squared);
  HO2m->SetVanDerVaalsRadius(0.25 * nm);

  molTable->CreateConfiguration("Oxy", G4Oxygen::Definition());
  molTable->GetConfiguration("Oxy")->SetDiffusionCoefficient(D_oxy * (m2 / s));
  molTable->GetConfiguration("Oxy")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration* Om =
    molTable->CreateConfiguration("Om",  // just a tag to store and retrieve from
                                                            // G4MoleculeTable
                                                      G4Om,
                                                     -1,  // charge
                                                     D_om * (m2 / s));
  Om->SetMass(15.99829 * g / Avogadro * c_squared);
  Om->SetVanDerVaalsRadius(0.25 * nm);

  molTable->CreateConfiguration("O2", G4O2::Definition());
  molTable->GetConfiguration("O2")->SetDiffusionCoefficient(D_o2 * (m2 / s));
  molTable->GetConfiguration("O2")->SetVanDerVaalsRadius(0.17 * nm);

  G4MolecularConfiguration* O2m =
    molTable->CreateConfiguration("O2m",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4O2::Definition(),
                                                     -1,  // charge
                                                     D_o2m * (m2 / s));
  O2m->SetMass(31.99602 * g / Avogadro * c_squared);
  O2m->SetVanDerVaalsRadius(0.22 * nm);

  molTable->CreateConfiguration("O3", G4O3::Definition());
  molTable->GetConfiguration("O3")->SetDiffusionCoefficient(D_o3 * (m2 / s));
  molTable->GetConfiguration("O3")->SetVanDerVaalsRadius(0.20 * nm);

  G4MolecularConfiguration* O3m =
    molTable->CreateConfiguration("O3m",  // just a tag to store and retrieve
                                                             // from G4MoleculeTable
                                                     G4O3::Definition(),
                                                     -1,  // charge
                                                     D_o3m * (m2 / s));
  O3m->SetMass(47.99375 * g / Avogadro * c_squared);
  O3m->SetVanDerVaalsRadius(0.20 * nm);

  molTable->CreateConfiguration("H2O(B)",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                   G4H2O::Definition(),
                                                   0,  // charge
                                                   0 * (m2 / s));

  molTable->CreateConfiguration("H3Op(B)",  // just a tag to store and retrieve
                                                               // from G4MoleculeTable
                                                   G4H3O::Definition(),
                                                   1,  // charge
                                                   0 * (m2 / s));

  molTable->CreateConfiguration("OHm(B)",  // just a tag to store and retrieve
                                                              // from G4MoleculeTable
                                                    G4OHm,
                                                   -1,  // charge
                                                   0 * (m2 / s));

  molTable->CreateConfiguration("NoneM", G4FakeMolecule::Definition());     //Method ConstructMolecule: Implementation of temperature-dependent diffusion coefficients of totally 15 species.

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ChemDissociationChannels_option1::ConstructDissociationChannels()
{
  //-----------------------------------
  // Get the molecular configuration
  auto molTable = G4MoleculeTable::Instance();
  G4MolecularConfiguration* OH = molTable->GetConfiguration("°OH");
  G4MolecularConfiguration* OHm = molTable->GetConfiguration("OHm");
  G4MolecularConfiguration* e_aq = molTable->GetConfiguration("e_aq");
  G4MolecularConfiguration* H2 = molTable->GetConfiguration("H2");
  G4MolecularConfiguration* H3O = molTable->GetConfiguration("H3Op");
  G4MolecularConfiguration* H = molTable->GetConfiguration("H");
  G4MolecularConfiguration* O = molTable->GetConfiguration("Oxy");

  //-------------------------------------
  // Define the decay channels
  G4MoleculeDefinition* water = G4H2O::Definition();
  G4MolecularDissociationChannel* decCh1;
  G4MolecularDissociationChannel* decCh2;
  G4MolecularDissociationChannel* decCh3;
  G4MolecularDissociationChannel* decCh4;
  G4MolecularDissociationChannel* decCh5;

  G4ElectronOccupancy* occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

  //////////////////////////////////////////////////////////
  //            EXCITATIONS                               //
  //////////////////////////////////////////////////////////
  G4DNAWaterExcitationStructure waterExcitation;
  //--------------------------------------------------------
  //---------------Excitation on the fifth layer------------

  decCh1 = new G4MolecularDissociationChannel("A^1B_1_Relax");
  decCh2 = new G4MolecularDissociationChannel("A^1B_1_DissociDecay");
  // Decay 1 : OH + H
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
  decCh1->SetProbability(0.35);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::NoDisplacement);

  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.65);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  occ->RemoveElectron(4, 1);  // this is the transition form ground state to
  occ->AddElectron(5, 1);  // the first unoccupied orbital: A^1B_1

  water->NewConfigurationWithElectronOccupancy("A^1B_1", *occ);
  water->AddDecayChannel("A^1B_1", decCh1);
  water->AddDecayChannel("A^1B_1", decCh2);

  //--------------------------------------------------------
  //---------------Excitation on the fourth layer-----------
  decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relax_Channel");
  decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociDecay");
  decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIoni_Channel");
  decCh4 = new G4MolecularDissociationChannel("A^1B_1_DissociDecay");
  decCh5 = new G4MolecularDissociationChannel("B^1A_1_DissociDecay2");

  // Decay 1 : energy
  decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
  decCh1->SetProbability(0.175);

  // Decay 2 : 2OH + H_2
  decCh2->AddProduct(H2);
  decCh2->AddProduct(OH);
  decCh2->AddProduct(OH);
  decCh2->SetProbability(0.0325);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  // Decay 3 : OH + H_3Op + e_aq
  decCh3->AddProduct(OH);
  decCh3->AddProduct(H3O);
  decCh3->AddProduct(e_aq);
  decCh3->SetProbability(0.50);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay 4 :  H + OH
  decCh4->AddProduct(H);
  decCh4->AddProduct(OH);
  decCh4->SetProbability(0.2535);
  decCh4->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  // Decay 5 : 2H + O
  decCh5->AddProduct(O);
  decCh5->AddProduct(H);
  decCh5->AddProduct(H);
  decCh5->SetProbability(0.039);
  decCh5->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay2);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3);  // this is the transition form ground state to
  occ->AddElectron(5, 1);  // the first unoccupied orbital: B^1A_1

  water->NewConfigurationWithElectronOccupancy("B^1A_1", *occ);
  water->AddDecayChannel("B^1A_1", decCh1);
  water->AddDecayChannel("B^1A_1", decCh2);
  water->AddDecayChannel("B^1A_1", decCh3);
  water->AddDecayChannel("B^1A_1", decCh4);
  water->AddDecayChannel("B^1A_1", decCh5);

  //-------------------------------------------------------
  //-------------------Excitation of 3rd layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci3rdLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci3rdLayer_Relax_Channel");

  // Decay channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(2));
  decCh2->SetProbability(0.5);

  // Electronic configuration of this decay
  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  occ->AddElectron(5, 1);

  // Configure the water molecule
  water->NewConfigurationWithElectronOccupancy("Exci3rdLayer", *occ);
  water->AddDecayChannel("Exci3rdLayer", decCh1);
  water->AddDecayChannel("Exci3rdLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 2nd layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci2ndLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci2ndLayer_Relax_Channel");

  // Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);

  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(3));
  decCh2->SetProbability(0.5);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  occ->AddElectron(5, 1);

  water->NewConfigurationWithElectronOccupancy("Exci2ndLayer", *occ);
  water->AddDecayChannel("Exci2ndLayer", decCh1);
  water->AddDecayChannel("Exci2ndLayer", decCh2);

  //-------------------------------------------------------
  //-------------------Excitation of 1st layer-----------------
  decCh1 = new G4MolecularDissociationChannel("Exci1stLayer_AutoIoni_Channel");
  decCh2 = new G4MolecularDissociationChannel("Exci1stLayer_Relax_Channel");

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  occ->AddElectron(5, 1);

  // Decay Channel 1 : : OH + H_3Op + e_aq
  decCh1->AddProduct(OH);
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(e_aq);
  decCh1->SetProbability(0.5);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

  // Decay channel 2 : energy
  decCh2->SetEnergy(waterExcitation.ExcitationEnergy(4));
  decCh2->SetProbability(0.5);

  water->NewConfigurationWithElectronOccupancy("Exci1stLayer", *occ);
  water->AddDecayChannel("Exci1stLayer", decCh1);
  water->AddDecayChannel("Exci1stLayer", decCh2);

  /////////////////////////////////////////////////////////
  //                  IONISATION                         //
  /////////////////////////////////////////////////////////
  //--------------------------------------------------------
  //------------------- Ionisation -------------------------

  decCh1 = new G4MolecularDissociationChannel("Ioni_Channel");

  // Decay Channel 1 : : OH + H_3Op
  decCh1->AddProduct(H3O);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(4, 1);
  // this is a ionized h2O with a hole in its last orbital
  water->NewConfigurationWithElectronOccupancy("Ioni5", *occ);
  water->AddDecayChannel("Ioni5", decCh1);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(3, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni4", *occ);
  water->AddDecayChannel("Ioni4", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(2, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni3", *occ);
  water->AddDecayChannel("Ioni3", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(1, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni2", *occ);
  water->AddDecayChannel("Ioni2", new G4MolecularDissociationChannel(*decCh1));

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->RemoveElectron(0, 1);
  water->NewConfigurationWithElectronOccupancy("Ioni1", *occ);
  water->AddDecayChannel("Ioni1", new G4MolecularDissociationChannel(*decCh1));

  //////////////////////////////////////////////////////////
  //            Dissociative Attachment                   //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("DissociAttachment_ch1");

  // Decay 1 : OHm + H
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OHm);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(1);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::DissociativeAttachment);

  *occ = *(water->GetGroundStateElectronOccupancy());
  occ->AddElectron(5, 1);  // H_2O^-

  water->NewConfigurationWithElectronOccupancy("DissociAttachment_ch1", *occ);
  water->AddDecayChannel("DissociAttachment_ch1", decCh1);

  //////////////////////////////////////////////////////////
  //            Electron-hole recombination               //
  //////////////////////////////////////////////////////////
  decCh1 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay1");
  decCh2 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay2");
  decCh3 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay3");
  decCh4 = new G4MolecularDissociationChannel("H2Ovib_DissociDecay4");

  // Decay 1 : 2OH + H_2
  decCh1->AddProduct(H2);
  decCh1->AddProduct(OH);
  decCh1->AddProduct(OH);
  decCh1->SetProbability(0.1365);
  decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

  // Decay 2 : OH + H
  decCh2->AddProduct(OH);
  decCh2->AddProduct(H);
  decCh2->SetProbability(0.3575);
  decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

  // Decay 3 : 2H + O(3p)
  decCh3->AddProduct(O);
  decCh3->AddProduct(H);
  decCh3->AddProduct(H);
  decCh3->SetProbability(0.156);
  decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay2);

  // Decay 4 : relaxation
  decCh4->SetProbability(0.35);

  const auto pH2Ovib = G4H2O::Definition()->NewConfiguration("H2Ovib");
  assert(pH2Ovib != nullptr);

  water->AddDecayChannel(pH2Ovib, decCh1);
  water->AddDecayChannel(pH2Ovib, decCh2);
  water->AddDecayChannel(pH2Ovib, decCh3);
  water->AddDecayChannel(pH2Ovib, decCh4);

  delete occ;
}
