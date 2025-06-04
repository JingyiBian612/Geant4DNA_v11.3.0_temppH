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
/*
 * G4EmDNAChemistry_option3.cc
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
 *  Updated: Hoang Tran : add SBS and IRT_syn models
*/

//These next two lines are included by Enger lab and are not part of original version
#define _USE_MATH_DEFINES
#include <cmath>

#include "G4EmDNAChemistry_option3.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4ProcessManager.hh"
// *** Processes and models for Geant4-DNA

#include "G4DNAElectronSolvation.hh"

#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNAMolecularIRTModel.hh"
#include "G4DNAIndependentReactionTimeModel.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNAIRT.hh"
#include "G4DNAElectronHoleRecombination.hh"
// particles
#include "G4Electron.hh"
#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4PhysicsListHelper.hh"

/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ProcessTable.hh"
#include "G4MolecularConfiguration.hh"
/****/
// factory
#include "G4PhysicsConstructorFactory.hh"
#include "G4ChemDissociationChannels_option1.hh"
//Parameter
#include "G4EmParameters.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistry_option3);

G4EmDNAChemistry_option3::G4EmDNAChemistry_option3(G4double T, G4double P) ://Implemenation of the chemistry constructor that takes temperature and pH as inputs.
    G4VUserChemistryList(true), 
    fDissociationChannel(T),Temperature(T), pH(P)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
  G4DNAChemistryManager::Instance()->setTemp(T);
  G4DNAChemistryManager::Instance()->setpH(pH);  // Call methods "setTemp(Temperature)" and "setpH(pH)" which initializes the values in the G4DNAChemistryManager class.
}

G4EmDNAChemistry_option3::G4EmDNAChemistry_option3() :
    G4VUserChemistryList(true)
{
  G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructMolecule()
{
  fDissociationChannel.ConstructMolecule(Temperature); //pass temperature to ConstructMolecule in G4ChemDissociationChannels_option1 for temperature dependent diffusion coefficient 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructDissociationChannels()
{
  G4ChemDissociationChannels_option1::ConstructDissociationChannels();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructReactionTablePhTemp(G4DNAMolecularReactionTable*
                                              theReactionTable, double t, double p)     //t: Temperature in Kelvin
{

  G4double T = t - 273.15; //Temperature in Celsius
  G4double inpH = p;
  G4double pH_neutral = 6.59e-5*pow(T,2) - 1.99e-2*T + 7.46; //neutral pH changes with temperature in Celsius 
  G4double inpH2 = (2*pH_neutral) - p;
  G4double H2O_concentration = 55.50 + 6.075e-3*T - 4.110e-4*pow(T,2) + 1.496e-6*pow(T,3) - 2.619e-9*pow(T,4); // water concentration in mol/L
  G4double Hp_concentration = pow(10,-inpH);
  G4double OHm_concentration = pow(10,-inpH2);

 auto model = G4EmParameters::Instance()->GetTimeStepModel();
  //-----------------------------------
  //Get the molecular configuration
  G4MolecularConfiguration* OH =
   G4MoleculeTable::Instance()->GetConfiguration("°OH");
  G4MolecularConfiguration* OHm =
   G4MoleculeTable::Instance()->GetConfiguration("OHm");
  G4MolecularConfiguration* e_aq =
   G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  G4MolecularConfiguration* H2 =
   G4MoleculeTable::Instance()->GetConfiguration("H2");
  G4MolecularConfiguration* H3Op =
   G4MoleculeTable::Instance()->GetConfiguration("H3Op");
  G4MolecularConfiguration* H =
   G4MoleculeTable::Instance()->GetConfiguration("H");
  G4MolecularConfiguration* H2O2 =
   G4MoleculeTable::Instance()->GetConfiguration("H2O2");
  G4MolecularConfiguration* HO2 =
   G4MoleculeTable::Instance()->GetConfiguration("HO2°");
  G4MolecularConfiguration* HO2m =
   G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration* O =
   G4MoleculeTable::Instance()->GetConfiguration("Oxy");
  G4MolecularConfiguration* Om =
   G4MoleculeTable::Instance()->GetConfiguration("Om");
  G4MolecularConfiguration* O2 =
   G4MoleculeTable::Instance()->GetConfiguration("O2");
  G4MolecularConfiguration* O2m =
   G4MoleculeTable::Instance()->GetConfiguration("O2m");
  G4MolecularConfiguration* O3 =
   G4MoleculeTable::Instance()->GetConfiguration("O3");
  G4MolecularConfiguration* O3m =
   G4MoleculeTable::Instance()->GetConfiguration("O3m");

  G4MolecularConfiguration* H2OB =
   G4MoleculeTable::Instance()->GetConfiguration("H2O(B)");
  G4MolecularConfiguration* H3OpB =
   G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  G4MolecularConfiguration* OHmB =
   G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");

  G4MolecularConfiguration* None =
   G4MoleculeTable::Instance()->GetConfiguration("NoneM");



  // Type I //
  //------------------------------------------------------------------
  G4double R = 8.314462618;   //gas constant (JK^-1)
  // *H + *H -> H2
  G4double R1_1 = 2.70e12*pow(M_E, -1867.5/t);   //Elliot and Bartels (2009)
  G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(
      //0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  R1_1 * (1e-3 * m3 / (mole * s)), H, H);   //replace 25C reaction rate with temperature-dependent polynomials, t in Kelvin.
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  reactionData->setTemp(t);       //Call setTemp method to set temperature in  reactionData class
  //------------------------------------------------------------------
  // e_aq + H* + H2O -> H2 + OH-
  G4double R1_2 = 1.14e13*pow(M_E,-1795.7/t);      //Elliot and Bartels (2009)
  reactionData = new G4DNAMolecularReactionData(
      //2.50e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
  R1_2 * (1e-3 * m3 / (mole * s)), e_aq, H);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
 //------------------------------------------------------------------
  // H + O(3p) -> OH
  G4double R1_3 = 2.03e10*pow(M_E, -12.6/(R*t));   //Burns (1981)
  reactionData = new G4DNAMolecularReactionData(
      //2.02e10 * (1e-3 * m3 / (mole * s)), H, O);
      R1_3 * (1e-3 * m3 / (mole * s)), H, O);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O- -> OH-
  G4double R1_4 = R1_3;  
  reactionData = new G4DNAMolecularReactionData(
      //2.00e10 * (1e-3 * m3 / (mole * s)), H, Om);
      R1_4 * (1e-3 * m3 / (mole * s)), H, Om);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O(3p) -> HO2
  G4double R1_5 = 2.03e10*pow(M_E, -12.6/(R*t));   
  reactionData = new G4DNAMolecularReactionData(
      //2.02e10 * (1e-3 * m3 / (mole * s)), OH, O);
      R1_5 * (1e-3 * m3 / (mole * s)), OH, O);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O(3p) -> O2
  G4double R1_6 = 2.03e10*pow(M_E,-12.6/(R*t));
  reactionData = new G4DNAMolecularReactionData(
      //2.02e10 * (1e-3 * m3 / (mole * s)), HO2, O);
      R1_6 * (1e-3 * m3 / (mole * s)), HO2, O);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + O(3p) -> O2
  G4double R1_7 = 2.21e10*pow(M_E,-12.6/(R*t));
  reactionData = new G4DNAMolecularReactionData(
      //2.20e10 * (1e-3 * m3 / (mole * s)), O, O);
      R1_7 * (1e-3 * m3 / (mole * s)), O, O);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);


// Type III //
  //------------------------------------------------------------------
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  G4double logR3_1 = 12.281-3.786e2/t - 6.673e4/(pow(t,2)) - 1.075e7/(pow(t,3));
  G4double R3_1 = pow(10,logR3_1);
 // G4double R3_1 = pow(10, 12.281-3.786e2/t-6.673e4/(pow(t,2))-1.075e7/(pow(t,3))); //for t<150 
  reactionData = new G4DNAMolecularReactionData(
      //0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
      R3_1 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  G4double LogR3_2 = 20.934 - 1.236e4/t + 6.364e6/pow(t,2) - 1.475e9/pow(t,3) + 1.237e11/pow(t,4);
  G4double R3_2 = pow(10,LogR3_2);
 // G4double R3_2 = pow(10, 20.934-1.236e4/t + 6.364e6/(pow(t,2))-1.475e9/(pow(t,3))+1.237e11/(pow(t,4)));
  reactionData = new G4DNAMolecularReactionData(
      //1.13e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
      R3_2 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O3- -> OH + O2
  G4double LogR3_3 = 16.410 - 4.888e3/t + 1.622e6/pow(t,2) - 2.004e8/pow(t,3);
  G4double R3_3 = pow(10,LogR3_3); 
//  G4double R3_3 = pow(10, 16.410-4.888e3/t+1.622e6/(pow(t,2)) - 2.004e8/(pow(t,3)));
  reactionData = new G4DNAMolecularReactionData(
      //9.0e10 * (1e-3 * m3 / (mole * s)), H3Op, O3m);
      R3_3 * (1e-3 * m3 / (mole * s)), H3Op, O3m);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

// Type II //

  //------------------------------------------------------------------
  // *OH + *H -> H2O
  G4double R2_1 = 4.26e11*pow(M_E, -1091.9/t);
  reactionData = new G4DNAMolecularReactionData(
     // 1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  R2_1 * (1e-3 * m3 / (mole * s)), OH, H);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O2 -> OH
  G4double R2_2 = 1.79e11*pow(M_E,-2533.6/t);
  reactionData = new G4DNAMolecularReactionData(
     // 3.50e7 * (1e-3 * m3 / (mole * s)), H, H2O2);
   R2_2 * (1e-3 * m3 / (mole * s)), H, H2O2);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH- -> eaq-
  G4double LogR2_3 = 22.970 -1.971e4/t + 1.137e7/pow(t,2) - 2.991e9/pow(t,3) + 2.803e11/pow(t,4);
  G4double R2_3 = pow(10,LogR2_3);
//  G4double R2_3 = pow(10, 22.970 -1.971e4/t + 1.137e7/(pow(t,2)) - 2.991e9/(pow(t,3)) + 2.803e11/(pow(t,4)));
  reactionData = new G4DNAMolecularReactionData(
     // 2.51e7 * (1e-3 * m3 / (mole * s)), H, OHm);
   R2_3 * (1e-3 * m3 / (mole * s)), H, OHm);
  reactionData->AddProduct(e_aq);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2 -> HO2
  G4double LogR2_4 = 10.704 + 2.840e2/t - 1.369e5/pow(t,2);
  G4double R2_4 = pow(10,LogR2_4);
//  G4double R2_4 = pow(10, 10.704+2.840e2/t-1.369e5/(pow(t,2)));
  reactionData = new G4DNAMolecularReactionData(
     // 2.10e10 * (1e-3 * m3 / (mole * s)), H, O2);
   R2_4 * (1e-3 * m3 / (mole * s)), H, O2);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + HO2 -> H2O2
  G4double R2_5 = 5.17e12*pow(M_E, -1824.2/t);
  reactionData = new G4DNAMolecularReactionData(
    //  1.00e10 * (1e-3 * m3 / (mole * s)), H, HO2);
    R2_5 * (1e-3 * m3 / (mole * s)), H, HO2);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2- -> HO2-
  G4double R2_6 = 5.17e12*pow(M_E,-1824.2/t);
  reactionData = new G4DNAMolecularReactionData(
      //1.00e10 * (1e-3 * m3 / (mole * s)), H, O2m);
   R2_6 * (1e-3 * m3 / (mole * s)), H, O2m);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *OH -> H2O2
  G4double LogR2_7  =8.054 + 2.193e3/t - 7.395e5/pow(t,2) + 6.870e7/pow(t,3);
  G4double R2_7 = pow(10,LogR2_7);
 // G4double R2_7 = pow(10, 8.054+2.193e3/t-7.395e5/(pow(t,2))+6.870e7/(pow(t,3)));
  reactionData = new G4DNAMolecularReactionData(
     // 0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
    R2_7 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2O2 -> HO2
  G4double R2_8 = 7.68e9*pow(M_E, -1661.4/t);
  reactionData = new G4DNAMolecularReactionData(
     // 2.88e7 * (1e-3 * m3 / (mole * s)), OH, H2O2);
   R2_8 * (1e-3 * m3 / (mole * s)), OH, H2O2);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2 -> H
  G4double LogR2_9 = -11.556 + 3.2546e4/t - 1.8623e7/pow(t,2) + 4.5543e9/pow(t,3) - 4.1364e11/pow(t,4);
  G4double R2_9 = pow(10,LogR2_9);
//  G4double R2_9 = pow(10, -11.556+3.2546e4/t-1.8623e7/(pow(t,2))+4.5543e9/(pow(t,3))-4.1364e11/(pow(t,4)));
  reactionData = new G4DNAMolecularReactionData(
     // 3.28e7 * (1e-3 * m3 / (mole * s)), OH, H2);
   R2_9 * (1e-3 * m3 / (mole * s)), OH, H2);
  reactionData->AddProduct(H);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + *OH -> OH-
  G4double LogR2_10 = 13.123 - 1.023e3/t + 7.634e4/pow(t,2);          
  G4double R2_10 = pow(10,LogR2_10);
//  G4double R2_10 = pow(10, 13.123-1.023e3/t+7.634e4/(pow(t,2)));
  reactionData = new G4DNAMolecularReactionData(
     // 2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  R2_10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH- -> O-
  G4double LogR2_11 = 13.339 - 2.220e3/t + 7.333e5/pow(t,2) - 1.065e8/pow(t,3);
  G4double R2_11 = pow(10,LogR2_11);
//  G4double R2_11 = pow(10, 13.339 -2.220e3/t +7.333e5/(pow(t,2))-1.065e8/(pow(t,3)));
  reactionData = new G4DNAMolecularReactionData(
    //  6.30e9 * (1e-3 * m3 / (mole * s)), OH, OHm);
  R2_11 * (1e-3 * m3 / (mole * s)), OH, OHm);
  reactionData->AddProduct(Om);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2 -> O2
  G4double R2_12 = 1.29e11*pow(M_E,-799.2/t);
  reactionData = new G4DNAMolecularReactionData(
     // 7.90e9 * (1e-3 * m3 / (mole * s)), OH, HO2);
  R2_12 * (1e-3 * m3 / (mole * s)), OH, HO2);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O2- -> O2 + OH-
  G4double R2_13 = 8.77e11*pow(M_E,-1306.2/t);
  reactionData = new G4DNAMolecularReactionData(
    //  1.07e10 * (1e-3 * m3 / (mole * s)), OH, O2m);
  R2_13 * (1e-3 * m3 / (mole * s)), OH, O2m);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2- -> HO2 + OH-
  G4double R2_14 = 4.5e12*pow(M_E,-1877.71360/t);
  reactionData = new G4DNAMolecularReactionData(
     // 8.32e9 * (1e-3 * m3 / (mole * s)), OH, HO2m);
  R2_14 * (1e-3 * m3 / (mole * s)), OH, HO2m);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O- -> HO2-
  G4double R2_15 = 1.45e13*pow(M_E, -2928.5/t);
  reactionData = new G4DNAMolecularReactionData(
    //  1.00e9 * (1e-3 * m3 / (mole * s)), OH, Om);
  R2_15 * (1e-3 * m3 / (mole * s)), OH, Om);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O3- -> O2- + HO2
  G4double R2_16 = R2_13;
  reactionData = new G4DNAMolecularReactionData(
    //  8.50e9 * (1e-3 * m3 / (mole * s)), OH, O3m);
  R2_16 * (1e-3 * m3 / (mole * s)), OH, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H2O2 -> OH- + *OH
  G4double R2_17 = 7.70e12*pow(M_E, -1889.6/t);
  reactionData = new G4DNAMolecularReactionData(
     // 1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  R2_17 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH- -> HO2-
  G4double LogR2_18 = 13.339 - 2.220e3/t + 7.333e5/pow(t,2) - 1.065e8/pow(t,3); 
  G4double R2_18 = pow(10,LogR2_18);  
//  G4double R2_18 = pow(10, 13.339-2.220e3/t+7.333e5/(pow(t,2)) - 1.065e8/(pow(t,3)));
  reactionData = new G4DNAMolecularReactionData(
     // 4.71e8 * (1e-3 * m3 / (mole * s)), H2O2, OHm);
  R2_18 * (1e-3 * m3 / (mole * s)), H2O2, OHm);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O(3p) -> HO2 + OH
  G4double R2_19 = 2.99e11*pow(M_E, -1876.36438/t);
  reactionData = new G4DNAMolecularReactionData(
    //  1.60e9 * (1e-3 * m3 / (mole * s)), H2O2, O);
  R2_19 * (1e-3 * m3 / (mole * s)), H2O2, O);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O- -> HO2 + OH-
  G4double R2_20 = 2.99176401e11*pow(M_E,-1.87636438e3/t);
  reactionData = new G4DNAMolecularReactionData(
     // 5.55e8 * (1e-3 * m3 / (mole * s)), H2O2, Om);
  R2_20 * (1e-3 * m3 / (mole * s)), H2O2, Om);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O(3p) -> H + OH
  G4double R2_21 = 4.837e3*pow(M_E, -34.6/(R*t));
  reactionData = new G4DNAMolecularReactionData(
     // 4.77e3 * (1e-3 * m3 / (mole * s)), H2, O);
  R2_21 * (1e-3 * m3 / (mole * s)), H2, O);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O- -> H + OH-
  G4double R2_22 = 2.32e10*pow(M_E, -1550.5/t);
  reactionData = new G4DNAMolecularReactionData(
   //   1.21e8 * (1e-3 * m3 / (mole * s)), H2, Om);
   R2_22 * (1e-3 * m3 / (mole * s)), H2, Om);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + O2 -> O2-
  G4double R2_23 = 2.52e12*pow(M_E, -1401.5/t);
  reactionData = new G4DNAMolecularReactionData(
    //  1.74e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  R2_23 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq + HO2 -> HO2-
  G4double R2_24 = 2.46e12*pow(M_E, -1563.6/t);
  reactionData = new G4DNAMolecularReactionData(
      //1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, HO2);
  R2_24 * (1e-3 * m3 / (mole * s)), e_aq, HO2);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + HO2 -> O2-
  G4double LogR2_25 = 13.339 - 2.220e3/t + 7.333e5/pow(t,2) - 1.065e8/pow(t,3);
  G4double R2_25 = pow(10,LogR2_25);
//  G4double R2_25 = pow(10, 13.339-2.220e3/t+7.333e5/(pow(t,2)) - 1.065e8/(pow(t,3)));
  reactionData = new G4DNAMolecularReactionData(
      //6.30e9 * (1e-3 * m3 / (mole * s)), OHm, HO2);
  R2_25 * (1e-3 * m3 / (mole * s)), OHm, HO2);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + O(3p) -> HO2-
  G4double R2_26 = R2_15;
  reactionData = new G4DNAMolecularReactionData(
     // 4.20e8 * (1e-3 * m3 / (mole * s)), OHm, O);
  R2_26 * (1e-3 * m3 / (mole * s)), OHm, O);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O(3p) -> O3
  G4double R2_27 = 3.41e11*pow(M_E,-1344.9/t);
  reactionData = new G4DNAMolecularReactionData(
     // 4.00e9 * (1e-3 * m3 / (mole * s)), O2, O);
   R2_27 * (1e-3 * m3 / (mole * s)), O2, O);
  reactionData->AddProduct(O3);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O- -> O3-
  G4double R2_28 = 3.41e11*pow(M_E,-1344.9/t);
  reactionData = new G4DNAMolecularReactionData(
    //  3.70e9 * (1e-3 * m3 / (mole * s)), O2, Om);
   R2_28 * (1e-3 * m3 / (mole * s)), O2, Om);
  reactionData->AddProduct(O3m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + HO2 -> H2O2 + O2
  G4double R2_29 = 2.78e9*pow(M_E, -2416.4/t);
  reactionData = new G4DNAMolecularReactionData(
     // 9.80e5 * (1e-3 * m3 / (mole * s)), HO2, HO2);
   R2_29 * (1e-3 * m3 / (mole * s)), HO2, HO2);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O2- -> HO2- + O2
  G4double R2_30 = 2.63e9*pow(M_E,-974.3/t);
  reactionData = new G4DNAMolecularReactionData(
     // 9.70e7 * (1e-3 * m3 / (mole * s)), HO2, O2m);
    R2_30 * (1e-3 * m3 / (mole * s)), HO2, O2m);
  reactionData->AddProduct(HO2m);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O(3p) -> O2- + OH
  G4double R2_31 = R1_6; 
  reactionData = new G4DNAMolecularReactionData(
     // 5.30e9 * (1e-3 * m3 / (mole * s)), HO2m, O);
   R2_31 * (1e-3 * m3 / (mole * s)), HO2m, O);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);

// Type IV //
  //------------------------------------------------------------------
  // e_aq + H3O+ -> H* + H2O
  G4double LogR4_1 = 39.127 - 3.888e4/t + 2.054e7/pow(t,2) - 4.899e9/pow(t,3) + 4.376e11/pow(t,4);
  G4double R4_1 = pow(10,LogR4_1); 
 // G4double R4_1 = pow(10, 39.127-3.888e4/t+2.054e7/(pow(t,2)) -4.899e9/(pow(t,3))+4.376e11/(pow(t,4)));
  reactionData = new G4DNAMolecularReactionData(
     // 2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
   R4_1 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(H);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O2- -> H2O2 + OH- + OH-
  G4double R4_2 = 2.46e12*pow(M_E,-1563.6/t);
  reactionData = new G4DNAMolecularReactionData(
   //   1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, O2m);
  R4_2 * (1e-3 * m3 / (mole * s)), e_aq, O2m);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + HO2- -> O- + OH-
  G4double R4_3 = 1.75e12*pow(M_E,-1852.7977/t);
  reactionData = new G4DNAMolecularReactionData(
    //  3.51e9 * (1e-3 * m3 / (mole * s)), e_aq, HO2m);
  R4_3 * (1e-3 * m3 / (mole * s)), e_aq, HO2m);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O- -> OH- + OH-
  G4double R4_4 = 5.6e11*pow(M_E,-951.64/t);
  reactionData = new G4DNAMolecularReactionData(
    //  2.31e10 * (1e-3 * m3 / (mole * s)), e_aq, Om);
   R4_4 * (1e-3 * m3 / (mole * s)), e_aq, Om);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O2- -> HO2
  G4double LogR4_5 = 16.410 - 4.888e3/t + 1.622e6/pow(t,2) - 2.004e8/pow(t,3);
  G4double R4_5 = pow(10,LogR4_5); 
 // G4double R4_5 = pow(10, 16.410-4.888e3/t+1.622e6/(pow(t,2))-2.004e8/(pow(t,3)));
  reactionData = new G4DNAMolecularReactionData(
    //  4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, O2m);
   R4_5* (1e-3 * m3 / (mole * s)), H3Op, O2m);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + HO2- -> H2O2
  G4double LogR4_6 = 16.410 - 4.888e3/t + 1.622e6/pow(t,2) - 2.004e8/pow(t,3);
  G4double R4_6 = pow(10,LogR4_6); 
  reactionData = new G4DNAMolecularReactionData(
     // 5.00e10 * (1e-3 * m3 / (mole * s)), H3Op, HO2m);
   R4_6 * (1e-3 * m3 / (mole * s)), H3Op, HO2m);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O- -> OH
  G4double LogR4_7 = 16.410 - 4.888e3/t + 1.622e6/pow(t,2) - 2.004e8/pow(t,3);
  G4double R4_7 = pow(10,LogR4_7);
  reactionData = new G4DNAMolecularReactionData(
    // 4.78e10  * (1e-3 * m3 / (mole * s)), H3Op, Om);
  R4_7 * (1e-3 * m3 / (mole * s)), H3Op, Om);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + O- -> O2 + OH- + OH-
  G4double R4_8 = R2_28; 
  //G4double R4_8 = 3.41e11*pow(M_E,-1344.9/t);
  reactionData = new G4DNAMolecularReactionData(
   //   6.00e8 * (1e-3 * m3 / (mole * s)), O2m, Om);
   R4_8 * (1e-3 * m3 / (mole * s)), O2m, Om);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O- -> O2- + OH-
  G4double R4_9 = 1.45e13*pow(M_E,-2928.5/t);
  reactionData = new G4DNAMolecularReactionData(
    //  3.50e8 * (1e-3 * m3 / (mole * s)), HO2m, Om);
   R4_9 * (1e-3 * m3 / (mole * s)), HO2m, Om);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O- -> H2O2 + OH- + OH-
  G4double R4_10 = R1_7; 
  reactionData = new G4DNAMolecularReactionData(
    //  1.00e8 * (1e-3 * m3 / (mole * s)), Om, Om);
  R4_10 * (1e-3 * m3 / (mole * s)), Om, Om);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O3- -> O2- + O2-
  G4double R4_11 = R2_28; 
  reactionData = new G4DNAMolecularReactionData(
   //   7.00e8 * (1e-3 * m3 / (mole * s)), Om, O3m);
  R4_11 * (1e-3 * m3 / (mole * s)), Om, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);

// Type VI
// pH is reflected in the system through type 6 reactions, whose reaction rates (Ksca) depend on their observed reaction rate (Kobs) and concentrations of OH- and H3O+  
  // First order reaction
  //------------------------------------------------------------------
  // O3- -> O- + O2
  G4double R6_1 = 3.20e11*pow(M_E,-5552.1/t);
  reactionData = new G4DNAMolecularReactionData(
    //  2.66e3 / s, O3m,None);
    R6_1 / s, O3m,None);
  reactionData->AddProduct(H3Op);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);

  // Scavenging reactions

  //------------------------------------------------------------------ 
  // HO2 + H2O -> H3O+ + O2-
  G4double pKHO2= 4.943-6.230e-3*T+4.125e-5*pow(T,2)-8.182e-9*pow(T,3);
  G4double KHO2=pow(10,-pKHO2);
  G4double R28f = R4_5*KHO2;
  G4double R6_2 = R28f*H2O_concentration;
  reactionData = new G4DNAMolecularReactionData(
    //  7.15e5 / s, HO2,H2OB);
   R6_2 / s, HO2,H2OB);
  reactionData->AddProduct(H3Op);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O -> eaq- + H3O+ 5.94 / s
  G4double pKH = 10.551 - 4.430e-2*T+1.902e-4*pow(T,2)-4.661e-7*pow(T,3)+5.980e-10*pow(T,4); 
  G4double KH = pow(10,-pKH);
  G4double R30f = R4_1*KH;
  G4double R6_3 = R30f*H2O_concentration; 
  reactionData = new G4DNAMolecularReactionData(
    //  5.94e0 / s, H,H2OB);
   R6_3 / s, H,H2OB);
  reactionData->AddProduct(e_aq);
  reactionData->AddProduct(H3Op);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H2O -> H + OH- 15.8 / s
  G4double pKH2O=16.690-4.262e-2*T+2.071e-4*pow(T,2)-5.594e-7*pow(T,3)+7.161e-10*pow(T,4);
  G4double KH2O=pow(10,-pKH2O);
  G4double R31b = (R2_3*KH2O)/KH;
  G4double R6_4 = R31b*H2O_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 1.58e1 / s, e_aq,H2OB);
   R6_4 / s, e_aq,H2OB);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H2O -> HO2 + OH- 0.15 / s
  G4double R6_5 = 2.266e10*pow(M_E, -2897.95134/t)*H2O_concentration;
  reactionData = new G4DNAMolecularReactionData(
    //  1.50e-1 / s, O2m,H2OB);
    R6_5 / s, O2m,H2OB);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H2O -> H2O2 + OH- 1.36e6 / s
  G4double pKH2O2 = 12.383-3.020e-2*T+1.700e-4*pow(T,2)-5.151e-7*pow(T,3)+6.960e-10*pow(T,4);
  G4double KH2O2 = pow(10, -pKH2O2);
  G4double R6_6 = ((R2_18*KH2O)/KH2O2)*H2O_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 1.36e6 / s, HO2m,H2OB);
   R6_6 / s, HO2m,H2OB);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + H2O -> OH + OH 1.90e3 / s
  G4double R6_7 = (18.4048*pow(M_E, -41.1/(R*t)))*H2O_concentration;     //Burns (1981) (Ea = 41.1)
  reactionData = new G4DNAMolecularReactionData(
    //  1.00e3 / s, O,H2OB);
   R6_7 / s, O,H2OB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H2O -> OH + OH- 1.36e6 / s
  G4double R6_8 = R6_6;
  reactionData = new G4DNAMolecularReactionData(
    //  1.36e6 / s, Om,H2OB);
    R6_8 / s, Om,H2OB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H3O+(B) -> H + H2O 2.09e3 / s
  G4double R6_9 = R4_1 * Hp_concentration;         //B: background  
  reactionData = new G4DNAMolecularReactionData(
   //   2.09e3 / s, e_aq,H3OpB);
    R6_9 / s, e_aq,H3OpB);
  reactionData->AddProduct(H);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s
  G4double R6_10 = R4_5 * Hp_concentration;
  reactionData = new G4DNAMolecularReactionData(
    //  4.73e3 / s, O2m,H3OpB);
    R6_10 / s, O2m,H3OpB);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + H3O+(B) -> 2H2O 1.11e4 / s
  G4double R6_11 = R3_2 * Hp_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 1.12e4 / s, OHm,H3OpB);
    R6_11 / s, OHm,H3OpB);
  theReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // H3O+ + OH-(B) -> 2H2O 1.11e4 / s
  // opposite description of OH- + H3O+(B) -> 2H2O
  G4double R6_12 = R3_2 * OHm_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 1.12e4 / s, H3Op,OHmB);
   R6_12 / s, H3Op,OHmB);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H3O+(B) -> H2O2 + H2O 4.98e3 / s
  G4double R6_13 = R4_6 * Hp_concentration;
  reactionData = new G4DNAMolecularReactionData(
    //  4.95e3 / s, HO2m,H3OpB);
   R6_13 / s, HO2m,H3OpB);
  reactionData->AddProduct(H2O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H3O+(B) -> OH + H2O 4.73e3 / s
  G4double R6_14 = R4_7 * Hp_concentration;
  reactionData = new G4DNAMolecularReactionData(
    //  4.73e3 / s, Om,H3OpB);
    R6_14 / s, Om,H3OpB);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O3- + H3O+(B) -> OH + O2 + H2O 8.91e3 / s
  G4double R6_15 = R3_3 * Hp_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 8.91e3 / s, O3m,H3OpB);
    R6_15 / s, O3m,H3OpB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH-(B) -> H2O + eaq- 2.49e3 / s
  G4double R6_16 = R2_3 * OHm_concentration;
  reactionData = new G4DNAMolecularReactionData(
    //  2.48e0 / s, H,OHmB);
   R6_16 / s, H,OHmB);
  reactionData->AddProduct(e_aq);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH-(B) -> O- + H2O 6.24e2 / s
  G4double R6_17 = R2_11 * OHm_concentration;
  reactionData = new G4DNAMolecularReactionData(
   //   6.24e2 / s, OH,OHmB);
    R6_17 / s, OH,OHmB);
  reactionData->AddProduct(Om);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s
  G4double R6_18 = R2_18 * OHm_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 4.66e1 / s, H2O2,OHmB);
    R6_18 / s, H2O2,OHmB);
  reactionData->AddProduct(HO2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + OH-(B) -> O2- + H2O 6.24e2 / s
  G4double R6_19 = R2_25 * OHm_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 6.24e2 / s, HO2,OHmB); 
    R6_19 / s, HO2,OHmB);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + OH-(B) -> HO2- 4.16e1 / s
  G4double R6_20 = R2_26 * OHm_concentration;
  reactionData = new G4DNAMolecularReactionData(
     // 4.16e1 / s, O,OHmB);
    R6_20 / s, O,OHmB);
  reactionData->AddProduct(HO2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------

}




void G4EmDNAChemistry_option3::ConstructReactionTable(G4DNAMolecularReactionTable*
                                              theReactionTable)
{
  auto model = G4EmParameters::Instance()->GetTimeStepModel();
  //-----------------------------------
  //Get the molecular configuration
  G4MolecularConfiguration* OH =
   G4MoleculeTable::Instance()->GetConfiguration("°OH");
  G4MolecularConfiguration* OHm =
   G4MoleculeTable::Instance()->GetConfiguration("OHm");
  G4MolecularConfiguration* e_aq =
   G4MoleculeTable::Instance()->GetConfiguration("e_aq");
  G4MolecularConfiguration* H2 =
   G4MoleculeTable::Instance()->GetConfiguration("H2");
  G4MolecularConfiguration* H3Op =
   G4MoleculeTable::Instance()->GetConfiguration("H3Op");
  G4MolecularConfiguration* H =
   G4MoleculeTable::Instance()->GetConfiguration("H");
  G4MolecularConfiguration* H2O2 =
   G4MoleculeTable::Instance()->GetConfiguration("H2O2");
  G4MolecularConfiguration* HO2 =
   G4MoleculeTable::Instance()->GetConfiguration("HO2°");
  G4MolecularConfiguration* HO2m =
   G4MoleculeTable::Instance()->GetConfiguration("HO2m");
  G4MolecularConfiguration* O =
   G4MoleculeTable::Instance()->GetConfiguration("Oxy");
  G4MolecularConfiguration* Om =
   G4MoleculeTable::Instance()->GetConfiguration("Om");
  G4MolecularConfiguration* O2 =
   G4MoleculeTable::Instance()->GetConfiguration("O2");
  G4MolecularConfiguration* O2m =
   G4MoleculeTable::Instance()->GetConfiguration("O2m");
  G4MolecularConfiguration* O3 =
   G4MoleculeTable::Instance()->GetConfiguration("O3");
  G4MolecularConfiguration* O3m =
   G4MoleculeTable::Instance()->GetConfiguration("O3m");

  G4MolecularConfiguration* H2OB =
   G4MoleculeTable::Instance()->GetConfiguration("H2O(B)");
  G4MolecularConfiguration* H3OpB =
   G4MoleculeTable::Instance()->GetConfiguration("H3Op(B)");
  G4MolecularConfiguration* OHmB =
   G4MoleculeTable::Instance()->GetConfiguration("OHm(B)");

  G4MolecularConfiguration* None =
   G4MoleculeTable::Instance()->GetConfiguration("NoneM");

  // Type I //
  //------------------------------------------------------------------
  // *H + *H -> H2
  G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(
      0.503e10 * (1e-3 * m3 / (mole * s)), H, H);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H* + H2O -> H2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      2.50e10 * (1e-3 * m3 / (mole * s)), e_aq, H);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);

  // H + O(3p) -> OH
  reactionData = new G4DNAMolecularReactionData(
      2.02e10 * (1e-3 * m3 / (mole * s)), H, O);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O- -> OH-
  reactionData = new G4DNAMolecularReactionData(
      2.00e10 * (1e-3 * m3 / (mole * s)), H, Om);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O(3p) -> HO2
  reactionData = new G4DNAMolecularReactionData(
      2.02e10 * (1e-3 * m3 / (mole * s)), OH, O);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O(3p) -> O2
  reactionData = new G4DNAMolecularReactionData(
      2.02e10 * (1e-3 * m3 / (mole * s)), HO2, O);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + O(3p) -> O2
  reactionData = new G4DNAMolecularReactionData(
      2.20e10 * (1e-3 * m3 / (mole * s)), O, O);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

  // Type III //
  //------------------------------------------------------------------
  // e_aq + e_aq + 2H2O -> H2 + 2OH-
  reactionData = new G4DNAMolecularReactionData(
      0.636e10 * (1e-3 * m3 / (mole * s)), e_aq, e_aq);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(H2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + OH- -> 2H2O
  reactionData = new G4DNAMolecularReactionData(
      1.13e11 * (1e-3 * m3 / (mole * s)), H3Op, OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O3- -> OH + O2
  reactionData = new G4DNAMolecularReactionData(
      9.0e10 * (1e-3 * m3 / (mole * s)), H3Op, O3m);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);

  // Type II //

  //------------------------------------------------------------------
  // *OH + *H -> H2O
  reactionData = new G4DNAMolecularReactionData(
      1.55e10 * (1e-3 * m3 / (mole * s)), OH, H);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O2 -> OH
  reactionData = new G4DNAMolecularReactionData(
      3.50e7 * (1e-3 * m3 / (mole * s)), H, H2O2);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH- -> eaq-
  reactionData = new G4DNAMolecularReactionData(
      2.51e7 * (1e-3 * m3 / (mole * s)), H, OHm);
  reactionData->AddProduct(e_aq);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2 -> HO2
  reactionData = new G4DNAMolecularReactionData(
      2.10e10 * (1e-3 * m3 / (mole * s)), H, O2);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + HO2 -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      1.00e10 * (1e-3 * m3 / (mole * s)), H, HO2);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + O2- -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      1.00e10 * (1e-3 * m3 / (mole * s)), H, O2m);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // *OH + *OH -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      0.55e10 * (1e-3 * m3 / (mole * s)), OH, OH);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2O2 -> HO2
  reactionData = new G4DNAMolecularReactionData(
      2.88e7 * (1e-3 * m3 / (mole * s)), OH, H2O2);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + H2 -> H
  reactionData = new G4DNAMolecularReactionData(
      3.28e7 * (1e-3 * m3 / (mole * s)), OH, H2);
  reactionData->AddProduct(H);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + *OH -> OH-
  reactionData = new G4DNAMolecularReactionData(
      2.95e10 * (1e-3 * m3 / (mole * s)), e_aq, OH);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH- -> O-
  reactionData = new G4DNAMolecularReactionData(
      6.30e9 * (1e-3 * m3 / (mole * s)), OH, OHm);
  reactionData->AddProduct(Om);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2 -> O2
  reactionData = new G4DNAMolecularReactionData(
      7.90e9 * (1e-3 * m3 / (mole * s)), OH, HO2);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O2- -> O2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.07e10 * (1e-3 * m3 / (mole * s)), OH, O2m);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + HO2- -> HO2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      8.32e9 * (1e-3 * m3 / (mole * s)), OH, HO2m);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O- -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      1.00e9 * (1e-3 * m3 / (mole * s)), OH, Om);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + O3- -> O2- + HO2
  reactionData = new G4DNAMolecularReactionData(
      8.50e9 * (1e-3 * m3 / (mole * s)), OH, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + H2O2 -> OH- + *OH
  reactionData = new G4DNAMolecularReactionData(
      1.10e10 * (1e-3 * m3 / (mole * s)), e_aq, H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH- -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      4.71e8 * (1e-3 * m3 / (mole * s)), H2O2, OHm);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O(3p) -> HO2 + OH
  reactionData = new G4DNAMolecularReactionData(
      1.60e9 * (1e-3 * m3 / (mole * s)), H2O2, O);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + O- -> HO2 + OH-
  reactionData = new G4DNAMolecularReactionData(
      5.55e8 * (1e-3 * m3 / (mole * s)), H2O2, Om);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O(3p) -> H + OH
  reactionData = new G4DNAMolecularReactionData(
      4.77e3 * (1e-3 * m3 / (mole * s)), H2, O);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2 + O- -> H + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.21e8 * (1e-3 * m3 / (mole * s)), H2, Om);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + O2 -> O2-
  reactionData = new G4DNAMolecularReactionData(
      1.74e10 * (1e-3 * m3 / (mole * s)), e_aq, O2);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq + HO2 -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, HO2);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + HO2 -> O2-
  reactionData = new G4DNAMolecularReactionData(
      6.30e9 * (1e-3 * m3 / (mole * s)), OHm, HO2);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + O(3p) -> HO2-
  reactionData = new G4DNAMolecularReactionData(
      4.20e8 * (1e-3 * m3 / (mole * s)), OHm, O);
  reactionData->AddProduct(HO2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O(3p) -> O3
  reactionData = new G4DNAMolecularReactionData(
      4.00e9 * (1e-3 * m3 / (mole * s)), O2, O);
  reactionData->AddProduct(O3);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2 + O- -> O3-
  reactionData = new G4DNAMolecularReactionData(
      3.70e9 * (1e-3 * m3 / (mole * s)), O2, Om);
  reactionData->AddProduct(O3m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + HO2 -> H2O2 + O2
  reactionData = new G4DNAMolecularReactionData(
      9.80e5 * (1e-3 * m3 / (mole * s)), HO2, HO2);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + O2- -> HO2- + O2
  reactionData = new G4DNAMolecularReactionData(
      9.70e7 * (1e-3 * m3 / (mole * s)), HO2, O2m);
  reactionData->AddProduct(HO2m);
  reactionData->AddProduct(O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O(3p) -> O2- + OH
  reactionData = new G4DNAMolecularReactionData(
      5.30e9 * (1e-3 * m3 / (mole * s)), HO2m, O);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);

  // Type IV //
  //------------------------------------------------------------------
  // e_aq + H3O+ -> H* + H2O
  reactionData = new G4DNAMolecularReactionData(
      2.11e10 * (1e-3 * m3 / (mole * s)), e_aq, H3Op);
  reactionData->AddProduct(H);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O2- -> H2O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.29e10 * (1e-3 * m3 / (mole * s)), e_aq, O2m);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + HO2- -> O- + OH-
  reactionData = new G4DNAMolecularReactionData(
      3.51e9 * (1e-3 * m3 / (mole * s)), e_aq, HO2m);
  reactionData->AddProduct(Om);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // e_aq + O- -> OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      2.31e10 * (1e-3 * m3 / (mole * s)), e_aq, Om);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O2- -> HO2
  reactionData = new G4DNAMolecularReactionData(
      4.78e10 * (1e-3 * m3 / (mole * s)), H3Op, O2m);
  reactionData->AddProduct(HO2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + HO2- -> H2O2
  reactionData = new G4DNAMolecularReactionData(
      5.00e10 * (1e-3 * m3 / (mole * s)), H3Op, HO2m);
  reactionData->AddProduct(H2O2);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H3O+ + O- -> OH
  reactionData = new G4DNAMolecularReactionData(
     4.78e10  * (1e-3 * m3 / (mole * s)), H3Op, Om);
  reactionData->AddProduct(OH);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + O- -> O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      6.00e8 * (1e-3 * m3 / (mole * s)), O2m, Om);
  reactionData->AddProduct(O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + O- -> O2- + OH-
  reactionData = new G4DNAMolecularReactionData(
      3.50e8 * (1e-3 * m3 / (mole * s)), HO2m, Om);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O- -> H2O2 + OH- + OH-
  reactionData = new G4DNAMolecularReactionData(
      1.00e8 * (1e-3 * m3 / (mole * s)), Om, Om);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  reactionData->AddProduct(OHm);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + O3- -> O2- + O2-
  reactionData = new G4DNAMolecularReactionData(
      7.00e8 * (1e-3 * m3 / (mole * s)), Om, O3m);
  reactionData->AddProduct(O2m);
  reactionData->AddProduct(O2m);
  SetReactionType(reactionData,model);//partially diffusion-controlled
  theReactionTable->SetReaction(reactionData);

  // Type VI
  // First order reaction
  //------------------------------------------------------------------
  // O3- -> O- + O2
  reactionData = new G4DNAMolecularReactionData(
      2.66e3 / s, O3m,None);
  reactionData->AddProduct(H3Op);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);

  // Scavenging reactions

  //------------------------------------------------------------------
  // HO2 + H2O -> H3O+ + O2-
  reactionData = new G4DNAMolecularReactionData(
      7.15e5 / s, HO2,H2OB);
  reactionData->AddProduct(H3Op);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + H2O -> eaq- + H3O+ 5.94 / s
  reactionData = new G4DNAMolecularReactionData(
      5.94e0 / s, H,H2OB);
  reactionData->AddProduct(e_aq);
  reactionData->AddProduct(H3Op);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H2O -> H + OH- 15.8 / s
  reactionData = new G4DNAMolecularReactionData(
      1.58e1 / s, e_aq,H2OB);
  reactionData->AddProduct(H);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H2O -> HO2 + OH- 0.15 / s
  reactionData = new G4DNAMolecularReactionData(
      1.50e-1 / s, O2m,H2OB);
  reactionData->AddProduct(HO2);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H2O -> H2O2 + OH- 1.36e6 / s
  reactionData = new G4DNAMolecularReactionData(
      1.36e6 / s, HO2m,H2OB);
  reactionData->AddProduct(H2O2);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + H2O -> OH + OH 1.90e3 / s
  reactionData = new G4DNAMolecularReactionData(
      1.00e3 / s, O,H2OB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H2O -> OH + OH- 1.36e6 / s
  reactionData = new G4DNAMolecularReactionData(
      1.36e6 / s, Om,H2OB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(OHm);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // eaq- + H3O+(B) -> H + H2O 2.09e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.09e3 / s, e_aq,H3OpB);
  reactionData->AddProduct(H);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O2- + H3O+(B) -> HO2 + H2O 4.73e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.73e3 / s, O2m,H3OpB);
  reactionData->AddProduct(HO2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH- + H3O+(B) -> 2H2O 1.11e4 / s
  reactionData = new G4DNAMolecularReactionData(
      1.12e4 / s, OHm,H3OpB);
  theReactionTable->SetReaction(reactionData);

  //------------------------------------------------------------------
  // H3O+ + OH-(B) -> 2H2O 1.11e4 / s
  // opposite description of OH- + H3O+(B) -> 2H2O
  reactionData = new G4DNAMolecularReactionData(
      1.12e4 / s, H3Op,OHmB);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2- + H3O+(B) -> H2O2 + H2O 4.98e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.95e3 / s, HO2m,H3OpB);
  reactionData->AddProduct(H2O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O- + H3O+(B) -> OH + H2O 4.73e3 / s
  reactionData = new G4DNAMolecularReactionData(
      4.73e3 / s, Om,H3OpB);
  reactionData->AddProduct(OH);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O3- + H3O+(B) -> OH + O2 + H2O 8.91e3 / s
  reactionData = new G4DNAMolecularReactionData(
      8.91e3 / s, O3m,H3OpB);
  reactionData->AddProduct(OH);
  reactionData->AddProduct(O2);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H + OH-(B) -> H2O + eaq- 2.49e3 / s
  reactionData = new G4DNAMolecularReactionData(
      2.48e0 / s, H,OHmB);
  reactionData->AddProduct(e_aq);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // OH + OH-(B) -> O- + H2O 6.24e2 / s
  reactionData = new G4DNAMolecularReactionData(
      6.24e2 / s, OH,OHmB);
  reactionData->AddProduct(Om);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // H2O2 + OH-(B) -> HO2- + H2O 4.66e2 / s
  reactionData = new G4DNAMolecularReactionData(
      4.66e1 / s, H2O2,OHmB);
  reactionData->AddProduct(HO2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // HO2 + OH-(B) -> O2- + H2O 6.24e2 / s
  reactionData = new G4DNAMolecularReactionData(
      6.24e2 / s, HO2,OHmB);
  reactionData->AddProduct(O2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------
  // O(3p) + OH-(B) -> HO2- 4.16e1 / s
  reactionData = new G4DNAMolecularReactionData(
      4.16e1 / s, O,OHmB);
  reactionData->AddProduct(HO2m);
  theReactionTable->SetReaction(reactionData);
  //------------------------------------------------------------------

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructProcess()
{
  auto ChemModel = G4EmParameters::Instance()->GetTimeStepModel();
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //===============================================================
  // Extend vibrational to low energy
  // Anyway, solvation of electrons is taken into account from 7.4 eV
  // So below this threshold, for now, no accurate modeling is done
  //
  G4VProcess* process =
      G4ProcessTable::GetProcessTable()->
        FindProcess("e-_G4DNAVibExcitation", "e-");

  if (process)
  {
    G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*) process;
    G4VEmModel* model = vibExcitation->EmModel();
    G4DNASancheExcitationModel* sancheExcitationMod =
        dynamic_cast<G4DNASancheExcitationModel*>(model);
    if(sancheExcitationMod)
    {
      sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
    }
  }

  //===============================================================
  // *** Electron Solvatation ***
  //
  process =
  G4ProcessTable::GetProcessTable()->
  FindProcess("e-_G4DNAElectronSolvation", "e-");
  
  if (process == 0)
  {
    ph->RegisterProcess(
        new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"),
        G4Electron::Definition());
  }


  //===============================================================
  // Define processes for molecules
  //
  G4MoleculeTable* theMoleculeTable = G4MoleculeTable::Instance();
  G4MoleculeDefinitionIterator iterator =
      theMoleculeTable->GetDefintionIterator();
  iterator.reset();
  while (iterator())
  {
    G4MoleculeDefinition* moleculeDef = iterator.value();

    if (moleculeDef != G4H2O::Definition())
    {
      if(ChemModel != G4ChemTimeStepModel::IRT)
      {
        auto* brown = new G4DNABrownianTransportation();
        ph->RegisterProcess(brown, moleculeDef);
      }
    }
    else
    {
      moleculeDef->GetProcessManager()
                      ->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
      auto* dissociationProcess =
          new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
      dissociationProcess->SetDisplacer(
          moleculeDef, new G4DNAWaterDissociationDisplacer);
      dissociationProcess->SetVerboseLevel(3);

      moleculeDef->GetProcessManager()
                ->AddRestProcess(dissociationProcess, 1);
    }
    /*
     * Warning : end of particles and processes are needed by
     * EM Physics builders
     */
  }

  G4DNAChemistryManager::Instance()->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3::ConstructTimeStepModel(G4DNAMolecularReactionTable*
                                              /*reactionTable*/)
{
  auto model = G4EmParameters::Instance()->GetTimeStepModel();
  if(model == G4ChemTimeStepModel::IRT)
  {
    RegisterTimeStepModel(new G4DNAMolecularIRTModel(), 0);
  }else if(model == G4ChemTimeStepModel::SBS)
  {
    RegisterTimeStepModel(new G4DNAMolecularStepByStepModel(), 0);
  }else if(model == G4ChemTimeStepModel::IRT_syn)
  {
    RegisterTimeStepModel(new G4DNAIndependentReactionTimeModel(), 0);
  }
}

void G4EmDNAChemistry_option3::SetReactionType(G4DNAMolecularReactionData* pData,
                                               G4ChemTimeStepModel model)
{
  if(model != G4ChemTimeStepModel::SBS) { pData->SetReactionType(1); }
}
