# Geant4DNA_v11.3.0_temppH

This project is a modified version of Geant4-DNA v11.3.0  adapted to simulate water radiolysis under different temperature and pH conditions.

## Acknowledgement

This repository is a modification of the original open-source project [Geant4DNA](http://geant4-dna.org)

The base code and structure are developed by the Geant4DNA collaboration and licensed under [original license name](http://cern.ch/geant4/license ).

## Features

- Updated G4EmDNAChemistry_option3 with the addition of temperature- and pH- dependency
- 
- Temperature and pH values can be set in chem6.cc


## Detailed description of the modification

/geant4_11.3.0_source/examples/extended/medical/dna/chem6/chem6.cc

-Creation of pH and temperature variables (Line 71-73)
                  G4double pH = 7;    
                  G4double T = 25;  //Celsius
                  T = T + 273.5;   //Kelvin          

-Instances of classes PhysicsList and DetectorConstruction are called with modified constructors that take as input T and pH. (Line 87-88)

                   runManager->SetUserInitialization(new PhysicsList(T,pH));
                   runManager->SetUserInitialization(new DetectorConstruction(T));

/geant4_11.3.0_source/examples/extended/medical/dna/chem6/beam.in

Comment out:
Line 23 :# /chem/reaction/UI
Line 25 - 36


/geant4_11.3.0_source/examples/extended/medical/dna/chem6/include/PhysicsList.hh

-#include "G4ChemDissociationChannels_option1.hh"   (Line 50)

-Definition of  constructor that takes temperature and pH as input. (Line 57)

           explicit PhysicsList(G4double T, G4double P);
  
-Definition of variables Temperature and pH. (Line 69-70)

          G4double Temperature;
          G4double pH;   //Definition of variables Temperature and pH

/geant4_11.3.0_source/examples/extended/medical/dna/chem6/src/PhysicsList.cc

-Implementation of constructor that takes temperature and pH as input and initializes values of variables temperature and pH. (Line 62)

PhysicsList::PhysicsList(G4double T, G4double P) : G4VModularPhysicsList(), Temperature(T), pH(P)

-Initialization of chemistry_option3 with contructor that takes temperature and pH as input, (Line 163)


fEmDNAChemistryList = std::unique_ptr<G4EmDNAChemistry_option3>(new G4EmDNAChemistry_option3(Temperature, pH)); 

/geant4_11.3.0_source/examples/extended/medical/dna/chem6/include/DetectorConstruction.hh

-Definition of constructor that takes temperature as input (Line 57)

             DetectorConstruction(G4double T);

-Definition of variable Temperature. (Line 65)

             G4double Temperature;

/geant4_11.3.0_source/examples/extended/medical/dna/chem6/include/DetectorConstruction.cc

-Implementation of a constructor that takes temperature as input and initializes variable temperature. (Line 63)

 DetectorConstruction::DetectorConstruction(G4double T) : G4VUserDetectorConstruction(), Temperature(T) {}

-Temperature-dependent density polynomial to create water material at different temperatures. (Line 75-77)

  double t = Temperature - 273.15;  //t in Celsius
  double rho = 0.999 + 1.094e-4*t - 7.397e-6*pow(t,2) + 2.693e-8*pow(t,3) - 4.714e-11*pow(t,4);
  double density = rho * g/cm3; 

G4Material* water =man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED","G4_WATER",density);

/geant4_11.3.0_source/source/physics_lists/constructors/electromagnetic/include/G4EmDNAChemistry_option3.hh

- #include "G4ChemDissociationChannels_option1.hh" (Line 41)

-class G4EmDNAChemistry_option3 : public G4VUserChemistryList, public G4VPhysicsConstructor, public G4ChemDissociationChannels_option1 //Inherit from G4EmChemDissociationChannels_option1 class (Line 48)

-Definition of constructor that takes temperature and pH as input (Line 54)

           G4EmDNAChemistry_option3(G4double T, G4double P); 

-Definition of void ConstructReactionTablePhTemp (Line 67)

 void ConstructReactionTablePhTemp(G4DNAMolecularReactionTable* reactionTable, double T, double pH) override;

-Definition of fDissociationChannels to help pass temperature to ConstructMolecule method (Line 70)

G4ChemDissociationChannels_option1 fDissociationChannel;

-Definition of variables Temperature and pH. (Line 71-72)

    G4double Temperature;
    G4double pH;


/geant4_11.3.0_source/source/physics_lists/constructors/electromagnetic/include/G4EmDNAChemistry_option3.cc

-#define _USE_MATH_DEFINES
 #include <cmath>   (Line 36-37)

-Implentation of constructor that takes temperature and pH as input, initializes variables of temperature and pH (Line 78-80)

G4EmDNAChemistry_option3::G4EmDNAChemistry_option3(G4double T, G4double P) :
 G4VUserChemistryList(true), fDissociationChannel(T), pH(P),Temperature(T)

 -Calls methods “setTemp(Temperature)” and “setpH(ph)” which initializes the aforementioned values in the G4DNAChemistryManger class. (Line 83-84)

  G4DNAChemistryManager::Instance()->setTemp(T);
  G4DNAChemistryManager::Instance()->setpH(pH);

-Pass temperature to ConstructMolecule method in G4ChemDissociationChannels_option1 for temperature dependent diffusion coefficient. (Line 94-98)

void G4EmDNAChemistry_option3::ConstructMolecule()
{
  fDissociationChannel.ConstructMolecule(Temperature);}

-Implementation of ConstructReactionTablePhTemp(G4DNAMolecularReactionTable* reactionTable, double T, double pH), which has temperature -dependent polynomials for observed reaction rates (Kobs) of all reactions, and concentrations of water, hydroxil ion (OH-) and hydronium ion (H3O+)  which vary depending by the input pH. pH is reflected in the system through type 6 reactions, whose reaction rates (Ksca) depend on their observed reaction rate (Kobs) and concentrations of OH- and H3O+. (Line 109-844)

-Call setTemp method to set temperature in  reactionData class (Line 177)
       reactionData->setTemp(t);   

/geant4_11.3.0_source/source/physics_lists/constructors/electromagnetic/include/G4ChemDissociationChannels_option1.hh

 - G4ChemDissociationChannels_option1(); (Line38)

 - Definition of the method that takes temperature as input (Line41)
static void ConstructMolecule(double Temperature);  

-Definition of Temperature(Line44)
  double Temperature;

-Constructor to initialize Temperature. (Line42)
  G4ChemDissociationChannels_option1(double initialTemperature);  

/geant4_11.3.0_source/source/physics_lists/constructors/electromagnetic/include/G4ChemistryDissociationChannels_option1.cc

-#include <cmath> //Multiple powers (Line 31)
 #include "G4EmDNAChemistry_option3.hh" (Line 33)
 #include "G4DNAChemistryManager.hh" (Line 34)

-Constructor to initialize Temperature.(Line 56-57)
G4ChemDissociationChannels_option1::G4ChemDissociationChannels_option1(double initialTemperature)  : Temperature(initialTemperature) {}

-To keep the constructor without parameters G4ChemDissociationChannels_option().(Line 59-63)
G4ChemDissociationChannels_option1::G4ChemDissociationChannels_option1()
    : Temperature(300.0)
{ G4cout << "Default constructor called. Temperature = " << Temperature << G4endl;
}                                                                                                

-Instantiation of all molecules with their appropriate temperature -dependent diffusion coefficient value. This is all within the implementation of the method constructMolecule (Line 59-249)

/geant4_11.3.0_source/source/processes/electromagnetic/dna/utils/include/G4VUserChemistryList.hh

-Definition of virtual void ConstructReactionTablePhTemp (Line 91)

  virtual void ConstructReactionTablePhTemp(G4DNAMolecularReactionTable* /*reactiontable*/, double /*T*/, double /*P*/) {} 

/geant4_11.3.0_source/source/processes/electromagnetic/dna/utils/include/G4DNAChemistryManager.hh

-Define two methods void setTemp(G4double) and setpH(G4double). (Line 116-117)
    void setTemp(G4double);
    void setpH(G4double);

-Define variables pH and Temperature. (Line 223-224)
    G4double Temperature;
    G4double pH;  

/geant4_11.3.0_source/source/processes/electromagnetic/dna/utils/src/G4DNAChemistryManager.cc

-Constructor initializes Temperature and pH with 0 (but their value is set by the user when the methods setTemp and setPh is called in G4EmDNAChemistry_option3 constructor). (Line 116-117)
    Temperature = 0;
    pH = 0; 

-Implementation of methods setTemp and setPh. (Line 196-205)
void G4DNAChemistryManager::setTemp(G4double T)
{
    Temperature = T;
}
void G4DNAChemistryManager::setpH(G4double P)
{
    pH = P;
}     

-Call to method ConstructReactionTablePhTemp. (Line 426)

fpUserChemistryList->ConstructReactionTablePhTemp(G4DNAMolecularReactionTable::GetReactionTable(),Temperature,pH);


/geant4_11.3.0_source/source/processes/electromagnetic/dna/utils/include/G4DNAMolecularReactionTable.hh

-Define setter and getter for temperature.  (Line 148-149)
    void setTemp(G4double T);
    G4double getTemp() const; 

-Definition inline static G4double Temp = 273.15;  (Line 156)

/geant4_11.3.0_source/source/processes/electromagnetic/dna/utils/src/G4DNAMolecularReactionTable.cc

-Implementation for setter and getter for temperature. (Line 114-121)

void G4DNAMolecularReactionData::setTemp(G4double T){
    Temp = T;
}                                            
G4double G4DNAMolecularReactionData::getTemp() const{
    return Temp;
}      

-Initialization of onsager radius with temperature dependent values inside
void G4DNAMolecularReactionData::ComputeEffectiveRadius (Line 141-145)

 G4double t = G4DNAMolecularReactionData::getTemp(); 
 G4double epsilon_H2O = 5321*pow(t,-1) + 233.76 - 0.9297*t +0.1417e-2*pow(t,2) - 0.8292e-6*pow(t,3);
fReactionRadius = fEffectiveReactionRadius;
fOnsagerRadius = (fpReactant1->GetCharge() * fpReactant2->GetCharge())/(4*pi*epsilon0*k_Boltzmann) / (t * epsilon_H2O);

