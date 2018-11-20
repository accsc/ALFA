/*
 *
 *	ALFA.cpp
 *
 * 	This program is free software; you can redistribute it and/or modify
 * 	it under the terms of the GNU General Public License as published by
 * 	the Free Software Foundation version 2 of the License.
 *
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 	GNU General Public License for more details.
 *
 */

/*
 *	
 *      v5.0 - June/2016.
 *             Improved library of torsionals (TorLib2 - Guba et al. 2016. JCIM
 *             Library generator mode added for multi-MOL2 inputs
 *             Bug fixing (nasty segfaults, etc.)
 *
 *	v4.3 - Published version 2014. Centro de Biologia Molecular Severo Ochoa.
 *	       Improved Empirical rules.
 *
 *	v4.0 - Initial OpenBabel port
*/
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <openbabel/mol.h>
#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/typer.h>
#include <openbabel/rand.h>
#include <openbabel/forcefield.h>
#include <openbabel/atom.h>
#include <openbabel/parsmart.h>
#include <openbabel/obiter.h>
#include <openbabel/conformersearch.h>
#include <openbabel/math/align.h>
#include <openbabel/canon.h>

#include "Constants.h"
#include "Rules.h"
#include "RotatableBond.h"
#include "TorsionRule.h"

#include "USR_lib.h"

using namespace OpenBabel;
using namespace std;

struct conformer_struct {
  unsigned long long conf;
  unsigned long long id;
  //vector<unsigned int> id;
  double coord[MAX_ATOM];
  double energy;
  double rmsd;
  double moments[NUM_MOMENTS];
  double meanDist;
};typedef struct conformer_struct tConformer;

//! Functions headers
void getParameters(int piArgc, char **piArgv);
void help(void);
unsigned long long string2ulonglong(string piInput);
vector < OBMol > readInputMolecule(string piInputFile, vector<unsigned int>& poFileID2molID, vector<unsigned int>& poMolID2FileID);
void getFileVsMolRelationIDs(OBMol& piMol, vector<unsigned int>& poFileID2molID, vector<unsigned int>& poMolID2fileID);
OBMol readMolecule(string piInputFile);
void loadTorsionRules(vector<TorsionRule> *torsionRulesList);
string simplifyWhitesAndTabs(char * inStr);
TorsionRule readTorsionRule(string piLine);
string simplifyWhitesAndTabs(char * inStr);
vector<RotatableBond> getRotatableBonds(OBMol& piMol, vector<TorsionRule> piTorsionRules);
void addOriginalAngles(OBMol& piMol, vector<RotatableBond>& piRotatableBonds);
// void cleanEquivalentAngles(OBMol& piMol, vector<RotatableBond>& piRotatableBonds);
void cleanEquivalentAngles(OBMol& piMol, vector<RotatableBond>& piRotatableBonds);
vector<unsigned int> getAmberTypes(OBMol& piMol);
vector< vector<unsigned int> > getMinimalDistances(OBMol& piMol);
vector< vector<unsigned int> > getTestAtoms4vdw(OBMol& piMol, const vector< vector<unsigned int> >& piMinimalDistances, vector<RotatableBond> rotBonds);
vector<unsigned long long> getSizeByIndexPosition(vector<RotatableBond>& piRotatableBonds);
list<tConformer> generateConformers(OBMol& piMol, vector<RotatableBond>& piRotatableBonds, const vector< vector<unsigned int> >& piMinimalDistances, const vector< vector<unsigned int> >& piTestAtoms4vdw, const vector<unsigned int>& piAmberTypes);
unsigned int stericClashFilter(const vector<unsigned int>& piAtomIdxs, const double* piCoordinates, const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw, vector< vector<float> >& poDistances);
double getVDWEnergy(const vector<unsigned int>& piAtomIdxs, const double* piCoordinates, const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw, const vector<vector<unsigned int> >& piMinimalDistances,vector< vector<float> > & piDistances);
bool isRigid(OBMol& piMol, vector<RotatableBond>& piRotatableBonds);
vector<unsigned int> getNextConformerMCSA(const vector<unsigned int>& piLastConformer, vector<RotatableBond>& piRotatableBonds);
vector<unsigned int> getIDFromIdConformer(unsigned long long piIdConformer, vector<RotatableBond>& piRotatableBonds);
vector<unsigned int> getNextConformer(const vector<unsigned int>& piLastConformer, vector<RotatableBond>& piRotatableBonds);
unsigned long long getConformerIndex(vector<unsigned int> piConformer);
double getScore(const vector<unsigned int>& piAtomIdxs, OBMol& piMol, double* piCoordinates, const vector< vector<unsigned int> >& piMinimalDistances,const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw);
list< tConformer > getRefinedConformerList(vector < list < tConformer > > piConfs, vector < OBMol > piMol, const vector< vector<unsigned int> >& piMinimalDistances,const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw);
int getSillasBotes(OBMol& piMol, int& poNumSillas);
void writeMol2Output(OBMol& piMol, list<tConformer>& piConformers);
string ulonglong2string(unsigned long long piInput);
string double2string(double piInput);
void writePDBOutput(OBMol& piMol, list<tConformer>& piConformers);
void writeXMLOutput(list<tConformer>& piConformers, vector< vector<unsigned int> > piAmberTypes, vector< vector<RotatableBond> > piRotatableBonds);
list< tConformer > getLibraryConformerList(list < tConformer > piConfs, OBMol piMol, const vector< vector<unsigned int> >& piMinimalDistances,const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw);


//! Command line parameters
string addRulesFile;
double cutOff;
unsigned int howManySelect;
string inputFile;
unsigned long long maxCombinations;
string outputFile;
string outputType;
string referenceFile;
string useRulesFile;
unsigned int clustering;
unsigned int EnergySelection;
unsigned int outMini;

unsigned long long minSteps;
string minMethod;
double minEconv;

unsigned long long useInputInRMSD;

//! Useful global variables
double bestGeneratedRMSD = 1000.0;
double bestGeneratedRMSDEnergy = 1000.0;
double bestSelectedRMSD = 1000.0;
double bestSelectedRMSDEnergy = 1000.0;
vector<unsigned int> fileID2molID;
unsigned int generatedConformers = 0;
vector<unsigned int> molID2fileID;
unsigned long long possibleConformers  = 1; // Because we have the input conformer
OBMol referenceMolecule;
vector<unsigned long long> sizeByIndexPosition;
unsigned int numConfig = 0;

double dist_cutoff = 0.80f;


int main(int argc, char** argv){
	
  //! Variables
  vector<unsigned int> amberTypesTMP;
  vector< vector<unsigned int> > amberTypes;
  vector< list<tConformer> > conformers;
  list<tConformer> conformersRefined;
  vector< vector<unsigned int> > minimalDistances;
  vector < OBMol > molecule;
  OBMol moleculeTMP1;
  OBMol moleculeTMP2;
  vector<RotatableBond> rotatableBondsTMP;
  vector< vector<RotatableBond> > rotatableBonds;
  string smiles;
  vector< vector<unsigned int> > testAtoms4vdw;
  OBStopwatch timer;
  vector<TorsionRule> torsionRules;
  srand( time(NULL) );
	
  //! Start timer
  timer.Start();

  //! Get command line parameters
  getParameters(argc, argv);


  cout << endl << endl << LOGO1 ;
  cout << endl << LOGO2 ;
  cout << endl << LOGO3 ;
  cout << endl << LOGO4 ;
  cout << endl << LOGO5 ;
  cout << endl << LOGO6 ;
  cout << endl << LOGO7 << endl << endl;

  cout << "2012-2016 Universidad de Alcala. Alvaro Cortes and Federico Gago" << endl;
  cout << "2008-2012 Centro de Biologia Molecular Severo Ochoa. Javier Klett and Ruben Gil" << endl;

  cout << "- Version: " << VERSION << endl << endl;
  cout << "- Parameters" << endl;
  cout << "\taddRulesFile: " << addRulesFile << endl;
  cout << "\tcutOff: " << cutOff << endl;
  cout << "\thowManySelect: " << howManySelect << endl;
  cout << "\tinputFile: " << inputFile << endl;
  cout << "\tmaxCombinations: " << maxCombinations << endl;
  cout << "\toutputFile: " << outputFile << endl;
  cout << "\toutputType: " << outputType << endl;
  cout << "\tuseRulesFile: " << useRulesFile << endl << endl;

  //! Read molecule
  molecule = readInputMolecule(inputFile, fileID2molID, molID2fileID); 

  cout << "- Total input molecules: " << molecule.size() << endl << endl;

  //! Test reference file
  if (referenceFile != "") {
    referenceMolecule = readMolecule(referenceFile);
    moleculeTMP1 = molecule.at(0);
    moleculeTMP2 = referenceMolecule;
    moleculeTMP1.DeleteHydrogens();
    moleculeTMP2.DeleteHydrogens();
    
    if (moleculeTMP1.NumAtoms() != moleculeTMP2.NumAtoms()) {
      cout << "ERROR: Different number of heavy-atoms in the input and reference molecule." << endl;
      cout << "#Atoms Input: " << moleculeTMP1.NumAtoms() << " #Atoms Ref: " << moleculeTMP2.NumAtoms() << endl;
      exit(-1);
    }
  }
  
  //! Load torsion rules
  loadTorsionRules(&torsionRules);
#ifdef _DEBUG
  cout << "Number of loaded torsion rules: " << torsionRules.size() << endl;
#endif

//   if(clustering != ""){
  if(clustering == 1)
  {
    for(unsigned int k=0;k<molecule.size();k++)
    {
      numConfig++;
      cout << "Processing input molecule " << numConfig << endl;

      //! Assign rotatable bonds
      rotatableBondsTMP = getRotatableBonds(molecule.at(k), torsionRules);
      rotatableBonds.push_back(rotatableBondsTMP);

      if( rotatableBondsTMP.size() >= 12)
	dist_cutoff = 0.95f;
      else if( rotatableBondsTMP.size() == 11)
        dist_cutoff = 0.90f;
      else if( rotatableBondsTMP.size() == 10)
        dist_cutoff = 0.90f;
      else if( rotatableBondsTMP.size() == 9)
        dist_cutoff = 0.85f;
      else if( rotatableBondsTMP.size() == 8)
        dist_cutoff = 0.85f;
      else if( rotatableBondsTMP.size() == 7)
        dist_cutoff = 0.80f;
      else if( rotatableBondsTMP.size() <= 6)
        dist_cutoff = 0.75f;
      else 
        dist_cutoff = 0.85f;

      //! Add original angles
      if (useInputInRMSD == 1) addOriginalAngles(molecule.at(k), rotatableBonds.at(k));

      //! Clean equivalent angles
      cleanEquivalentAngles(molecule.at(k), rotatableBonds.at(k));

      //! Assign AMBER types
      amberTypesTMP = getAmberTypes(molecule.at(k));
      amberTypes.push_back(amberTypesTMP);

      //! Get minimal distances
      minimalDistances = getMinimalDistances(molecule.at(k));

      //! Get test atoms for vdw
      testAtoms4vdw = getTestAtoms4vdw(molecule.at(k), minimalDistances, rotatableBonds.at(k));

      //! Get size by index position
      sizeByIndexPosition = getSizeByIndexPosition(rotatableBonds.at(k));
      vector<unsigned int> atomIdxs(molecule.at(k).NumAtoms(), 0);
      unsigned int i = 0;
      for( OBMolAtomIter atomsIter(molecule.at(k)); atomsIter; ++atomsIter){
	atomIdxs.at(i) = atomsIter->GetIdx();
	i++;
      }
      //! Get conformer list for current input molecule
      conformers.push_back(generateConformers(molecule.at(k), rotatableBonds.at(k), minimalDistances, testAtoms4vdw, amberTypes.at(k)));
    }
    conformersRefined = getRefinedConformerList(conformers, molecule,  minimalDistances, amberTypes.at(0), testAtoms4vdw); // be aware that we are passing amberTypes.at(0) this may leed to uncorrect vdw ene if amberTypes changes for different input configurations.
    //getRefinedConformerList(conformers, molecule,  minimalDistances, amberTypes.at(0), testAtoms4vdw); // be aware that we are passing amberTypes.at(0) this may leed to uncorrect vdw ene if amberTypes changes for different input configurations.
  }else{

    for(unsigned int k=0;k<molecule.size();k++)
    {
      numConfig++;
      cout << "****** Processing input molecule " << numConfig << " ******" << endl;

	//! Assign rotatable bonds
	rotatableBondsTMP = getRotatableBonds(molecule.at(k), torsionRules);
	rotatableBonds.push_back(rotatableBondsTMP);

	//! Add original angles
	if (useInputInRMSD == 1) {
	  addOriginalAngles(molecule.at(k), rotatableBonds.at(k));
	}

	//! Clean equivalent angles
	cleanEquivalentAngles(molecule.at(k), rotatableBonds.at(k));

	//! Assign AMBER types
	amberTypesTMP = getAmberTypes(molecule.at(k));
	amberTypes.push_back(amberTypesTMP);

	//! Get minimal distances
	minimalDistances = getMinimalDistances(molecule.at(k));

	//! Get test atoms for vdw
	testAtoms4vdw = getTestAtoms4vdw(molecule.at(k), minimalDistances,rotatableBonds.at(k));

	//! Get size by index position
	sizeByIndexPosition = getSizeByIndexPosition(rotatableBonds.at(k));

	vector<unsigned int> atomIdxs(molecule.at(k).NumAtoms(), 0);
	unsigned int i = 0;
	for( OBMolAtomIter atomsIter(molecule.at(k)); atomsIter; ++atomsIter){
	  atomIdxs.at(i) = atomsIter->GetIdx();
	  i++;
	}

        if( rotatableBondsTMP.size() >= 12)
          dist_cutoff = 0.95f;
        else if( rotatableBondsTMP.size() == 11)
          dist_cutoff = 0.90f;
        else if( rotatableBondsTMP.size() == 10)
          dist_cutoff = 0.90f;
        else if( rotatableBondsTMP.size() == 9)
          dist_cutoff = 0.85f;
        else if( rotatableBondsTMP.size() == 8)
          dist_cutoff = 0.85f;
        else if( rotatableBondsTMP.size() == 7)
          dist_cutoff = 0.80f;
        else if( rotatableBondsTMP.size() <= 6)
          dist_cutoff = 0.75f;
        else
          dist_cutoff = 0.85f;

       
//        if( rotatableBonds.at(k).size() > 7)
//        {
//                cout << "Rotable bonds larger than 7. Skipping" << endl;
//                list<tConformer> empty;
//		conformers.push_back( empty );
 //               continue;
  //      }


	conformers.push_back(generateConformers(molecule[k], rotatableBonds.at(k), minimalDistances, testAtoms4vdw, amberTypes.at(k)));
        conformersRefined = getLibraryConformerList(conformers.at(k), molecule[k], minimalDistances, amberTypes.at(k), testAtoms4vdw);

	cout << "- Number of atoms: " << molecule.at(k).NumAtoms() << endl;
	cout << "- Rotatable bonds: " << rotatableBonds.at(k).size() << endl << endl;
	unsigned int index = 0;
	for (vector<RotatableBond>::iterator itRotBonds = rotatableBonds.at(k).begin(); itRotBonds != rotatableBonds.at(k).end(); itRotBonds++) {
	  index++;
	  cout << "\t" << index << " (" << itRotBonds->getType() << "); atoms: ";
	  cout << molID2fileID[itRotBonds->getAtoms()[0]->GetIdx()] << " ";
	  cout << molID2fileID[itRotBonds->getAtoms()[1]->GetIdx()] << " ";
	  cout << molID2fileID[itRotBonds->getAtoms()[2]->GetIdx()] << " ";
	  cout << molID2fileID[itRotBonds->getAtoms()[3]->GetIdx()];
	  cout << "; angles:";
	  for (unsigned int i = 0; i < itRotBonds->getAngles().size(); i++) {
	    cout << " " << itRotBonds->getAngles()[i];
	  }
	  cout << endl;
	}
	cout << endl;
	cout << "- Possible conformers: " << possibleConformers << endl;
	cout << "- Generated conformers: " << generatedConformers << endl;
        cout << "- Accepted conformers: " << conformersRefined.size() << endl << endl;
        writeMol2Output(molecule.at(k), conformersRefined);
    }
  }

  //! Write outupt
  if(clustering != 0){
	if (outputType == "mol2" || outputType == "mol2-xml") {
	  writeMol2Output(molecule.at(0), conformersRefined);
	}
	if (outputType == "pdb"){
	  writePDBOutput(molecule.at(0), conformersRefined);
	}
	if (outputType == "xml" || outputType == "mol2-xml") {
	  writeXMLOutput(conformersRefined, amberTypes, rotatableBonds);
	}
/*  }else{
	if (outputType == "mol2" || outputType == "mol2-xml") {
	  writeMol2Output(molecule.at(0), conformers.at(0));
	}
	if (outputType == "pdb"){
	  writePDBOutput(molecule.at(0), conformers.at(0));
	}
	if (outputType == "xml" || outputType == "mol2-xml") {
	  writeXMLOutput(conformers.at(0), amberTypes, rotatableBonds);
	}*/
  }

  double seconds = timer.Elapsed();

//  cout << "- Number of atoms: " << molecule.at(0).NumAtoms() << endl;
//  cout << "- Rotatable bonds: " << rotatableBonds.at(0).size() << endl << endl;
//  unsigned int index = 0;
//  for (vector<RotatableBond>::iterator itRotBonds = rotatableBonds.at(0).begin(); itRotBonds != rotatableBonds.at(0).end(); itRotBonds++) {
//    index++;
//    cout << "\t" << index << " (" << itRotBonds->getType() << "); atoms: ";
//    cout << molID2fileID[itRotBonds->getAtoms()[0]->GetIdx()] << " ";
//    cout << molID2fileID[itRotBonds->getAtoms()[1]->GetIdx()] << " ";
//    cout << molID2fileID[itRotBonds->getAtoms()[2]->GetIdx()] << " ";
//    cout << molID2fileID[itRotBonds->getAtoms()[3]->GetIdx()];
//    cout << "; angles:";
//    for (unsigned int i = 0; i < itRotBonds->getAngles().size(); i++) {
//      cout << " " << itRotBonds->getAngles()[i];
//    }
//    cout << endl;
//  }
//  cout << endl;
//  cout << "- Possible conformers: " << possibleConformers << endl;
//  cout << "- Generated conformers: " << generatedConformers << endl;
  if (referenceFile != "") cout << "- Best generated RMSD: " << bestGeneratedRMSD << " (energy = " << bestGeneratedRMSDEnergy << ")" << endl;
  if(clustering != 0){
	cout << "- Accepted conformers: " << conformersRefined.size() << endl;
//  }else{
//	cout << "- Accepted conformers: " << conformers.at(0).size() << endl;
  }
  if (referenceFile != "") cout << "- Best selected RMSD: " << bestSelectedRMSD << " (energy = " << bestSelectedRMSDEnergy << ")" << endl;
  //cout << "- Time per conformer: " << (seconds / generatedConformers) * 1000 << " milliseconds" << endl;
  cout << "- Total time: " << seconds << " seconds" << endl << endl;


  exit(0);
}

void getParameters(int piArgc, char **piArgv) {

  int i;

  if (piArgc == 1)help();

  //! Setting up default values
  howManySelect = 20;
  cutOff = 0.0;
  maxCombinations = 300000;
  addRulesFile  = "";
  referenceFile = "";
  useRulesFile  = "";
  referenceFile = "";
  outputType = "mol2";
  useInputInRMSD = 1;
  outMini = 0;
  clustering = 1;
  EnergySelection = 0;
  minSteps = 1000;
  minMethod = "conjugate";
  minEconv = 1.0e-6;
  
  for(i = 1; i < piArgc; i++){

    //! Getting input parameters
    if (strncmp(piArgv[i],"-help",5)==0){
		    help();
    }else if (strncmp(piArgv[i],"-inputFile",10) == 0){
	    inputFile = piArgv[i+1];
	    //printf("INPUTFILE >> %i %s\n",i,inputFile.c_str());
    }else if (strncmp(piArgv[i],"-outputFile",11) == 0){
	    outputFile = piArgv[i+1];
	    //printf("OUTPUTFILE >> %i %s\n",i,outputFile.c_str());
    }else if (strncmp(piArgv[i],"-howManySelect",14) == 0){
	    howManySelect = (unsigned int) atoi(piArgv[i+1]);
	    //printf("HOWMANYSELECT >> %i %i\n",i, howManySelect);
    }else if (strncmp(piArgv[i],"-cutOff", 7) == 0){
	    //printf("CUTOFF\n");
	    cutOff = atof(piArgv[i+1]);
    }else if (strncmp(piArgv[i],"-maxCombinations",16) == 0){
	    //printf("MAXCOMBINATIONS\n");
	    maxCombinations = string2ulonglong(piArgv[i+1]);
    }else if (strncmp(piArgv[i],"-useInputInRMSD",15) == 0){
	    //printf("USEINPUTINRMSD\n");
	    useInputInRMSD = (unsigned int)atof(piArgv[i+1]);
    }else if (strncmp(piArgv[i],"-outputType",11) == 0){
	    //printf("OUTPUTTYPE\n");
	    outputType = piArgv[i+1];
    }else if (strncmp(piArgv[i],"-addRulesFile",13) == 0){
	    addRulesFile = piArgv[i+1];
	    //printf("ADDRULESFILE\n");
    }else if (strncmp(piArgv[i],"-referenceFile",14) == 0){
	    //printf("REFERENCEFILE\n");
	    referenceFile = piArgv[i+1];
    }else if (strncmp(piArgv[i],"-useRulesFile",13) == 0){
	    //printf("USERULESFILE\n");
	    useRulesFile = piArgv[i+1];
    }else if (strncmp(piArgv[i],"-energySelection",11) == 0){
	    EnergySelection = 1;
    }else if (strncmp(piArgv[i],"-clustering",11) == 0){
	    clustering = 1;
    }else if(strncmp(piArgv[i],"-minimize",8)==0){
	    outMini = 1;
    }else if(strncmp(piArgv[i],"-minMethod",8)==0){
	    minMethod = piArgv[i+1];
    }else if(strncmp(piArgv[i],"-minSteps",8)==0){
	    minSteps =  string2ulonglong(piArgv[i+1]);
    }else if(strncmp(piArgv[i],"-h",2)==0){
	    help();
    }
  }

  if ( EnergySelection == 1){
    clustering = 0;
  }

  return;
}

void help(void){
  fprintf(stderr, "\nALFA - %s \n\n", VERSION.c_str());
  fprintf(stderr, "Another program to generate conformers based on rules\n");
  fprintf(stderr, "2012 - 2016 Universidad de Alcala. Alvaro Cortes <alvarocortesc@gmail.com> and Federico Gago.\n");
  fprintf(stderr, "2008 - 2012 Centro de Biologia Molecular Severo Ochoa. Javier Klett and Ruben Gil.\n\n");
  fprintf(stderr, "Complete parameter list:\n");
  fprintf(stderr, "   -addRulesFile : File containing torsion rules in order to add they (replacing if needed) to the default torsion rules\n");
  fprintf(stderr, "   -cutOff : Cut off for the energies in the final list (allows values only of minEnergy + cutOff). -- DEFAULT: 0.0 (no cutOff)\n");
  fprintf(stderr, "   -howManySelect : Maximum number of selected conformers. -- DEFAULT: 100\n");
  fprintf(stderr, "   -inputFile : Name of the mol2 input file\n");
  fprintf(stderr, "   -maxCombinations : Maximum number of generated conformers. -- DEFAULT: 300000\n");
  fprintf(stderr, "   -outputFile : Name of the output file (without extension)\n");
  fprintf(stderr, "   -outputType : Type of the output file. -- DEFAULT: mol2\n");
  fprintf(stderr, "   -referenceFile : Name of the reference file in order to perform RMSD calculations\n");
  fprintf(stderr, "   -useInputInRMSD : Say if input structure must be taken into account for RMSD calculations. Be careful, if you activate this option then ALFA rules\n");
  fprintf(stderr, "                            can include the input molecule angles. Allowed values: 0 (no), 1 (yes). -- DEFAULT: 1\n");
  fprintf(stderr, "   -useRulesFile : File containing torsion rules in order to use they (not taking into account default torsion rules\n");
  fprintf(stderr, "   -clustering : conformers are chosen by similirity criteria -- DEFAULT\n");
  fprintf(stderr, "   -energySelection : conformers are chosen by minimum energy criteria\n");
  fprintf(stderr, "   -minimize : Minimize output conformations\n");
  fprintf(stderr, "   -minMethod : Allowed values: conjugate / steepest -- DEFAULT conjugate\n");
  fprintf(stderr, "   -minSteps : Number of steeps in minimization-- DEFAULT 1000\n");
  fprintf(stderr, "   -minEconv : Energy convergence criteria.-- DEFAULT 1e-6\n");

  exit(1);
}

unsigned long long string2ulonglong(string piInput) {
  //if (!isNumber(piInput)) {
  //throw VSDBException(piInput + " is not a number");
  //}
  char* pEnd;
  return strtoull(piInput.c_str(), &pEnd, 10);
}

vector < OBMol > readInputMolecule(string piInputFile, vector<unsigned int>& poFileID2molID, vector<unsigned int>& poMolID2FileID) {
  vector < OBMol > mol;
  OBMol * molTMP;
  ifstream ifs(piInputFile.c_str());

  if(!ifs)
  {
    cout << "Error. Cannot open input file. Aborting..." << endl;
    exit(-1);
  }

	OBConversion conv;
	OBFormat* inFormat = conv.FormatFromExt(piInputFile.c_str());
	conv.SetInFormat(inFormat);

#ifdef _DEBUG
	cout << "Calculating conformers for: " << piInputFile << endl;
#endif
	
  if ( strcmp(inFormat->GetID(),"sy2") == 0 ) {
#ifdef _DEBUG
		cout << "MOL2 file" << endl;
#endif
		conv.SetInFormat("MOL2");
    for( molTMP = new OBMol ; conv.Read(molTMP,&ifs); molTMP = new OBMol ){
			molTMP->FindRingAtomsAndBonds();
			OBAromaticTyper aromTyper;
			aromTyper.AssignAromaticFlags(*molTMP);
			OBAtomTyper atomTyper;
			atomTyper.AssignHyb(*molTMP);
			mol.push_back(*molTMP);
		}
  } else if (strcmp(inFormat->GetID(),"ent") == 0) {
#ifdef _DEBUG
		cout << "PDB file" << endl;
#endif
		conv.SetInFormat("PDB");
    for( molTMP = new OBMol ; conv.Read(molTMP,&ifs); molTMP = new OBMol ){
	    molTMP->ConnectTheDots();
	    molTMP->FindRingAtomsAndBonds();
	    molTMP->PerceiveBondOrders();
	    OBAromaticTyper aromTyper;
	    aromTyper.AssignAromaticFlags(*molTMP);
	    OBAtomTyper atomTyper;
	    atomTyper.AssignHyb(*molTMP);
	    mol.push_back(*molTMP);
		}
		
  } else {
    cout << "ERROR: unknow type for file " << piInputFile << endl;
    exit(-1);
  }

	ifs.close();

#ifdef _DEBUG
		cout << "Num. Input conformations " << mol.size() << endl;
		cout << "Num. Atoms " << mol[0].NumAtoms() << endl;
#endif
  
  getFileVsMolRelationIDs(mol[0], poFileID2molID, poMolID2FileID);

  return mol;
}

void getFileVsMolRelationIDs(OBMol& piMol, vector<unsigned int>& poFileID2molID, vector<unsigned int>& poMolID2fileID) {
  poFileID2molID.reserve(piMol.NumAtoms() + 1);
  poFileID2molID.assign(0, piMol.NumAtoms() + 1);
  poMolID2fileID.reserve(piMol.NumAtoms() + 1);
  poMolID2fileID.assign(0, piMol.NumAtoms() + 1);
  OBMolAtomIter atom(piMol);
  unsigned int index;
  for(index = 1; atom; index++, ++atom ){
    poFileID2molID[index] = atom->GetIdx();
    poMolID2fileID[atom->GetIdx()] = index;
  }
}

OBMol readMolecule(string piInputFile) {
  
  OBMol mol;
  ifstream ifs(piInputFile.c_str());
  OBConversion conv;
	OBFormat* inFormat = conv.FormatFromExt(piInputFile.c_str());
	conv.SetInFormat(inFormat);	
	
#ifdef _DEBUG
	cout << "Referece file: " << piInputFile << endl;
#endif
  if (strcmp(inFormat->GetID(),"sy2") == 0) {
#ifdef _DEBUG
		cout << "MOL2 file" << endl;
#endif
		conv.SetInFormat("MOL2");
		conv.Read(&mol,&ifs);
		mol.FindRingAtomsAndBonds();
		OBAromaticTyper aromTyper;
		aromTyper.AssignAromaticFlags(mol);
		OBAtomTyper atomTyper;
		atomTyper.AssignHyb(mol);		
  } else if (strcmp(inFormat->GetID(),"ent") == 0) {
#ifdef _DEBUG
		cout << "PDB file" << endl;
#endif
		conv.SetInFormat("PDB");
		conv.Read(&mol,&ifs);		
		mol.ConnectTheDots();
		mol.FindRingAtomsAndBonds();
	  mol.PerceiveBondOrders();
		OBAromaticTyper aromTyper;
		aromTyper.AssignAromaticFlags(mol);
		OBAtomTyper atomTyper;
		atomTyper.AssignHyb(mol);
  } else {
    cout << "ERROR: unknow type for file " << piInputFile << endl;
    exit(-1);
  }
 	ifs.close();
  return mol;
}

void loadTorsionRules(vector<TorsionRule> *torsionRulesList) {
  
  ifstream ifs;
  string line; 

  if (useRulesFile != "") {
    ifs.open(useRulesFile.c_str());
    while (getline(ifs, line)) {
      line = simplifyWhitesAndTabs((char *)line.c_str());
      if ((line.find_first_of("#") == 0) || (line.length() == 0)) {
	continue; // Commentary or blank line
      }
      readTorsionRule(line);
      torsionRulesList->push_back(readTorsionRule(line));
    }
    ifs.close();
  } else {
    torsionRulesList->reserve(NDTR);
    for (unsigned int i = 0; i < NDTR; i++) {
      list<string> stLineDTR;
      torsionRulesList->push_back(readTorsionRule( simplifyWhitesAndTabs((char *)DTR[i].c_str())));
    }
  }
  if (addRulesFile != "") {
    ifs.open(addRulesFile.c_str());
    while (getline(ifs, line)) {
      line = simplifyWhitesAndTabs((char *)line.c_str());
      if ((line.find_first_of("#") == 0) || (line.length() == 0)) {
	continue;
      }
      torsionRulesList->push_back(readTorsionRule(line));
    }
    ifs.close();
  }
  return;
}

string simplifyWhitesAndTabs(char * inStr){
  string outStr;
  char * pch;
  pch = strtok (inStr," \t");
  while (pch != NULL){
    outStr.append(pch);
    outStr.append(" ");
    pch = strtok (NULL, " \t");
  }
  return outStr;

}

list<string> separateByWhites(char * inStr){
  list<string> outStr;
  char * pch;
  pch = strtok (inStr," \t");
  while (pch != NULL){
    outStr.push_back(pch);
    pch = strtok (NULL, " \t");
  }
  return outStr;
}

TorsionRule readTorsionRule(string piLine) {
  list<string> stLine;
  stLine = separateByWhites((char *)piLine.c_str());
  list<string>::iterator itList;
  unsigned int index;
  string name;
  string pattern;
  int point;
  vector<unsigned int> points;
  points.reserve(4);
  float angle;
  vector<double> angles;
  for (itList = stLine.begin(), index = 1; itList != stLine.end(); itList++, index++) {
    switch (index) {
    case 1: name = *itList; break;
    case 2: pattern = *itList; break;
    case 3:case 4:case 5:case 6:
      sscanf(itList->c_str(), "%d", &point );
      points.push_back(point);
      break;
    default:
      sscanf(itList->c_str(), "%f", &angle );
      angles.push_back((double)angle);
      break;
    }
  }
  TorsionRule torsionRule(name, pattern, points, angles);
  return torsionRule;
}

vector<RotatableBond> getRotatableBonds(OBMol& piMol, vector<TorsionRule> piTorsionRules) {
  vector<RotatableBond> rotatableBonds;
	//! Take rotatable bonds
  for( OBMolBondIter iterRotBonds(piMol); iterRotBonds; ++iterRotBonds ){
    if(iterRotBonds->IsRotor()){
      vector<double> angles;
      vector<OBAtom*> atoms;
      OBAtom *atom1, *atom2, *atom3, *atom4;
      //! main atoms
      atom2 = iterRotBonds->GetBeginAtom();
      atom3 = iterRotBonds->GetEndAtom();
      //! random references
      for( OBAtomAtomIter it(atom2); it; ++it ){
        if(it->GetIdx() != atom2->GetIdx()){ atom1 = &(*it); break;}
      }
      for( OBAtomAtomIter it(atom3); it; ++it ){
        if(it->GetIdx() != atom3->GetIdx() && it->GetIdx() != atom2->GetIdx()){ atom4 = &(*it); break; }
      }
      atoms.push_back(atom1); atoms.push_back(atom2); atoms.push_back(atom3); atoms.push_back(atom4);
        
      RotatableBond rotatableBond("UNKNOWN", atoms, angles);
      rotatableBonds.push_back(rotatableBond);
    }
  }

  //! Assign types, references and angles
  for (unsigned int i = 0; i < piTorsionRules.size(); i++) {
    TorsionRule torsionRule = piTorsionRules.at(i);

    OBSmartsPattern smartPatter;
    smartPatter.Init(torsionRule.getPattern().c_str());
    smartPatter.Match(piMol);
    
    if(smartPatter.NumMatches() == 0) continue;
    vector<vector<int> > maplist;
    maplist = smartPatter.GetMapList();
    
    for (unsigned int j = 0; j < rotatableBonds.size(); j++) {
      vector<vector<int> >::iterator i;
      for (i = maplist.begin();i != maplist.end();++i){
	vector<OBAtom*> points;
	for (unsigned int k = 0; k < 4; k++) points.push_back(piMol.GetAtom(i->at(k)));
	if ((rotatableBonds.at(j).getAtoms().at(1)->GetIdx() == points.at(1)->GetIdx()) && (rotatableBonds.at(j).getAtoms().at(2)->GetIdx() == points.at(2)->GetIdx())) {
	  rotatableBonds.at(j).updateAtom(0, points.at(0));
	  rotatableBonds.at(j).updateAtom(3, points.at(3));
	  rotatableBonds.at(j).setType(torsionRule.getName());
	  rotatableBonds.at(j).setAngles(torsionRule.getAngles());
	  break;
	}else if ((rotatableBonds.at(j).getAtoms().at(1)->GetIdx() == points.at(2)->GetIdx()) && (rotatableBonds.at(j).getAtoms().at(2)->GetIdx() == points.at(1)->GetIdx())) {
	  rotatableBonds.at(j).updateAtom(0, points.at(3));
	  rotatableBonds.at(j).updateAtom(3, points.at(0));
	  rotatableBonds.at(j).setType(torsionRule.getName());
	  rotatableBonds.at(j).setAngles(torsionRule.getAngles());
	  break;
	}
      }
    }
  }
  //! Clear not rotatable bonds (torsionals without angles)
  vector<RotatableBond> tmpRotatableBonds = rotatableBonds;
  rotatableBonds.clear();
  for (unsigned int i = 0; i < tmpRotatableBonds.size(); i++) {
    if (tmpRotatableBonds[i].getAngles().size() != 0) {
      rotatableBonds.push_back(tmpRotatableBonds[i]);
    }
  }
#ifdef _DEBUG_ANG
  for (unsigned int h = 0; h < rotatableBonds.size(); h++){
    cout << rotatableBonds.at(h).getType() << " --> ";
    cout << rotatableBonds.at(h).getAtoms().at(0)->GetIdx() << " ";
    cout << rotatableBonds.at(h).getAtoms().at(1)->GetIdx() << " ";
    cout << rotatableBonds.at(h).getAtoms().at(2)->GetIdx() << " ";
    cout << rotatableBonds.at(h).getAtoms().at(3)->GetIdx() << endl;
    for (unsigned int g = 0; g < rotatableBonds.at(h).getAngles().size(); g++){
      cout << rotatableBonds.at(h).getAngles().at(g) << " ";
    }
    cout << endl;
  }
#endif
      
  return rotatableBonds;
}

void addOriginalAngles(OBMol& piMol, vector<RotatableBond>& piRotatableBonds) {
  
  for (unsigned int i = 0; i < piRotatableBonds.size(); i++) {
    double originalAngle = piMol.GetTorsion(piRotatableBonds[i].getAtoms()[0], piRotatableBonds[i].getAtoms()[1], piRotatableBonds[i].getAtoms()[2], piRotatableBonds[i].getAtoms()[3]) ;
    unsigned int howManyAngles = piRotatableBonds[i].getAngles().size();
    for (int j = (howManyAngles - 1); j >= 0; j--) {
      int diffAngle = (int) rint(piRotatableBonds[i].getAngles()[j] - originalAngle);
      if (diffAngle > 180) {
	diffAngle = 360 - diffAngle;
      }
      if (abs(diffAngle) < 30) { //! Remove similar angle to the original
	piRotatableBonds[i].removeAngle(j);
      }
    }
    piRotatableBonds[i].addAngle(rint(originalAngle)); //! Add original angle
#ifdef _DEBUG_ANG    
    cout << (int)originalAngle << " " ;
#endif
  }
#ifdef _DEBUG_ANG
  cout << "<--ANGLES IN ORIG. MOL." << endl;
#endif
}

void cleanEquivalentAngles(OBMol& piMol, vector<RotatableBond>& piRotatableBonds) {

  OBAlign *rmsd = new OBAlign(false,true);
  OBMol referenceMol_1 = piMol;
  OBMol referenceMol_2 = piMol;
  double DELETED_VALUE = 1000;
  double *origCoord = new double[piMol.NumAtoms() * 3];
  double *origCoordTmp = piMol.GetCoordinates(); //! Save the original coordinates
  
  for(unsigned int i = 0; i < piMol.NumAtoms() * 3; i = i+3){
    origCoord[i  ] = origCoordTmp[i  ];
    origCoord[i+1] = origCoordTmp[i+1];
    origCoord[i+2] = origCoordTmp[i+2];
  }
  for (unsigned int i = 0; i < piRotatableBonds.size(); i++) {
    vector<double> newAngles;
    vector<double> angles = piRotatableBonds[i].getAngles();
    vector<OBAtom*> atoms = piRotatableBonds[i].getAtoms();
#ifdef _DEBUG_ANG
    cout << piRotatableBonds[i].getType() << " Atms: ";
    cout << " " << piRotatableBonds.at(0).getAtoms().at(0)->GetIdx();
    cout << " " << piRotatableBonds.at(0).getAtoms().at(1)->GetIdx();
    cout << " " << piRotatableBonds.at(0).getAtoms().at(2)->GetIdx();
    cout << " " << piRotatableBonds.at(0).getAtoms().at(3)->GetIdx() << endl;
#endif
    for (unsigned int j = 0; j < angles.size() - 1; j++) {
      if (angles[j] == DELETED_VALUE) continue;
      referenceMol_1.SetCoordinates(origCoord);
      OBAtom *a_1 = referenceMol_1.GetAtom(piRotatableBonds.at(0).getAtoms().at(0)->GetIdx());
      OBAtom *b_1 = referenceMol_1.GetAtom(piRotatableBonds.at(0).getAtoms().at(1)->GetIdx());
      OBAtom *c_1 = referenceMol_1.GetAtom(piRotatableBonds.at(0).getAtoms().at(2)->GetIdx());
      OBAtom *d_1 = referenceMol_1.GetAtom(piRotatableBonds.at(0).getAtoms().at(3)->GetIdx());
      referenceMol_1.SetTorsion(a_1,b_1,c_1,d_1, angles[j] * DEG_TO_RAD);
#ifdef _DEBUG_ANG
      cout  <<j << "  " <<angles[j] << "\tRef_1 "<< referenceMol_1.GetTorsion(a_1, b_1, c_1, d_1) << endl;
#endif
      for (unsigned int k = j + 1; k < angles.size(); k++) {
	if (angles[k] == DELETED_VALUE) continue;
	referenceMol_2.SetCoordinates(origCoord);
	OBAtom *a_2 = referenceMol_2.GetAtom(piRotatableBonds.at(0).getAtoms().at(0)->GetIdx());
	OBAtom *b_2 = referenceMol_2.GetAtom(piRotatableBonds.at(0).getAtoms().at(1)->GetIdx());
	OBAtom *c_2 = referenceMol_2.GetAtom(piRotatableBonds.at(0).getAtoms().at(2)->GetIdx());
	OBAtom *d_2 = referenceMol_2.GetAtom(piRotatableBonds.at(0).getAtoms().at(3)->GetIdx());
	referenceMol_2.SetTorsion(a_2,b_2,c_2,d_2, angles[k] * DEG_TO_RAD);
#ifdef _DEBUG_ANG
	cout << k<< "  " <<angles[k] << "\tRef_2 " << referenceMol_2.GetTorsion(a_2, b_2, c_2, d_2) << endl;
#endif
	rmsd->SetRefMol(referenceMol_1);
	rmsd->SetTargetMol(referenceMol_2);
	rmsd->Align();
	rmsd->UpdateCoords(&referenceMol_2);
	if (rmsd->GetRMSD () < DELTA) {
#ifdef _DEBUG_ANG
	  cout << "DELETED " << rmsd->GetRMSD() << angles[k] << endl;
#endif
	  angles[k] = DELETED_VALUE;
	}
      }
    }
    for (unsigned int j = 0; j < angles.size(); j++) {
      if (angles[j] != DELETED_VALUE) newAngles.push_back(angles[j]);
    }
    piRotatableBonds[i].setAngles(newAngles);
  }
  piMol.SetCoordinates(origCoord); //! Restore original coordinates, may be not needed
  return;
}

vector<unsigned int> getAmberTypes(OBMol& piMol) {

  vector<unsigned int> result;
  result.assign(piMol.NumAtoms() + 1, 0);
  for( OBMolAtomIter atomsIter(piMol); atomsIter; ++atomsIter ){
    if( atomsIter->IsHydrogen() ){
#ifdef _DEBUG_AMBTYPES			
			cout << "Atom " << atomsIter->GetIdx() << " is Hyd? " << atomsIter->IsHydrogen() << endl;
#endif			
			continue;
		}
    if (atomsIter->GetAtomicNum() == 6) { // carbon
      if (atomsIter->GetHyb() <= 2) { // sp2 o sp1
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Carbon sp2 ó sp1. AtomicNum = " << atomsIter->GetAtomicNum() << " Hybridation = " << atomsIter->GetHyb()<< endl;
#endif
				result[atomsIter->GetIdx()] = 2;
      } else if (atomsIter->GetHyb() >= 3) { //sp3
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Carbon sp3. AtomicNum = " << atomsIter->GetAtomicNum() << " Hybridation = " << atomsIter->GetHyb()<< endl;
#endif
				result[atomsIter->GetIdx()] = 1;
      }
#ifdef _DEBUG_AMBTYPES
      cout << "Vecinos de carbonos  -Z " ;
 #endif     
      //! Assing nprotonate state and nheteroatoms
      bool nprotonate = false;
      unsigned int nheteroatoms = 0;      
      for (OBAtomAtomIter neighborAtomsIter( *atomsIter ); neighborAtomsIter; ++neighborAtomsIter ) {
#ifdef _DEBUG_AMBTYPES
	cout << neighborAtomsIter->GetIdx() << " - " ;
#endif
	if(neighborAtomsIter->IsCarbon() ||  neighborAtomsIter->IsHydrogen() || neighborAtomsIter->IsPhosphorus() ) continue;
	if (neighborAtomsIter->IsNitrogen()) {
		if (neighborAtomsIter->GetValence() > 3) {
			nprotonate = true;
		} else {
			nheteroatoms++;
		}
	} else {
		nheteroatoms++;
	}
      }
#ifdef _DEBUG_AMBTYPES
      cout << endl;
#endif
      //! Assign H type
      unsigned int hType = 8; //! Generic
      if (atomsIter->GetHyb() <= 2) {
	if (nheteroatoms == 0) {
		hType = 9;
	} else if (nheteroatoms == 1) {
		hType = 20;
	} else if (nheteroatoms == 2) {
		hType = 21;
	}
      } else if (atomsIter->GetHyb() >= 3) {
	if (nprotonate) {
		hType = 22;
	} else if (nheteroatoms == 0) {
		hType = 8;
	} else if (nheteroatoms == 1) {
		hType = 17;
	} else if (nheteroatoms == 2) {
		hType = 18;
	} else if (nheteroatoms == 3) {
		hType = 19;
	}
      }
      for (OBAtomAtomIter neighborHAtomsIter( *atomsIter ); neighborHAtomsIter; ++neighborHAtomsIter ) {
       if (!neighborHAtomsIter->IsHydrogen()) continue;
      	result[ neighborHAtomsIter->GetIdx()] = hType;
#ifdef _DEBUG_AMBTYPES
	cout <<  "Hyd: " << neighborHAtomsIter->GetIdx() << " Tipo: " << hType << endl;
#endif
      }
    }	else if (atomsIter->IsOxygen() ) {
      result[atomsIter->GetIdx()] = 5; //! Default O
      if (atomsIter->GetValence() - atomsIter->GetHvyValence() >= 1) {
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Oxygen. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
				result[atomsIter->GetIdx()] = 4;
      } else if (atomsIter->GetHvyValence() == 2) {
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Oxygen sp?. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
				result[atomsIter->GetIdx()] = 6;
      }
      //! Assign H types
      for (OBAtomAtomIter neighborHAtomsIter( *atomsIter ); neighborHAtomsIter; ++neighborHAtomsIter ) {
				if (!neighborHAtomsIter->IsHydrogen()) continue;
				result[neighborHAtomsIter->GetIdx()] = 11;
			}
    }
    else if ( atomsIter->IsNitrogen() ) {
      result[atomsIter->GetIdx()] = 23; //! Default N
      if (atomsIter->GetValence() <= 2) {
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Nitrogen. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
				result[atomsIter->GetIdx()] = 3;
      }
      //! Assign H types
      for (OBAtomAtomIter neighborHAtomsIter( *atomsIter ); neighborHAtomsIter; ++neighborHAtomsIter ) {
				if (!neighborHAtomsIter->IsHydrogen()) continue;
					result[neighborHAtomsIter->GetIdx()] = 10;
      }
    } else if ( atomsIter->IsSulfur() ) {
      result[atomsIter->GetIdx()] = 7;
      //! Assign H types
      for (OBAtomAtomIter neighborHAtomsIter( *atomsIter ); neighborHAtomsIter; ++neighborHAtomsIter ) {
				if (!neighborHAtomsIter->IsHydrogen()) continue;
				result[neighborHAtomsIter->GetIdx()] = 10;
			}
    } else if (atomsIter->IsPhosphorus()) {
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Phosphorus. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
      result[atomsIter->GetIdx()] = 12;
    } else if (atomsIter->GetAtomicNum() == 53) { //Iodo
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Phosphorus. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
      result[atomsIter->GetIdx()] = 13;
    } else if (atomsIter->GetAtomicNum() == 35) { // Br
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Br. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
      result[atomsIter->GetIdx()] = 14;
    } else if (atomsIter->GetAtomicNum() == 17) {
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is Cl. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif			
      result[atomsIter->GetIdx()] = 15;
    } else if (atomsIter->GetAtomicNum() == 9) {
#ifdef _DEBUG_AMBTYPES			
				cout << "Atom " << atomsIter->GetIdx() << " is F. AtomicNum = " << atomsIter->GetAtomicNum() << endl;
#endif
      result[atomsIter->GetIdx()] = 16;
    }  
  }
  
  //! Test unassigned types and assign type - 1 (in order to match c arrays)
  for( OBMolAtomIter atomsIter(piMol); atomsIter; ++atomsIter ){
    if (result.at(atomsIter->GetIdx()) == 0) {
      cout << "Unknow type for atom " << molID2fileID[atomsIter->GetIdx()];
      cout << ", taking value ";
      if (atomsIter->IsHydrogen()) {
				result[atomsIter->GetIdx()] = 8;
				cout << result[atomsIter->GetIdx()] << " (generic H)" << endl;
      } else {
				result[atomsIter->GetIdx()] = 2;
				cout << result[atomsIter->GetIdx()] << " (generic C)" << endl;
      }
    }
    result[atomsIter->GetIdx()] = result[atomsIter->GetIdx()] - 1;  //// CAUTIONNNNN lo he quitado para probar  
    }
    //result.pop_back();
  return result;
}

vector< vector<unsigned int> > getMinimalDistances(OBMol& piMol) {

  vector< vector<unsigned int> > result(piMol.NumAtoms()+1, vector<unsigned int>(piMol.NumAtoms()+1, 0));
  vector<unsigned int> idxVector;
  idxVector.reserve(piMol.NumAtoms());
  for( OBMolAtomIter atom(piMol); atom; ++atom ){
    idxVector.push_back(atom->GetIdx());
  }
  for (unsigned int i = 0; i < idxVector.size() - 1; i++) {
    for (unsigned int j = i + 1; j < idxVector.size(); j++) {
      //if (piMol.GetAtom(HasAtomIdx(idxVector[i]))->IsConnected(piMol.GetAtom(HasAtomIdx(idxVector[j])))) {
      if (piMol.GetAtom(idxVector[i])->IsConnected(piMol.GetAtom(idxVector[j]))) {
				result[idxVector[i]][idxVector[j]] = 1;
				result[idxVector[j]][idxVector[i]] = 1;
			} else {
				result[idxVector[i]][idxVector[j]] = 1000; // Like infinite
				result[idxVector[j]][idxVector[i]] = 1000; // Like infinite
      }
    }
  }
  // Floyd-Warshall
  for (unsigned int i = 0; i < idxVector.size(); i++) {
    for (unsigned int j = 0; j < idxVector.size(); j++) {
      for (unsigned int k = 0; k < idxVector.size(); k++) {
				if ((j == i) || (k == i) || (j == k)) continue;
				unsigned int distancejk = result[idxVector[j]][idxVector[k]];
				unsigned int distancejik = result[idxVector[j]][idxVector[i]] + result[idxVector[i]][idxVector[k]];
				if (distancejik < distancejk) {
					result[idxVector[j]][idxVector[k]] = distancejik;
				}
      }
    }
  }
  return result;
}

vector< vector<unsigned int> > getTestAtoms4vdw(OBMol& piMol, const vector< vector<unsigned int> >& piMinimalDistances, vector<RotatableBond> rotBonds) {

  vector< vector<unsigned int> > result;
  vector<unsigned int> atoms4vdw;
  result.assign(piMol.NumAtoms()+1, atoms4vdw);
//  cout << "Number of atoms: " << piMol.NumAtoms() << " Number of bonds rota: " << rotBonds.size() << endl;
// 	OBMol *brokenMolecule = new(OBMol);
  OBMol brokenMolecule = piMol;
  
  for(unsigned int i = 0; i < rotBonds.size(); i++){
	  brokenMolecule.DeleteBond( brokenMolecule.GetBond( rotBonds.at(i).getAtoms().at(1)->GetIdx(), rotBonds.at(i).getAtoms().at(2)->GetIdx() ) );
  }
  
  vector< vector< int > > cfl;
  brokenMolecule.ContigFragList( cfl);
  unsigned int *blockIDs = new unsigned int[brokenMolecule.NumAtoms()+1];
//  cout << "Broken molecule size: " << brokenMolecule.NumAtoms() << " Fragments: " << cfl.size() <<endl;
  for (unsigned int i = 0; i < cfl.size(); i ++){
    for(unsigned int j = 0; j < cfl.at(i).size(); j++){
//      cout << brokenMolecule.GetAtom( cfl.at(i).at(j))->GetIdx() << endl;
      blockIDs[brokenMolecule.GetAtom( cfl.at(i).at(j))->GetIdx()] = i; // This might be = i... be aware of indexes in blockIDs variable.
    }
  }
 
  for( OBMolAtomIter iterAtoms(piMol); iterAtoms; ++iterAtoms ){
#ifdef _DEBUG_TESTATM4VDW
    cout << iterAtoms->GetIdx() << " ---> ";
#endif
    for( OBMolAtomIter iterAtoms2 = iterAtoms; iterAtoms2; ++iterAtoms2 ){
      if ((piMinimalDistances[iterAtoms->GetIdx()][iterAtoms2->GetIdx()] >= 3) && (*(blockIDs + iterAtoms->GetIdx()) != *(blockIDs + iterAtoms2->GetIdx()))) {
#ifdef _DEBUG_TESTATM4VDW
	cout << iterAtoms2->GetIdx() << " " ;
#endif
        result.at(iterAtoms->GetIdx()).push_back(iterAtoms2->GetIdx());
      }
    }
#ifdef _DEBUG_TESTATM4VDW
    cout << endl;
#endif

  }
  delete[] blockIDs;
  
  return result;
}

vector<unsigned long long> getSizeByIndexPosition(vector<RotatableBond>& piRotatableBonds) {
  vector<unsigned long long> result(piRotatableBonds.size(), 0);
  if (piRotatableBonds.size() > 0) {
    result[0] = 1;
    if (piRotatableBonds.size() > 1) {
      result[1] = piRotatableBonds[0].getAngles().size();
      for (unsigned int i = 1; i < piRotatableBonds.size(); i++) {
				result[i] = piRotatableBonds[i - 1].getAngles().size() * result[i - 1];
      }
    }
  }
  return result;
}

list<tConformer> generateConformers(OBMol& piMol, vector<RotatableBond>& piRotatableBonds, const vector< vector<unsigned int> >& piMinimalDistances, const vector< vector<unsigned int> >& piTestAtoms4vdw, const vector<unsigned int>& piAmberTypes){
  list<tConformer> result;
  vector<unsigned int> lastConformer(piRotatableBonds.size(), 0);
  vector<unsigned int> currentConformer(piRotatableBonds.size(), 0);
  vector<unsigned int> atomIdxs(piMol.NumAtoms(), 0);
  unsigned int maxAtomIdx = piMol.NumAtoms();
  double *origCoord = new double[piMol.NumAtoms() * 3];
  double *origCoordTmp; //! Save the original coordinates
  //double origCoord[MAX_ATOM];
  double *currCoord= new double[maxAtomIdx * 3];
  double *currCoordTmp;
  double *referenceCoord= new double[maxAtomIdx * 3];
  double *referenceCoordTmp;
  double vdwEnergy = 0;
  double rmsd = 0.0;
  double lastAcceptedEnergyMCSA = 1000;
  double temperatureMCSA = INITIAL_MCSA_TEMPERATURE;
  unsigned int maxGeneratedPerRoundMCSA = MAX_GENERATED_PER_MCSA_ROUND * piRotatableBonds.size();
  unsigned int maxAcceptedPerRoundMCSA = (unsigned int) floor(maxGeneratedPerRoundMCSA * MAX_ACCEPTED_PER_MCSA_ROUND);
  unsigned int generatedConformersInRoundMCSA = 0;
  unsigned int acceptedConformersInRoundMCSA = 1; //! The original conformer is the first accepted
  unsigned int acceptedConformersMCSA = 1; //! The original conformer is the first accepted
  unsigned int roundsMCSA = 1;
  vector<unsigned int> lastAcceptedConformerMCSA(piRotatableBonds.size(), 0);
  unsigned int howManySelectBAK = howManySelect;
  bool doMCSA = false;
  bool useAlternativeRMSD = false;
  unsigned int tooClose = 0;
//   OBMolAtomIter atomsIter;
  unsigned int i,k;
  unsigned int print;
  double x[maxAtomIdx];
  double y[maxAtomIdx];
  double z[maxAtomIdx];
  double q[maxAtomIdx];
  double meanEnergy;
  unsigned int cluster_conf =0;
  double d =0;
  double dMin = 0;
  double * momentsTmp2 = new double[NUM_MOMENTS];
  unsigned int casos = 0;
  list<tConformer>::iterator itConformers ;
  list<tConformer>::iterator itConformers2 ;
  vector< vector<float> > distances(atomIdxs.size(),vector<float>(atomIdxs.size(),0));
  unsigned int NUM_PREREFINED_TMP = 200;

#ifdef _DEBUG_TIMES
  OBStopwatch TMPTIME;
  TMPTIME.Start();
  cout << "Generate--> " << TMPTIME.Elapsed() << endl;
#endif

  OBAlign *rmsd_dist = new OBAlign(false,true);

  //! Check number of possible conformers 
  possibleConformers = 1;
  for (unsigned int i = 0; i < piRotatableBonds.size(); i++) {
    possibleConformers = possibleConformers * piRotatableBonds[i].getAngles().size();
  }

#ifdef _DEBUG
  cout << "#Posible-Confs: " << possibleConformers << " Max. allowed: "<< maxCombinations << endl;
#endif

  //! molId array
  i = 0;
  for ( OBMolAtomIter atomsIter(piMol); atomsIter; ++atomsIter) {
    atomIdxs[i] = atomsIter->GetIdx();
    i++;
  }
  
  origCoordTmp = piMol.GetCoordinates(); //! Save the original coordinates
  for(unsigned int i = 0; i < piMol.NumAtoms() * 3; i = i+3){
    origCoord[i  ] = origCoordTmp[i  ];
    origCoord[i+1] = origCoordTmp[i+1];
    origCoord[i+2] = origCoordTmp[i+2];
  }
    
  //! Check steric clashes
  tooClose = stericClashFilter(atomIdxs, origCoord, piAmberTypes, piTestAtoms4vdw, distances);

  if ( tooClose == 1 ){  //! Some atoms are too below 2 amstrongs
    cout << "Some atoms have esteric clashes in the original conformer" << endl;
    //cout << "EXITING!!" << endl;
    //return result;
  }

  //! Get VDW Energy
  vdwEnergy = getVDWEnergy(atomIdxs, origCoord, piAmberTypes, piTestAtoms4vdw, piMinimalDistances,distances);
 
  double * momentsTmp = new double[NUM_MOMENTS];

  if(clustering != 2){

    //! Centroid Calculation
    OBMolAtomIter atom(piMol);
    k=0;
    for(i=0; i< maxAtomIdx * 3; i=i+3){
      x[k] = origCoord[i];
      y[k] = origCoord[i+1];
      z[k] = origCoord[i+2];
      q[k] = atom->GetPartialCharge();
      ++atom;
      k++;
    }
 
    getCentroids(x,y,z,q,maxAtomIdx);
    momentsTmp = calculateMoments(x, y, z, q,maxAtomIdx,momentsTmp);
    cluster_conf++;
  }

  //! First, keep the original conformer
  print=0;
  //! RMSD
  if ((referenceFile != "") && (useInputInRMSD == 1)) {
    //! Sometimes the RMSD cannot be calculate (mabye because of a bug in RMSD of Openeye...)
    //! but the RMSD only taken into account coordinates can be used in these cases
    if (!useAlternativeRMSD) {	
      //OBAlign *rmsd_dist = new OBAlign(false,true);
      rmsd_dist->SetRefMol(piMol);
      rmsd_dist->SetTargetMol(referenceMolecule);
      rmsd_dist->Align();
      rmsd_dist->UpdateCoords(&referenceMolecule);
      rmsd = rmsd_dist->GetRMSD();
      if (rmsd == -1) {
	cout << "error computing RMSD 1" << endl;
	//! Not used alternative RMSD by the moment
	//! useAlternativeRMSD = true;
      }
    }
    if (rmsd < bestGeneratedRMSD) {
      bestGeneratedRMSD = rmsd;
      bestGeneratedRMSDEnergy = vdwEnergy;
    }
  }
  
  tConformer * origConformer = new tConformer;
  origConformer->conf = numConfig;
  origConformer->id = 0;
  for(unsigned int j=0;j< maxAtomIdx * 3; j++) origConformer->coord[j] = origCoord[j];
  origConformer->energy = vdwEnergy;
  origConformer->rmsd = rmsd;
  for(unsigned int j=0;j< NUM_MOMENTS; j++) origConformer->moments[j] = momentsTmp[j];
  origConformer->meanDist=0;

  result.clear();
  result.push_back(*origConformer);
 
  delete [] momentsTmp;
 
  meanEnergy = vdwEnergy;
 
  //! Apply input energy cut of
  if(clustering != 2){
    if(cutOff != 0.0){
      MAX_ENE = cutOff;
      cutOff = 0.0;
    };
  }
 
   
  if (possibleConformers > maxCombinations) {
    doMCSA = true;
    if(howManySelect> NUM_PREREFINED_CONF){
      NUM_PREREFINED_TMP = (unsigned int) floor(howManySelect * 1.5);
    }else{
      NUM_PREREFINED_TMP = (unsigned int) floor(NUM_PREREFINED_CONF * 1.5);
    }
  }
 
  //! Take the reference coord if needed
  if (referenceFile != "") {
    referenceCoordTmp = referenceMolecule.GetCoordinates();
    for(unsigned int i = 0; i < piMol.NumAtoms() * 3; i = i+3){
      referenceCoord[i  ] = referenceCoordTmp[i  ];
      referenceCoord[i+1] = referenceCoordTmp[i+1];
      referenceCoord[i+2] = referenceCoordTmp[i+2];
    }   
  }
 
  unsigned long long int numConf = 0;
  unsigned long long int cont1 =0;
  unsigned long long int cont2 =0;
  unsigned long long int step=0;
 
  if(possibleConformers < maxCombinations && clustering!=2){
     if(possibleConformers/NUM_PREREFINED_TMP > 1 ){
      step = (unsigned long long int)ceil(possibleConformers/NUM_PREREFINED_TMP);
    }else{
      step = 1;
    }
  }

  if (!isRigid(piMol, piRotatableBonds)) {
    for (unsigned long long i = 0; i < maxCombinations && i < possibleConformers; i++) {
      //! Select wich conformer to build
      if (doMCSA) { currentConformer = getNextConformerMCSA(lastAcceptedConformerMCSA, piRotatableBonds);
      } else if ( clustering != 2 ){
 	//! Select the constructors to build in a pseudo random way
 	if ( cont1*step + cont2 < possibleConformers ){
 	  numConf = cont1*step + cont2;
 	  cont1++;
 	}else{
 	  cont1 = 0;
 	  cont2++;
 	  numConf = cont1*step + cont2;
 	  cont1++;
 	}
 	currentConformer = getIDFromIdConformer(numConf,piRotatableBonds);
      } else { 	currentConformer = getNextConformer(lastConformer, piRotatableBonds);   }

      //! Get torsions for current conformer
      for (unsigned int j = 0; j < currentConformer.size(); j++) {
 	if (currentConformer[j] != lastConformer[j]) { //! Rotate
	  piMol.SetTorsion(piRotatableBonds[j].getAtoms()[0], piRotatableBonds[j].getAtoms()[1], piRotatableBonds[j].getAtoms()[2], piRotatableBonds[j].getAtoms()[3], piRotatableBonds[j].getAngles()[currentConformer[j] - 1] * DEG_TO_RAD);
 	}
      }
      currCoordTmp = piMol.GetCoordinates();
      for(unsigned int h = 0; h < piMol.NumAtoms() * 3; h = h+3){
	currCoord[h  ] = currCoordTmp[h  ];
	currCoord[h+1] = currCoordTmp[h+1];
	currCoord[h+2] = currCoordTmp[h+2];
      }
      generatedConformers++;
      lastConformer = currentConformer;

      //! Check esteric crashes
      tooClose = stericClashFilter(atomIdxs, currCoord,piAmberTypes,piTestAtoms4vdw, distances);
      if ( tooClose == 1 ){  //! Some atoms are too below 2 amstrongs
#ifdef _DEBUG_CONFS
 	cout << "Some atoms clash in conf: " << getConformerIndex(currentConformer)  << endl;
#endif
 	continue; //! Conformer not taken
      }
       
      //! Get vdw energy
      vdwEnergy = getVDWEnergy(atomIdxs, currCoord, piAmberTypes, piTestAtoms4vdw, piMinimalDistances,distances);
     
      //vdwEnergy = getScore(atomIdxs, piMol, currCoord, piMinimalDistances,piAmberTypes);
      //! RMSD
      if (referenceFile != "") {
 	//! El error puede estar al estar mal colocados los grupos alogenos en el fichero de entrada?????
 	rmsd_dist->SetRefMol(piMol);
	rmsd_dist->SetTargetMol(referenceMolecule);
	rmsd_dist->Align();
	rmsd_dist->UpdateCoords(&referenceMolecule);
	rmsd = rmsd_dist->GetRMSD();
 	if (rmsd == -1) {
	  cout << "error computing RMSD 2" << endl;
	//! Not used alternative RMSD by the moment
	//! useAlternativeRMSD = true;
 	}
 	if (rmsd < bestGeneratedRMSD) {
 	  bestGeneratedRMSD = rmsd;
 	  bestGeneratedRMSDEnergy = vdwEnergy;
 	}
      }
      
      if (doMCSA){
 	generatedConformersInRoundMCSA++;
 	//! Metropolis criterion
	// CHECKOUT randomrumbergeneration /////////////////////////////////////////////
	int random_num = ((double) rand() / (RAND_MAX));
	
 	double diffEnergy = vdwEnergy - lastAcceptedEnergyMCSA;
 	if ( (diffEnergy >= 0) && (random_num >= exp(diffEnergy - temperatureMCSA))) {
 	  continue; //! Conformer not taken
 	}
 	lastAcceptedEnergyMCSA = vdwEnergy;
 	lastAcceptedConformerMCSA = currentConformer;
 	acceptedConformersInRoundMCSA++;
	acceptedConformersMCSA++;
      }


      if(clustering != 2){
 	k = 0;
        //! Insert?
 	if( MAX_ENE <= vdwEnergy ){ //! Conformer with vdW clashes taken
	//if (( result.size() == NUM_PREREFINED_CONF ) && ( MAX_ENE <= vdwEnergy )){
 	  continue; //! Conformer not taken
 	}else{
 	  //! Compute centroids and moments for current molecule
 	  for(unsigned int j=0; j< maxAtomIdx * 3; j=j+3){
	    x[k] = currCoord[j];
	    y[k] = currCoord[j+1];
	    z[k] = currCoord[j+2];
	    k++;
	  }
	  getCentroids(x,y,z,q,maxAtomIdx);
	  momentsTmp2 = calculateMoments(x, y, z, q,maxAtomIdx, momentsTmp2);
 
	  //! Set current conformer into a tConformer struct
 	  tConformer * conformer = new tConformer;
 	  conformer->conf = numConfig;
 	  conformer->id = getConformerIndex(currentConformer);
	  
 	  for(unsigned int jjj = 0; jjj< maxAtomIdx * 3; jjj++) conformer->coord[jjj] = currCoord[jjj];
 	  conformer->energy = vdwEnergy;
 	  conformer->rmsd = rmsd;
 	  for(unsigned int jjj=0;jjj< NUM_MOMENTS; jjj++){
	    conformer->moments[jjj] = momentsTmp2[jjj];
	  }
	  
 	  bool inserted = false;
 	  int countConf=0;
 	  int replaceConf = 0;
 	  double energyTmp = 0;
 	  double meanDist = 0;
     	  int erase = 0;
 	  dMin = 0;
 	  meanDist = 0;
 	  //! Compute the mean distance of the current conformer to the rest of conformers results list and select the candidate to be replaces (the most similar one)
 	  for ( itConformers = result.begin(); itConformers != result.end(); itConformers++) {
	    d = shapeDistance(conformer->moments, itConformers->moments);
	    if(d > dMin){
	      dMin = d;
	      replaceConf = countConf;
	      energyTmp = itConformers->energy;
	    }
	    meanDist += d;
 	    countConf++;
	  }
 
	  if(countConf > 1 ){
	     conformer->meanDist = (meanDist-dMin)/(double)(countConf-1);
 	  }else if (countConf == 1 ){
	    conformer->meanDist = meanDist;
 	  }else if (countConf == 0) {
	    conformer->meanDist = 1;
 	  }
 
 	  //! If the list is not full and current conformer is disimiliar enought to the rest of conformers at result list, insert it!!
 	  if ((result.size() < NUM_PREREFINED_TMP)  && (dMin < dist_cutoff)) {
	    double distInsertConf = 0;
	    double distInsertConfTmp = 0;
	    int size = result.size();
  
	    if(size > 1){
	    //! Update mean distances on the list
	      for ( list<tConformer>::iterator confDist1 = result.begin(); confDist1 != result.end(); confDist1++) {
		//! Distance from the conformer on loop to the conformer to be replaced
		distInsertConfTmp = shapeDistance(conformer->moments, confDist1->moments);
		distInsertConf += distInsertConfTmp;
		//! Update mean distance by substracting the distance to the conformer to be replaced and adding the distance to the conformer to insert
		confDist1->meanDist = ( (confDist1->meanDist)*(size-1) +  distInsertConfTmp ) / (size);
	      }
	      conformer->meanDist = ( conformer->meanDist * (double)(size-1) + dMin ) / (double) size;
	    }
 
	    result.push_back(*conformer);

            //cout << i << " list is NOT full " << result.size() << " " << NUM_PREREFINED_TMP << endl; ///////////////

	    
	    inserted = true;
	    casos = 0;
	  //! If the list is not full but current conformer is very similar to one of the in the list OR the list is full
	  } else if ( ( (result.size() < NUM_PREREFINED_TMP)  && (dMin >= dist_cutoff) )  || (result.size() <= NUM_PREREFINED_TMP) ) {
 
	    //! Iterator pointing te candidate to be replaced
	    countConf=0;
	    erase = 0;
	    for(itConformers2 = result.begin(); itConformers2 != result.end(); itConformers2++){
	      if(countConf == replaceConf) break;
	      countConf ++;
	    }
 
	    //! TAFF WAY: Compute the mean dist of each conformer to the rest of conformers in the list EVERY TIME (VERY EXPENSIVE, BUT WORKS)
	    //int cc1 = 0;
	    //int cc2 = 0;
	    //for ( list<tConformer>::iterator confDist1 = result.begin(); confDist1 != result.end(); confDist1++) {
	      //confDist1->meanDist = 0;
	      //cc2 = 0;
	      //for ( list<tConformer>::iterator confDist2 = result.begin(); confDist2 != result.end(); confDist2++) {
		    //if( cc2 != cc1 ){
		      //confDist1->meanDist += shapeDistance(confDist1->moments, confDist2->moments);
		    //}
		    //cc2++;
	      //}
	      //confDist1->meanDist = confDist1->meanDist/(NUM_PREREFINED_CONF-1);
	      //cc1++;
	    //}
 
	    //! Insert conformer if its mean distance is larger than the mean distance of the candidate to be replaced
	    if( itConformers2->meanDist > conformer->meanDist ){
 
	      erase = 1;
 
	      double distReplaceConf = 0;
	      double distReplaceConfTmp = 0;
	      double distInsertConf = 0;
	      double distInsertConfTmp = 0;
	      int size = result.size();
 
	      if(size > 1){
		//! Update mean distances on the list
		for ( list<tConformer>::iterator confDist1 = result.begin(); confDist1 != result.end(); confDist1++) {
		  //! Distance from the conformer on loop to the conformer to be replaced
		  distReplaceConfTmp = shapeDistance(itConformers2->moments, confDist1->moments);
		  distReplaceConf += distReplaceConfTmp;
		  distInsertConfTmp = shapeDistance(conformer->moments, confDist1->moments);
		  distInsertConf += distInsertConfTmp;
		  //! Update mean distance by substracting the distance to the conformer to be replaced and adding the distance to the conformer to insert
		  confDist1->meanDist = ( (confDist1->meanDist)*(size-1) - distReplaceConfTmp +  distInsertConfTmp ) / (size-1);
		}
	      }
 
	      result.erase(itConformers2);
	      result.push_back(*conformer);	      
	      inserted = true;
	    }
	  }
	  delete [] conformer;
	}

      } else if ((result.size() == howManySelect) && (result.back().energy <= vdwEnergy)) {
	continue;
      } else {
 
	list<tConformer>::iterator itConformers;
	tConformer * conformer = new tConformer;
	conformer->conf = numConfig;
	conformer->id = getConformerIndex(currentConformer);
	for(unsigned int jjj =0;jjj< maxAtomIdx * 3; jjj++) conformer->coord[jjj] = currCoord[jjj];
	conformer->energy = vdwEnergy;
	conformer->rmsd = rmsd;
 
	bool inserted = false;
	for (itConformers = result.begin(); itConformers != result.end(); itConformers++) {
	  if (vdwEnergy < itConformers->energy) {
	    result.insert(itConformers, *conformer);
	    inserted = true;
	    break;
	  }
	}
 
	if (!inserted) {
	  result.insert(itConformers, *conformer); //! Insert last
	}
	if (result.size() > howManySelect) {
	  result.pop_back(); //! Remove last
	}
	delete [] conformer;
      }
      if (doMCSA) { //! MCSA exit conditions
	if ((generatedConformersInRoundMCSA == maxGeneratedPerRoundMCSA) ||(acceptedConformersInRoundMCSA == maxAcceptedPerRoundMCSA)) {
	  if (acceptedConformersInRoundMCSA == 0) {
	    break;
	  }
	  roundsMCSA++;
	  temperatureMCSA = temperatureMCSA * ANNEALING_SCHEDULE;
	  generatedConformersInRoundMCSA = 0;
	  acceptedConformersInRoundMCSA = 0;
	}
      }
    }
  } else {
    generatedConformers = 1; //! Because we have the original conformer
  }
#ifdef _DEBUG_TIMES
  cout << "CLUTERINT POST--> " << TMPTIME.Elapsed() << "   "<< result.size() << endl; //////////////////////////////
#endif

  if (doMCSA) { //! Clean repeated
    list<tConformer> resultMCSA;
    howManySelect = howManySelectBAK;
    if (howManySelect > 0) {
      resultMCSA.push_back(*result.begin());
      list<tConformer>::iterator itConformer = result.begin();
      itConformer++;
      while (itConformer != result.end() && resultMCSA.size() < howManySelect) {
	bool repeatedConformer = false;

	for (list<tConformer>::iterator itConformerFinal = resultMCSA.begin(); itConformerFinal != resultMCSA.end(); itConformerFinal++) {
	  if (itConformer->id == itConformerFinal->id) {
	    repeatedConformer = true;
	    break;
	  }
	}
	if (!repeatedConformer) {
	  resultMCSA.push_back(*itConformer);
	}
	//else {
	//delete [] itConformer->coord;
	//}
	itConformer++;
      }
    }
    result = resultMCSA;
  }
 
  //! Apply cutOff
  if (cutOff > 0.0) {
    double applyCutOff = result.front().energy + cutOff;
    while (result.back().energy > applyCutOff) {
      result.pop_back();
    }
  }
 
  //! Get best selected RMSD (if needed) !!!!!!!!!!!!!!! ESTO SE HACE 2 VECES: AQUÍ Y EN EL REFINADO
  if (referenceFile != "") {
    for (list<tConformer>::iterator itConformer = result.begin(); itConformer != result.end(); itConformer++) {
      if (itConformer->rmsd < bestSelectedRMSD && ((itConformer->id != 0) || (useInputInRMSD == 1))) {
	bestSelectedRMSD = itConformer->rmsd;
	bestSelectedRMSDEnergy = itConformer->energy;
      }
    }
  }

   piMol.SetCoordinates(origCoord); //! Restore the original coordinates

#ifdef _DEBUG_CONFS
  cout << "-------------*************************--------------------" << endl;
  unsigned int  countConf = 0;
  for(itConformers2 = result.begin(); itConformers2 != result.end(); itConformers2++){
    countConf++;
    cout << "Conf Num: " << countConf << " id: " << itConformers2->id << " dist " << itConformers2->meanDist << endl;
  }
#endif

  return result;
}

unsigned int stericClashFilter(const vector<unsigned int>& piAtomIdxs, const double* piCoordinates, const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw, vector< vector<float> >& poDistances){

  int idxAtomA,idxAtomB;
  for (unsigned int i = 0; i < piAtomIdxs.size(); i++) {
    idxAtomA = piAtomIdxs[i];
    for (unsigned int j = 0; j < piTestAtoms4vdw[idxAtomA].size(); j++) {
      idxAtomB = piTestAtoms4vdw.at(idxAtomA).at(j);
      poDistances[i][j] = (piCoordinates[3 * (idxAtomA-1)] - piCoordinates[3 * (idxAtomB-1)]) * (piCoordinates[3 * (idxAtomA-1)] - piCoordinates[3 * (idxAtomB-1)]) +
			  (piCoordinates[3 * (idxAtomA-1) + 1] - piCoordinates[3 * (idxAtomB-1) + 1]) * (piCoordinates[3 * (idxAtomA-1) + 1] - piCoordinates[3 * (idxAtomB-1) + 1]) +
			  (piCoordinates[3 * (idxAtomA-1) + 2] - piCoordinates[3 * (idxAtomB-1) + 2]) * (piCoordinates[3 * (idxAtomA-1) + 2] - piCoordinates[3 * (idxAtomB-1) + 2]);
#ifdef _DEBUG_ATOM_DIST
      cout << "Atoms: "<< idxAtomA <<" , "<< idxAtomB << "  =  "<< poDistances[i][j] << " , "<< RadiiMin[ piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]] << endl;
#endif
      //if(poDistances[i][j] < 2.25){
	
	if(poDistances[i][j] < RadiiMin[ piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]] ){
	  return(1);
	}
      //}
    }
  }
  return(0);
}

double getVDWEnergy(const vector<unsigned int>& piAtomIdxs, const double* piCoordinates, const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw, const vector<
		      vector<unsigned int> >& piMinimalDistances,vector< vector<float> > & piDistances) {

  double result = 0;
  unsigned int idxAtomA;
  unsigned int idxAtomB;

  double ra;
  double rb;
  double vdwTmp;
  double maxvdwTmp = 0;

  for (unsigned int i = 0; i < piAtomIdxs.size(); i++) {
    idxAtomA = piAtomIdxs[i];
    for (unsigned int j = 0; j < piTestAtoms4vdw[idxAtomA].size(); j++) {
      idxAtomB = piTestAtoms4vdw[idxAtomA][j];

      rb = piDistances[i][j] * piDistances[i][j] * piDistances[i][j];
      ra = rb * rb;

      //! soft VDW
      //vdwTmp = (VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/ra) - (VDWMATRIXB[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/rb);
      //! FULL LEONARD JONES POTENTIAL
      //vdwTmp = (VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/ra) - (VDWMATRIXB[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/rb);
      //! ONLY REPULSIVE TERM
      vdwTmp = (VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/ra);

      if (piMinimalDistances[idxAtomA][idxAtomB] == 3) { // 1-4 interaction
		vdwTmp = vdwTmp / 2.0; // correction
      }
      if( vdwTmp > maxvdwTmp)
         maxvdwTmp = vdwTmp;

      result += vdwTmp;

      //! PRINT VANDERWAALS ENERGY
#ifdef _DEBUG_ATOM_DIST
      cout << i << " - " << idxAtomA << " -> " << j << " - " << idxAtomB << " = " << piDistances[i][j] << endl;
      cout << "Amb. typs.->" << piAtomAmberTypes[idxAtomA] << " " << piAtomAmberTypes[idxAtomB] << " vdW Params.-> " << VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]] << " - " << VDWMATRIXB[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]] << " wdV ene-> " << vdwTmp << endl; 
#endif
    }
  }  
  return result;
}

bool isRigid(OBMol& piMol, vector<RotatableBond>& piRotatableBonds) {
  bool result = true;
  if (possibleConformers == 1) {
    for (unsigned int i = 0; i < piRotatableBonds.size(); i++) {
      double originalAngle = piMol.GetTorsion(piRotatableBonds[i].getAtoms()[0], piRotatableBonds[i].getAtoms()[1], piRotatableBonds[i].getAtoms()[2], piRotatableBonds[i].getAtoms()[3]);
      if (abs((int) rint(originalAngle - piRotatableBonds[i].getAngles()[0])) > 1.0) {
        result = false;
        break;
      }
    }
  } else {
    result = false;
  }
  return result;
}

vector<unsigned int> getNextConformerMCSA(const vector<unsigned int>& piLastConformer, vector<RotatableBond>& piRotatableBonds) {

  vector<unsigned int> result = piLastConformer;
  unsigned int torsionalToChange = 0;
  unsigned int selectedAngle = 0;
  unsigned int numOfTorsionals = piRotatableBonds.size();
  vector<bool> selectedTorsionalToChange(numOfTorsionals, false);
  unsigned int numOfChanges = 0;
  double random0to1 = 0.0;
  double random_num;
  
  if (result[0] == 0) { // Initialization
    for (unsigned int i = 0; i < piRotatableBonds.size(); i++) {
      unsigned int numOfAngles = piRotatableBonds[i].getAngles().size();
      random_num = ((double) rand() / (RAND_MAX));
      unsigned int initialAngle = 1 + (unsigned int) floor(random_num * numOfAngles);
      result[i] = initialAngle;
    }
  }
  random0to1 = ((double) rand() / (RAND_MAX));
  numOfChanges = CHANGES_IN_MCSA;
  if (random0to1 <= PROB_EXTRA_CHANGE) {
    numOfChanges++;
  }
  numOfChanges = numOfChanges <= numOfTorsionals ? numOfChanges : numOfTorsionals;

  for (unsigned int i = 0; i < numOfChanges; i++) {
    random_num = ((double) rand() / (RAND_MAX));
    torsionalToChange = (int) floor(random_num * numOfTorsionals);
    while (selectedTorsionalToChange[torsionalToChange]) {
      torsionalToChange++;
      if (torsionalToChange == numOfTorsionals) torsionalToChange = 0;
    }
    unsigned int numOfAngles = piRotatableBonds[torsionalToChange].getAngles().size();
    random_num = ((double) rand() / (RAND_MAX));
    selectedAngle = 1 + (unsigned int) floor(random_num * numOfAngles);
    if (selectedAngle == piLastConformer[torsionalToChange]) {
      selectedAngle++;
      if (selectedAngle > numOfAngles) selectedAngle = 1;
    }
    result[torsionalToChange] = selectedAngle;
    //cout << torsionalToChange << " " << selectedAngle << endl;
  }
  return result;
}

vector<unsigned int> getIDFromIdConformer(unsigned long long piIdConformer, vector<RotatableBond>& piRotatableBonds) {

  //vector<unsigned int> result(piRotatableBonds.size(),0);
  vector<unsigned int> result(piRotatableBonds.size(),0);
  vector<int> modules(piRotatableBonds.size(),0);
  int selectedAngle;

  //! Create me modules vector
  modules.at(0) = 1;
  for(unsigned int i = 1; i < piRotatableBonds.size(); i ++){
        modules.at(i) = piRotatableBonds.at(i-1).getAngles().size()*modules.at(i-1);
  }
  if (piIdConformer == 0) { // Original conformer
        for(unsigned int i = 0; i < piRotatableBonds.size(); i++){
          result.at(i) = piRotatableBonds.at(i).getAngles().size()-1;
        }
  } else {
    unsigned long long rest = piIdConformer - 1; // Because if we have n combinations then their numbers are 1..n, not 0..(n-1)
        result.clear();
    for (int i = (piRotatableBonds.size() - 1); i >= 0; i--) {
      selectedAngle = rest / modules.at(i);
      rest = rest % modules.at(i);
      result.insert(result.begin(), 1, selectedAngle + 1);
    }
  }

  //cout << "____>>>" << piIdConformer << "---" << result.at(0) <<result.at(1) <<result.at(2) <<result.at(3) <<result.at(4) <<result.at(5) <<result.at(6) <<result.at(7) << endl;

  return result;
}

vector<unsigned int> getNextConformer(const vector<unsigned int>& piLastConformer, vector<RotatableBond>& piRotatableBonds) {
  unsigned int carry = 1;
  unsigned int partialSum = 0;
  vector<unsigned int> result = piLastConformer;
  for (unsigned int i = 0; i < piRotatableBonds.size(); i++) {
    partialSum = piLastConformer[i] + carry;
    if (partialSum > piRotatableBonds[i].getAngles().size()) {
      partialSum = 1;
      carry = 1;
    } else {
      carry = 0;
    }
    if (partialSum == 0) partialSum = 1; // Initialization
    result[i] = partialSum;
  }
  return result;
}

unsigned long long getConformerIndex(vector<unsigned int> piConformer) {
  unsigned long long result = 1;
  for (unsigned int i = 0; i < piConformer.size(); i++) {
    result = result + ((piConformer[i] - 1) * sizeByIndexPosition[i]);
  }
  return result;
}

list< tConformer > getLibraryConformerList(list < tConformer > piConfs, OBMol piMol, const vector< vector<unsigned int> >& piMinimalDistances,const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw){

  vector<unsigned int> atomIdxs(piMol.NumAtoms()+1, 0);
  unsigned int i = 0;
  list< tConformer > resultTmp;

  for ( OBMolAtomIter atomsIter(piMol); atomsIter; ++atomsIter) {
    atomIdxs[i] = atomsIter->GetIdx();
    i++;
  }

  double min_ene=1000;

  for(unsigned int k=0;k<1;k++){
    int l =0;
    double ene = 0;
    list<tConformer>::iterator itConformers;
    for(itConformers = piConfs.begin(); itConformers != piConfs.end(); itConformers++){
      l++;
      ene = getScore(atomIdxs, piMol, itConformers->coord, piMinimalDistances, piAtomAmberTypes, piTestAtoms4vdw);
      itConformers->energy = ene;
      if (min_ene > ene) min_ene = ene;
#ifdef _DEBUG_CONFS_OUTINFO
      cout << "Molecule: " << k << " conf: " << l ;
      cout << " RMSD: " << itConformers->rmsd ;
      cout << " ene: " << ene;
      cout << " NumConfig. " << itConformers->conf ;
      cout << " id. " << itConformers->id;
      cout << " mdist "<<itConformers->meanDist << endl;
#endif
    }
  }

  list<tConformer>::iterator itConformers2;

  for(unsigned int l=0; l < 1; l++){
        for(list<tConformer>::iterator itConformers = piConfs.begin(); itConformers != piConfs.end(); itConformers++){

          if (min_ene + MAX_ENE_DIFF <= itConformers->energy ){
                continue;
          }else{

                bool inserted = false;
                int countConf=0;
                int replaceConf = 0;
                double energyTmp = 0;
                double meanDist = 0;
                int erase = 0;
                double dMin = 0;
                double d = 0;

                for ( itConformers2 = resultTmp.begin(); itConformers2 != resultTmp.end(); itConformers2++) {
                  d = shapeDistance(itConformers->moments, itConformers2->moments);
                  if(d > dMin){
                        dMin = d;
                        replaceConf = countConf;
                        energyTmp = itConformers->energy;
                  }
                  meanDist += d;
                  countConf++;
                }

                if(countConf > 1 ){
                  itConformers->meanDist = (meanDist-dMin)/(double)(countConf-1);
                }else if (countConf == 1 ){
                  itConformers->meanDist = meanDist;
                }else if (countConf == 0) {
                  itConformers->meanDist = 0;
                }

                if ((resultTmp.size() < howManySelect)  && (dMin < dist_cutoff)) {

                        double distInsertConf = 0;
                        double distInsertConfTmp = 0;
                        int size = resultTmp.size();

                        if(size > 1){
                          for ( list<tConformer>::iterator confDist1 = resultTmp.begin(); confDist1 != resultTmp.end(); confDist1++) {
                                distInsertConfTmp = shapeDistance(itConformers->moments, confDist1->moments);
                                distInsertConf += distInsertConfTmp;
                                confDist1->meanDist = ( (confDist1->meanDist)*(size-1) +  distInsertConfTmp ) / (size);
                          }
                          itConformers->meanDist = ( itConformers->meanDist * (double)(size-1) + dMin ) / (double) size;
                        }


                  resultTmp.push_back(*itConformers);
                  inserted = true;

                } else if ( (( (resultTmp.size() < howManySelect)  && (dMin >= dist_cutoff) )) || (resultTmp.size() <= howManySelect) ) {
                  countConf=0;
                  erase = 0;
                  for(itConformers2 = resultTmp.begin(); itConformers2 != resultTmp.end(); itConformers2++){
                        if(countConf == replaceConf) break;
                        countConf ++;
                  }


                  if( itConformers2->meanDist > itConformers->meanDist ){
                        erase = 1;

                        double distReplaceConf = 0;
                        double distReplaceConfTmp = 0;
                        double distInsertConf = 0;
                        double distInsertConfTmp = 0;
                        int size = resultTmp.size();

                        if(size > 1){
                          for ( list<tConformer>::iterator confDist1 = resultTmp.begin(); confDist1 != resultTmp.end(); confDist1++) {
                                distReplaceConfTmp = shapeDistance(itConformers2->moments, confDist1->moments);
                                distReplaceConf += distReplaceConfTmp;
                                distInsertConfTmp = shapeDistance(itConformers->moments, confDist1->moments);
                                distInsertConf += distInsertConfTmp;
                                confDist1->meanDist = ( (confDist1->meanDist)*(size-1) - distReplaceConfTmp +  distInsertConfTmp ) / (size-1);
                          }
                        }

                        resultTmp.erase(itConformers2);
                        resultTmp.push_back(*itConformers);
                        inserted = true;
                  }
                }
          }
        }
  }

  list<tConformer>::iterator itConformer1;
  list<tConformer>::iterator itConformerTmp = resultTmp.begin();
  list<tConformer> result;
  for(unsigned int n=0; n < resultTmp.size(); n++){
    for(itConformer1 = resultTmp.begin(); itConformer1 != resultTmp.end(); itConformer1++){
      if(itConformer1->meanDist < itConformerTmp->meanDist){
            itConformerTmp = itConformer1;
      }
    }
    result.push_back(*itConformerTmp);
    itConformerTmp->meanDist = 2.0;
    int idxAtomA = 2;
  }
  return result;
}

list< tConformer > getRefinedConformerList(vector < list < tConformer > > piConfs, vector < OBMol > piMol, const vector< vector<unsigned int> >& piMinimalDistances,const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw){
  
  vector<unsigned int> atomIdxs(piMol.at(0).NumAtoms()+1, 0);
  unsigned int i = 0;
  list< tConformer > resultTmp;

  for ( OBMolAtomIter atomsIter(piMol.at(0)); atomsIter; ++atomsIter) {
    atomIdxs[i] = atomsIter->GetIdx();
    i++;
  }	

  double min_ene=1000;

  for(unsigned int k=0;k<piConfs.size();k++){
    int l =0;
    double ene = 0;
    list<tConformer>::iterator itConformers;
    for(itConformers = piConfs.at(k).begin(); itConformers != piConfs.at(k).end(); itConformers++){
      l++;      
      ene = getScore(atomIdxs, piMol.at(k), itConformers->coord, piMinimalDistances, piAtomAmberTypes, piTestAtoms4vdw);
      itConformers->energy = ene;
      if (min_ene > ene) min_ene = ene;
#ifdef _DEBUG_CONFS_OUTINFO      
      cout << "Molecule: " << k << " conf: " << l ;
      cout << " RMSD: " << itConformers->rmsd ;
      cout << " ene: " << ene;
      cout << " NumConfig. " << itConformers->conf ;
      cout << " id. " << itConformers->id;
      cout << " mdist "<<itConformers->meanDist << endl;
#endif
    }
  }

  list<tConformer>::iterator itConformers2;

  for(unsigned int l=0; l < piMol.size(); l++){
	for(list<tConformer>::iterator itConformers = piConfs.at(l).begin(); itConformers != piConfs.at(l).end(); itConformers++){

	  if (min_ene + MAX_ENE_DIFF <= itConformers->energy ){
		continue; 
	  }else{

		bool inserted = false;
		int countConf=0;
		int replaceConf = 0;
		double energyTmp = 0;
		double meanDist = 0;
		int erase = 0;
		double dMin = 0;
		double d = 0;

		for ( itConformers2 = resultTmp.begin(); itConformers2 != resultTmp.end(); itConformers2++) {
		  d = shapeDistance(itConformers->moments, itConformers2->moments);
		  if(d > dMin){
			dMin = d;
			replaceConf = countConf;
			energyTmp = itConformers->energy;
		  }
		  meanDist += d;
		  countConf++;
		}

		if(countConf > 1 ){
		  itConformers->meanDist = (meanDist-dMin)/(double)(countConf-1);
		}else if (countConf == 1 ){
		  itConformers->meanDist = meanDist;
		}else if (countConf == 0) {
		  itConformers->meanDist = 0;
		}

		if ((resultTmp.size() < howManySelect)  && (dMin < dist_cutoff)) {

			double distInsertConf = 0;
			double distInsertConfTmp = 0;
			int size = resultTmp.size();

			if(size > 1){
			  for ( list<tConformer>::iterator confDist1 = resultTmp.begin(); confDist1 != resultTmp.end(); confDist1++) {
				distInsertConfTmp = shapeDistance(itConformers->moments, confDist1->moments);
				distInsertConf += distInsertConfTmp;
				confDist1->meanDist = ( (confDist1->meanDist)*(size-1) +  distInsertConfTmp ) / (size);
			  }
			  itConformers->meanDist = ( itConformers->meanDist * (double)(size-1) + dMin ) / (double) size;
			}


		  resultTmp.push_back(*itConformers);
		  inserted = true;

		} else if ( (( (resultTmp.size() < howManySelect)  && (dMin >= dist_cutoff) )) || (resultTmp.size() <= howManySelect) ) {
		  countConf=0;
		  erase = 0;
		  for(itConformers2 = resultTmp.begin(); itConformers2 != resultTmp.end(); itConformers2++){
			if(countConf == replaceConf) break;
			countConf ++;
		  }

		  if( itConformers2->meanDist > itConformers->meanDist ){
			erase = 1;

			double distReplaceConf = 0;
			double distReplaceConfTmp = 0;
			double distInsertConf = 0;
			double distInsertConfTmp = 0;
			int size = resultTmp.size();

			if(size > 1){
			  for ( list<tConformer>::iterator confDist1 = resultTmp.begin(); confDist1 != resultTmp.end(); confDist1++) {
				distReplaceConfTmp = shapeDistance(itConformers2->moments, confDist1->moments);
				distReplaceConf += distReplaceConfTmp;
				distInsertConfTmp = shapeDistance(itConformers->moments, confDist1->moments);
				distInsertConf += distInsertConfTmp;
				confDist1->meanDist = ( (confDist1->meanDist)*(size-1) - distReplaceConfTmp +  distInsertConfTmp ) / (size-1);
			  }
			}

			resultTmp.erase(itConformers2);
			resultTmp.push_back(*itConformers);
			inserted = true;
		  }
		}
	  }

	  int countConf = 0;
	  for(itConformers2 = resultTmp.begin(); itConformers2 != resultTmp.end(); itConformers2++){
		countConf++;
#ifdef _DEBUG_CONFS		
		cout << "REFINING: Conf Num: " << countConf << " id: " << itConformers2->id << " dist " << itConformers2->meanDist << endl;
#endif
	  }  
	}
  }

  bestSelectedRMSD = 1000;

  if (referenceFile != "") {
    for (list<tConformer>::iterator itConformer = resultTmp.begin(); itConformer != resultTmp.end(); itConformer++) {
      if (itConformer->rmsd < bestSelectedRMSD && ((itConformer->id != 0) || (useInputInRMSD == 1)) ) {
		bestSelectedRMSD = itConformer->rmsd;
		bestSelectedRMSDEnergy = itConformer->energy;
      }
    }
  }

  double *origCoord = new double[piMol.at(0).NumAtoms() * 3];
  double *origCoordTmp; //! Save the original coordinates

  //! Store Original coordinates
  OBMolAtomIter atom(piMol.at(0));
  origCoordTmp = piMol.at(0).GetCoordinates(); //! Save the original coordinates
  for(unsigned int i = 0; i < piMol.at(0).NumAtoms() * 3; i = i+3){
    origCoord[i  ] = origCoordTmp[i  ];
    origCoord[i+1] = origCoordTmp[i+1];
    origCoord[i+2] = origCoordTmp[i+2];
  }
  
  double minEne;
  double *newCoord = new double[piMol.at(0).NumAtoms() * 3];
  //double *newCoordTmp; //! Save the original coordinates
  unsigned int kk = 0;

  for(list<tConformer>::iterator itConformer = resultTmp.begin(); itConformer != resultTmp.end(); itConformer++){
    kk++;
    piMol.at(0).SetCoordinates(itConformer->coord);
    string ff = "GAFF";
    //string ff = "MMFF94";
    OBForceField* pFF = OBForceField::FindForceField(ff);
    pFF->Setup(piMol.at(0));
//! Perform Minimization
    if(outMini != 0){
#ifdef _DEBUG_MINI
      cout << kk << " Ene-ini: " << pFF->Energy() << endl;
#endif
      if(strncmp(minMethod.c_str(), "conjugate", 9) == 0){
	pFF->ConjugateGradients(minSteps, minEconv);
#ifdef _DEBUG_MINI
	cout << kk << " Ene-ConjugateGradients: " << pFF->Energy() << endl;
#endif
      }else if( strncmp(minMethod.c_str(), "steepest", 8) == 0){
	pFF->SteepestDescent(minSteps, minEconv);
#ifdef _DEBUG_MINI
	cout << kk << " Ene-SteepestDescent: " << pFF->Energy() << endl;
#endif

      }
    }

    itConformer->energy = pFF->Energy();
    if (minEne > itConformer->energy) minEne = itConformer->energy;
    
    pFF->UpdateCoordinates(piMol.at(0));
    newCoord = piMol.at(0).GetCoordinates(); //! Save the original coordinates
    for(unsigned int i = 0; i < piMol.at(0).NumAtoms() * 3; i = i+3){
      itConformer->coord[i  ] = newCoord[i  ];
      itConformer->coord[i+1] = newCoord[i+1];
      itConformer->coord[i+2] = newCoord[i+2];
    }   
  }

  list<tConformer>::iterator itConformer1;
  list<tConformer>::iterator itConformerTmp = resultTmp.begin();
  list<tConformer> result;
  for(unsigned int n=0; n < resultTmp.size(); n++){
    for(itConformer1 = resultTmp.begin(); itConformer1 != resultTmp.end(); itConformer1++){
      if(itConformer1->meanDist < itConformerTmp->meanDist){
	    itConformerTmp = itConformer1;
      }
      //! sort by energy
      //if(itConformer1->energy < itConformerTmp->energy){
	    //itConformerTmp = itConformer1;
      //}
    }
    result.push_back(*itConformerTmp);
    itConformerTmp->meanDist = 2.0;
    int idxAtomA = 2;
  }

  piMol.at(0).SetCoordinates(origCoord);

  return result;
}

double getScore(const vector<unsigned int>& piAtomIdxs, OBMol& piMol, double* piCoordinates, const vector< vector<unsigned int> >& piMinimalDistances,const vector<unsigned int>& piAtomAmberTypes, const vector< vector<unsigned int> >& piTestAtoms4vdw) {

  double vdwEne = 0,coulomb=0;
  unsigned int idxAtomA;
  unsigned int idxAtomB;
  double ra;
  double rb;
  double vdwTmp;
  int numSillas = 0;

  double distances;

  for (unsigned int i = 0; i < piAtomIdxs.size(); i++) {
    idxAtomA = piAtomIdxs[i];

    for (unsigned int j = 0; j < piTestAtoms4vdw[idxAtomA].size(); j++) {
      idxAtomB = piTestAtoms4vdw[idxAtomA][j];

      distances = (piCoordinates[3 * (idxAtomA-1)    ] - piCoordinates[3 * (idxAtomB-1)    ]) * (piCoordinates[3 * (idxAtomA-1)    ] - piCoordinates[3 * (idxAtomB-1)    ]) +
		  (piCoordinates[3 * (idxAtomA-1) + 1] - piCoordinates[3 * (idxAtomB-1) + 1]) * (piCoordinates[3 * (idxAtomA-1) + 1] - piCoordinates[3 * (idxAtomB-1) + 1]) +
		  (piCoordinates[3 * (idxAtomA-1) + 2] - piCoordinates[3 * (idxAtomB-1) + 2]) * (piCoordinates[3 * (idxAtomA-1) + 2] - piCoordinates[3 * (idxAtomB-1) + 2]);

      rb = distances * distances * distances;
      ra = rb * rb;
      
      //! soft VDW
      //vdwTmp = (VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/ra) - (VDWMATRIXB[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/rb);
      //! FULL LEONARD JONES POTENTIAL
      //vdwTmp = (VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/ra) - (VDWMATRIXB[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/rb);
      //! ONLY REPULSIVE TERM
      vdwTmp = (VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]/ra);

      if (piMinimalDistances[idxAtomA][idxAtomB] == 3) { // 1-4 interaction
	vdwTmp = vdwTmp / 2.0; // correction
      } else if (piMinimalDistances[idxAtomA][idxAtomB] <= 2){
	vdwTmp = 0;
      }

      vdwEne += vdwTmp;

      if (piMinimalDistances[idxAtomA][idxAtomB] > 3){
	coulomb += (332*piMol.GetAtom(idxAtomA)->GetPartialCharge() * piMol.GetAtom(idxAtomB)->GetPartialCharge()) / distances;
      }

      //! Debugging: PRINT VANDERWAALS ENERGY
      //cout << "idx" << i << "," << j << ": " << coulomb << " " << atomsIter1->GetPartialCharge() <<  " " << atomsIter2->GetPartialCharge()  << "  " << distances << endl;
      //cout << "VdW( " << idxAtomA << " , " << idxAtomB << " ) = " << vdwTmp << endl;  //////////////////////
      //cout << "dist2: "<< distances <<" A= "<< VDWMATRIXA[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]<<" B= " <<VDWMATRIXB[piAtomAmberTypes[idxAtomA]][piAtomAmberTypes[idxAtomB]]<< " Atom type A B: "<< piAtomAmberTypes[idxAtomA] <<" , " <<piAtomAmberTypes[idxAtomB] << "\n" << endl;

    }
  }

  getSillasBotes(piMol, numSillas);
  if(isnan(vdwEne)) vdwEne = 0;
  if(isnan(coulomb)) coulomb = 0;
  //double ene = vdwEne + coulomb + (double)numSillas * SILLA_SCORE;
  double ene = vdwEne  + (double)numSillas * SILLA_SCORE;

  //! Debuging scores
  //cout << "VDW: " << vdwEne << "--- COUL: " << coulomb << "--- silla: "<< numSillas*SILLA_SCORE << endl;

  return ene;
}

int getSillasBotes(OBMol& piMol, int& poNumSillas){
  

  //! Check whether rings are are aromatic or not in the imput molecules output, in case not check for sillas

  bool isSillaBoteTMP = false;
  bool isSillaTMP	= false;
  bool notAromatic = false;
  unsigned int* parts = new unsigned int[piMol.NumAtoms()];
  unsigned int* rings = new unsigned int[piMol.NumAtoms()];

  poNumSillas = 0;

  vector<OBRing*> vr;
  vr = piMol.GetSSSR();
  unsigned int countRings = vr.size();


//! iterating rings and getting atoms ide  
  vector<OBRing*>::iterator ringIter;
  vector<int>::iterator atmid;
  vector<OBRing*> *rlist = (vector<OBRing*>*)piMol.GetData("RingList");

  //! I look for the atoms involved in each ring
  for (ringIter = vr.begin();ringIter != vr.end();++ringIter){
      vector<int> sillaBoteAtoms;
      if( !(*ringIter)->IsAromatic() && (*ringIter)->Size() == 6){ 
	isSillaBoteTMP = true;
	unsigned int j = 0;
	for(atmid = (*ringIter)->_path.begin(); atmid != (*ringIter)->_path.end(); ++atmid){
	  //cout << *atmid << " -- " ;
	  isSillaBoteTMP = true;
	  sillaBoteAtoms.push_back(*atmid);
	  j++;
	}
	//cout << endl;
      }else {isSillaBoteTMP = false; }

      if(isSillaBoteTMP){ notAromatic = true; }else{ notAromatic = false; }

      if(notAromatic){

	//! Start looking for neighbours from the firs atom belonging to the ring
	unsigned int currIdx = 0;
	unsigned int nextIdx = 0;
	unsigned int prevIdx = 2000;
	unsigned int n = 0;

	//! Debuggin torsion angles
	//cout << "Tor1: " << piMol.GetTorsion(sillaBoteAtoms.at(0),sillaBoteAtoms.at(1),sillaBoteAtoms.at(2),sillaBoteAtoms.at(3))<< " : " << sillaBoteAtoms.at(0) << "-" << sillaBoteAtoms.at(1) << "-" << sillaBoteAtoms.at(2) << "-" <<	sillaBoteAtoms.at(3) << endl;
	//cout << "Tor2: " << piMol.GetTorsion(sillaBoteAtoms.at(1),sillaBoteAtoms.at(2),sillaBoteAtoms.at(3),sillaBoteAtoms.at(4))<< " : " << sillaBoteAtoms.at(1) << "-" << sillaBoteAtoms.at(2) << "-" << sillaBoteAtoms.at(3) << "-" <<	sillaBoteAtoms.at(4) << endl;
	//cout << "Tor3: " << piMol.GetTorsion(sillaBoteAtoms.at(2),sillaBoteAtoms.at(3),sillaBoteAtoms.at(4),sillaBoteAtoms.at(5))<< " : " << sillaBoteAtoms.at(2) << "-" << sillaBoteAtoms.at(3) << "-" << sillaBoteAtoms.at(4) << "-" <<	sillaBoteAtoms.at(5) << endl;
	//cout << "Tor4: " << piMol.GetTorsion(sillaBoteAtoms.at(3),sillaBoteAtoms.at(4),sillaBoteAtoms.at(5),sillaBoteAtoms.at(0))<< " : " << sillaBoteAtoms.at(3) << "-" << sillaBoteAtoms.at(4) << "-" << sillaBoteAtoms.at(5) << "-" <<	sillaBoteAtoms.at(0) << endl;
	//cout << "Tor5: " << piMol.GetTorsion(sillaBoteAtoms.at(4),sillaBoteAtoms.at(5),sillaBoteAtoms.at(0),sillaBoteAtoms.at(1))<< " : " << sillaBoteAtoms.at(4) << "-" << sillaBoteAtoms.at(5) << "-" << sillaBoteAtoms.at(0) << "-" <<	sillaBoteAtoms.at(1) << endl;
	//cout << "Tor6: " << piMol.GetTorsion(sillaBoteAtoms.at(5),sillaBoteAtoms.at(0),sillaBoteAtoms.at(1),sillaBoteAtoms.at(2))<< " : " << sillaBoteAtoms.at(5) << "-" << sillaBoteAtoms.at(0) << "-" << sillaBoteAtoms.at(1) << "-" <<	sillaBoteAtoms.at(2) << endl;

	double torAngs[6];

	torAngs[0]=piMol.GetTorsion(sillaBoteAtoms.at(0),sillaBoteAtoms.at(1),sillaBoteAtoms.at(2),sillaBoteAtoms.at(3));
	torAngs[1]=piMol.GetTorsion(sillaBoteAtoms.at(1),sillaBoteAtoms.at(2),sillaBoteAtoms.at(3),sillaBoteAtoms.at(4));
	torAngs[2]=piMol.GetTorsion(sillaBoteAtoms.at(2),sillaBoteAtoms.at(3),sillaBoteAtoms.at(4),sillaBoteAtoms.at(5));
	torAngs[3]=piMol.GetTorsion(sillaBoteAtoms.at(3),sillaBoteAtoms.at(4),sillaBoteAtoms.at(5),sillaBoteAtoms.at(0));
	torAngs[4]=piMol.GetTorsion(sillaBoteAtoms.at(4),sillaBoteAtoms.at(5),sillaBoteAtoms.at(0),sillaBoteAtoms.at(1));
	torAngs[5]=piMol.GetTorsion(sillaBoteAtoms.at(5),sillaBoteAtoms.at(0),sillaBoteAtoms.at(1),sillaBoteAtoms.at(2));

	for(n=1;n<6;n++){

	  if ( (torAngs[n-1] > 0 && torAngs[n] > 0 ) || (torAngs[n-1] < 0 && torAngs[n] < 0) ){
	    isSillaTMP = false;
	    break;
	  }
	  isSillaTMP = true;
	}
	if (isSillaTMP){
	      poNumSillas++;
	}
      }
    sillaBoteAtoms.clear();
  }

  delete[] parts;
  delete[] rings;

  return(0);
}

void writeMol2Output(OBMol& piMol, list<tConformer>& piConformers) {

  ofstream ofs;
  OBConversion * conv = new OBConversion();
  conv->SetOutFormat("MOL2");
  string outputFileTmp = outputFile + ".mol2";
  ofs.open(outputFileTmp.c_str(),ios::out | ios::app);

  double *origCoord = new double[piMol.NumAtoms() * 3];
  double *origCoordTmp = piMol.GetCoordinates(); //! Save the original coordinates
  for(unsigned int i = 0; i < piMol.NumAtoms() * 3; i = i+3){
    origCoord[i  ] = origCoordTmp[i  ];
    origCoord[i+1] = origCoordTmp[i+1];
    origCoord[i+2] = origCoordTmp[i+2];
  }


  for (list<tConformer>::iterator itConformer = piConformers.begin(); itConformer != piConformers.end(); itConformer++) {
    piMol.SetCoordinates (itConformer->coord);
    if (referenceFile != "") {
      string s = ulonglong2string(itConformer->conf) + "_" + ulonglong2string(itConformer->id) + ", energy = " + double2string(itConformer->energy) +
		 ", RMSD = " + double2string(itConformer->rmsd) + ", meanDist = " +  double2string(itConformer->meanDist);
      piMol.SetTitle(s);
    } else {
      string s = ulonglong2string(itConformer->conf) + "_" + ulonglong2string(itConformer->id) + ", energy = " + double2string(itConformer->energy);
      piMol.SetTitle(s);
    }
    piMol.Center();
    conv->Write(&piMol,&ofs);
  }
  
  piMol.SetCoordinates(origCoord); // Restore original coordinates
  ofs.close();
  free(conv);
  free(origCoord);
}

string ulonglong2string(unsigned long long piInput) {
  ostringstream outs;
  outs << piInput;
  return outs.str();
}

string double2string(double piInput) {
  ostringstream outs;
  outs << piInput; ///// Cuidadín a ver cuantos decimales mete...
  return outs.str();
}

void writePDBOutput(OBMol& piMol, list<tConformer>& piConformers) {

  ofstream ofs;
  OBConversion * conv = new OBConversion();
  conv->SetOutFormat("PDB");
  string outputFileTmp = outputFile + ".pdb";
  ofs.open(outputFileTmp.c_str());
  
  double *origCoord = new double[piMol.NumAtoms() * 3];
  double *origCoordTmp = piMol.GetCoordinates(); //! Save the original coordinates
  for(unsigned int i = 0; i < piMol.NumAtoms() * 3; i = i+3){
    origCoord[i  ] = origCoordTmp[i  ];
    origCoord[i+1] = origCoordTmp[i+1];
    origCoord[i+2] = origCoordTmp[i+2];
  }
  
  for (list<tConformer>::iterator itConformer = piConformers.begin(); itConformer != piConformers.end(); itConformer++) {
    piMol.SetCoordinates(itConformer->coord);

    if (referenceFile != "") {
      string s = ulonglong2string(itConformer->conf) + "_" + ulonglong2string(itConformer->id) + ", energy = " + double2string(itConformer->energy) +
		 ", RMSD = " + double2string(itConformer->rmsd) + ", meanDist = " +  double2string(itConformer->meanDist);
      piMol.SetTitle(s);

    } else {
      string s = ulonglong2string(itConformer->conf) + "_" + ulonglong2string(itConformer->id) + ", energy = " + double2string(itConformer->energy);
      piMol.SetTitle(s);
    }
    piMol.Center();
    conv->Write(&piMol,&ofs);
  }
  
  piMol.SetCoordinates(origCoord); // Restore original coordinates
  ofs.close();
  free(conv);
  free(origCoord);
}

void writeXMLOutput(list<tConformer>& piConformers, vector< vector<unsigned int> > piAmberTypes, vector< vector<RotatableBond> > piRotatableBonds){
  ofstream ofs;
  ofs.open((outputFile + ".xml").c_str());
  ofs << "<ALFA-output>" << endl;

  for (unsigned int j = 0; j < piAmberTypes.size(); j++) {
	ofs << "<atom-types idCorina = \"" << (j+1) << "\" >" << endl;
	for (unsigned int i = 0; i < piAmberTypes.at(j).size(); i++) {
	  ofs << piAmberTypes[j][i] + 1 << " "; // Restore de original amber type
	}
	ofs << endl;
	ofs << "</atom-types>" << endl;
  }
  ofs << "<configurations>" << endl;
  for (unsigned int k = 0; k < piRotatableBonds.size(); k++) {
	ofs << "<torsionals idCorina = \"" << (k+1) << "\" numTorsionals=\"" << piRotatableBonds[k].size() << "\">" << endl;
	for (unsigned int i = 0; i < piRotatableBonds[k].size(); i++) {

	  ofs << "<torsional idCorina = \"" << (k+1) << "\" atoms = \"";
	  ofs << molID2fileID[piRotatableBonds[k][i].getAtoms().at(0)->GetIdx()] << " ";
	  ofs << molID2fileID[piRotatableBonds[k][i].getAtoms().at(1)->GetIdx()] << " ";
	  ofs << molID2fileID[piRotatableBonds[k][i].getAtoms().at(2)->GetIdx()] << " ";
	  ofs << molID2fileID[piRotatableBonds[k][i].getAtoms().at(3)->GetIdx()] << "\"";
	  ofs << " type = \"" << piRotatableBonds[k][i].getType() << "\"";
	  ofs << " numAngles = \"" << piRotatableBonds[k][i].getAngles().size() << "\"";
	  ofs << " angles = \"";
	  for (unsigned int j = 0; j < piRotatableBonds[k][i].getAngles().size(); j++) {
		if (j != 0) {
		  ofs << " ";
		}
		ofs << piRotatableBonds[k][i].getAngles().at(j);
	  }
	  ofs << "\"/>" << endl;
	}
	ofs << "</torsionals>" << endl;
  }
  ofs << "</configurations>" << endl;

  ofs << "<conformers>" << endl;
  for (list<tConformer>::iterator itConformers = piConformers.begin(); itConformers != piConformers.end(); itConformers++) {
    ofs << "<conformer id-conformer = \"" << itConformers->conf << "_" << itConformers->id << "\" vdw-energy = \"";
    ofs << setw(9) << setprecision(4) << right << fixed << itConformers->energy;
    ofs << "\"/>" << endl;
  }
  ofs << "</conformers>" << endl;
  ofs << "</ALFA-output>" << endl;
  ofs.close();
}
