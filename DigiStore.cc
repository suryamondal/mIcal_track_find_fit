#include <iostream>
#include <fstream>
#include <iomanip>
#include "TMath.h"   

#include "DigiStore.h"

using namespace std;

ClassImp(DigiOpt1)

TClonesArray *DigiStore::fgDigiOpt1 = 0;

DigiStore::DigiStore() {
  // Create an Hit object.
  if (!fgDigiOpt1) fgDigiOpt1 = new TClonesArray("DigiOpt1", 1000);
  fDigiOpt1 = fgDigiOpt1;
  NDigiOpt1=0;
}

//______________________________________________________________________________
DigiStore::~DigiStore() {
  
}

//______________________________________________________________________________
DigiOpt1 * DigiStore::AddDigi1Store() {
  // Add a new track to the list of tracks for this event.
  TClonesArray &digi = *fDigiOpt1;
  DigiOpt1 *digipos = new(digi[NDigiOpt1++]) DigiOpt1();
  return digipos;
}

void DigiStore::ClearDigi() {
  NDigiOpt1 = 0; // !!! reset counter !!!
  if(fDigiOpt1) {
    // fDigiOpt1->Clear("C"); // will also call Track::Clear
    fDigiOpt1->Delete(); // will also call Track::Clear
  }
  ENum = 0;			// Event Number
  REnum = 0;			// Event Number
  CEnum = 0;
}

