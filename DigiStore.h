#ifndef ROOT_DigiStore
#define ROOT_DigiStore

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include "DigiOpt1.h"

using namespace std;

/* 

 */


class DigiStore : public TObject {
  
public:

  TTimeStamp EventTime;
  UInt_t ENum;			// Event Number
  UInt_t REnum;			// Event Number
  UInt_t CEnum;
  
  int NDigiOpt1; // Number of DigiStore not sets!
  
  TClonesArray  *fDigiOpt1;    //->array with all DigiStore

  static TClonesArray *fgDigiOpt1;

public:

  DigiStore();
  virtual ~DigiStore();
  void ClearDigi();// clears previous track objects..
  
  DigiOpt1  *AddDigi1Store();
  TClonesArray *GetDigiOpt1() const {return fDigiOpt1;}
  
  ClassDef(DigiStore,1);  //Event structure 
};

#endif
 
