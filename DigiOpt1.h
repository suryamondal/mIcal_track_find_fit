#ifndef ROOT_digi_opt1
#define ROOT_digi_opt1

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include <vector>
using namespace std;
class DigiOpt1 : public TObject {
 public:
  // unsigned int stripid;
  // unsigned int digitime;
  // unsigned int digienr;

  UShort_t rpcId;
  //      2 bit for INO module  
  //      3 bit for X-row
  //      3 bit for Y row
  //      8 bit for Z-layer
  
  ULong64_t strips[2];
  
  vector<UInt_t> tdc;
  // first two element is tdc_ref_l and tdc_ref_t
  // next element is <tdc time> <l or t: 1bit> <tdc id: 4bit>
  
  
  DigiOpt1(){};
  // DigiOpt1(Float_t random){random=0;};
  
  virtual ~DigiOpt1(){};
  ClassDef(DigiOpt1,1);
};

#endif
