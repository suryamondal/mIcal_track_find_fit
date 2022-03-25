#ifndef ROOT_HealthStore
#define ROOT_HealthStore

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include <vector>

//#include <pair>
using namespace std;

/* 

 */


class HealthStore : public TObject {
 public:
  
  TTimeStamp               HealthTime;
  
  vector<UShort_t>         rpcId;
  vector<UInt_t>           TPHdata;
  vector<UShort_t>         HVdata[2];	// for X and Y side in Volt
  vector<UShort_t>         AMPdata[2];	// for X and Y side in nA
  vector<vector<UInt_t>>   NOISEdata;	// one vector<UInt_t> per layer
  
public:
  HealthStore() {};
  virtual ~HealthStore() {};
  ClassDef(HealthStore,1);  //Event structure 
};

#endif
 
