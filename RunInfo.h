#ifndef ROOT_run_info
#define ROOT_run_info

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TString.h"
#include <vector>

using namespace std;

class RunInfo : public TObject {
  
 public:
  
  UShort_t   MagnetCurrent;	  /* in Amp */
  TTimeStamp StartTime, StopTime; /*  */
  TString    TriggerInfo;	  /*  */
  TString    Creator;		  /*  */
  TTimeStamp CreatedOn;		  /* UTC TimeStamp of FileCreation */
  TString    Comment;		  /*  */
  
  RunInfo(){};
  /* RunInfo(Float_t random){random=0;}; */
  
  virtual ~RunInfo(){};
  ClassDef(RunInfo,1);
};

#endif
