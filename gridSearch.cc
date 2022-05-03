
#include "gridSearch.h"

gridSearch::gridSearch(Int_t maxpar) {

  fMaxPar    = maxpar;

  fParNames  = new TString[fMaxPar];
  fParStart  = new Double_t[fMaxPar];
  fParStep   = new Double_t[fMaxPar];
  fParLimLow = new Double_t[fMaxPar];
  fParLimUp  = new Double_t[fMaxPar];

  fFCN       = 0;
  fIter      = 1;
}

gridSearch::~gridSearch() {
  delete [] fParNames;
  delete [] fParStart;
  delete [] fParStep;
  delete [] fParLimLow;
  delete [] fParLimUp;
}

void gridSearch::SetIter(Int_t iter) {
  fIter = iter;
}

void gridSearch::SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t)) {
  fFCN = fcn;
}

void gridSearch::SetParams(Int_t kpar, TString parname, Double_t vstart, Double_t vstep, Double_t vlow, Double_t vup) {

  if(kpar<fMaxPar) {
    fParNames[kpar]  = parname;
    fParStart[kpar]  = vstart;
    fParStep[kpar]   = vstep;
    fParLimLow[kpar] = vlow;
    fParLimUp[kpar]  = vup;

    fParVal[kpar]    = vstart;
  }
}

void gridSearch::GetParams(Int_t kpar, TString &parname, Double_t &val, Double_t &err) const {
  if(kpar<fMaxPar) {
    parname = fParNames[kpar];
    val     = fParVal[kpar];
    err     = 0;
  }
}

double gridSearch::GetChi2() {
  return fChi2;
}

void gridSearch::Eval(Int_t npar, Double_t *gin, Double_t &fval, Double_t *par, Int_t flag) {
  if(fFCN) (*fFCN)(npar, gin, fval, par, flag);
}

void gridSearch::DoMinimization() {
  
  double prevchi2 = 10000;
  double fnew;
  
  double* prevPars = new double[fMaxPar];
  double* tempPars = new double[fMaxPar];
  double* gin      = new double[fMaxPar];
  for(int ijpar=0;ijpar<fMaxPar;ijpar++) {
    prevPars[ijpar] = fParStart[ijpar];
    tempPars[ijpar] = fParStart[ijpar];}

  for(int iter=0;iter<fIter;iter++) {
    
    for(int ijpar=0;ijpar<fMaxPar;ijpar++) {
      
      bool fForward = false;
      bool fBackward = false;
      int stepCnt = 0;
      while(1) {
	
	stepCnt++;
	if(fForward && fBackward) break;
	  
	int fDir = 0;
	if(!fForward) {fDir = 1;} else {fDir = -1;}
	
	tempPars[ijpar] += fParStep[ijpar]*fDir/(iter + 1);
	
	int tmpnpar = 0;
	Eval(tmpnpar, gin, fnew, tempPars, tmpnpar);
	if(fnew<prevchi2) {
	  for(int ijpar=0;ijpar<fMaxPar;ijpar++) {
	    prevPars[ijpar] = tempPars[ijpar];}
	  prevchi2 = fnew;
	} else {
	  if(!fForward) {fForward = true;} else if(!fBackward) {fBackward = true;}
	}
      }
    } // for(int ijpar=0;ijpar<5;ijpar++) {
    prevchi2 = 10000;
  } // for(int iter=0;iter<iter;iter++) {
  
  for(int ijpar=0;ijpar<fMaxPar;ijpar++) {
    fParVal[ijpar] = prevPars[ijpar];}
  fChi2 = fnew;

  delete [] prevPars;
  delete [] tempPars;
  delete [] gin;
    
}
