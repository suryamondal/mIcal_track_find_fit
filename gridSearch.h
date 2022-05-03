#ifndef gridSearch__H
#define gridSearch__H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <bitset>

#include "TMath.h"
#include "TString.h"

class gridSearch {

 public:
  gridSearch(Int_t maxpar);
  ~gridSearch();

  void   SetIter(Int_t iter);
  void   SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t));
  void   SetParams(Int_t kpar, TString parname, Double_t vstart, Double_t vstep, Double_t vlow, Double_t vup);
  void   GetParams(Int_t kpar, TString &parname, Double_t &val, Double_t &err) const;
  void   DoMinimization();
  void   DoMinimizationFullRange(Int_t range);
  double GetChi2();
  
 private:

  Int_t        fIter;		/* Iteration */
  Int_t        fMaxPar;		/* Maximum number of parameters */
  TString*     fParNames;	/* Parameter names */
  Double_t*    fParStart;	/* Parameter start values */
  Double_t*    fParStep;	/* Parameter step */
  Double_t*    fParLimLow;	/* Parameter lower limit */
  Double_t*    fParLimUp;	/* Parameter upper limit */

  void         (*fFCN)(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t flag);

  Double_t*    fParVal; 	/* parameter start values */
  Double_t     fChi2;

  void   Eval(Int_t npar, Double_t *gin, Double_t &fval, Double_t *par, Int_t flag);

};

#endif
