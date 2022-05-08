



// #define isIter

#define isSimData
// #define isSpecial		// sim data taken from reco file

// #define isTriggerRateOnly

#define isTrgCheck

// #define isDebug
// #define isSpclDebug
// #define isStrpMulti

// #define removeEarlierHit	// see top description

#define isMiniICAL
#define isTimePosFit
#define isTDCstrpCorr

// #define isRungeKutta
#define isEloss

#define isAnalysis
// #define isCorrection		// reconstrution happens if not define

// #define is5of8			// trigger for efficiency or correction

// #define isLifetime		// TMinuit part is disabled if this flag is defined

#define is2dpos



#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <bitset>

#include "TTimeStamp.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TObject.h"
#include "TRandom.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMinuit.h"
#include "TF1.h"


#ifndef isSimData
#include "RunInfo.h"
#include "DigiStore.h"
#include "CauStore.h"
#include "HealthStore.h"
#else  // #ifndef isSimData
#ifndef isSpecial
#include "T2.h"
#endif	// #ifndef isSpecial
#endif	// #ifndef isSimData


// #ifdef isTDCstrpCorr
// #include "mIcal_correction_20210325.txt"
// #endif	// #ifdef isTDCstrpCorr


#include "mystuff.h"
#include "data.h"
#include "circle.h"
#include "Utilities.cpp"
// #include "CircleFitByTaubin.cpp"
// #include "CircleFitByPratt.cpp"
// #include "CircleFitByKasa.cpp"
// #include "CircleFitByHyper.cpp"
#include "CircleFitByChernovLesort.cpp"
// #include "CircleFitByChernovHoussam.cpp"
// #include "CircleFitByLevenbergMarquardtFull.cpp"
// #include "CircleFitByLevenbergMarquardtReduced.cpp"
// #include <time.h>



using namespace std;

/*
  
  
*/


// #ifdef isIter
// const int itermodule = 0;
// const int iterxrow = 0;
// const int iteryrow = 0;
// const int iterlayer = 0;
// #endif	// #ifdef isIter

#ifdef isCorrection
#ifdef is5of8
const int ntrigLayers1 = 5;
const int trigLayers1[5] = {1,2,3,4,5};
#else  // #ifdef is5of8
const int ntrigLayers1 = 4;
const int trigLayers1[4] = {6,7,8,9};
#endif	// #ifdef is5of8
const int MinLayHit = ntrigLayers1; 
#else	// #ifdef isCorrection
const int ntrigLayers1 = 4;
const int trigLayers1[4] = {6,7,8,9};
#endif	// #ifdef isCorrection


const int        nside         =   2;
const int        nlayer        =  10;
const int        nxrow         =   1;
const int        nyrow         =   1;
const int        nmodule       =   1;
const int        nstrip        =  64;
const int        nTDC          =   8;
const double     tdc_least     =   0.1;	 // in ns
const double     strpwidth     =   0.03; // in m

const int        blockM        =   nstrip/4;

const int        MaxTotHit     =   nstrip/2;
const int        MaxTotPos     =   3;
const int        MaxMulti      =   3;
const int        nmxhits       =   6; // important variable
#ifndef isCorrection
const int        MinLayHit     =   5;
#endif	// #ifndef isCorrection
const int        MinClusterSep =   0;
const double     maxPosDev     =   1.25;
const double     maxPosDev_Effi=   2.5;
const int        strpTimeFit   =   4;


const double     muMass        =   1.88353e-28; /* in Kg */
const double     muCharge      =   1.60218e-19;	/* in Coulombs */

const double     maxTrkLen     =   4.;		/* in meter */
const double     gevtojoule    =   1.60218e-10; /* GeV to Joules */
const double     minPartE      =   0.03;	/* in GeV */

const double     uniformField  =   1.5;	   /* in Tesla */
const double     stepSize      =   0.001;  /* in m     */
const double     timeStepFit   =   3.336e-12; /* in Seconds */

const double     thetaSpread   =  15.*TMath::Pi()/180.;	/* */
const double     phiSpread     =  30.*TMath::Pi()/180.;	/* */
const double     posSpread     =   3.*strpwidth;
// const double     momSpread     =  3.; // in GeV
const double     momSpread     =  0.5; // in GeV
const double     minMom        =  0.4; // in GeV

const double     maxtime       =  22.e3;      // in ns
const double     spdl          =   5.;	      // ns/m
const double     cval          =   0.29979;   // light speed in m/ns
const double     cval1         =   0.29979e9; /* light speed in m/s */

const double     airDensity    = 0.0012; /* in g/cm3 */
const double     gapDensity    = 0.4;	 /* in g/cm3 */
const double     ironDensity   = 7.86;	 /* in g/cm3 */

const int        circlePt      =  nlayer;
// const double     circleLen     =  circlePt*(0.056+0.045);

const int        multiDiv1     =  15;
const int        multiDiv2     =   5;


#ifndef isSimData
// #include "mIcal_correction_20210325.txt"
// #include "mIcal_correction_20210702.txt"
#include "mIcal_correction_TDCstrp_20210702.txt"
#include "mIcal_correction_PosTime_20210702.txt"
#else  // #ifndef isSimData
#include "mIcal_correction_sim.txt"
#endif	// #ifndef isSimData



TGraph *eLossCurve;
TH2D *xyvsbxin, *xyvsbyin;



// const double poserrsq_ref[6] = {0.08,0.08,0.08,0.08,0.08,0.08};
// const double timeErr_ref[6] = {1.,1.,1.,1.,1.,1.};

// const double poserrsq_ref[6] =
//   {
//    0.103728,	// m0 xr0 yr0 z4 x mul1
//    0.0794917,	// m0 xr0 yr0 z4 x mul2
//    0.131974,	// m0 xr0 yr0 z4 x mul3
//    0.449816,	// m0 xr0 yr0 z4 x mul4
//    0.849692,	// m0 xr0 yr0 z4 x mul5
//    1.23988	// m0 xr0 yr0 z4 x mul6
//   };
// const double timeErr_ref[6] =
//   {
//    1.2939,	// m0 xr0 yr0 z4 x mul1
//    1.43214,	// m0 xr0 yr0 z4 x mul2
//    1.07453,	// m0 xr0 yr0 z4 x mul3
//    1.416,	// m0 xr0 yr0 z4 x mul4
//    1.81177,	// m0 xr0 yr0 z4 x mul5
//    1.69344	// m0 xr0 yr0 z4 x mul6
//   };

// const double poserrsq[1][1][1][10][2][6] = { // nmod,xrow,yrow,zlay,side,nmxhits
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08,
//    0.08,0.08,0.08,0.08,0.08,0.08
// };

// const double timeerrsq[1][1][1][10][2][6] = { // nmod,xrow,yrow,zlay,side,nmxhits
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1,
//    1,1,1,1,1,1
// };


TVector3 ttxx, ttyy;
double calPointDist(TVector3 xx, TVector3 yy) {
  double dist = sqrt(pow(xx.X()-yy.X(),2.) +
		     pow(xx.Y()-yy.Y(),2.) +
		     pow(xx.Z()-yy.Z(),2.));
  return dist;
};


#ifdef isMiniICAL
const double   gapThickness      = 0.008; // 
const double   ironThickness     = 0.056; // 
const double   airGap            = 0.045; //
const double   rpcZShift         =-0.005; // 
const double   rpcXdistance      = 2.;	  // 
const double   rpcYdistance      = 2.1;	  // 
const double   rpcXOffset        = 0.;	  // 
const double   rpcYOffset        = 0.;	  // 
const double   moduleDistance    = 0.;	  // 

TVector3 getRPCpos(Int_t module,
		   Int_t xrow,
		   Int_t yrow,
		   Int_t zlay) {
  double zpos = (zlay + 1)*(ironThickness + airGap) - airGap/2. + rpcZShift; // z offset of rpc in from the middle ~0.5cm
  double ypos = rpcYdistance*yrow + rpcYOffset;
  double xpos = rpcXdistance*xrow + rpcXOffset;
  TVector3 pos(xpos,ypos,zpos);
  return pos;			// in meter
};


TVector3 GetMagneticField(double x, double y, double z) {
  /* Getting B field by Changing the coordinates to match histogram coordinate **/
  
  // cout << " xpos " << x << " ypos " << y << " bin " << xyvsbxin->FindBin((x + 1.)*1000.,(y + 1.)*1000.) << endl;
#ifdef isSimData
  double Bx = -0.85*xyvsbxin->GetBinContent(xyvsbxin->FindBin((x + 1.)*1000.,(y + 1.)*1000.));
  double By = -0.85*xyvsbyin->GetBinContent(xyvsbyin->FindBin((x + 1.)*1000.,(y + 1.)*1000.));
#else
  double Bx = -0.85*xyvsbxin->GetBinContent(xyvsbxin->FindBin((x + 1.)*1000.,(y + 1.)*1000.));
  double By = -0.85*xyvsbyin->GetBinContent(xyvsbyin->FindBin((x + 1.)*1000.,(y + 1.)*1000.));
#endif
  // cout << " xpos " << x << " ypos " << y << " Bx " << Bx << " By " << By << endl;
  
  // double Bx = 0.;
  // double By = -uniformField;
  
  double maxzpos = (nlayer + 1)*ironThickness + nlayer*airGap;
  
  TVector3 magField(0.,0.,0.);
  if(z>=0. && z<=maxzpos) {
    if(fmod(z+airGap/2.,ironThickness+airGap)>airGap/2. && fmod(z+airGap/2.,ironThickness+airGap)<ironThickness+airGap/2.) {
      magField.SetXYZ(Bx,By,0.);
    }
  }
  return magField;
};	  // TVector3 GetMagneticField(double x, double y, double z) {


double GetMaterialDensity(double x, double y, double z) {
  double maxzpos = (nlayer + 1)*ironThickness + nlayer*airGap;
  double density = 0;
  if(z>=0. && z<=maxzpos) {
    if(fmod(z+airGap/2.,ironThickness+airGap)>airGap/2. && fmod(z+airGap/2.,ironThickness+airGap)<ironThickness+airGap/2.) {
      density = ironDensity;
    } else {
      density = gapDensity;
    }
  } else {
    density = airDensity;
  }
  return density;
};	  // double GetMaterialDensity(double x, double y, double z) {

#endif





const double H_Theta_BinW =  5.*TMath::Pi()/180.;
const double   H_Phi_BinW = 10.*TMath::Pi()/180.;
const double     H_R_BinW =  4.;
const int           ndf3D =  4;
const int maxHoughDist    =  nlayer;
const int MinHoughEntry   =  2;
const int MinHoughPt      =  1;
const int MaxHoughPt      = 10000;
const int MaxMeanHoughPt  = 4.833;
const int MaxStdHoughPt   = 0.15;


void getRPCId(UShort_t  rpcId,
	      Int_t    &module,
	      Int_t    &xrow,
	      Int_t    &yrow,
	      Int_t    &zlay) {
  // rpcId : 2 bit for INO module  
  //         3 bit for X-row
  //         3 bit for Y row
  //         8 bit for Z-layer
  zlay   = (rpcId    )&0xFF;//0b11111111;
  yrow   = (rpcId>> 8)&0b111;
  xrow   = (rpcId>>11)&0b111;
  module = (rpcId>>14)&0b11;
};

UShort_t constructRPCId(Int_t module,
			Int_t xrow,
			Int_t yrow,
			Int_t zlay) {
  UShort_t rpcId = module;
  rpcId        <<= 3;
  rpcId         +=  xrow;
  rpcId        <<= 3;
  rpcId         +=  yrow;
  rpcId        <<= 8;
  rpcId         +=  zlay;
  return rpcId;
};

void getTDC(UInt_t  tdc,
	    Bool_t &isTrail,	// leading or trailing
	    Bool_t &side,	// x or y
	    Int_t  &tdcno,	// tdc no 0-7
	    Int_t  &tdcval) {
  tdcno     = (tdc   )&0b1111;
  isTrail   = (tdc>>4)&0b1;
  tdcval    = (tdc>>5);
  side      = (tdcno>=nTDC)?1:0;
  if(tdcno>=nTDC) {tdcno-=nTDC;}
};

void getCAU(ULong64_t  tdc,
	    Bool_t    &isTrail,	// leading or trailing
	    UShort_t  &rpcId,
	    Int_t     &tdcval) {
  rpcId     = (tdc    )&0xFF;
  isTrail   = (tdc>>16)&0b1;
  tdcval    = (tdc>>17);
};

struct CauInfo {
  TTimeStamp       Cautime;
  vector<UShort_t> rpcId;
  vector<Double_t> val;
};

struct rawstrp {
  Int_t            strp;
  vector<Double_t> tdc[2];	// leading & trailing
};

struct rawlayer {
  int module, xrow, yrow, zlay;
  vector<double> rawTDC[nside][nTDC][2]; // raw TDC hits for debug purpose
  vector<rawstrp> hit[nside];	/* Hits are in reverse order */
  double tdc_ref[2];		// leading and trailing
  vector<vector<rawstrp>> cluster[nside]; /* Hits are in order */
};

struct PosInfo {
  Double_t pos;
  int multi;
};

struct PosInSpace {
  rawlayer     rawhitInfo;
  TVector3     RawHitPos;	// Middle of the Cluster
  TVector3     RawHitPos1;	// Fitting strip times
  TVector3     PosCorrection;
  
  // Double_t     oriPosXY[nside];
  // Double_t     posZ;
  // Double_t     posXY[nside];
  // Int_t        multiXY[nside];
  // Double_t     poserrXY[nside];
  // Double_t     timeXY[nside];
};

struct TrackInfo {
  TVector3 iMom;
  int charge;			// -1 -> mu-, +1 -> mu+
  // double ipos[3];
  TVector3         ipos;
  vector<TVector3> xyzpos;
  vector<UShort_t> xyzId;
  vector<TVector3> posCorr;
  vector<TVector3> posmulti;
  vector<TVector2> xyerr;
  vector<TVector3> xyztime;
  // vector<TVector3> timeCorr;
  vector<TVector2> xyterr;
  
  double opos[3];
  TVector3 fitMom;
  int ndf, ndfout;
  int layfirst, laylast;
  double chi2, xychi2[nside], whichMax, whichMin, chi2new;
  double trkAsymmetry;
  double energyLoss;		// in GeV
  double trackLength;		// in GeV
  vector<double> xyext[nside];	// nhits
  vector<double> minDist;	// nhits
  vector<TVector3> momext;
  vector<double> mulScAngle[nside]; // scattering angle in iron
};


struct GroupInfo {
  /* Contains only grouped hits, check allrawlay for all raw hits */
  vector<rawlayer>   alllayer;
  vector<PosInSpace> PosForTrack;
  vector<vector<PosInSpace>> HoughCluster;
  int maxHoughPt;
  double meanHough, stdHough;
};


struct HoughBinInfo {
  vector<int> posId;
  double xbin[ndf3D];
  int xbinI[ndf3D];
  double binC;
};

void getHoughCluster(vector<HoughBinInfo> &H_Bins1,
		     vector<HoughBinInfo> &outBins) {
  
  // cout << " entry " << H_Bins1.size() << endl;
  // vector<HoughBinInfo> H_Bins1 = H_Bins;
  // outBins.push_back(H_Bins1.front());
  // H_Bins1.erase(H_Bins1.begin());
  
  // for(int ki=0;ki<int(H_Bins1.size());ki++) {
  //   if(int(H_Bins1[ki].posId.size())>0) {
  //     outBins.push_back(H_Bins1.front());
  //     break;
  //   }
  //   H_Bins1.erase(H_Bins1.begin()+ki);
  // }
  
  for(int ki=0;ki<int(H_Bins1.size());ki++) {
    // if(int(H_Bins1[ki].posId.size())==0) {continue;}
    // for(int px=0;px<ndf3D;px++) {
    //   // cout << " test " << outBins.back().xbinI[px]
    //   // 	   << " " << H_Bins1[ki].xbinI[px] << endl;
    //   if(abs(outBins.back().xbinI[px] -
    // 	     H_Bins1[ki].xbinI[px])<=1) {
    // 	testb += 1;}
    // }

    for(int ij=0;ij<int(outBins.size());ij++) {

      int testb[2]  = {0};

      for(int pxx=0;pxx<2;pxx++) {
	if(abs(outBins[ij].xbinI[0] -
	       H_Bins1[ki].xbinI[0])<=1+pxx) { // R
	  testb[pxx] += 1;}
	if((abs(outBins[ij].xbinI[1] -
		H_Bins1[ki].xbinI[1])<=1+pxx) ||
	   (abs(outBins[ij].xbinI[1] -
		H_Bins1[ki].xbinI[1])>=int(2.*TMath::Pi()/H_Phi_BinW)-1-2*pxx)) { // Phi
	  testb[pxx] += 1;}
	if(abs(outBins[ij].xbinI[2] -
	       H_Bins1[ki].xbinI[2])<=1+pxx) { // Theta
	  testb[pxx] += 1;}
	if((abs(outBins[ij].xbinI[3] -
		H_Bins1[ki].xbinI[3])<=1+pxx) ||
	   (abs(outBins[ij].xbinI[3] -
		H_Bins1[ki].xbinI[3])>=int(2.*TMath::Pi()/H_Phi_BinW)-1-2*pxx)) { // Phi
	  testb[pxx] += 1;}
      }	// for(int pxx=0;pxx<2;pxx) {
      
      // cout << " some " << ki
      // 	   << " " << testb[0]
      // 	   << " " << testb[1]
      // 	   << endl;
      if((testb[0]==ndf3D) ||
	 (testb[0]==ndf3D-1 && testb[1]==ndf3D)) {
	// outBins.back().binC += H_Bins1[ki].binC;
	// H_Bins1[ki].binC  = 0;
	// for(int px=0;px<int(H_Bins1[ki].posId.size());px++) {
	// 	outBins.back().posId.push_back(H_Bins1[ki].posId[px]);
	// }
	// H_Bins1[ki].posId.clear();
	// H_Bins1.erase(H_Bins1.begin()+ki);
	outBins.push_back(H_Bins1[ki]);
	H_Bins1.erase(H_Bins1.begin()+ki);
	// cout << " entry1 " << H_Bins1.size() << endl;
	if(int(H_Bins1.size())>0) {
	  getHoughCluster(H_Bins1,outBins);
	}
	// outBins.insert(outBins.begin(),outBins.back());
	// outBins.pop_back();
	ki = 0;
	break;
      }
    } // for(int ij=0;ij<int(outBins.size());ij++) {
  } // for(int ki=jk+1;ki<int(H_Bins1.size());ki++) {

};



void calHitPos(PosInSpace &pos) {
  
  int nm = pos.rawhitInfo.module;
  int nx = pos.rawhitInfo.xrow;
  int ny = pos.rawhitInfo.yrow;
  int nl = pos.rawhitInfo.zlay;
  
  pos.RawHitPos  = getRPCpos(nm, nx, ny, nl)*(1./strpwidth);
  pos.RawHitPos1 = getRPCpos(nm, nx, ny, nl)*(1./strpwidth);
  
  int multi[nside];
  double calps[nside];
  double calps1[nside];
  for(int nj=0;nj<nside;nj++) {
    multi[nj] = pos.rawhitInfo.cluster[nj][0].size();
    calps[nj] = (pos.rawhitInfo.cluster[nj][0][0].strp
			      + multi[nj]*0.5);
    calps1[nj] = calps[nj];
    // cout << " lay " << pos.rawhitInfo.zlay << " " << nj
    // 	 << " " << calps[nj] << endl;
  } // for(int nj=0;nj<nside;nj++) {
  for(int nj=0;nj<nside;nj++) {
#ifndef isSimData
    vector<double> xpos, ytm;
    xpos.clear(); ytm.clear();
    int blkx = int(calps1[0]*blockM/nstrip);
    int blky = int(calps1[1]*blockM/nstrip);
    for(int cnt=0;cnt<multi[nj];cnt++) {
      if(int(pos.rawhitInfo.cluster[nj][0][cnt].tdc[0].size())) {
	xpos.push_back(pos.rawhitInfo.cluster[nj][0][cnt].strp + 0.5);
	// ytm.push_back(pos.rawhitInfo.cluster[nj][0][cnt].tdc[0][0]
	// 	      - pos.rawhitInfo.tdc_ref[0]);
	ytm.push_back(pos.rawhitInfo.cluster[nj][0][cnt].tdc[0][0]-
		      pos.rawhitInfo.tdc_ref[0]-
		      calps1[int(!nj)]*strpwidth*spdl-
		      blocktimeoff[nm][nx][ny][nl][nj][blkx][blky]-
		      multimeoff[nm][nx][ny][nl][nj][multi[nj]-1]
		      );
	// cout << " " << cnt
	//      << " " << xpos.back() << " " << ytm.back() << endl;
      }
    }	// for(int cnt=0;cnt<multi[nj];cnt++) {
    if(int(xpos.size())>=strpTimeFit) {
      TGraph gr1(int(xpos.size()),&xpos[0],&ytm[0]);
      int status = gr1.Fit(TString::Format("pol%i",(int(xpos.size())==3?2:3)),"SQ");
      // int status = gr1.Fit(TString::Format("strfn2","SQ");
      // int status = gr1.Fit(TString::Format("pol%i",int(xpos.size())-1),"SQ");
      if(status==0) {
	TF1 *f1 = (TF1 *)gr1.GetFunction(TString::Format("pol%i",(int(xpos.size())==3?2:3)));
	// TF1 *f1 = (TF1 *)gr1.GetFunction("strfn2");
	// TF1 *f1 = (TF1 *)gr1.GetFunction(TString::Format("pol%i",int(xpos.size())-1));
	double X0 = f1->GetMinimumX(xpos.front(),xpos.back());
	// cout << " lay " << pos.rawhitInfo.zlay
	//      << " " << nj << " " << X0 << endl;
	if(fabs(X0-calps1[nj])<multi[nj]*0.5-1.) {
	  calps1[nj] = X0;}
      }	// if(status==0) {
    }
#endif	// #ifndef isSimData
  } // for(int nj=0;nj<nside;nj++) {
  
  TVector3 posx0(calps[0],calps[1],0);
  pos.RawHitPos  += posx0;
  TVector3 posx1(calps1[0],calps1[1],0);
  pos.RawHitPos1 += posx1;
  
  int blkx = int(posx0.X()*blockM/nstrip);
  int blky = int(posx0.Y()*blockM/nstrip);
  pos.PosCorrection.SetX(blockposoff[nm][nx][ny][nl][0][blkx][blky]);
  pos.PosCorrection.SetY(blockposoff[nm][nx][ny][nl][1][blkx][blky]);
  pos.PosCorrection.SetZ(0.);
  
};


void LinearVectorFit(bool              isTime, // time iter
		     vector<TVector3>  pos,
		     vector<TVector2>  poserr,
		     vector<bool>      occulay,
		     TVector2         &slope,
		     TVector2         &inter,
		     TVector2         &chi2,
		     vector<TVector3> &ext,
		     vector<TVector3> &exterr) {
  
  double szxy[nside] = {0};
  double   sz[nside] = {0};
  double  sxy[nside] = {0};
  double   sn[nside] = {0};
  double  sz2[nside] = {0};
  
  double     slp[nside] = {-10000,-10000};
  double  tmpslp[nside] = {-10000,-10000};
  double intersect[nside] = {-10000,-10000};
  double    errcst[nside] = {-10000,-10000};
  double    errcov[nside] = {-10000,-10000};
  double    errlin[nside] = {-10000,-10000};
  
  for(int ij=0;ij<int(pos.size());ij++) {
    if(int(occulay.size()) && !occulay[ij]) {continue;}
    // cout << " ij " << ij << endl;
    double xyzval[3] = {pos[ij].X(),
			pos[ij].Y(),
			pos[ij].Z()};
    double xyerr[2]  = {poserr[ij].X(),
			poserr[ij].Y()};
    for(int nj=0;nj<nside;nj++) {
      szxy[nj] += xyzval[2]*xyzval[nj]/xyerr[nj];
      sz[nj]   += xyzval[2]/xyerr[nj];
      sz2[nj]  += xyzval[2]*xyzval[2]/xyerr[nj];
      sxy[nj]  += xyzval[nj]/xyerr[nj];
      sn[nj]   += 1/xyerr[nj];
    }   // for(int nj=0;nj<nside;nj++) {
  } // for(int ij=0;ij<int(pos.size());ij++){
  
  for(int nj=0;nj<nside;nj++) {
    if(sn[nj]>0. && sz2[nj]*sn[nj] - sz[nj]*sz[nj] !=0.) { 
      slp[nj] = (szxy[nj]*sn[nj] -
		   sz[nj]*sxy[nj])/(sz2[nj]*sn[nj] - sz[nj]*sz[nj]);
      tmpslp[nj] = slp[nj]; 
      if(isTime) { //time offset correction
        // if(fabs((cval*1.e-9)*slope+1)<3.30) { 
	tmpslp[nj] = -1./cval;
	// }
      }
      intersect[nj] = sxy[nj]/sn[nj] - tmpslp[nj]*sz[nj]/sn[nj];

      double determ = (sn[nj]*sz2[nj] - sz[nj]*sz[nj]);
      errcst[nj] = sz2[nj]/determ;
      errcov[nj] = -sz[nj]/determ;
      errlin[nj] = sn[nj]/determ;
    }
  } // for(int nj=0;nj<nside;nj++) {
  slope.SetX(tmpslp[0]);
  slope.SetY(tmpslp[1]);
  inter.SetX(intersect[0]);
  inter.SetY(intersect[1]);
  
  // theta = atan(sqrt(pow(tmpslp[0],2.)+pow(tmpslp[1],2.)));
  // phi = atan2(tmpslp[1],tmpslp[0]);
  
  double sumx = 0, sumy = 0;
  ext.clear(); exterr.clear();
  TVector3 xxt;
  TVector3 xxtt;
  for(int ij=0;ij<int(pos.size());ij++){
    xxt.SetX(tmpslp[0]*pos[ij].Z()+intersect[0]);
    xxt.SetY(tmpslp[1]*pos[ij].Z()+intersect[1]);
    xxt.SetZ(pos[ij].Z());
    ext.push_back(xxt);
    xxtt.SetX(errcst[0] + 2*errcov[0]*pos[ij].Z()+
	     errlin[0]*pos[ij].Z()*pos[ij].Z());
    xxtt.SetY(errcst[1] + 2*errcov[1]*pos[ij].Z()+
	     errlin[1]*pos[ij].Z()*pos[ij].Z());
    exterr.push_back(xxtt);
    // cout << " " << int(exterr.size())
    // 	 << " " << 1./exterr.back().X()
    // 	 << " " << 1./exterr.back().Y() << endl;
    if(int(occulay.size())==0 || occulay[ij]) {
      sumx += pow(xxt.X()-pos[ij].X(), 2.)/poserr[ij].X(); 
      sumy += pow(xxt.Y()-pos[ij].Y(), 2.)/poserr[ij].Y(); 
    }
  } // for(int ij=0;ij<int(pos.size());ij++){
  chi2.SetX(sumx);
  chi2.SetY(sumy);
};


/*
  This retuns the vector for minimum distance between two lines
  represented by (pos0, dir0) and (pos1, dir1).
  Note: Here dir0 and dir1 are unit vectors.
  
  This function can also be used to get minimum distance between
  a line and a point by making dir0=dir1.
*/
TVector3 getMinDist(TVector3 pos0,
		    TVector3 dir0,
		    TVector3 pos1,
		    TVector3 dir1) {
  
  double l0 = dir0.X();
  double m0 = dir0.Y();
  double n0 = dir0.Y();
  double l1 = dir1.X();
  double m1 = dir1.Y();
  double n1 = dir1.Z();
  
  double alpha =  (l0*(pos0.X() - pos1.X()) +
		   m0*(pos0.Y() - pos1.Y()) +
		   n0*(pos0.Z() - pos1.Z()));
	      
  double beta  =  (l1*(pos0.X() - pos1.X()) +
		   m1*(pos0.Y() - pos1.Y()) +
		   n1*(pos0.Z() - pos1.Z()));
  
  double gamma =  l0*l1 + m0*m1 + n0*n1;
  
  double r0 = (-alpha +  beta*gamma)/(1. - gamma*gamma);
  double r1 = ( beta  - alpha*gamma)/(1. - gamma*gamma);
  
  double px0 = pos0.X() + l0*r0;
  double py0 = pos0.Y() + m0*r0;
  double pz0 = pos0.Z() + n0*r0;
	      
  double qx1 = pos1.X() + l1*r1;
  double qy1 = pos1.Y() + m1*r1;
  double qz1 = pos1.Z() + n1*r1;
  
  TVector3 minDistV(qx1-px0,qy1-py0,qz1-pz0);
  return minDistV;
  
};				// TVector3 getMinDist(



#ifdef isRungeKutta

void ComputeRightHandSide(double* P,double* dPdS) {
  double P3vec[3]= {P[0], P[1], P[2]}; // Time is P[7]
  TVector3 MagField;
  MagField = GetMagneticField(P[0], P[1], P[2]);
  double particleCharge = P[6];
  double m_mom   = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5]); 
  double m_imom  = 1./m_mom;
  double m_cof   = particleCharge/*muCharge*cval*/*m_imom;                  
  dPdS[0] = P[3]*m_imom; // dx /ds
  dPdS[1] = P[4]*m_imom; // dy /ds
  dPdS[2] = P[5]*m_imom; // dz /ds
  dPdS[3] = m_cof*(P[4]*MagField.Z()-P[5]*MagField.Y()) ; // dPx/ds
  dPdS[4] = m_cof*(P[5]*MagField.X()-P[3]*MagField.Z()) ; // dPy/ds
  dPdS[5] = m_cof*(P[3]*MagField.Y()-P[4]*MagField.X()) ; // dPz/ds
  // cout<<"dpds"<<"   "<<dPdS[0]<<"   "<<dPdS[1]<<"   "<<dPdS[2]<<"   "<<dPdS[3]<<"   "<<dPdS[4]<<"   "<<dPdS[5]<<"   "<<MagField.X()<<"   "<<MagField.Y()<<"    "<<MagField.Z()<<endl;

} // void ComputeRightHandSide(double* P,double* dPdS) {


void Stepper(double* P, double *dPdS, double Step, double *Po, double *Err) {
  double R[3] = {P[0],   P[1] ,    P[2]};
  double A[3] = {dPdS[0], dPdS[1], dPdS[2]};
  double m_iPoint[3];
  double m_mPoint[3];
  double m_fPoint[3];
  TVector3 MagField;
  double particleCharge = P[6];
  double m_mom   = sqrt(P[3]*P[3]+P[4]*P[4]+P[5]*P[5]); 
  double m_imom  = 1./m_mom;
  double m_cof   = particleCharge/*muCharge*cval*/*m_imom;                  
  double m_lastField[3];
  m_iPoint[0]=R[0]; m_iPoint[1]=R[1]; m_iPoint[2]=R[2];
  
  const double one_sixth= 1./6.;
  double S  =      Step;
  double S5 =  0.5*Step;
  double S4 = 0.25*Step;
  double S6 =      Step*one_sixth;   // Step / 6.;
  
  // Point 1
  //
  double K1[3] = { m_imom*dPdS[3], m_imom*dPdS[4], m_imom*dPdS[5] };
  
  // Point2
  //
  double p[3] = {R[0]+S5*(A[0]+S4*K1[0]),
		 R[1]+S5*(A[1]+S4*K1[1]),
		 R[2]+S5*(A[2]+S4*K1[2])};
  MagField = GetMagneticField(p[0],p[1],p[2]);
  m_lastField[0] = MagField.X();
  m_lastField[1] = MagField.Y();
  m_lastField[2] = MagField.Z();
  // cout<<"check magfield"<<"   "<<p[0]<<" "<<p[1]<<"  "<<p[2]<<"   "<<m_lastField[0]<<"   "<< m_lastField[1]<<"   "<<m_lastField[2]<<"   "<<endl;
  double A2[3] = {A[0]+S5*K1[0],A[1]+S5*K1[1],A[2]+S5*K1[2]};
  double K2[3] = {(A2[1]*m_lastField[2]-A2[2]*m_lastField[1])*m_cof,
		    (A2[2]*m_lastField[0]-A2[0]*m_lastField[2])*m_cof,
		    (A2[0]*m_lastField[1]-A2[1]*m_lastField[0])*m_cof};
  
  m_mPoint[0]=p[0]; m_mPoint[1]=p[1]; m_mPoint[2]=p[2];
  
  // Point 3 with the same magnetic field
  //
  double A3[3] = {A[0]+S5*K2[0],A[1]+S5*K2[1],A[2]+S5*K2[2]};
  double K3[3] = {(A3[1]*m_lastField[2]-A3[2]*m_lastField[1])*m_cof,
		    (A3[2]*m_lastField[0]-A3[0]*m_lastField[2])*m_cof,
		    (A3[0]*m_lastField[1]-A3[1]*m_lastField[0])*m_cof};
  
  // Point 4
  //
  p[0] = R[0]+S*(A[0]+S5*K3[0]);
  p[1] = R[1]+S*(A[1]+S5*K3[1]);
  p[2] = R[2]+S*(A[2]+S5*K3[2]);             
  
  MagField = GetMagneticField(p[0],p[1],p[2]);
  m_lastField[0] = MagField.X();
  m_lastField[1] = MagField.Y();
  m_lastField[2] = MagField.Z();
  
  double A4[3] = {A[0]+S*K3[0],A[1]+S*K3[1],A[2]+S*K3[2]};
  double K4[3] = {(A4[1]*m_lastField[2]-A4[2]*m_lastField[1])*m_cof,
		    (A4[2]*m_lastField[0]-A4[0]*m_lastField[2])*m_cof,
		    (A4[0]*m_lastField[1]-A4[1]*m_lastField[0])*m_cof};
  
  // New position
  //
  Po[0] = P[0]+S*(A[0]+S6*(K1[0]+K2[0]+K3[0]));
  Po[1] = P[1]+S*(A[1]+S6*(K1[1]+K2[1]+K3[1]));
  Po[2] = P[2]+S*(A[2]+S6*(K1[2]+K2[2]+K3[2]));
  
  m_fPoint[0]=Po[0]; m_fPoint[1]=Po[1]; m_fPoint[2]=Po[2];
 
  // New direction
  //
  Po[3] = A[0]+S6*(K1[0]+K4[0]+2.*(K2[0]+K3[0]));
  Po[4] = A[1]+S6*(K1[1]+K4[1]+2.*(K2[1]+K3[1]));
  Po[5] = A[2]+S6*(K1[2]+K4[2]+2.*(K2[2]+K3[2]));
  
  // Errors
  //
  Err[3] = S*fabs(K1[0]-K2[0]-K3[0]+K4[0]);
  Err[4] = S*fabs(K1[1]-K2[1]-K3[1]+K4[1]);
  Err[5] = S*fabs(K1[2]-K2[2]-K3[2]+K4[2]);
  Err[0] = S*Err[3]                       ;
  Err[1] = S*Err[4]                       ;
  Err[2] = S*Err[5]                       ;
  Err[3]*= m_mom                          ;
  Err[4]*= m_mom                          ;
  Err[5]*= m_mom                          ;
  
  // Normalize momentum
  //
  double normF = m_mom/sqrt(Po[3]*Po[3]+Po[4]*Po[4]+Po[5]*Po[5]);
  Po [3]*=normF; Po[4]*=normF; Po[5]*=normF; 
  
  //cout<<"Stepper"<<"    "<<P[0]<<"   "<<P[1]<<"   "<<P[2]<<"   "<<Po[0]<<"   "<<Po[1]<<"   "<<Po[2]<<"   "<<endl;
  // Pass Energy, time unchanged -- time is not integrated !!
  //Po[6]=P[6]; Po[7]=P[7];

} // void Stepper (double* P, double *dPdS, double Step, double *Po, double *Err) {


void PropagateTrack(TrackInfo &inPoints) {
  
  double maxzpos = (nlayer + 1)*ironThickness + nlayer*airGap;
  
  const int totalPts = inPoints.xyzpos.size();
  inPoints.ndf = totalPts;
  inPoints.layfirst = nlayer;
  inPoints.laylast = -1;
  inPoints.chi2 = 0.;
  inPoints.minDist.clear();
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.minDist.push_back(0.);
  }
  inPoints.xyext[0].clear();inPoints.xyext[1].clear();
  inPoints.momext.clear();
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.xyext[0].push_back(-100.);
    inPoints.xyext[1].push_back(-100.);
    inPoints.momext.push_back(TVector3(-100,1,1));
  }
  inPoints.energyLoss = 0.;
  
  double in_Points[7] = {inPoints.ipos[0],
			 inPoints.ipos[1],
			 inPoints.ipos[2],
			 inPoints.iMom.X(),
			 inPoints.iMom.Y(),
			 inPoints.iMom.Z(),
			 double(inPoints.charge)};
  
  // cout << " ipos " << in_Points[0]///strptom
  //      << " " << in_Points[1]///strptom
  //      << " " << in_Points[2]///strptom
  //      << " mom " << in_Points[3]
  //      << " " << in_Points[4]
  //      << " " << in_Points[5]
  //      << " q " << in_Points[6] << endl;
  
  TVector3 MomDir;
  double temp_trk_wt[totalPts];
  double temp_xext[totalPts];
  double temp_yext[totalPts];
  double temp_minDist[totalPts];
  TVector3 temp_momext[nlayer];
  double totalLength;
  double eLoss = 0.;		// in GeV

  for(int ij=0;ij<totalPts;ij++) {
    temp_xext[ij] = -1000.;
    temp_yext[ij] = -1000.;
    temp_minDist[ij] = -1000.;
    temp_momext[ij].SetMagThetaPhi(-1000.,1,1);
    temp_trk_wt[ij] = 1.;
  }
  
  double dxdt, dydt, dzdt, dvxdt, dvydt, dvzt;
  
  double dPdS[6] = {0.};
  double out_Points[6] = {0.};
  double PError[6] = {0.};
  // double maxarc  = 4.; //maxium arc length of 4m.
  double isteparc = 0.;
  // double stepSize = 0.003;// 3mm
    
  int extFlag = nlayer-1;
  double totalTime = 0.;

  int nbreakpoint = 0;

  while(1) {
    
    nbreakpoint++;
    if(nbreakpoint==6000) break;
    
    // isteparc = sqrt(pow(out_Points[0]-in_Points[0],2)
    // 		    +pow(out_Points[1]-in_Points[1],2)
    // 		    +pow(out_Points[2]-in_Points[2],2));
    if(0
       // || isteparc>maxTrkLen
       || in_Points[2]<0
       || in_Points[2]>maxzpos
       || in_Points[0]<0
       || in_Points[0]>nstrip*strpwidth
       || in_Points[1]<0
       || in_Points[1]>nstrip*strpwidth) { break;}

    // cout << "  in " << in_Points[0]
    // 	 << " " << in_Points[1]
    // 	 << " " << in_Points[2] << endl;
          
    ComputeRightHandSide(in_Points, dPdS);
    Stepper(in_Points, dPdS, stepSize, out_Points, PError);
    
    // cout << " out " << out_Points[0]
    // 	 << " " << out_Points[1]
    // 	 << " " << out_Points[2]
    // 	 << endl << endl;
    
    MomDir.SetXYZ(out_Points[3],out_Points[4],out_Points[5]);

    ttxx.SetXYZ(out_Points[0],out_Points[1],out_Points[2]);
    ttyy.SetXYZ( in_Points[0], in_Points[1], in_Points[2]);
    double tStepDistance = calPointDist(ttxx,ttyy);
    // double tStepDistance = sqrt(pow(out_Points[0]-in_Points[0],2.)
    // 				+pow(out_Points[1]-in_Points[1],2.)
    // 				+pow(out_Points[2]-in_Points[2],2.));
    totalTime += tStepDistance;
    // cout << " 5 totalTime " << totalTime << endl;
    // cout << " tStepDistance " << 1000.*tStepDistance << endl;
    if(totalTime>maxTrkLen) {
      // cout << " 5 totalTime " << totalTime << endl;
      // inPoints.chi2 = 1000.;
      break;
    }
    
    // cout << " vals " << in_Points[0]/strptom
    // 	 << " " << in_Points[1]/strptom
    // 	 << " " << in_Points[2]/strptom << endl;
    
    // if(in_Points[0]<0. || in_Points[0]>stripwidth*nstrip
    //    || in_Points[1]<0. || in_Points[1]>nstrip*stripwidth
    //    || in_Points[2]<0. || in_Points[2]>zz[nlayer-1]+airGap/2.+ironThickness) {break;}

    for(int ij=0;ij<totalPts;ij++) {
      double lowz  = TMath::Min(in_Points[2], out_Points[2]);
      double highz = TMath::Max(in_Points[2], out_Points[2]);
      if(lowz  <= inPoints.xyzpos[ij].Z() &&
	 highz >= inPoints.xyzpos[ij].Z()) {
	/** interpolate x and y in the layer */
	temp_xext[ij] = in_Points[0] +
	  (out_Points[0] - in_Points[0]) * ((inPoints.xyzpos[ij].Z() - in_Points[2]) /
					    (out_Points[2]           - in_Points[2]));
	temp_yext[ij] = in_Points[1] +
	  (out_Points[1] - in_Points[1]) * ((inPoints.xyzpos[ij].Z() - in_Points[2]) /
					    (out_Points[2]           - in_Points[2]));
	temp_minDist[ij] =
	  sqrt(pow(temp_xext[ij] - inPoints.xyzpos[ij].X(),2.) +
	       pow(temp_yext[ij] - inPoints.xyzpos[ij].Y(),2.));
	temp_momext[ij] = MomDir;
    	// temp_trk_wt[ij] = TMath::Exp(-totalTime);
    	// temp_trk_wt[ij] = 1./sqrt(1. + totalTime/(ironThickness + airGap));
    	temp_trk_wt[ij] = 1.;
    	// temp_trk_wt[ij] = 1./pow(TMath::E(),sqrt(totalTime));
	// cout << " ij " << ij
	//      << " xext " << temp_xext[ij]/strpwidth
	//      << " xext " << temp_yext[ij]/strpwidth
	//      << endl;
	break;
      }
    } // for(int ij=0;ij<totalPts;ij++) {
        
    double tmpMag = MomDir.Mag();
#ifdef isEloss
    double tmpELoss = (100.*tStepDistance)*(eLossCurve->Eval(1000.*tmpMag)*GetMaterialDensity(out_Points[0],out_Points[1],out_Points[2])*0.001);
    eLoss += tmpELoss;
    tmpMag -= tmpELoss;
#endif	// #ifdef isEloss
    if(tmpMag<minPartE) {
      // cout << " tmpMag " << tmpMag << " lay " << extFlag << endl;
      break;
    }
    // cout << " tmpMag " << tmpMag << " eloss " << eLoss << endl;
    // cout << " tmpELoss " << tmpELoss << " stepLoss " << eLossCurve->Eval(1000.*tmpMag) << " tmpMag " << tmpMag << " eloss " << eLoss << endl;
    
    MomDir.SetMag(tmpMag);
    
    for(int ij=0;ij<3;ij++) {
      in_Points[ij] = out_Points[ij];
    }
    in_Points[3] = MomDir.X();
    in_Points[4] = MomDir.Y();
    in_Points[5] = MomDir.Z();
    
    // TVector3 MagField;
    // cout << " " << in_Points[0]
    // 	 << " " << in_Points[1]
    // 	 << " " << in_Points[2]
    // 	 << " " << MomDir.Mag()
    // 	 << " " << MomDir.Theta()*180./TMath::Pi()
    //   // << " " << zz[nlayer-1]+airGap/2.+ironThickness
    // 	 << " " << GetMagneticField(in_Points[0],in_Points[1],in_Points[2]).Y()
    // 	 << endl;
    
    // cout << " Final  " << endl;
    // if(extFlag<=0) { break;}
  }	// while(1) {
  
  totalLength = totalTime;
  
  double chi2Each = 0.;
  double chi2EachMax = 0.;
  double xychi2Each[nside] = {0};
  double sumW = 0;
  double xysumW[nside] = {0};
  double leastDist = 100., mostDist = 0;
  for(int ij=0;ij<totalPts;ij++) {
    if(temp_xext[ij]>-1000.) {
      double ttDistEr =
	pow(temp_xext[ij]-inPoints.xyzpos[ij].X(),2.) /
	inPoints.xyerr[ij].X();
      xychi2Each[0] += ttDistEr * temp_trk_wt[ij];
    }
    if(temp_yext[ij]>-1000.) {
      double ttDistEr =
	pow(temp_yext[ij]-inPoints.xyzpos[ij].Y(),2.) /
	inPoints.xyerr[ij].Y();
      xychi2Each[1] += ttDistEr * temp_trk_wt[ij];
    }
    
    if(temp_xext[ij]>-1000. && temp_yext[ij]>-1000.) {
      double ttDistEr =
	(pow(temp_xext[ij]-inPoints.xyzpos[ij].X(),2.)/inPoints.xyerr[ij].X() +
	 pow(temp_yext[ij]-inPoints.xyzpos[ij].Y(),2.)/inPoints.xyerr[ij].Y());
      // cout<<" ij "<<ij
      // 	  <<" xdev "<<temp_xext[ij]-inPoints.xyzpos[ij].X()
      // 	  <<" ydev "<<temp_yext[ij]-inPoints.xyzpos[ij].Y()
      // 	  <<" ttDistEr "<<ttDistEr<<" temp_trk_wt "<<temp_trk_wt[ij]<<endl;

      chi2Each += ttDistEr * temp_trk_wt[ij];
      chi2EachMax = TMath::Max(chi2EachMax,ttDistEr);

      sumW += 1./(inPoints.xyerr[ij].X()+inPoints.xyerr[ij].Y());
      if(leastDist>ttDistEr) {leastDist = ttDistEr;}
      if(mostDist<ttDistEr) {mostDist = ttDistEr;}
    }
    // cout << " ij " << ij
    // 	 << " x " << inPoints.hits[ij].posX/strptom
    // 	 << " y " << inPoints.hits[ij].posY/strptom
    // 	 << " z " << inPoints.hits[ij].posZi
    // 	 << " " << temp_minDist[ij]/strptom << endl;
    // cout << " ij " << ij
    // 	 << " z " << inPoints.hits[ij].posZi
    // 	 << " " << temp_minDist[ij]
    //   // << " " << pow(temp_minDist[ij],2.)
    // 	 << " errx " << sqrt(inPoints.hits[ij].poserrXY[0])
    // 	 << " erry " << sqrt(inPoints.hits[ij].poserrXY[1])
    // 	 << " chi2 " << chi2Each
    // 	 << endl;
  }
  // cout << " \t\tfinChi2 " << chi2Each << endl;
  
  // if(sumW>0) {
  //   chi2Each /= sumW;
  // }
  // if(xychi2Each[0]>0) {
  //   xychi2Each[0] /= xysumW[0];
  //   xychi2Each[1] /= xysumW[1];
  // }
  
  
  /********* Writting Output *********/
  
  inPoints.energyLoss = eLoss;
  inPoints.chi2 = chi2Each==0.?100000.:chi2Each;
  inPoints.whichMax = chi2EachMax==0.?1000.:chi2EachMax;
  inPoints.xychi2[0] = xychi2Each[0]==0.?100000.:xychi2Each[0];
  inPoints.xychi2[1] = xychi2Each[1]==0.?100000.:xychi2Each[1];
  inPoints.trkAsymmetry = (mostDist-leastDist==-100)?100000.:(mostDist-leastDist)*inPoints.chi2;
  
  inPoints.trackLength = totalLength;
  
  inPoints.ndfout = 0;
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.minDist[ij] = temp_minDist[ij];
    if(temp_minDist[ij]<maxPosDev*strpwidth) {inPoints.ndfout++;}
    // cout << " dist " << ij << " " << temp_minDist[ij]/strpwidth << endl;
    // cout << ij << endl;
  }
  // cout << inPoints.ndfout << endl;
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.xyext[0][ij] = temp_xext[ij];
    inPoints.xyext[1][ij] = temp_yext[ij];
    inPoints.momext[ij].SetMagThetaPhi(temp_momext[ij].Mag(),temp_momext[ij].Theta(),temp_momext[ij].Phi());
    if(temp_momext[ij].Mag()<10.) {
      if(ij>inPoints.laylast) inPoints.laylast = ij;
      if(ij<inPoints.layfirst) inPoints.layfirst = ij;
    }
    // cout << " lay " << ij << " X " << temp_xext[ij]/strpwidth << " Y " << temp_yext[ij]/strpwidth << " mom " << temp_momext[ij].Mag()*cval/gevtojoule << " theta " << temp_momext[ij].Theta()*180/TMath::Pi() << " phi " << temp_momext[ij].Phi()*180/TMath::Pi() << endl;
  }
  // cout << " first " << inPoints.layfirst << " last " << inPoints.laylast << endl;
  // inPoints.fitMom = MomDir;
  // cout << " mom " << inPoints.fitMom.Mag()*cval/gevtojoule << endl;
  
};			 // void propagateTrack(TrackInfo &inPoints) {



#else  // #ifdef isRungeKutta



void PropagateTrack(TrackInfo &inPoints) {
  
  double maxzpos = (nlayer + 1)*ironThickness + nlayer*airGap;

  const int totalPts = inPoints.xyzpos.size();
  inPoints.ndf = totalPts;
  inPoints.layfirst = nlayer;
  inPoints.laylast = -1;
  inPoints.chi2 = 0.;
  inPoints.minDist.clear();
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.minDist.push_back(0.);
  }
  inPoints.xyext[0].clear();inPoints.xyext[1].clear();
  inPoints.momext.clear();
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.xyext[0].push_back(-100.);
    inPoints.xyext[1].push_back(-100.);
    inPoints.momext.push_back(TVector3(-100,1,1));
  }
  inPoints.energyLoss = 0.;
  

  
  TVector3 MomDir;
  MomDir.SetMagThetaPhi(inPoints.iMom.Mag()*gevtojoule/cval1,inPoints.iMom.Theta(),inPoints.iMom.Phi());
  double X0 = inPoints.ipos[0];
  double Y0 = inPoints.ipos[1];
  double Z0 = inPoints.ipos[2];//zz[laylast] + airGap/2. + ironThickness;
  double fitCharge = muCharge*inPoints.charge;
  
  // cout << " ipos " << X0/strpwidth << " " << Y0/strpwidth << " " << Z0/strpwidth << " mom " << inPoints.iMom.Mag() << " " << inPoints.iMom.Theta()*180./TMath::Pi() << " " << inPoints.iMom.Phi()*180./TMath::Pi() << " q " << inPoints.charge << endl;

  double temp_trk_wt[totalPts];
  double temp_xext[totalPts];
  double temp_yext[totalPts];
  double temp_minDist[totalPts];
  TVector3 temp_momext[totalPts];
  double totalLength;
  double eLoss = 0.;		// in GeV
  
  double xVal[2], yVal[2], zVal[2];
  double vxVal[2], vyVal[2], vzVal[2];
  double axVal[2], ayVal[2], azVal[2];
  TVector3 MagField;
  double gammaF;
    
  for(int ij=0;ij<totalPts;ij++) {
    temp_xext[ij] = -1000.;
    temp_yext[ij] = -1000.;
    temp_minDist[ij] = -1000.;
    temp_momext[ij].SetMagThetaPhi(-1000.,1,1);
    temp_trk_wt[ij] = 1.;
  }
  
  xVal[0] = X0;
  yVal[0] = Y0;
  zVal[0] = Z0;
  // MomDir.SetMagThetaPhi(Magn*gevtojoule/cval1,iTheta,iPhi);
  gammaF = sqrt(1. + pow(MomDir.Mag()/(muMass*cval1),2.));
  MagField = GetMagneticField(xVal[0],yVal[0],zVal[0]);
  vxVal[0] = MomDir.X()/(gammaF*muMass);
  vyVal[0] = MomDir.Y()/(gammaF*muMass);
  vzVal[0] = MomDir.Z()/(gammaF*muMass);
  axVal[0] = (fitCharge/(gammaF*muMass))*(vyVal[0]*MagField.Z()-vzVal[0]*MagField.Y());
  ayVal[0] = (fitCharge/(gammaF*muMass))*(vzVal[0]*MagField.X()-vxVal[0]*MagField.Z());
  azVal[0] = (fitCharge/(gammaF*muMass))*(vxVal[0]*MagField.Y()-vyVal[0]*MagField.X());
  
  int extFlag = nlayer-1;
  double totalTime = 0.;
  
  while(1) {
    
    // xVal[1] = xVal[0]+vxVal[0]*timeStepFit*(direction==false?1.:-1.);
    // yVal[1] = yVal[0]+vyVal[0]*timeStepFit*(direction==false?1.:-1.);
    // zVal[1] = zVal[0]+vzVal[0]*timeStepFit*(direction==false?1.:-1.);
    xVal[1] = xVal[0]+vxVal[0]*timeStepFit;
    yVal[1] = yVal[0]+vyVal[0]*timeStepFit;
    zVal[1] = zVal[0]+vzVal[0]*timeStepFit;
    
    // double tStepDistance = sqrt(pow(xVal[1]-xVal[0],2.)+pow(yVal[1]-yVal[0],2.)+pow(zVal[1]-zVal[0],2.));
    ttxx.SetXYZ(xVal[1],yVal[1],zVal[1]);
    ttyy.SetXYZ(xVal[0],yVal[0],zVal[0]);
    double tStepDistance = calPointDist(ttxx,ttyy);
    totalTime += tStepDistance;
    // cout << " 5 totalTime " << totalTime << endl;
    // cout << " tStepDistance " << 1000.*tStepDistance << endl;
    if(totalTime>maxTrkLen) {
      // cout << " 5 totalTime " << totalTime << endl;
      // inPoints.chi2 = 1000.;
      break;
    }
    
    // cout << " vals0 " << xVal[0]/strpwidth << " " << yVal[0]/strpwidth << " " << zVal[0]/strpwidth << endl;
    // cout << " vals1 " << xVal[1]/strpwidth << " " << yVal[1]/strpwidth << " " << zVal[1]/strpwidth << endl;
    // cout << " vvals0 " << vxVal[0] << " " << vyVal[0] << " " << vzVal[0] << endl;
    
    // if(xVal[0]<0. || xVal[0]>stripwidth*nstrip
    //    || yVal[0]<0. || yVal[0]>nstrip*stripwidth
    //    || zVal[0]<0. || zVal[0]>zz[nlayer-1]+airGap/2.+ironThickness) {break;}
    
    if(xVal[0]<0. || xVal[0]>nstrip*strpwidth
       || yVal[0]<0. || yVal[0]>nstrip*strpwidth
       || zVal[0]<0. || zVal[0]>maxzpos) {break;}
    
            
    // vxVal[1] = vxVal[0]+axVal[0]*timeStepFit*(direction==false?1.:-1.);
    // vyVal[1] = vyVal[0]+ayVal[0]*timeStepFit*(direction==false?1.:-1.);
    // vzVal[1] = vzVal[0]+azVal[0]*timeStepFit*(direction==false?1.:-1.);
    vxVal[1] = vxVal[0]+axVal[0]*timeStepFit;
    vyVal[1] = vyVal[0]+ayVal[0]*timeStepFit;
    vzVal[1] = vzVal[0]+azVal[0]*timeStepFit;
    
    double tmpMag = MomDir.Mag()*cval1/gevtojoule;
#ifdef isEloss
    double tmpELoss = (100.*tStepDistance)*(eLossCurve->Eval(1000.*tmpMag)*GetMaterialDensity(xVal[1],yVal[1],zVal[1])*0.001);
    // double tmpELoss = (100.*tStepDistance)*(2.*GetMaterialDensity(xVal[1],yVal[1],zVal[1])*0.001);
    // eLoss += tmpELoss*(direction==false?1.:-1.);
    // tmpMag -= tmpELoss*(direction==false?1.:-1.);
    eLoss += tmpELoss;
    tmpMag -= tmpELoss;
#endif	// #ifdef isEloss
    if(tmpMag<minPartE) {
      // cout << " tmpMag " << tmpMag << " lay " << extFlag << endl;
      break;
    }
    // cout << " tmpMag " << tmpMag << " eloss " << eLoss << endl;
    // cout << " tmpELoss " << tmpELoss << " stepLoss " << eLossCurve->Eval(1000.*tmpMag) << " tmpMag " << tmpMag << " eloss " << eLoss << endl;
    
    MomDir.SetXYZ(gammaF*muMass*vxVal[1],gammaF*muMass*vyVal[1],gammaF*muMass*vzVal[1]);
    // double tmpMag = pow(10.,Magn);//-totalTime*eLoss;
    // if(tmpMag<minPartE) break;
    MomDir.SetMag(tmpMag*gevtojoule/cval1);
    
    gammaF = sqrt(1. + pow(MomDir.Mag()/(muMass*cval1),2.));
    
    // cout << " gamma " << gammaF << endl;
    
    vxVal[1] = MomDir.X()/(gammaF*muMass);
    vyVal[1] = MomDir.Y()/(gammaF*muMass);
    vzVal[1] = MomDir.Z()/(gammaF*muMass);
    
    MagField = GetMagneticField(xVal[1],yVal[1],zVal[1]);
    axVal[1] = (fitCharge/(gammaF*muMass))*(vyVal[1]*MagField.Z()-vzVal[1]*MagField.Y());
    ayVal[1] = (fitCharge/(gammaF*muMass))*(vzVal[1]*MagField.X()-vxVal[1]*MagField.Z());
    azVal[1] = (fitCharge/(gammaF*muMass))*(vxVal[1]*MagField.Y()-vyVal[1]*MagField.X());
    
    for(int ij=0;ij<totalPts;ij++) {
      double lowz  = TMath::Min(zVal[0], zVal[1]);
      double highz = TMath::Max(zVal[0], zVal[1]);
      if(lowz  <= inPoints.xyzpos[ij].Z() &&
	 highz >= inPoints.xyzpos[ij].Z()) {
	/** interpolate x and y in the layer */
	temp_xext[ij] = xVal[0] +
	  (xVal[1] - xVal[0]) * ((inPoints.xyzpos[ij].Z() - zVal[0]) /
				 (zVal[1]                 - zVal[0]));
        temp_yext[ij] = yVal[0] +
	  (yVal[1] - yVal[0]) * ((inPoints.xyzpos[ij].Z() - zVal[0]) /
				 (zVal[1]                 - zVal[0]));
	temp_minDist[ij] =
	  sqrt(pow(temp_xext[ij] - inPoints.xyzpos[ij].X(),2.) +
	       pow(temp_yext[ij] - inPoints.xyzpos[ij].Y(),2.));
	temp_momext[ij] = MomDir;
    	// temp_trk_wt[ij] = TMath::Exp(-totalTime);
    	// temp_trk_wt[ij] = 1./sqrt(1. + totalTime/(ironThickness + airGap));
    	temp_trk_wt[ij] = 1.;
    	// temp_trk_wt[ij] = 1./pow(TMath::E(),sqrt(totalTime));
	// cout << " ij " << ij
	//      << " xext " << temp_xext[ij]/strpwidth
	//      << " xext " << temp_yext[ij]/strpwidth
	//      << endl;
	break;
      }
    } // for(int ij=0;ij<totalPts;ij++) {
    
    xVal[0] = xVal[1];
    yVal[0] = yVal[1];
    zVal[0] = zVal[1];
    vxVal[0] = vxVal[1];
    vyVal[0] = vyVal[1];
    vzVal[0] = vzVal[1];
    axVal[0] = axVal[1];
    ayVal[0] = ayVal[1];
    azVal[0] = azVal[1];
    
    // cout << " " << xVal[0]
    // 	 << " " << yVal[0]
    // 	 << " " << zVal[0]
    //   // << " " << MomDir.Mag()
    //   // << " " << MomDir.Theta()*180./Pi()
    //   // << " " << zz[nlayer-1]+airGap/2.+ironThickness
    //   // << " " << MagField.Y()
    // 	 << endl;
    
    // cout << " Final  " << endl;
    // if(extFlag<=0) { break;}
  }	// while(1) {
  
  totalLength = totalTime;
  
  double chi2Each = 0.;
  double chi2EachNew = 0.;
  double chi2EachMax = 0.;
  double chi2EachMin = 1000.;
  double xychi2Each[nside] = {0};
  double sumW = 0;
  double xysumW[nside] = {0};
  double leastDist = 100., mostDist = 0;
  for(int ij=0;ij<totalPts;ij++) {
    if(temp_xext[ij]>-1000.) {
      double ttDistEr =
	pow(temp_xext[ij]-inPoints.xyzpos[ij].X(),2.) /
	inPoints.xyerr[ij].X();
      xychi2Each[0] += ttDistEr * temp_trk_wt[ij];
    }
    if(temp_yext[ij]>-1000.) {
      double ttDistEr =
	pow(temp_yext[ij]-inPoints.xyzpos[ij].Y(),2.) /
	inPoints.xyerr[ij].Y();
      xychi2Each[1] += ttDistEr * temp_trk_wt[ij];
    }

    if(temp_xext[ij]>-1000. && temp_yext[ij]>-1000.) {
      double ttDistEr =
	(pow(temp_xext[ij]-inPoints.xyzpos[ij].X(),2.)/inPoints.xyerr[ij].X() +
	 pow(temp_yext[ij]-inPoints.xyzpos[ij].Y(),2.)/inPoints.xyerr[ij].Y());
      // cout<<" ij "<<ij
      // 	  <<" xdev "<<temp_xext[ij]-inPoints.xyzpos[ij].X()
      // 	  <<" ydev "<<temp_yext[ij]-inPoints.xyzpos[ij].Y()
      // 	  <<" ttDistEr "<<ttDistEr<<" temp_trk_wt "<<temp_trk_wt[ij]<<endl;

      chi2EachMax = TMath::Max(chi2EachMax,ttDistEr);
      chi2EachMin = TMath::Max(chi2EachMin,ttDistEr);

      chi2Each += ttDistEr * temp_trk_wt[ij];
      chi2EachNew = chi2EachNew==0?ttDistEr:chi2EachNew*ttDistEr;
      
      sumW += 1./(inPoints.xyerr[ij].X()+inPoints.xyerr[ij].Y());
      if(leastDist>ttDistEr) {leastDist = ttDistEr;}
      if(mostDist<ttDistEr) {mostDist = ttDistEr;}
    }
    // cout << " ij " << ij
    // 	 << " x " << inPoints.hits[ij].posX/strpwidth
    // 	 << " y " << inPoints.hits[ij].posY/strpwidth
    // 	 << " z " << inPoints.hits[ij].posZi
    // 	 << " " << temp_minDist[ij]/strpwidth << endl;
    // cout << " ij " << ij
    // 	 << " z " << inPoints.hits[ij].posZi
    // 	 << " " << temp_minDist[ij]
    //   // << " " << pow(temp_minDist[ij],2.)
    // 	 << " errx " << sqrt(inPoints.hits[ij].poserrXY[0])
    // 	 << " erry " << sqrt(inPoints.hits[ij].poserrXY[1])
    // 	 << " chi2 " << chi2Each
    // 	 << endl;
  }
  // cout << " \t\tfinChi2 " << chi2Each << endl;

  // if(sumW>0) {
  //   chi2Each /= sumW;
  // }
  // if(xychi2Each[0]>0) {
  //   xychi2Each[0] /= xysumW[0];
  //   xychi2Each[1] /= xysumW[1];
  // }
  
  /********* Writting Output *********/
  
  inPoints.energyLoss = eLoss;
  inPoints.chi2 = chi2Each==0.?1000.:chi2Each;
  inPoints.chi2new = chi2EachNew==0.?1000.:chi2EachNew;
  inPoints.whichMin = chi2EachMin==0.?1000.:chi2EachMin;
  inPoints.whichMax = chi2EachMax==0.?1000.:chi2EachMax;
  inPoints.xychi2[0] = xychi2Each[0]==0.?1000.:xychi2Each[0];
  inPoints.xychi2[1] = xychi2Each[1]==0.?1000.:xychi2Each[1];
  inPoints.trkAsymmetry = (mostDist-leastDist==-100)?1000.:(mostDist-leastDist)*inPoints.chi2;
  
  inPoints.trackLength = totalLength;
  
  inPoints.ndfout = 0;
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.minDist[ij] = temp_minDist[ij];
    if(temp_minDist[ij]<maxPosDev*strpwidth) {inPoints.ndfout++;}
    // cout << " dist " << ij << " " << temp_minDist[ij]/strpwidth << endl;
    // cout << ij << endl;
  }
  // cout << inPoints.ndfout << endl;
  for(int ij=0;ij<totalPts;ij++) {
    inPoints.xyext[0][ij] = temp_xext[ij];
    inPoints.xyext[1][ij] = temp_yext[ij];
    inPoints.momext[ij].SetMagThetaPhi(temp_momext[ij].Mag()*cval1/gevtojoule,temp_momext[ij].Theta(),temp_momext[ij].Phi());
    if(temp_momext[ij].Mag()*cval1/gevtojoule<10.) {
      if(ij>inPoints.laylast) inPoints.laylast = ij;
      if(ij<inPoints.layfirst) inPoints.layfirst = ij;
    }
    // cout << " lay " << ij << " X " << temp_xext[ij]/strpwidth << " Y " << temp_yext[ij]/strpwidth << " mom " << temp_momext[ij].Mag()*cval1/gevtojoule << " theta " << temp_momext[ij].Theta()*180/TMath::Pi() << " phi " << temp_momext[ij].Phi()*180/TMath::Pi() << endl;
  }
  // cout << " first " << inPoints.layfirst << " last " << inPoints.laylast << endl;
  // inPoints.fitMom = MomDir;
  // cout << " mom " << inPoints.fitMom.Mag()*cval1/gevtojoule << endl;
  
};			 // void propagateTrack(TrackInfo &inPoints) {


#endif	// #ifdef isRungeKutta


TrackInfo inPoints1;

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

  inPoints1.ipos[0] = par[0];
  inPoints1.ipos[1] = par[1];
  // inPoints1.iMom.SetXYZ(par[2],par[3],par[4]);
  inPoints1.iMom.SetMagThetaPhi(par[2],par[3],par[4]);

  // cout << " " << par[0]
  //      << " " << par[1]
  //      << " " << par[2]
  //      << " " << par[3]
  //      << " " << par[4]
  //      // << " " << f
  //      << endl;

  PropagateTrack(inPoints1);

  f = TMath::Max(inPoints1.xychi2[0], inPoints1.xychi2[1]);
  
  // cout << " " << par[0]
  //      << " " << par[1]
  //      << " " << par[2]
  //      << " " << par[3]
  //      << " " << par[4]
  //      << " " << f
  //      << endl;

};












int main(int argc, char** argv) {

  /* 
     argv[0] : main
     argv[1] : filename
     argv[2] : start event
     argv[3] : number of event from begining
     argv[4] : output file number
     argv[5] : itermodule
     argv[6] : iterxrow
     argv[7] : iteryrow
     argv[8] : iterlayer
     argv[9] :
     argv[10] :
 */

  
#ifdef isIter
  int itermodule = stoi(argv[5]);
  int iterxrow = stoi(argv[6]);
  int iteryrow = stoi(argv[7]);
  int iterlayer = stoi(argv[8]);
#endif	// #ifdef isIter
  

  const char *sideMark[nside] = {"x","y"};

  TF1 *strfn2 = new TF1("strfn2","[0]+[1]*x*x+[2]*x*x*x*x",0,nstrip);
  
  eLossCurve = new TGraph("muon_energy_loss_in_fe.txt");
  
  /** Magnetic Field Histograms **/
  TFile *fileB = new TFile("B_mical_hist.root", "read");
  xyvsbxin = (TH2D*)((TH2D*)fileB->Get("xyvsbxin"))->Clone("xyvsbxin");
  xyvsbyin = (TH2D*)((TH2D*)fileB->Get("xyvsbyin"))->Clone("xyvsbyin");
  // fileB->Close();
  
// #ifndef isCorrection
//   unsigned int triggerinfo_ref = 0;
//   for(int ij=0;ij<ntrigLayers;ij++) {
//     triggerinfo_ref<<=1;
//     triggerinfo_ref+=1;}
// #else  // #ifndef isCorrection
//   unsigned int triggerinfo_ref1 = 0;
//   for(int ij=0;ij<ntrigLayers1;ij++) {
//     triggerinfo_ref1<<=1;
//     triggerinfo_ref1+=1;}
// #endif	// #ifndef isCorrection
  
  unsigned int triggerinfo_ref1 = 0;
  for(int ij=0;ij<ntrigLayers1;ij++) {
    triggerinfo_ref1<<=1;
    triggerinfo_ref1+=1;}



  Long64_t start_s = clock();
  
#ifndef isSimData
  CauStore  *cau   = 0;
  DigiStore *event = 0;
  
  CauInfo tmpCau;
  vector<CauInfo> allCAUevents;
  double cauCorrections[nmodule][nxrow][nyrow][nlayer] = {0};
#endif	// #ifndef isSimData
  
  rawstrp            tmphit;
  // vector<rawstrp>    hit[nside][nlayer];
  vector<rawstrp>    tmpCluster;
  rawlayer           tmplayer;
  vector<rawlayer>   allrawlay;
    
  GroupInfo          tmpevt;
  vector<GroupInfo>  TimeGroup;
  
  PosInSpace         tmpPos;
  
  vector<PosInSpace> tmpHoughGroup;
   
  Double_t           tdc_ref[nlayer][2]; // leading and trailing
  Double_t           ttExecTimeVal;

  vector<TVector3> xyzpos, xyzpos1;
  vector<UShort_t> xyzId;
  vector<TVector3> posCorr;
  vector<TVector3> posmulti;
  vector<TVector3> totLayHit;
  vector<bool>     occulay;
  vector<TVector2> xyerr;
  vector<TVector3> xyzerr;
  vector<TVector3> xyext;
  vector<TVector3> xyexterr;
  TVector2         xychi2;
  TVector2         slope;
  TVector2         inter;
  double thetad, phid;
  
  vector<TVector3> xyztime;
  vector<TVector3> timeCorr;
  vector<TVector2> xyterr;
  vector<bool>     tocculay;
  TVector2         tslope;
  TVector2         tinter;
  vector<TVector3> xytext;
  vector<TVector3> xytexterr;
  TVector2         xytchi2;
  

  
  char name[300];
  
  // vector<PosInfo>    hitpos[nside][nlayer]; // Hit pos and Multiplicity
  //, PosAfterHough, PosAfterTime, PosAfterHou
  
  TH1D *strptimeraw = new TH1D("strptimeraw","strptimeraw",
			       int(maxtime/10.),0.,maxtime);
  TH1D *hbinIns = new TH1D("hbinIns","hbinIns",
			   50,0.5,50.5);
  
  TH1D *xybeta[nside];
  for(int nj=0;nj<nside;nj++) {
    sprintf(name,"xybeta_%s",sideMark[nj]);
    xybeta[nj] = new TH1D(name,name,320,-5.,3.);
  } // for(int nj=0;nj<nside;nj++) {

  TH1D *hxychi2[nside];
  for(int nj=0;nj<nside;nj++) {
    sprintf(name,"xychi2_%s",sideMark[nj]);
    hxychi2[nj] = new TH1D(name,name,500,0.,50.);
  } // for(int nj=0;nj<nside;nj++) {

  
  TH1D *circleS[nside];
  for(int nj=0;nj<nside;nj++) {
    sprintf(name,"circleS_%s",sideMark[nj]);
    circleS[nj] = new TH1D(name,name,200,0.,5.);
  } // for(int nj=0;nj<nside;nj++) {
  

#ifdef isCorrection
#ifndef isSimData

#ifndef isIter
#ifdef isStrpMulti 
  TH2D *strpCorrel2D[nmodule][nxrow][nyrow][nlayer][nside];
  TH3D *strpCorrel3D[nmodule][nxrow][nyrow][nlayer][nside];
  // TH1D *laymul[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
  TH2D *laymul[nmodule][nxrow][nyrow][nlayer][nside];
  // TH1D *layblockmul[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM][nmxhits];
  TH2D *layblockmul[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  TH3D *layblocktotmul[nmodule][nxrow][nyrow][nlayer][nside];
  // TH1D *layblocktotmul[nmodule][nxrow][nyrow][nlayer][nside][nstrip][nstrip];
  // TH1D *layblocktotmul[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  // TH1D *noiseStrpMul[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
  TH2D *noiseStrpMul2D[nmodule][nxrow][nyrow][nlayer][nside];
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
  	for(int nl=0;nl<nlayer;nl++) {
  	  for(int nj=0;nj<nside;nj++) {
	    sprintf(name,"strpCorrel2D_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    strpCorrel2D[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,nstrip,-0.5,nstrip-0.5,nstrip,-0.5,nstrip-0.5);
	    sprintf(name,"strpCorrel3D_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    int tbinno = (nmxhits%2==0?nmxhits+1:nmxhits);
	    strpCorrel3D[nm][nx][ny][nl][nj] =
	      new TH3D(name,name,nstrip,-0.5,nstrip-0.5,nstrip,-0.5,nstrip-0.5,
		       tbinno,-tbinno*0.5,tbinno*0.5);
	    sprintf(name,"laymul_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,30,-0.5,0.5,nmxhits+1,-0.5,nmxhits+0.5);
  	    for(int nsx=0;nsx<blockM;nsx++) {
  	      for(int nsy=0;nsy<blockM;nsy++) {
		sprintf(name,"layblockmul_m%i_xr%i_yr%i_%s%i_%i_%i",
			nm,nx,ny,sideMark[nj],nl,nsx,nsy);
		layblockmul[nm][nx][ny][nl][nj][nsx][nsy] =
		  new TH2D(name,name,30,-0.5,0.5,nmxhits+1,-0.5,nmxhits+0.5);
	      }}
	    sprintf(name,"layblocktotmul_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    layblocktotmul[nm][nx][ny][nl][nj] =
	      new TH3D(name,name,
		       nstrip,-0.5,nstrip-0.5,
		       nstrip,-0.5,nstrip-0.5,
		       nmxhits+1,-0.5,nmxhits+0.5);
	    // for(int nsx=0;nsx<nstrip;nsx++) {
  	    //   for(int nsy=0;nsy<nstrip;nsy++) {
  	    // 	sprintf(name,"layblocktotmul_m%i_xr%i_yr%i_%s%i_%i_%i",
  	    // 		nm,nx,ny,sideMark[nj],nl,nsx,nsy);
  	    // 	layblocktotmul[nm][nx][ny][nl][nj][nsx][nsy] =
  	    // 	  new TH1D(name,name,nmxhits,0.5,nmxhits+0.5);
  	    //   }}
	    // for(int nsx=0;nsx<blockM;nsx++) {
  	    //   for(int nsy=0;nsy<blockM;nsy++) {
  	    // 	sprintf(name,"layblocktotmul_m%i_xr%i_yr%i_%s%i_%i_%i",
  	    // 		nm,nx,ny,sideMark[nj],nl,nsx,nsy);
  	    // 	layblocktotmul[nm][nx][ny][nl][nj][nsx][nsy] =
  	    // 	  new TH1D(name,name,nmxhits,0.5,nmxhits+0.5);
  	    //   }}
  	    // for(int nsx=0;nsx<blockM;nsx++) {
  	    //   for(int nsy=0;nsy<blockM;nsy++) {
  	    // 	sprintf(name,"noiseStrpMul_m%i_xr%i_yr%i_%s%i_%i_%i",
  	    // 		nm,nx,ny,sideMark[nj],nl,nsx,nsy);
  	    // 	noiseStrpMul[nm][nx][ny][nl][nj][nsx][nsy] =
  	    // 	  new TH1D(name,name,nmxhits,0.5,nmxhits+0.5);
  	    //   }}
	    sprintf(name,"noiseStrpMul_vs_dist_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    noiseStrpMul2D[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,nmxhits,0.5,nmxhits+0.5,nstrip,0.,nstrip);
  	  }}}}}
#endif	// #ifdef isStrpMulti
#endif	// #ifndef isIter
  
  // TH1D *rawtdcStrp[nside][nlayer][nstrip];
  // for(int nj=0;nj<nside;nj++) {
  //   for(int nl=0;nl<nlayer;nl++) {
  //     for(int ns=0;ns<nstrip;ns++) {
  // 	sprintf(name,"rawtdcStrp_%s%i_s%i",sideMark[nj],nl,ns);
  // 	rawtdcStrp[nj][nl][ns] = new TH1D(name,name,240,80,140);
  //     }}}
  
  /* 64+64 strp in each tprofile */
  TProfile *rawtdcStrpProf[nmodule][nxrow][nyrow][nlayer];
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  sprintf(name,"rawtdcStrpProf_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  rawtdcStrpProf[nm][nx][ny][nl] =
	    new TProfile(name,name,2*nstrip,-0.5,2*nstrip-0.5,0.,200);
	}}}}
  
  TH1D *layPosExtErr[nmodule][nxrow][nyrow][nlayer];
  TH1D *layPosR[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
#ifdef isIter
  TH1D *blockPosR[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
#endif	// #ifdef isIter
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  sprintf(name,"layPosExtErr_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  layPosExtErr[nm][nx][ny][nl] =
	    new TH1D(name,name,300,0,0.3);
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      sprintf(name,"layPosR_m%i_xr%i_yr%i_%s%i_mul%i",
		      nm,nx,ny,sideMark[nj],nl,nh+1);
	      layPosR[nm][nx][ny][nl][nj][nh] =
		new TH1D(name,name,600,-10,10);
	    }
#ifdef isIter
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		sprintf(name,"blockPosR_m%i_xr%i_yr%i_%s%i_%i_%i",
			nm,nx,ny,sideMark[nj],nl,nsx,nsy);
		blockPosR[nm][nx][ny][nl][nj][nsx][nsy] =
		  new TH1D(name,name,300,-5.,5.);
	      }}
#endif	// #ifdef isIter
	  }}}}}
  
  TH1D *layTimeExtErr[nmodule][nxrow][nyrow][nlayer];
  TH1D *layTimeR[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
#ifdef isIter
  TH1D *blockTimeR[nmodule][nxrow][nyrow][nlayer][nside][blockM][blockM];
#endif	// #ifdef isIter
  // TH1D *strpTimeR[nmodule][nxrow][nyrow][nlayer][nside][nstrip][nstrip];
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  sprintf(name,"layTimeExtErr_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  layTimeExtErr[nm][nx][ny][nl] =
	    new TH1D(name,name,500,0,5.);
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      sprintf(name,"layTimeR_m%i_xr%i_yr%i_%s%i_mul%i",
		      nm,nx,ny,sideMark[nj],nl,nh+1);
	      layTimeR[nm][nx][ny][nl][nj][nh] =
		new TH1D(name,name,200,-20,20);
	    }
#ifdef isIter
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		sprintf(name,"blockTimeR_m%i_xr%i_yr%i_%s%i_%i_%i",
			nm,nx,ny,sideMark[nj],nl,nsx,nsy);
		blockTimeR[nm][nx][ny][nl][nj][nsx][nsy] =
		  new TH1D(name,name,200,-20.,20.);
	      }}
#endif  // #ifdef isIter
	  }}}}}

#endif	// #ifndef isSimData

  TH1D *h_ndfout = new TH1D("nhits_finder","nhits_finder",nlayer+1,-0.5,nlayer+0.5);  
  TH2D *triggereffi_evt[nmodule][nxrow][nyrow][nlayer][nside];
  /* both side signal absent */
  TH2D *inefficiency_cor[nmodule][nxrow][nyrow][nlayer];
  /* both side hit-strip absent */
  TH2D *inefficiency_cor_str[nmodule][nxrow][nyrow][nlayer];
  /* at least one side signal present, but both side hit-strip absent */
  TH2D *inefficiency_cor_strS[nmodule][nxrow][nyrow][nlayer];
  /* other side signal present, this side signal absent */
  TH2D *inefficiency_unc[nmodule][nxrow][nyrow][nlayer][nside];
  /* this side signal absent */
  TH2D *inefficiency_unc_tot[nmodule][nxrow][nyrow][nlayer][nside];
  // /* this side signal present, this side hit-strip absent */
  // TH2D *inefficiency_unc_strS_tot[nmodule][nxrow][nyrow][nlayer][nside];
  /* other side signal present, this side hit-strip absent */
  TH2D *inefficiency_unc_str[nmodule][nxrow][nyrow][nlayer][nside];
  /* this side hit-strip absent */
  TH2D *inefficiency_unc_str_tot[nmodule][nxrow][nyrow][nlayer][nside];
  TH2D *inefficiency_tot[nmodule][nxrow][nyrow][nlayer];
  // TH2D *inefficiency_tot30[nmodule][nxrow][nyrow][nlayer];
  // TH2D *inefficiency_tot3[nmodule][nxrow][nyrow][nlayer];
  TH2D *inefficiency_tot1[nmodule][nxrow][nyrow][nlayer][nside];
  
  /* Multiplicity Plots */
  TH2S *laymul_2Dstrp[nmodule][nxrow][nyrow][nlayer][nside][nstrip];
  TH2S *laymul_2DcTotFine[nmodule][nxrow][nyrow][nlayer][nside];
  TH2S *laymul_2DuTotFine[nmodule][nxrow][nyrow][nlayer][nside];
  TH2S *laymul_2DmTotFine[nmodule][nxrow][nyrow][nlayer][nside];

#ifdef is2dpos
  // TH2D *laymul_2D[nmodule][nxrow][nyrow][nlayer][nside][nmxhits];
  // TH2D *laymul_2Dm1[nmodule][nxrow][nyrow][nlayer][nside];
  TH2D *laymul_2DcTot[nmodule][nxrow][nyrow][nlayer][nside];
  TH2D *laymul_2DuTot[nmodule][nxrow][nyrow][nlayer][nside];
  TH2S *laymul_2DmTot[nmodule][nxrow][nyrow][nlayer][nside];
  TH3S *laymul_2Dpos[nmodule][nxrow][nyrow][nlayer][nside];
#endif	// #ifdef is2dpos

  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  sprintf(name,"inefficiency_tot_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  inefficiency_tot[nm][nx][ny][nl] =
	    new TH2D(name,name,
		     nstrip,0,nstrip,
		     nstrip,0,nstrip);
	  // sprintf(name,"inefficiency_tot30_m%i_xr%i_yr%i_l%i",
	  // 	  nm,nx,ny,nl);
	  // inefficiency_tot30[nm][nx][ny][nl] =
	  //   new TH2D(name,name,
	  // 	     30*nstrip,0,nstrip,
	  // 	     30*nstrip,0,nstrip);
	  // sprintf(name,"inefficiency_tot3_m%i_xr%i_yr%i_l%i",
	  // 	  nm,nx,ny,nl);
	  // inefficiency_tot3[nm][nx][ny][nl] =
	  //   new TH2D(name,name,
	  // 	     nstrip,0,nstrip,
	  // 	     nstrip,0,nstrip);
	  sprintf(name,"inefficiency_cor_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  inefficiency_cor[nm][nx][ny][nl] =
	    new TH2D(name,name,
		     nstrip,0,nstrip,
		     nstrip,0,nstrip);
	  sprintf(name,"inefficiency_cor_str_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  inefficiency_cor_str[nm][nx][ny][nl] =
	    new TH2D(name,name,
		     nstrip,0,nstrip,
		     nstrip,0,nstrip);
	  sprintf(name,"inefficiency_cor_strS_m%i_xr%i_yr%i_l%i",
		  nm,nx,ny,nl);
	  inefficiency_cor_strS[nm][nx][ny][nl] =
	    new TH2D(name,name,
		     nstrip,0,nstrip,
		     nstrip,0,nstrip);
	  
	  int sometmpval = (nmxhits+1)%2==0?(nmxhits+1)+1:(nmxhits+1);
	  for(int nj=0;nj<nside;nj++) {
	    
	    for(int nh=0;nh<nstrip;nh++) {
	      int tbinxx = (nj==0?(nmxhits+3)*multiDiv1:(nstrip+5)*multiDiv2);
	      int tbinyy = (nj==1?(nmxhits+3)*multiDiv1:(nstrip+5)*multiDiv2);
	      double tmnx = (nj==0?nh-0.5*(nmxhits+1)+0.5-2.:-4.);
	      double tmny = (nj==1?nh-0.5*(nmxhits+1)+0.5-2.:-4.);
	      double tmxx = (nj==0?nh+0.5*(nmxhits+1)+0.5:nstrip+1.);
	      double tmxy = (nj==1?nh+0.5*(nmxhits+1)+0.5:nstrip+1.);
	      sprintf(name,"laymul_2Dstrp_m%i_xr%i_yr%i_%s%i_s%i",
		      nm,nx,ny,sideMark[nj],nl,nh);
	      laymul_2Dstrp[nm][nx][ny][nl][nj][nh] =
		new TH2S(name,name,
			 tbinxx,tmnx,tmxx,
			 tbinyy,tmny,tmxy);
	    } // for(int nh=0;nh<nmxhits;nh++) {
	    
	    // for(int nh=0;nh<nmxhits;nh++) {
	    //   sprintf(name,"laymul_2D_m%i_xr%i_yr%i_%s%i_mul%i",
	    // 	      nm,nx,ny,sideMark[nj],nl,nh+1);
	    //   laymul_2D[nm][nx][ny][nl][nj][nh] =
	    // 	new TH2D(name,name,
	    // 		 nstrip*30,0,nstrip,
	    // 		 nstrip*30,0,nstrip);
	    // } // for(int nh=0;nh<nmxhits;nh++) {
	    // sprintf(name,"laymul_2Dm1_m%i_xr%i_yr%i_%s%i",
	    // 	    nm,nx,ny,sideMark[nj],nl);
	    // laymul_2Dm1[nm][nx][ny][nl][nj] =
	    //   new TH2D(name,name,
	    // 	       (nj==0?multiDiv1:multiDiv2)*nstrip,0,nstrip,
	    // 	       (nj==1?multiDiv1:multiDiv2)*nstrip,0,nstrip
	    // 	       // 30*nstrip,0,nstrip,
	    // 	       // 30*nstrip,0,nstrip
	    // 	       );
	    // int tbinxx = (nj==0?multiDiv1*nstrip+1:multiDiv2*nstrip);
	    // int tbinyy = (nj==1?multiDiv1*nstrip+1:multiDiv2*nstrip);
	    // double tmnx = (nj==0?-0.5*nstrip/tbinxx:0.);
	    // double tmny = (nj==1?-0.5*nstrip/tbinyy:0.);
	    // double tmxx = (nj==0?nstrip+0.5*nstrip/tbinxx:nstrip);
	    // double tmxy = (nj==1?nstrip+0.5*nstrip/tbinyy:nstrip);
	    sprintf(name,"laymul_2DuTotFine_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2DuTotFine[nm][nx][ny][nl][nj] =
	      new TH2S(name,name,
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+5),-4.,nstrip+1.,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+5),-4.,nstrip+1.);
	    sprintf(name,"laymul_2DcTotFine_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2DcTotFine[nm][nx][ny][nl][nj] =
	      new TH2S(name,name,
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+5),-4.,nstrip+1.,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+5),-4.,nstrip+1.);
	    sprintf(name,"laymul_2DmTotFine_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2DmTotFine[nm][nx][ny][nl][nj] =
	      new TH2S(name,name,
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+5),-4.,nstrip+1.,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+5),-4.,nstrip+1.);
#ifdef is2dpos
	    sprintf(name,"laymul_2DmTot_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2DmTot[nm][nx][ny][nl][nj] =
	      new TH2S(name,name,
		       // tbinxx,tmnx,tmxx,
		       // tbinyy,tmny,tmxy
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip
		       // 30*nstrip,0,nstrip,
		       // 30*nstrip,0,nstrip
		       );
	    sprintf(name,"laymul_2DcorTot_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2DcTot[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
		       // tbinxx,tmnx,tmxx,
		       // tbinyy,tmny,tmxy
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip
		       // 30*nstrip,0,nstrip,
		       // 30*nstrip,0,nstrip
		       );
	    sprintf(name,"laymul_2DuncTot_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2DuTot[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
		       // tbinxx,tmnx,tmxx,
		       // tbinyy,tmny,tmxy
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip
		       // 30*nstrip,0,nstrip,
		       // 30*nstrip,0,nstrip
		       );
	    sprintf(name,"laymul_2Dpos_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    laymul_2Dpos[nm][nx][ny][nl][nj] =
	      new TH3S(name,name,
		       // tbinxx,tmnx,tmxx,
		       // tbinyy,tmny,tmxy,
		       (nj==0?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip,
		       (nj==1?multiDiv1:multiDiv2)*(nstrip+1),-1.,nstrip,
		       // 30*nstrip,0,nstrip,
		       // 30*nstrip,0,nstrip,
		       sometmpval,-sometmpval*0.5,sometmpval*0.5
		       // nstrip,-0.5,nstrip-0.5
		       );
#endif	// #ifdef is2dpos
	    sprintf(name,"inefficiency_tot1_m%i_xr%i_yr%i_%s%i",
	    	    nm,nx,ny,sideMark[nj],nl);
	    inefficiency_tot1[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
	    	       nstrip,0,nstrip,
	    	       nstrip,0,nstrip);
	    sprintf(name,"triggereffi_evt_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    triggereffi_evt[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
		       nstrip,0,nstrip,
		       nstrip,0,nstrip);
	    sprintf(name,"inefficiency_unc_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    inefficiency_unc[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
		       nstrip,0,nstrip,
		       nstrip,0,nstrip);
	    sprintf(name,"inefficiency_unc_tot_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    inefficiency_unc_tot[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
		       nstrip,0,nstrip,
		       nstrip,0,nstrip);
	    sprintf(name,"inefficiency_unc_str_m%i_xr%i_yr%i_%s%i",
		    nm,nx,ny,sideMark[nj],nl);
	    inefficiency_unc_str[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
		       nstrip,0,nstrip,
		       nstrip,0,nstrip);
	    // sprintf(name,"inefficiency_unc_strS_tot_m%i_xr%i_yr%i_%s%i",
	    // 	    nm,nx,ny,sideMark[nj],nl);
	    // inefficiency_unc_strS_tot[nm][nx][ny][nl][nj] =
	    //   new TH2D(name,name,
	    // 	       nstrip,0,nstrip,
	    // 	       nstrip,0,nstrip);
	    sprintf(name,"inefficiency_unc_str_tot_m%i_xr%i_yr%i_%s%i",
	    	    nm,nx,ny,sideMark[nj],nl);
	    inefficiency_unc_str_tot[nm][nx][ny][nl][nj] =
	      new TH2D(name,name,
	    	       nstrip,0,nstrip,
	    	       nstrip,0,nstrip);
	  }}}}}
#endif	// #ifdef isCorrection
// #endif	// #ifndef isSimData
  

  
  char outfil[300] = {}; 
  char outfilx[300] = {};
  int len = strlen(argv[1]);
  strncpy(outfil, argv[1], len-5); // root file
  // sprintf(outfilx, "../rredata/temp/%s_o_%i.root", outfil, stoi(argv[4]));
  // sprintf(outfilx, "..//miniICAL_data_GMA/temp/%s_o_%i.root", outfil, stoi(argv[4]));
  sprintf(outfilx, "./temp/%s_20211205aqq_o_%05i.root", outfil, stoi(argv[4]));
  TFile *f1 = new TFile(outfilx,"RECREATE");
  
// #ifndef isSimData
#ifdef isCorrection
  // sprintf(outfilx, "../rredata/temp/%s_20210607ac_corr_%05i.root", outfil, stoi(argv[4]));
  // sprintf(outfilx, "/media/jim/INO2_mical_SSD/miniICAL_data_GMA/temp/%s_corr_%i.root", outfil, stoi(argv[4]));
#ifndef isIter
#ifdef is5of8
  sprintf(outfilx, "./temp/%s_20211205aqq_trg5of8_corr_%05i.root", outfil, stoi(argv[4]));
#else  // #ifdef is5of8
  sprintf(outfilx, "./temp/%s_20211205aqq_trg6789_corr_%05i.root", outfil, stoi(argv[4]));
#endif	// #ifdef is5of8
#else  // #ifndef isIter
  sprintf(outfilx, "./temp/%s_20211205aqq_corr_%01i_%02i_%02i_%03i_%05i.root",
	  outfil,itermodule,iterxrow,iteryrow,iterlayer,stoi(argv[4]));
#endif	// #ifndef isIter
  TFile *f2 = new TFile(outfilx,"RECREATE");
#endif	// #ifdef isCorrection
// #endif	// #ifndef isSimData
  
    
  /*
    C : a character string terminated by the 0 character
    B : an 8 bit signed integer (Char_t)
    b : an 8 bit unsigned integer (UChar_t)
    S : a 16 bit signed integer (Short_t)
    s : a 16 bit unsigned integer (UShort_t)
    I : a 32 bit signed integer (Int_t)
    i : a 32 bit unsigned integer (UInt_t)
    F : a 32 bit floating point (Float_t)
    f : a 24 bit floating point with truncated mantissa (Float16_t)
    D : a 64 bit floating point (Double_t)
    d : a 24 bit truncated floating point (Double32_t)
    L : a 64 bit signed integer (Long64_t)
    l : a 64 bit unsigned integer (ULong64_t)
    G : a long signed integer, stored as 64 bit (Long_t)
    g : a long unsigned integer, stored as 64 bit (ULong_t)
    O : [the letter o, not a zero] a boolean (Bool_t)
  */

  
  Double_t maxHough;
  Double_t meanHough;
  Double_t stdHough;
  
  TTree *outtreeh = new TTree("Outtreeh","Outtreeh");
  outtreeh->Branch("maxHough",&maxHough,"maxHough/D");
  outtreeh->Branch("meanHough",&meanHough,"meanHough/D");
  outtreeh->Branch("stdHough",&stdHough,"stdHough/D");

  Double_t evetimeFill, evesepFill, recotimeFill, recosepFill;

  Double_t momInFill, thetainFill, phiinFill;
  // Bool_t isTriggerFill;
  Int_t nhits_finder;
  Double_t momout, momerr, theout, theerr, phiout, phierr,
    dEdx, dEdxerr, eloss;
  Double_t trkLen, chi2n;
  Double_t csFill,crFill,cmomFill;
  
  Bool_t isForwardFill;
  UShort_t muDecayLayIDFill, eLayIDFill;
  Double_t decayPosInLayFill[2], ePosInLayFill[2], eTimeFill[2]; // x and y
  UShort_t eMultiFill[2];
  Double_t muExitAngleFill, eExitAngleFill, ironInteractionFill;
  
  TTree *TriggerTree = new TTree("TriggerTree","TriggerTree");
#ifndef isSimData
  TriggerTree->Branch("evetime",&evetimeFill,"evetime/D");
  TriggerTree->Branch("evesep",&evesepFill,"evesep/D");
  // TriggerTree->Branch("recotime",&recotimeFill,"recotime/D");
  // TriggerTree->Branch("recosep",&recosepFill,"recosep/D");
#else
  TriggerTree->Branch("momin",&momInFill,"momin/D");
  TriggerTree->Branch("thein",&thetainFill,"thein/D");
  TriggerTree->Branch("phiin",&phiinFill,"phiin/D");
#endif	// #ifndef isSimData
  
#ifdef isSimData
  TTree *SimTree = new TTree("SimTree","SimTree");
  SimTree->Branch("momin",&momInFill,"momin/D");
  SimTree->Branch("thein",&thetainFill,"thein/D");
  SimTree->Branch("phiin",&phiinFill,"phiin/D");
  // SimTree->Branch("isTrigger",&isTriggerFill,"isTrigger/O");
#endif	// #ifdef isSimData
  
  TTree *momTree = new TTree("MomTree","MomTree");
#ifndef isSimData
  momTree->Branch("recotime",&recotimeFill,"recotime/D");
  momTree->Branch("recosep",&recosepFill,"recosep/D");
#endif	// #ifndef isSimData
#ifdef isSimData
  momTree->Branch("momin",&momInFill,"momin/D");
  momTree->Branch("thein",&thetainFill,"thein/D");
  momTree->Branch("phiin",&phiinFill,"phiin/D");
#endif  // #ifdef isSimData
  momTree->Branch("crFill",&crFill,"crFill/D");
  momTree->Branch("csFill",&csFill,"csFill/D");
  momTree->Branch("cmomFill",&cmomFill,"cmomFill/D");
  momTree->Branch("trkLen",&trkLen,"trkLen/D");
  momTree->Branch("nhits_finder",&nhits_finder,"nhits_finder/I");
#ifndef isLifetime
  momTree->Branch("chi2",&chi2n,"chi2/D");
  momTree->Branch("momout",&momout,"momout/D");
  momTree->Branch("momerr",&momerr,"momerr/D");
  momTree->Branch("theout",&theout,"theout/D");
  momTree->Branch("theerr",&theerr,"theerr/D");
  momTree->Branch("phiout",&phiout,"phiout/D");
  momTree->Branch("phierr",&phierr,"phierr/D");
  // momTree->Branch("dEdx",&dEdx,"dEdx/D");
  // momTree->Branch("dEdxerr",&dEdxerr,"dEdxerr/D");
  momTree->Branch("eloss",&eloss,"eloss/D");
#endif	// #ifndef isLifetime
#ifdef isLifetime
  momTree->Branch("muDeacyLayID",&muDecayLayIDFill,"muDecayLayID/s");
  momTree->Branch("deacyPosInLay",decayPosInLayFill,"decayPosInLay[2]/D");
  momTree->Branch("muExitAngle",&muExitAngleFill,"muExitAngle/D");
  momTree->Branch("eExitAngle",&eExitAngleFill,"eExitAngle/D");
  // momTree->Branch("isForward",&isForwardFill,"isForward/O");
  momTree->Branch("eLayID",&eLayIDFill,"eLayID/s");
  momTree->Branch("ePosInLay",ePosInLayFill,"ePosInLay[2]/D");
  momTree->Branch("eMulti",eMultiFill,"eMulti[2]/s");
  momTree->Branch("eTime",&eTimeFill,"eTime/D");
  momTree->Branch("ironInteraction",&ironInteractionFill,"ironInteraction/D");
#endif	// #ifdef isLifetime
  
  
  
  char datafile[300] = {};
  strncpy(datafile,argv[1],300);
  Long64_t nentrymx = stoi(argv[3]);
  Long64_t nentrymn = stoi(argv[2]);
  char infile[300]  = {};

#ifdef isSimData
  // sprintf(infile, "../mIcal_Track_Fit_20201224/digiData/%s", datafile);
  // sprintf(infile, "./input_files/%s", datafile);
  // sprintf(infile, "/media/surya/Surya_5/DaqMadurai/maduraiData_mICAL/gmaData/%s", datafile);
  sprintf(infile, "/var/nfscondor/surya/sim/%s", datafile);
#else  // #ifdef isSimData
  // sprintf(infile, "./input_files/%s", datafile);
  // sprintf(infile, "/media/surya/Surya_5/DaqMadurai/maduraiData_mICAL/gmaData/%s", datafile);
  sprintf(infile, "/var/nfscondor/surya/%s", datafile);
#endif
  TFile *fileIn = new TFile(infile, "read");
  
  if(!fileIn->IsZombie()) {
    


    
// #ifndef isSimData
//     /******  Run Info   ******/
    
//     RunInfo *runinfo = (RunInfo*)fileIn->Get("run");
//     cout << " created by: " << runinfo->Creator << endl;
//     cout << " on: " << runinfo->CreatedOn << endl;
//     cout << " trigger: " << runinfo->TriggerInfo << endl;
// #endif








#ifndef isSimData
    /******  CAU Data    ******/
    
    TTree *cau_tree = (TTree*)fileIn->Get("CAUtree");
    cau_tree->SetBranchAddress("CauData", &cau);
    Long64_t Cnentry = cau_tree->GetEntries();
    cout << " cau entry " << Cnentry << endl;

    allCAUevents.clear();
    
    for(Long64_t iev=0;iev<Cnentry;iev++) {

      tmpCau.rpcId.clear();
      tmpCau.val.clear();
      
      fileIn->cd();
      cau_tree->GetEntry(iev);

      
      // CauInfo tmpCau;
      // vector<CauInfo> allCAUevents;
      
      // cout << " " << cau->Cautime
      // 	   << " " << cau->cau_tdc_ref1[0]//.size()
      // 	   << " " << cau->cau_tdc[0]//.size()
      // 	   << endl;
      // cout << " " << cau->Cautime << endl;
      
      tmpCau.Cautime = cau->Cautime;
      for(int ij=0;ij<int(cau->cau_tdc_ref1.size());ij++) {
	int tdcCnt; UShort_t rpcId; Bool_t lead;
	getCAU(cau->cau_tdc_ref1[ij],lead,rpcId,tdcCnt);
	// cout << " " << rpcId << " " << lead << " " << tdcCnt << endl;
	if(lead==0) {		// using only leading
	  if(int(tmpCau.rpcId.size())==0
	     || (int(tmpCau.rpcId.size())
		 && tmpCau.rpcId.back()!=rpcId)) {
	    tmpCau.rpcId.push_back(rpcId);
	    tmpCau.val.push_back(tdcCnt*tdc_least);
	  }
	}
      }	// for(int ij=0;ij<int(cau->cau_tdc_ref1.size());ij++) {
      allCAUevents.push_back(tmpCau);
    } // for(Long64_t iev=0;iev<Cnentry;iev++) {

// #ifdef isDebug
//     for(Long64_t iev=0;iev<int(allCAUevents.size());iev++) {
//       cout << " " << allCAUevents[iev].Cautime << endl;
//       for(int ij=0;ij<int(allCAUevents[iev].rpcId.size());ij++) {
// 	cout << "   " << allCAUevents[iev].rpcId[ij]
// 	     << " " << allCAUevents[iev].val[ij]
// 	     << endl;
//       }
//     } // for(Long64_t iev=0;iev<Cnentry;iev++) {
// #endif  // #ifdef isDebug
    
#endif	// #ifndef isSimData


    




    /******  Event Data    ******/
    
#ifdef isSimData
    TTree *event_tree = (TTree*)fileIn->Get("T2");
#ifndef isSpecial
    T2 *event = new T2(event_tree);
#else  // #ifndef isSpecial
    // Double_t        momInFill;
    // Double_t        thetainFill;
    // Double_t        phiinFill;
    UInt_t          ndigihtFill;
    UInt_t          stripidFill[200];   //[ndigihtFill]
    event_tree->SetBranchAddress("momInFill", &momInFill);
    event_tree->SetBranchAddress("thetainFill", &thetainFill);
    event_tree->SetBranchAddress("phiinFill", &phiinFill);
    event_tree->SetBranchAddress("ndigihtFill", &ndigihtFill);
    event_tree->SetBranchAddress("stripidFill", stripidFill);
#endif	// #ifndef isSpecial
#else  // #ifdef isSimData
    TTree *event_tree = (TTree*)fileIn->Get("RPCtree");
    event_tree->SetBranchAddress("EventData", &event);
#endif	// #ifdef isSimData
    
    Long64_t nentry = event_tree->GetEntries();
    nentry = min(nentry,nentrymx);
    cout << datafile <<" has "<<nentry<<" events "<<endl;
    

    for(Long64_t iev=nentrymn;iev<nentry;iev++) {
      
      // cout << " " << iev << endl;
      
#ifdef isDebug
      cout << endl << endl << " " << iev << endl;
#endif  // #ifdef isDebug
      
      if(iev%1000==0) {
	Long64_t stop_s = clock();
	cout << " iev " << iev
	     << " time " << (stop_s-start_s)/Double_t(CLOCKS_PER_SEC)
	     << endl;
      }
      
      ttExecTimeVal = 100000000;
      // for(int nj=0;nj<nside;nj++) {
      // 	for(int ij=0;ij<nlayer;ij++) {
      // 	  hit[nj][ij].clear();
      // 	}
      // }
      allrawlay.clear();
      
      fileIn->cd();
      event_tree->GetEntry(iev);
      
#ifndef isSimData
      TTimeStamp eventTime = event->EventTime;
      evesepFill = eventTime.AsDouble() - evetimeFill;
      evetimeFill = eventTime.AsDouble();
      // TriggerTree->Fill();
      // recosepFill = eventTime.AsDouble() - recotimeFill;
      // recotimeFill = eventTime.AsDouble();
#ifdef isDebug
      cout << " time " << eventTime << endl;
#endif  // #ifdef isDebug
#endif	// #ifndef isSimData

      


#ifdef isSimData
#ifndef isSpecial
      
      if((event->ndigiht)>200) continue; // hard-coded max-particle in T2.h
      if((event->ngent)>10) continue; // hard-coded max-particle in T2.h
      
      momInFill = event->momin[0];
      thetainFill = event->thein[0];
      phiinFill = event->phiin[0];
      
      SimTree->Fill();
      
      for(int ki=0;ki<event->ndigiht;ki++) {
	UInt_t strpID = event->stripid[ki];
	  
	int nInSide = (strpID>>31)==0?0:1; // x or y
	  
	strpID>>=8;
	int nInX = strpID%128;
	strpID>>=7;
	int nInCH = strpID%8;
	strpID>>=3;
	int nInMO = strpID%8;
	strpID>>=3;
	int nInLA = strpID%256;
	strpID>>=8;
	int nInDT = strpID%4;

	Int_t module = 0;
	Int_t xrow   = 0;
	Int_t yrow   = 0;
	Int_t zlay   = nInLA;
	
	int ispresent = -1;
	for(int nl=0;nl<int(allrawlay.size());nl++) {
	  if(allrawlay[nl].module==module &&
	     allrawlay[nl].xrow==xrow &&
	     allrawlay[nl].yrow==yrow &&
	     allrawlay[nl].zlay==zlay) {ispresent = nl; break;}
	} // for(int nl=0;nl<int();nl++) {
	
	tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
	tmphit.strp = nInX;
	tmphit.tdc[0].push_back(event->digitime[ki]*tdc_least);
	if(ispresent>=0) {
	  allrawlay[ispresent].hit[nInSide].push_back(tmphit);
	} else {
	  tmplayer.module = module;
	  tmplayer.xrow = xrow;
	  tmplayer.yrow = yrow;
	  tmplayer.zlay = zlay;
	  tmplayer.tdc_ref[0] = 0;
	  tmplayer.tdc_ref[1] = 0;
	  tmplayer.hit[0].clear();
	  tmplayer.hit[1].clear();
	  tmplayer.hit[nInSide].push_back(tmphit);
	  allrawlay.push_back(tmplayer);
	}
	
	     } // for(int ki=0;ki<event->ndigiht;ki++) {
      
#else	// #ifndef isSpecial
      
	     if(ndigihtFill>200) continue; // hard-coded
      
	     for(int ki=0;ki<ndigihtFill;ki++) {
	       UInt_t strpID = stripidFill[ki];
	  
	       int nInSide = (strpID>>31)==0?0:1; // x or y
	  
	       strpID>>=8;
	       int nInX = strpID%128;
	       strpID>>=7;
	       int nInCH = strpID%8;
	       strpID>>=3;
	       int nInMO = strpID%8;
	       strpID>>=3;
	       int nInLA = strpID%256;
	       strpID>>=8;
	       int nInDT = strpID%4;

	       Int_t module = 0;
	       Int_t xrow   = 0;
	       Int_t yrow   = 0;
	       Int_t zlay   = nInLA;
	
	       int ispresent = -1;
	       for(int nl=0;nl<int(allrawlay.size());nl++) {
		 if(allrawlay[nl].module==module &&
		    allrawlay[nl].xrow==xrow &&
		    allrawlay[nl].yrow==yrow &&
		    allrawlay[nl].zlay==zlay) {ispresent = nl; break;}
	       } // for(int nl=0;nl<int();nl++) {
	
	       tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
	       tmphit.strp = nInX;
	       tmphit.tdc[0].push_back(100.-nInLA*0.5);
	       if(ispresent>=0) {
		 allrawlay[ispresent].hit[nInSide].push_back(tmphit);
	       } else {
		 tmplayer.module = module;
		 tmplayer.xrow = xrow;
		 tmplayer.yrow = yrow;
		 tmplayer.zlay = zlay;
		 tmplayer.tdc_ref[0] = 0;
		 tmplayer.tdc_ref[1] = 0;
		 tmplayer.hit[0].clear();
		 tmplayer.hit[1].clear();
		 tmplayer.hit[nInSide].push_back(tmphit);
		 allrawlay.push_back(tmplayer);
	       }
	
	     } // for(int ki=0;ki<ndigiht;ki++) {
            
#endif  // #ifndef isSpecial
      
	     /* sorting hits */
	     for(int nl=0;nl<int(allrawlay.size());nl++) {
	       for(int nj=0;nj<nside;nj++) {
		 rawstrp key;
		 for(int ij=1;ij<int(allrawlay[nl].hit[nj].size());ij++) {
		   key = allrawlay[nl].hit[nj][ij];
		   int kj = ij-1;
		   while((kj>=0) &&
			 (allrawlay[nl].hit[nj][kj].strp<key.strp)) {
		     allrawlay[nl].hit[nj][kj+1] = allrawlay[nl].hit[nj][kj];
		     kj--;
		   }
		   allrawlay[nl].hit[nj][kj+1] = key;
		 }
	       } // for(int nj=0;nj<nside;nj++) {
	     }	// for(int nl=0;nl<int(allrawlay.size());nl++) {
	     
	     /* sorting layers */
	     // for(int nl=0;nl<int(allrawlay.size());nl++) {
	     //   for(int nj=0;nj<nside;nj++) {
	     {
	       rawlayer key;
	       for(int ij=1;ij<int(allrawlay.size());ij++) {
		 key = allrawlay[ij];
		 int kj = ij-1;
		 while((kj>=0) &&
		       (allrawlay[kj].zlay>key.zlay)) {
		   allrawlay[kj+1] = allrawlay[kj];
		   kj--;
		 }
		 allrawlay[kj+1] = key;
	       }
	     }
	     //   } // for(int nj=0;nj<nside;nj++) {
	     // }	// for(int nl=0;nl<int(allrawlay.size());nl++) {
		 
	     
	     // 	     /* Trigger Check */
	     // 	     unsigned int triggerinfo[nside] = {0};
	     // 	     for(int nj=0;nj<nside;nj++) {
	     // 	       for(int ij=0;ij<ntrigLayers;ij++) {
	     // 		 triggerinfo[nj]<<=1;
	     // 		 for(int nl=0;nl<int(allrawlay.size());nl++) {
	     // 		   if(allrawlay[nl].zlay==trigLayers[ij] &&
	     // 		      int(allrawlay[nl].hit[nj].size())>0) {
	     // 		     triggerinfo[nj]+=1;}
	     // 		 } // for(int nl=0;nl<int(allrawlay.size());nl++) {
	     // 	       }   // 
	     // 	     }	    // for(int nj=0;nj<nside;nj++) {
	     // #ifdef isDebug
	     // 	     cout << " \ttrigger" << bitset<10>(triggerinfo[0])
	     // 		  << " " << bitset<10>(triggerinfo[1])
	     // 		  << " " << bitset<10>(triggerinfo[0]|triggerinfo[1])
	     // 		  << endl;
	     // #endif  // #ifdef isDebug
	     // 	     // if(triggerinfo[0]!=0b1111 && triggerinfo[1]!=0b1111) {continue;}
	     // 	     if(triggerinfo[0]!=triggerinfo_ref && triggerinfo[1]!=triggerinfo_ref) {continue;}
	     // #ifdef isDebug
	     // 	     cout << " \t Not Continued"<< endl;
	     // #endif  // #ifdef isDebug
	     
	     
	     
#else  // #ifdef isSimData
	     
	     
	     
	     TClonesArray *tmpTObj = (TClonesArray *)event->GetDigiOpt1();
	     // DigiOpt1 *tmpTObj = (DigiOpt1 *)event->GetDigiOpt1();
	     // int objentry = event->GetDigiOpt1()->GetEntries();
	     int objentry = tmpTObj->GetEntries();
#ifdef isDebug
	     cout << " event entries : " << objentry << endl;
#endif  // #ifdef isDebug
      
	     for(int iobj=0;iobj<objentry;iobj++) {
	
	       DigiOpt1 *tmpDigiOpt1 = (DigiOpt1 *)tmpTObj->At(iobj);
#ifdef isDebug
	       cout << " " << bitset<16>(tmpDigiOpt1->rpcId) << endl
		    << " " << bitset<64>(tmpDigiOpt1->strips[0])
		    << " " << bitset<64>(tmpDigiOpt1->strips[1])
		    << " " <<        int(tmpDigiOpt1->tdc.size())
		    << endl;
#endif  // #ifdef isDebug
	
	       Int_t module, xrow, yrow, zlay;
	       getRPCId(tmpDigiOpt1->rpcId,module,xrow,yrow,zlay);
#ifdef isDebug
	       cout << " rpcId " << module
		    << " " << xrow
		    << " " << yrow
		    << " " << zlay
		    << endl;
#endif  // #ifdef isDebug
	       tmplayer.module = module;
	       tmplayer.xrow = xrow;
	       tmplayer.yrow = yrow;
	       tmplayer.zlay = zlay;

	       for(int nj=0;nj<nside;nj++) {
		 for(int nt=0;nt<nTDC;nt++) {
		   for(int nld=0;nld<2;nld++) {
		     tmplayer.rawTDC[nj][nt][nld].clear();
		   }}}
	       
	       tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
	       for(int nj=0;nj<nside;nj++) {
		 tmplayer.hit[nj].clear();
		 for(int kl=nstrip-1; kl>=0; kl--) {
		   if(((tmpDigiOpt1->strips[nj])>>kl)&0b1) {
		     tmphit.strp = kl;
		     // hit[nj][zlay].push_back(tmphit);
		     tmplayer.hit[nj].push_back(tmphit);
		   }
		 } // for(int kl=nstrip-1; kl>=0; kl--) {
	       } // for(int nj=0;nj<nside;nj++) {
	
	       tmplayer.tdc_ref[0] = tdc_least*tmpDigiOpt1->tdc[0];
	       tmplayer.tdc_ref[1] = tdc_least*tmpDigiOpt1->tdc[1];
	
	       tdc_ref[zlay][0] = tdc_least*tmpDigiOpt1->tdc[0];
	       tdc_ref[zlay][1] = tdc_least*tmpDigiOpt1->tdc[1];
	
#ifdef isDebug
	       cout << " tdc_ref " << tdc_ref[zlay][0]
		    << " "         << tdc_ref[zlay][1] << endl;
#endif  // #ifdef isDebug
	       for(int kl=2;kl<int(tmpDigiOpt1->tdc.size());kl++) {
		 Int_t tdcno, tdcval; Bool_t side; Bool_t lead;
		 getTDC(tmpDigiOpt1->tdc[kl],lead,side,tdcno,tdcval);
#ifdef isDebug
		 cout << " tdc " << lead
		      << " " << side
		      << " " << tdcno
		      << " " << tdcval*tdc_least
		      << endl;
#endif  // #ifdef isDebug
		 
		 tmplayer.rawTDC[side][tdcno][lead].push_back(tdcval*tdc_least);
		 
		 for(int ki=0;ki<int(tmplayer.hit[side].size());ki++) {
		   if(tmplayer.hit[side][ki].strp%nTDC==tdcno) {

		     double tdccorr = 0;
#ifdef isTDCstrpCorr
		     tdccorr = strpTDCcorr[tmplayer.module][tmplayer.xrow][tmplayer.yrow][tmplayer.zlay][side][tmplayer.hit[side][ki].strp];
		     // cout << " tdccorr " << tdccorr << endl;
#endif	// #ifdef isTDCstrpCorr
		     
		     tmplayer.hit[side][ki].tdc[lead].push_back(tdcval*tdc_least - tdccorr);
		     if(lead==0) {
		       Double_t tdctime = tdcval*tdc_least - tmplayer.tdc_ref[0] - tdccorr;
		       if(ttExecTimeVal>tdctime) {
			 ttExecTimeVal = tdctime;}
		     }
		   }
		 } // for(int ki=0;ki<int(tmplayer.hit[side][zlay].size());ki++) {
	       } // for(int kl=2;kl<int(tmpDigiOpt1->tdc.size());kl++) {
	       allrawlay.push_back(tmplayer);
	     }	// for(int iobj=0;iobj<objentry;iobj++) {

#endif	// #ifdef isSimData
      
	     /* Making all tdc_l +ve w.r.t. tdc_ref_l */
	     for(int nj=0;nj<nside;nj++) {
	       for(int nl=0;nl<int(allrawlay.size());nl++) {
		 for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
		   for(int ntd=0;ntd<2;ntd++) {
		     for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[ntd].size());nht++) {
#ifndef isSimData
		       allrawlay[nl].hit[nj][nh].tdc[ntd][nht] -= (ttExecTimeVal-100.);
#else  // #ifndef isSimData
		       allrawlay[nl].hit[nj][nh].tdc[ntd][nht] -= (-100.);
#endif	// #ifndef isSimData
		     }
		   } // for(int ntd=0;ntd<2;ntd++) {
		 }
	       } // for(int nl=0;nl<nlayer;nl++) {
	     }	// for(int nj=0;nj<nside;nj++) {
      
#ifdef isDebug
	     cout << " min time " << ttExecTimeVal << endl;
	     for(int nj=0;nj<nside;nj++) {
	       cout << " " << sideMark[nj] << endl;
	       for(int nl=0;nl<int(allrawlay.size());nl++) {
		 cout << "  l" << allrawlay[nl].zlay << " ";
		 for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
		   cout << " " << allrawlay[nl].hit[nj][nh].strp;
		   cout << " (";
		   for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[0].size());nht++) {
		     cout << " " << allrawlay[nl].hit[nj][nh].tdc[0][nht] - allrawlay[nl].tdc_ref[0];
		   }
		   cout << " ) ";
		 }
		 cout << endl;
	       }
	     }
#endif  // #ifdef isDebug
      
      
      
	     /*****  Grouping Hits as per time info *****/
	     /*****  Using only leading time stamp  *****/
      
	     TimeGroup.clear();
      
	     strptimeraw->Scale(0);
	     for(int nj=0;nj<nside;nj++) {
	       for(int nl=0;nl<int(allrawlay.size());nl++) {
		 for(int nh=0;nh<int(allrawlay[nl].hit[nj].size());nh++) {
		   for(int nht=0;nht<int(allrawlay[nl].hit[nj][nh].tdc[0].size());nht++) {
		     double tdctime = allrawlay[nl].hit[nj][nh].tdc[0][nht]-allrawlay[nl].tdc_ref[0];
		     strptimeraw->Fill(tdctime);
		   }
		 }
	       } // for(int nl=0;nl<int(allrawlay.size());nl++) {
	     }	// for(int nj=0;nj<nside;nj++) {
      
	     vector<rawlayer>   tallrawlay1 = allrawlay;
	     strptimeraw->SetBinContent(strptimeraw->GetNbinsX(),0);
	     int bnx = strptimeraw->GetNbinsX()-1;
	     // while(bnx<strptimeraw->GetNbinsX()) {
	     while(bnx>0) {
	       if((strptimeraw->GetBinContent(bnx+1))>0) {
		 double bEnd = (strptimeraw->GetBinLowEdge(bnx+1) +
				strptimeraw->GetBinWidth(bnx+1));
		 while((strptimeraw->GetBinContent(bnx+1))>0) {
		   bnx--;
		 }
		 double bStart = strptimeraw->GetBinLowEdge(bnx+2);
		 // cout << " " << bStart << " " << bEnd << endl;
		 tmpevt.alllayer.clear();
		 for(int tnl=0;tnl<int(tallrawlay1.size());tnl++) {
		   tmplayer = tallrawlay1[tnl];
		   for(int nj=0;nj<nside;nj++) {
		     tmplayer.hit[nj].clear(); // tmplayer.hit[nj].clear();
		     for(int nh=0;nh<int(tallrawlay1[tnl].hit[nj].size());nh++
			 ) {
		       tmphit.strp = tallrawlay1[tnl].hit[nj][nh].strp;
		       tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
		       for(int nht=0;nht<int(tallrawlay1[tnl].hit[nj][nh].tdc[0].size());nht++) {
			 double tdctime = tallrawlay1[tnl].hit[nj][nh].tdc[0][nht]-tallrawlay1[tnl].tdc_ref[0];
			 if(tdctime>=bStart && tdctime<=bEnd) {
			   tmphit.tdc[0].push_back(tallrawlay1[tnl].hit[nj][nh].tdc[0][nht]);
			   if(nht<int(tallrawlay1[tnl].hit[nj][nh].tdc[1].size())) {
			     tmphit.tdc[1].push_back(tallrawlay1[tnl].hit[nj][nh].tdc[1][nht]);
			   }
			 }
		       } // for(int nht=0;nht<int(hit[nj][tnl][nh].tdc[0].size());nht++)
		       if(int(tmphit.tdc[0].size())>0) {
			 tmplayer.hit[nj].push_back(tmphit);
#ifdef removeEarlierHit
			 /* Removing revious group's hits: Check 20210701 for details */
			 tallrawlay1[tnl].hit[nj].erase(tallrawlay1[tnl].hit[nj].begin()+nh);
			 nh--;
#endif		// #ifdef removeEarlierHit
		       } // if(int(tmphit.tdc[0].size())>0) {
		     }	// for(int nh=0;nh<int(hit[nj][tnl].size());nh++) {
		   }   // for(int nj=0;nj<nside;nj++) {
		   if(int(tmplayer.hit[0].size())>0
		      || int(tmplayer.hit[1].size())>0) {
		     tmpevt.alllayer.push_back(tmplayer);
		   }
		 } // for(int nl=0;nl<int(tallrawlay1.size());nl++) {
		 if(int(tmpevt.alllayer.size())>0) {
		   TimeGroup.push_back(tmpevt);}
	       } else {bnx--;}
	     }	// while(1) {
	     for(int bnx1=0;bnx1<int(TimeGroup.size());bnx1++) {
	       TimeGroup.insert(TimeGroup.begin()+bnx1,TimeGroup.back());
	       TimeGroup.pop_back();
	     }
	     
	     if(int(TimeGroup.size())<=0) {continue;}
      
#ifdef isDebug
	     cout << " TimeGroup " << " " << TimeGroup.size() << endl;
	     for(int isj=0;isj<int(TimeGroup.size());isj++) {
	       cout << " TimeGroup " << isj << endl;
	       for(int nj=0;nj<nside;nj++) {
		 cout << "   " << sideMark[nj] << endl;
		 for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
		   tmplayer = TimeGroup[isj].alllayer[nl];
		   cout << "    l" << tmplayer.zlay << " ";
		   for(int nh=0;nh<int(tmplayer.hit[nj].size());nh++) {
		     cout << " " << tmplayer.hit[nj][nh].strp;
		     cout << " (";
		     for(int nht=0;nht<int(tmplayer.hit[nj][nh].tdc[0].size());nht++) {
		       cout << " " << tmplayer.hit[nj][nh].tdc[0][nht] - tmplayer.tdc_ref[0];
		     }
		     cout << " ) ";
		   }
		   cout << endl;
		 }
	       }
	     }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {
#endif  // #ifdef isDebug
	     
	     // cout << " iev " << iev << endl;
	     
// #ifndef isCorrection
#ifdef isTrgCheck
	       
	     /* Trigger Check */
	     unsigned int triggerinfo[nside] = {0};
	     for(int nj=0;nj<nside;nj++) {
	       for(int ij=0;ij<ntrigLayers1;ij++) {
		 triggerinfo[nj]<<=1;
		 for(int nl=0;nl<int(TimeGroup[0].alllayer.size());nl++) {
		   if(TimeGroup[0].alllayer[nl].zlay==trigLayers1[ij] &&
		      int(TimeGroup[0].alllayer[nl].hit[nj].size())>0) {
		     triggerinfo[nj]+=1;}
		 }
	       } // 
	     }	 // for(int nj=0;nj<nside;nj++) {
#ifdef isDebug
	     cout << " \ttrigger " << bitset<10>(triggerinfo[0])
		  << " " << bitset<10>(triggerinfo[1])
		  << " " << bitset<10>(triggerinfo[0]|triggerinfo[1])
		  << endl;
#endif  // #ifdef isDebug
	     // if(triggerinfo[0]!=0b1111 && triggerinfo[1]!=0b1111) {continue;}
	     if(triggerinfo[0]!=triggerinfo_ref1 && triggerinfo[1]!=triggerinfo_ref1) {continue;}
#ifdef isDebug
	     cout << " \t Not Continued"<< endl;
#endif  // #ifdef isDebug
	     
#endif	// #ifdef isTrgCheck
// #endif	// #ifndef isCorrection
	     
	     TriggerTree->Fill();
#ifdef isTriggerRateOnly
	     continue;
#endif	// #ifdef isTriggerRateOnly
	     

	     
	     /*****  Forming Cluster of Strips in each Layer *****/
      
	     for(int isj=0;isj<int(TimeGroup.size());isj++) {
	       
	       for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
		 tmplayer = TimeGroup[isj].alllayer[nl];
	  
		 for(int nj=0;nj<nside;nj++) {
		   // TimeGroup[isj].alllayer[nl].cluster[nj].clear();
		   tmplayer.cluster[nj].clear();
		   int totHit = tmplayer.hit[nj].size();
		   // if(totHit>MaxTotHit) {continue;}
	    
		   tmpCluster.clear();
		   int tempX = 0;
		   for(int jk=0;jk<totHit;jk++) {
		     int lowX = tmplayer.hit[nj].back().strp;
		     tmpCluster.push_back(tmplayer.hit[nj].back());
		     tmplayer.hit[nj].insert(tmplayer.hit[nj].begin(),
					     tmplayer.hit[nj].back());
		     tmplayer.hit[nj].pop_back();
		     tempX++;
		     if((lowX+1!=tmplayer.hit[nj].back().strp)
			||(jk+1==totHit)) {
		       // if(tempX<=MaxMulti) {
		       //   tmplayer.cluster[nj].push_back(tmpCluster);
		       // }
		       tmplayer.cluster[nj].push_back(tmpCluster);
		       tempX=0;
		       tmpCluster.clear();
		     }
		   } // for(int jk=0;jk<totHit;jk++) {
		   // if(tmplayer.cluster[nj].size()>MaxTotPos) {tmplayer.cluster[nj].clear();}
	    
		   /*** checking gap between clusters ***/
		   for(int jk=0;jk<int(tmplayer.cluster[nj].size())-1;jk++) {
		     for(int ki=jk+1;ki<int(tmplayer.cluster[nj].size());) {
		       // cout << " some1 " << jk << " " << ki << endl;
		       int startPos = tmplayer.cluster[nj][jk].back().strp;
		       int endPos   = tmplayer.cluster[nj][ki].front().strp;
		       if(endPos-startPos<=MinClusterSep+1) {
			 for(int kii=startPos+1;kii<endPos;kii++) {
			   // filling gap strp without tdc entry
			   tmphit.tdc[0].clear(); tmphit.tdc[1].clear();
			   tmphit.strp = kii;
			   tmplayer.cluster[nj][jk].push_back(tmphit);
			 }
			 for(int kii=0;kii<int(tmplayer.cluster[nj][ki].size());kii++) {
			   tmplayer.cluster[nj][jk].push_back(tmplayer.cluster[nj][ki][kii]);
			 }
			 tmplayer.cluster[nj].erase(tmplayer.cluster[nj].begin()+ki);
		       } else {ki++;}
		       // cout << " something " << iev
		       //      << " " << tmplayer.cluster[nj][jk].back().strp
		       //      << " " << tmplayer.cluster[nj][ki].front().strp
		       //      << endl;
		     }	// for(int ki=jk+1;ki<int(tmplayer.cluster[nj].size());) {
		     // if(int(tmplayer.cluster[nj][jk].size())>MaxMulti) {
		     // 	tmplayer.cluster[nj].erase(tmplayer.cluster[nj].begin()+jk);
		     // } else {jk++;}
		   }	// for(int jk=0;jk<int(tmplayer.cluster[nj].size())-1;jk++) {
		 } // for(int nj=0;nj<nside;nj++) {
		 TimeGroup[isj].alllayer[nl] = tmplayer;

#ifndef isSimData	  
#ifdef isCorrection
		 if(int(tmplayer.cluster[0].size())==1 &&
		    int(tmplayer.cluster[1].size())==1 &&
		    int(tmplayer.cluster[0][0].size())==1 &&
		    int(tmplayer.cluster[1][0].size())==1) {
		   double ttttime = (tmplayer.cluster[0][0][0].tdc[0][0]-tmplayer.tdc_ref[0]) - tmplayer.cluster[1][0][0].strp*strpwidth*spdl;
		   // cout << " ttttime x " << ttttime << endl;
		   // rawtdcStrp[0][tmplayer.zlay][tmplayer.cluster[0][0][0].strp]->Fill(ttttime);
		   rawtdcStrpProf[tmplayer.module][tmplayer.xrow][tmplayer.yrow][tmplayer.zlay]->Fill(tmplayer.cluster[0][0][0].strp,ttttime);

		   ttttime = (tmplayer.cluster[1][0][0].tdc[0][0]-tmplayer.tdc_ref[0]) - tmplayer.cluster[0][0][0].strp*strpwidth*spdl;
		   // cout << " ttttime y " << ttttime << endl;
		   // rawtdcStrp[1][tmplayer.zlay][tmplayer.cluster[1][0][0].strp]->Fill(ttttime);
		   rawtdcStrpProf[tmplayer.module][tmplayer.xrow][tmplayer.yrow][tmplayer.zlay]->Fill(tmplayer.cluster[1][0][0].strp+nstrip,ttttime);
		 }
#endif	// #ifdef isCorrection
#endif	// #ifndef isSimData
	  
	       } // for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	
#ifdef isDebug
	       cout << " Cluster Pos " << isj << endl;
	       for(int nj=0;nj<nside;nj++) {
		 cout << " " << sideMark[nj] << endl;
		 for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
		   tmplayer = TimeGroup[isj].alllayer[nl];
		   cout << "  l" << tmplayer.zlay;
		   for(int jk=0;jk<int(tmplayer.cluster[nj].size());jk++) {
		     cout << " (";
		     for(int ki=0;ki<int(tmplayer.cluster[nj][jk].size());ki++) {
		       cout << " " << tmplayer.cluster[nj][jk][ki].strp;
		     }
		     cout << " )";
		   }
		   cout << endl;
		 }
	       } // for(int nj=0;nj<nside;nj++) {
#endif	// #ifdef isDebug
	
	       // #ifdef isSpclDebug
	       // 	// cout << " TimeGroup Pos " << isj << endl;
	       // 	for(int nj=0;nj<nside;nj++) {
	       // 	  // cout << " " << sideMark[nj] << endl;
	       // 	  for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	       // 	    tmplayer = TimeGroup[isj].alllayer[nl];
	       // 	    // cout << "  l" << tmplayer.zlay;
	       // 	    for(int jk=0;jk<int(tmplayer.cluster[nj].size());jk++) {
	       // 	      if(int(tmplayer.cluster[nj][jk].size())<5) {continue;}
	       // 	      cout << " (";
	       // 	      for(int ki=0;ki<int(tmplayer.cluster[nj][jk].size());ki++) {
	       // 		cout << " iev " << iev
	       // 		     << " " << sideMark[nj]
	       // 		     << " l" << tmplayer.zlay
	       // 		     << " " << tmplayer.cluster[nj][jk][ki].strp
	       // 		     << endl;
	       // 	      }
	       // 	      cout << " )" << endl;
	       // 	    }
	       // 	    // cout << endl;
	       // 	  }
	       // 	} // for(int nj=0;nj<nside;nj++) {
	       // #endif	// #ifdef isDebug
	
	     }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {


      
      

	     /*****  Forming the vector of PosInSpace *****/
	     
	     for(int isj=0;isj<int(TimeGroup.size());isj++) {
	       TimeGroup[isj].PosForTrack.clear();
	       for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
		 // tmplayer = TimeGroup[isj].alllayer[nl];
		 tmplayer = TimeGroup[isj].alllayer[nl];
		 tmpPos.rawhitInfo = TimeGroup[isj].alllayer[nl];
		 
		 if(int(tmplayer.hit[0].size())>MaxTotHit ||
		    int(tmplayer.hit[1].size())>MaxTotHit ||
		    int(tmplayer.cluster[0].size())>MaxTotPos ||
		    int(tmplayer.cluster[1].size())>MaxTotPos) {continue;}
	  
		 for(int xn=0;xn<int(tmplayer.cluster[0].size());xn++) {
		   for(int yn=0;yn<int(tmplayer.cluster[1].size());yn++) {
		     // cout << " in " << isj << " " << tmplayer.zlay
		     // 	   << " " << xn << " " << yn << endl;
		     if(int(tmplayer.cluster[0][xn].size())>MaxMulti ||
			int(tmplayer.cluster[1][yn].size())>MaxMulti) {continue;}
		     tmpPos.rawhitInfo.cluster[0].clear();
		     tmpPos.rawhitInfo.cluster[0].push_back(tmplayer.cluster[0][xn]);
		     tmpPos.rawhitInfo.cluster[1].clear();
		     tmpPos.rawhitInfo.cluster[1].push_back(tmplayer.cluster[1][yn]);	      
		     // cout << " out " << isj << " " << tmplayer.zlay
		     // 	   << " " << xn << " " << yn << endl;
		     // tmpPos.rawhitInfo.hit[0] = tmplayer.cluster[0][xn];
		     // tmpPos.rawhitInfo.hit[1] = tmplayer.cluster[1][yn];
		     calHitPos(tmpPos);
		     TimeGroup[isj].PosForTrack.push_back(tmpPos);
		     // cout << " size " << TimeGroup[isj].PosForTrack.size() << endl;
		   } // for(int yn=0;yn<int(tmplayer.cluster[1].size());yn++) {
		 } // for(int xn=0;xn<int(tmplayer.cluster[0].size());xn++) {
	       } // for(int nl=0;nl<int(TimeGroup[isj].alllayer.size());nl++) {
	       
#ifdef isDebug
	       cout << " TimeGroup Cluster List " << isj << endl;
	       for(int nj=0;nj<nside;nj++) {
		 cout << " " << sideMark[nj] << endl;
		 for(int nl=0;nl<int(TimeGroup[isj].PosForTrack.size());nl++) {
		   tmpPos = TimeGroup[isj].PosForTrack[nl];
		   cout << "  l" << tmpPos.rawhitInfo.zlay;
		   for(int ki=0;ki<int(tmpPos.rawhitInfo.cluster[nj].back().size());ki++) {
		     cout << " " << tmpPos.rawhitInfo.cluster[nj][0][ki].strp;
		     cout << " (";
		     for(int nht=0;nht<int(tmpPos.rawhitInfo.cluster[nj][0][ki].tdc[0].size());nht++) {
		       cout << " " << tmpPos.rawhitInfo.cluster[nj][0][ki].tdc[0][nht] - tmpPos.rawhitInfo.tdc_ref[0];
		     }
		     cout << " )";
		   }
		   cout << endl;
		 }
	       } // for(int nj=0;nj<nside;nj++) {
#endif	// #ifdef isDebug
	
	     }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {



#ifdef isAnalysis
	     
      
	     /*****  Analysis and Stuff *****/

	     vector<HoughBinInfo> H_Bins_Arr, H_Bins_Cal;
	     HoughBinInfo tttbinC;
      
	     for(int isj=0;isj<int(TimeGroup.size());isj++) {
	       // cout << " TimeGroup " << isj << endl;
	       if(int(TimeGroup[isj].PosForTrack.size())<MinLayHit) {continue;}
	       TimeGroup[isj].HoughCluster.clear();
	
	       // cout << " TimeGroup " << isj << endl;
	
	       H_Bins_Arr.clear();
	
	       int totCnt = TimeGroup[isj].PosForTrack.size();
	       for(int ki=0;ki<totCnt-1;ki++) {
		 for(int jk=ki+1;jk<totCnt;jk++) {
	    
		   int zlay0 = TimeGroup[isj].PosForTrack[ki].rawhitInfo.zlay;
		   int zlay1 = TimeGroup[isj].PosForTrack[jk].rawhitInfo.zlay;
	    
		   // if(fabs(zlay0-zlay1)!=1) {continue;}
		   // if(fabs(zlay0-zlay1)==0) {continue;}
		   if(fabs(zlay0-zlay1)<1) {continue;}
		   if(fabs(zlay0-zlay1)>maxHoughDist) {continue;}
	    
		   TVector3 pos0 = TimeGroup[isj].PosForTrack[ki].RawHitPos;
		   TVector3 pos1 = TimeGroup[isj].PosForTrack[jk].RawHitPos;
		   TVector3 tmpUnit = (pos0 - pos1).Unit();
		   TVector3 OOrigin(0,0,0);
	    
		   TVector3 minDistV = getMinDist(OOrigin,
						  tmpUnit,
						  pos0,
						  tmpUnit);
#ifdef isDebug
		   cout << " " << zlay0 << " " << zlay1
			<< " pos0 " << pos0.X()
			<< " " << pos0.Y()
			<< " " << pos0.Z()
			<< " pos1 " << pos1.X()
			<< " " << pos1.Y()
			<< " " << pos1.Z()
		     // << " unit " << tmpUnit.X()
		     // << " " << tmpUnit.Y()
		     // << " " << tmpUnit.Z()
			<< " dir " << (pos0 - pos1).Theta()*180./TMath::Pi()
			<< " " << (pos0 - pos1).Phi()*180./TMath::Pi()
		     // << " " << (pos0 - pos1).Z()
			<< endl;
		   cout << " " << minDistV.Mag()
			<< " " << minDistV.Phi()
			<< " " << tmpUnit.Theta()
			<< " " << tmpUnit.Phi()
			<< endl;
#endif	// #ifdef isDebug
	    
		   tttbinC.posId.clear();
		   tttbinC.posId.push_back(ki);
		   tttbinC.posId.push_back(jk);
	      
		   tttbinC.binC = 1;
	    
		   tttbinC.xbin[0] = (int(minDistV.Mag()/H_R_BinW) +
				      0.5)*H_R_BinW;

		   tttbinC.xbin[1] = minDistV.Phi();
		   if(tttbinC.xbin[1]<0) {tttbinC.xbin[1] += 2.*TMath::Pi();}
		   tttbinC.xbin[1] = (int(tttbinC.xbin[1]/H_Phi_BinW) +
				      0.5)*H_Phi_BinW;
	    
		   tttbinC.xbin[2] = (int(tmpUnit.Theta()/H_Theta_BinW) +
				      0.5)*H_Theta_BinW;

		   tttbinC.xbin[3] = tmpUnit.Phi();
		   if(tttbinC.xbin[3]<0) {tttbinC.xbin[3] += 2.*TMath::Pi();}
		   tttbinC.xbin[3] = (int(tttbinC.xbin[3]/H_Phi_BinW) +
				      0.5)*H_Phi_BinW;
	    
		   tttbinC.xbinI[0] = tttbinC.xbin[0]/H_R_BinW;
		   tttbinC.xbinI[1] = tttbinC.xbin[1]/H_Phi_BinW;
		   tttbinC.xbinI[2] = tttbinC.xbin[2]/H_Theta_BinW;
		   tttbinC.xbinI[3] = tttbinC.xbin[3]/H_Phi_BinW;
	    
#ifdef isDebug
		   for(int px=0;px<ndf3D;px++) {
		     cout << " \t" << px
			  << " " << tttbinC.xbinI[px]
			  << " " << tttbinC.xbin[px]//*180./TMath::Pi()
			  << endl;
		   }
#endif	// #ifdef isDebug
	    
		   H_Bins_Arr.push_back(tttbinC);
	    	    
		 } // for(int jk=ki+1;jk<totCnt;jk++) {
	       } // for(int ki=0;ki<totCnt-1;ki++) {
	
#ifdef isDebug
	       cout << " size0 " << int(H_Bins_Arr.size()) << endl;
#endif	// #ifdef isDebug
	
	/*** Adding Entries in same bins ***/
	       for(int jk=0;jk<int(H_Bins_Arr.size())-1;jk++) {
		 // cout << jk << endl;
		 for(int ki=jk+1;ki<int(H_Bins_Arr.size());) {
		   bool testb = 1;
		   for(int px=0;px<ndf3D;px++) {
		     if(H_Bins_Arr[jk].xbinI[px]!=H_Bins_Arr[ki].xbinI[px]) {
		       testb = 0;}
		   }
		   if(testb) {
		     H_Bins_Arr[jk].binC+=H_Bins_Arr[ki].binC;
		     for(int px=0;px<int(H_Bins_Arr[ki].posId.size());px++) {
		       H_Bins_Arr[jk].posId.push_back(H_Bins_Arr[ki].posId[px]);
		     }
		     H_Bins_Arr.erase(H_Bins_Arr.begin()+ki);
		   } else {ki++;}
		 } // for(int ki=jk+1;ki<H_Bins_Arr.size();) {
	       }   // for(int jk=0;jk<H_Bins_Arr.size()-1;jk++) {

	       int maxHoughE = 0;
	       hbinIns->Scale(0);
	       for(int jk=0;jk<int(H_Bins_Arr.size());jk++) {
		 if(H_Bins_Arr[jk].binC>1) {
		   hbinIns->Fill(H_Bins_Arr[jk].binC);}
		 if(maxHoughE<H_Bins_Arr[jk].binC) {
		   maxHoughE = H_Bins_Arr[jk].binC;
		 }
	       } // for(int jk=0;jk<int(H_Bins_Arr.size())-1;jk++) {
	
	       maxHough = maxHoughE;
	       meanHough = hbinIns->GetMean();
	       stdHough = hbinIns->GetStdDev();
	       // outtreeh->Fill();
	
	       TimeGroup[isj].maxHoughPt = maxHough;
	       TimeGroup[isj].meanHough = meanHough;
	       TimeGroup[isj].stdHough  = stdHough;
	
#ifdef isDebug
	       cout << " size1 " << int(H_Bins_Arr.size()) << endl;
	       for(int jk=0;jk<int(H_Bins_Arr.size());jk++) {
		 cout << " " << jk << " " << H_Bins_Arr[jk].binC << " (";
		 for(int px=0;px<int(H_Bins_Arr[jk].posId.size());px++) {
		   cout << " " << H_Bins_Arr[jk].posId[px];
		 }
		 cout << " )";
		 cout << " (";
		 for(int px=0;px<ndf3D;px++) {
		   cout << " " << H_Bins_Arr[jk].xbinI[px];
		 } 
		 cout << " )" << endl;
	       } // for(int jk=0;jk<int(H_Bins_Arr.size())-1;jk++) {
#endif	// #ifdef isDebug
	
	/*** removing entries with less points ***/
	       H_Bins_Cal.clear();
	       for(int jk=0;jk<int(H_Bins_Arr.size());jk++) {
		 for(int px=0;px<ndf3D;px++) {
		   tttbinC.xbin[px] = H_Bins_Arr[jk].xbin[px];
		   tttbinC.xbinI[px] = H_Bins_Arr[jk].xbinI[px];
		 }
		 tttbinC.binC = H_Bins_Arr[jk].binC;
		 tttbinC.posId = H_Bins_Arr[jk].posId;
		 if(tttbinC.binC>=MinHoughEntry) {H_Bins_Cal.push_back(tttbinC);}
	       }   // for(int jk=0;jk<H_Bins_Arr.size();jk++) {
	
#ifdef isDebug
	       cout << " size2 " << int(H_Bins_Cal.size()) << endl;
	       for(int jk=0;jk<int(H_Bins_Cal.size());jk++) {
		 cout << " " << H_Bins_Cal[jk].binC << " (";
		 for(int px=0;px<int(H_Bins_Cal[jk].posId.size());px++) {
		   cout << " " << H_Bins_Cal[jk].posId[px];
		 }
		 cout << " )";
		 cout << " (";
		 for(int px=0;px<ndf3D;px++) {
		   cout << " " << H_Bins_Cal[jk].xbinI[px];
		 } 
		 cout << " )" << endl;
	       } // for(int jk=0;jk<int(H_Bins_Cal.size())-1;jk++) {
#endif	// #ifdef isDebug
	
	       vector<HoughBinInfo> H_Bins_Cal1 = H_Bins_Cal;
	




	
	
	       /*** Clustering the hough points ***/
	
	       vector<HoughBinInfo> outBins;
	       outBins.clear();
	       vector<HoughBinInfo> outBinsTmp;

	       while(int(H_Bins_Cal1.size())) {
	  
		 outBinsTmp.clear();
	  
#ifdef isDebug
		 cout << " in size3 " << int(outBinsTmp.size())
		      << " " << int(H_Bins_Cal1.size()) << endl;
#endif	// #ifdef isDebug
	  
		 // removing elements with zero pos entry
		 for(int ij=0;ij<int(H_Bins_Cal1.size());ij++) {
		   if(int(H_Bins_Cal1[ij].posId.size())==0) {
		     H_Bins_Cal1.erase(H_Bins_Cal1.begin()+ij);
		   } else {ij++;}
		 } // for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
	  
		 if(int(H_Bins_Cal1.size())==0) {continue;}
	  
		 outBinsTmp.push_back(H_Bins_Cal1.front());
		 H_Bins_Cal1.erase(H_Bins_Cal1.begin());
	  
		 /* Checking nearby Hough Points */
		 getHoughCluster(H_Bins_Cal1,outBinsTmp);
		 
#ifdef isDebug
		 cout << " out size3 " << int(outBinsTmp.size())
		      << " " << int(H_Bins_Cal1.size()) << endl;
#endif	// #ifdef isDebug
	  
		 outBins.push_back(outBinsTmp.front());
		 outBinsTmp.erase(outBinsTmp.begin());
		 for(int ki=0;ki<int(outBinsTmp.size());ki++) {
		   for(int px=0;px<int(outBinsTmp[ki].posId.size());px++) {
		     outBins.back().posId.push_back(outBinsTmp[ki].posId[px]);
		   }
		 }
#ifdef isDebug
		 cout << " (";
		 for(int px=0;px<int(outBins.back().posId.size());px++) {
		   cout << " " << outBins.back().posId[px];
		 }
		 cout << " )";
		 cout << " (";
		 for(int px=0;px<ndf3D;px++) {
		   cout << " " << outBins.back().xbinI[px];
		 } 
		 cout << " )" << endl;
#endif	// #ifdef isDebug
	  
	       } // while(int(H_Bins_Cal1.size())) {
	
	       // if(int(outBins.size())>1) {
	       //   cout << " " << iev << " " << int(outBins.size()) << endl;
	       // }
	
	       for(int jk=0;jk<int(outBins.size());jk++) {
	  
		 for(int ij=0;ij<int(outBins[jk].posId.size())-1;ij++) {
		   for(int ki=ij+1;ki<int(outBins[jk].posId.size());) {
		     if((outBins[jk].posId[ij]==outBins[jk].posId[ki])) {
		       outBins[jk].posId.erase(outBins[jk].posId.begin()+ki);
		     } else {ki++;}
		   }
		 } // for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
	  
		 int key;
		 for(int ij=1;ij<int(outBins[jk].posId.size());ij++) {
		   key = outBins[jk].posId[ij];
		   int kj = ij-1;
		   while((kj>=0) && (outBins[jk].posId[kj]>key)) {
		     outBins[jk].posId[kj+1] = outBins[jk].posId[kj];
		     kj--;
		   }
		   outBins[jk].posId[kj+1] = key;
		 }
#ifdef isDebug
		 for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		   cout << " in " << ij
			<< " " << outBins[jk].posId[ij]
			<< endl;
		 }
		 for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		   cout << " out " << ij
			<< " " << outBins[jk].posId[ij]
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].rawhitInfo.zlay
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.X()
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.Y()
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.Z()
			<< endl;
		 }
#endif	// #ifdef isDebug

	       }   // for(int jk=0;jk<outBins.size();jk++) {

	
	       /* Hough Groups with highest to lowest */
	       if(1) {
		 HoughBinInfo key;
		 for(int ij=1;ij<int(outBins.size());ij++) {
		   key = outBins[ij];
		   int kj = ij-1;
		   while((kj>=0)&&
			 (outBins[kj].posId.size()<key.posId.size())) {
		     outBins[kj+1] = outBins[kj];
		     kj--;
		   }
		   outBins[kj+1] = key;
		 }
	       }
#ifdef isDebug
	       for(int jk=0;jk<outBins.size();jk++) {
		 for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		   cout << " out1 " << ij
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].rawhitInfo.zlay
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.X()
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.Y()
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.Z()
			<< endl;
		 }
	       }   // for(int jk=0;jk<outBins.size();jk++) {
#endif	// #ifdef isDebug
	
	
	
	/*** checking overlap ***/
	
	       for(int jk=0;jk<int(outBins.size())-1;jk++) {
		 for(int ki=jk+1;ki<int(outBins.size());) {
		   int overl = 0;
		   for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		     for(int iji=0;iji<int(outBins[ki].posId.size());iji++) {
		       if(outBins[jk].posId[ij]==outBins[ki].posId[iji]) {
			 overl++;}
		     }
		   }
#ifdef isDebug
		   cout << " " << jk << " " << ki
			<< " " << int(outBins[jk].posId.size())
			<< " " << int(outBins[ki].posId.size())
			<< " " << overl << endl;
#endif	// #ifdef isDebug
		   if(overl>=int(outBins[ki].posId.size())*0.5) {
#ifdef isDebug
		     cout << " \t delete " << iev << endl;
#endif	// #ifdef isDebug
		     outBins.erase(outBins.begin()+ki);
		   } else {ki++;}
		 } // for(int ki=jk+1;ki<outBins.size();ki++) {
	       } // for(int jk=0;jk<outBins.size();jk++) {
#ifdef isDebug
	       for(int jk=0;jk<outBins.size();jk++) {
		 for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		   cout << " out2 " << ij
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].rawhitInfo.zlay
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.X()
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.Y()
			<< " " << TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]].RawHitPos.Z()
			<< endl;
		 }
	       }   // for(int jk=0;jk<outBins.size();jk++) {
#endif	// #ifdef isDebug
	
	
	/*** Hough Group Info in vector of PosInSpace ***/
	       for(int jk=0;jk<outBins.size();jk++) {
		 tmpHoughGroup.clear();
		 for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		   tmpHoughGroup.push_back(TimeGroup[isj].PosForTrack[outBins[jk].posId[ij]]);
		 } // for(int ij=0;ij<int(outBins[jk].posId.size());ij++) {
		 TimeGroup[isj].HoughCluster.push_back(tmpHoughGroup);
	       } // for(int jk=0;jk<outBins.size();jk++) {
#ifdef isDebug
	       for(int jk=0;jk<int(TimeGroup[isj].HoughCluster.size());jk++) {
		 for(int ij=0;ij<int(TimeGroup[isj].HoughCluster[jk].size());ij++) {
		   cout << " out3 " << jk << " " << ij
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].rawhitInfo.zlay
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].RawHitPos.X()
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].RawHitPos.Y()
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].RawHitPos.Z()
			<< endl;
		 }
	       }
#endif	// #ifdef isDebug
	
	/*** remove hits if in same RPC ***/
	       for(int jk=0;jk<int(TimeGroup[isj].HoughCluster.size());jk++) {
		 for(int ij=0;ij<int(TimeGroup[isj].HoughCluster[jk].size())-1;) {
		   int testb = 0;
		   int testinta = TimeGroup[isj].HoughCluster[jk][ij].rawhitInfo.module;
		   int testintb = TimeGroup[isj].HoughCluster[jk][ij].rawhitInfo.zlay;
		   int testintc = TimeGroup[isj].HoughCluster[jk][ij].rawhitInfo.yrow;
		   int testintd = TimeGroup[isj].HoughCluster[jk][ij].rawhitInfo.xrow;
	    
		   for(int iji=ij+1;iji<int(TimeGroup[isj].HoughCluster[jk].size());) {
		     int testinta1 = TimeGroup[isj].HoughCluster[jk][iji].rawhitInfo.module;
		     int testintb1 = TimeGroup[isj].HoughCluster[jk][iji].rawhitInfo.zlay;
		     int testintc1 = TimeGroup[isj].HoughCluster[jk][iji].rawhitInfo.yrow;
		     int testintd1 = TimeGroup[isj].HoughCluster[jk][iji].rawhitInfo.xrow;
		     if(testinta==testinta1 &&
			testintb==testintb1 &&
			testintc==testintc1 &&
			testintd==testintd1) {
		       TimeGroup[isj].HoughCluster[jk].erase(TimeGroup[isj].HoughCluster[jk].begin()+iji);
		       testb++;
		     } else {iji++;}
		   }
		   if(testb) {
		     // cout << " " << iev << endl;
		     TimeGroup[isj].HoughCluster[jk].erase(TimeGroup[isj].HoughCluster[jk].begin()+ij);
		   } else {ij++;}
		 }
	       }
	
	       /* sorting layerwise : lowest first*/
	       for(int jk=0;jk<int(TimeGroup[isj].HoughCluster.size());jk++) {
		 PosInSpace key;
		 for(int ij=1;ij<int(TimeGroup[isj].HoughCluster[jk].size());ij++) {
		   key = TimeGroup[isj].HoughCluster[jk][ij];
		   int kj = ij-1;
		   while((kj>=0) &&
			 (TimeGroup[isj].HoughCluster[jk][kj].rawhitInfo.zlay>key.rawhitInfo.zlay)) {
		     TimeGroup[isj].HoughCluster[jk][kj+1] = TimeGroup[isj].HoughCluster[jk][kj];
		     kj--;
		   }
		   TimeGroup[isj].HoughCluster[jk][kj+1] = key;
		 }
	       } // for(int jk=0;jk<int(TimeGroup[isj].HoughCluster.size());jk++) {
	
#ifdef isDebug
	       for(int jk=0;jk<int(TimeGroup[isj].HoughCluster.size());jk++) {
		 for(int ij=0;ij<int(TimeGroup[isj].HoughCluster[jk].size());ij++) {
		   cout << " out4 " << jk << " " << ij
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].rawhitInfo.zlay
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].RawHitPos.X()
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].RawHitPos.Y()
			<< " " << TimeGroup[isj].HoughCluster[jk][ij].RawHitPos.Z()
			<< endl;
		 }
	       }
#endif	// #ifdef isDebug
	
	       if(int(TimeGroup[isj].HoughCluster.size())==1) {
		 maxHough = TimeGroup[isj].maxHoughPt;
		 meanHough = TimeGroup[isj].meanHough;
		 stdHough = TimeGroup[isj].stdHough;
		 outtreeh->Fill();
	       }
	
	       // cout << " something " << isj << endl;
	
	     }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {
      
	     // cout << " " << iev << " " << int(TimeGroup.size()) << endl;


      
      
	     for(int isj=0;isj<int(TimeGroup.size());isj++) {
	
	       if(int(TimeGroup[isj].HoughCluster.size())!=1 ||
		  int(TimeGroup[isj].HoughCluster[0].size())<MinLayHit ||
		  TimeGroup[isj].maxHoughPt<MinHoughPt ||
		  TimeGroup[isj].maxHoughPt>MaxHoughPt) {continue;}
	       // cout << " " << TimeGroup[isj].maxHoughPt << endl;
	
	       xyzpos.clear();
	       xyzpos1.clear();
	       xyzId.clear();
	       totLayHit.clear();
	       posCorr.clear();
	       posmulti.clear();
	       xyerr.clear();
	       occulay.clear();
	       
	       vector<PosInSpace> tCluster = TimeGroup[isj].HoughCluster[0];
	       vector<rawlayer>  talllayer = TimeGroup[isj].alllayer;
	       int tbotLayer = tCluster.front().rawhitInfo.zlay-1;
	       int ttopLayer =  tCluster.back().rawhitInfo.zlay+1;
	       if(tbotLayer<0) {tbotLayer=0;}
	       if(ttopLayer>nlayer-1) {ttopLayer=nlayer-1;}
	
	       for(int nl=tbotLayer;nl<=ttopLayer;nl++) {
	  
		 // cout << " nl0 " << nl << endl;
	  
		 TVector3 tpos(0,0,0);
		 TVector3 tpos1(0,0,0);
		 TVector2 terr(0.08,0.08);
		 TVector3 tmul(0,0,0);
		 TVector3 ttothit(0,0,0);
		 TVector3 tposCorr(0,0,0);
		 UShort_t tid;
	  
		 int tm, tx, ty, tzl;
		 if(int(tCluster.size())) {
		   tm = tCluster.front().rawhitInfo.module;
		   tx = tCluster.front().rawhitInfo.xrow;
		   ty = tCluster.front().rawhitInfo.yrow;
		   // tpos = getRPCpos(nl)*(1./strpwidth);
		   // tid = constructRPCId(tCluster.front().rawhitInfo.module,
		   // 			 tCluster.front().rawhitInfo.xrow,
		   // 			 tCluster.front().rawhitInfo.yrow,
		   // 			 nl);
		 } else {
		   getRPCId(xyzId.back(), tm, tx, ty, tzl);
		 }
		 tid = constructRPCId(tm, tx, ty, nl);
		 tpos = tpos1 = getRPCpos(tm, tx, ty, nl)*(1./strpwidth);
		 int toccu = false;
		 for(int ij=0;ij<int(talllayer.size());ij++) {
		   if(talllayer[ij].module==tm &&
		      talllayer[ij].xrow==tx &&
		      talllayer[ij].yrow==ty &&
		      talllayer[ij].zlay==nl) {
		     ttothit.SetX(talllayer[ij].hit[0].size());
		     ttothit.SetY(talllayer[ij].hit[1].size());
		     talllayer.erase(talllayer.begin()+ij);
		     break;
		   }
		 } // for(int ij=0;ij<int(talllayer.size());ij++) {
		 
		 for(int ij=0;ij<int(tCluster.size());ij++) {
		   int tzl = tCluster[ij].rawhitInfo.zlay;
		   if(tzl==nl) {
		     tm = tCluster[ij].rawhitInfo.module;
		     tx = tCluster[ij].rawhitInfo.xrow;
		     ty = tCluster[ij].rawhitInfo.yrow;
		     tpos  = tCluster[ij].RawHitPos  - tCluster[ij].PosCorrection;
		     tpos1 = tCluster[ij].RawHitPos1 - tCluster[ij].PosCorrection;
		     // cout << " " << tpos.Y() << " " << tpos1.Y() << endl;
		     tid = constructRPCId(tm, tx, ty, nl);
		     tposCorr = tCluster[ij].PosCorrection;
		     tmul.SetX(tCluster[ij].rawhitInfo.cluster[0][0].size());
		     tmul.SetY(tCluster[ij].rawhitInfo.cluster[1][0].size());
		     // terr.SetX(poserrsq_ref[int(tmul.X())-1]);
		     // terr.SetY(poserrsq_ref[int(tmul.Y())-1]);
		     terr.SetX(poserrsq[tm][tx][ty][nl][0][int(tmul.X())-1]);
		     terr.SetY(poserrsq[tm][tx][ty][nl][1][int(tmul.Y())-1]);
		     toccu = true;
		     tCluster.erase(tCluster.begin()+ij);
		     break;
		   }
		 } // for(int ij=0;ij<int(tCluster.size());ij++) {
		 xyzpos.push_back(tpos);
		 xyzpos1.push_back(tpos1);
		 xyzId.push_back(tid);
		 posCorr.push_back(tposCorr);
		 posmulti.push_back(tmul);
		 totLayHit.push_back(ttothit);
		 xyerr.push_back(terr);
		 occulay.push_back(toccu);
		 
		 // cout << " nl" << nl
		 //      << " occu " << occulay.back()
		 //      << " poscorr " << posCorr.back().X()
		 //      << " " << posCorr.back().Y()
		 //      << " pos " << xyzpos.back().X()
		 //      << " " << xyzpos.back().Y()
		 //      << " pos1 " << xyzpos1.back().X()
		 //      << " " << xyzpos1.back().Y()
		 //      << " m " << posmulti.back().X()
		 //      << " " << posmulti.back().Y()
		 //      << endl;
		 
		 // cout << " " << nl
		 //      << " " << occulay.back()
		 //      << endl;
	       } // for(int ij=0;ij<int(TimeGroup[isj].HoughCluster[0].size());ij++) {
	
	       // cout<<" xyzpos1.size() "<<int(xyzpos1.size())<<endl;
	       
// #ifndef isSimData
#ifdef isCorrection
	       
	       // h_ndfout->Fill(int(xyzpos1.size()));
	       
	       
	       /**** Efficiency and Position ****/
	       
	       for(int nl=0;nl<int(xyzpos1.size());nl++) {
		 
#ifdef isIter
		 if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=xyzId[nl]) {continue;}
#endif	// #ifdef isIter
		 
		 int tm, tx, ty, tzl;
		 getRPCId(xyzId[nl], tm, tx, ty, tzl);
		 // cout<<" nl "<<nl<<" tzl "<<tzl<<endl;
		 
		 // /* Trigger Check */
		 // unsigned int triggerinfo1[nside] = {0};
		 // for(int nj=0;nj<nside;nj++) {
		 //   for(int ij=0;ij<ntrigLayers1;ij++) {
		 //     triggerinfo1[nj]<<=1;
		 //     for(int nll=0;nll<int(xyzpos1.size());nll++) {
		 //       int tm1, tx1, ty1, tzl1;
		 //       getRPCId(xyzId[nll], tm1, tx1, ty1, tzl1);
		 //       if(tzl1==trigLayers1[ij] && totLayHit[nll][nj]>0) {
		 // 	 triggerinfo1[nj]+=1;}
		 //     }} // for(int ij=0;ij<ntrigLayers1;ij++) {
		 // }   // for(int nj=0;nj<nside;nj++) {
		 
		 // /* Trigger Check */
		 // int startTrglay;// = tzl - ntrigLayers1;
		 // if(tzl>=ntrigLayers1) {startTrglay = tzl - ntrigLayers1;}
		 // else {startTrglay = tzl + 1;}
		 // unsigned int triggerinfo1[nside] = {0};
		 // for(int nj=0;nj<nside;nj++) {
		 //   for(int ij=startTrglay;ij<startTrglay+ntrigLayers1;ij++) {
		 //     // cout<<" ij "<<ij<<endl;
		 //     triggerinfo1[nj]<<=1;
		 //     for(int nll=0;nll<int(xyzpos1.size());nll++) {
		 //       int tm1, tx1, ty1, tzl1;
		 //       getRPCId(xyzId[nll], tm1, tx1, ty1, tzl1);
		 //       if(tzl1==ij && totLayHit[nll][nj]>0) {
		 // 	 triggerinfo1[nj]+=1;}
		 //     }
		 //   } // 
		 // }   // for(int nj=0;nj<nside;nj++) {

		 // int triggerinfo1t[nside] = {0};
		 // for(int nj=0;nj<nside;nj++) {
		 //   for(int ij=max(0,tzl-MinLayHit-1);
		 //       ij<min(nlayer-1,max(0,tzl-MinLayHit-2)+MinLayHit*2+2);ij++) {
		 //     for(int nll=0;nll<int(xyzpos1.size());nll++) {
		 //       int tm1, tx1, ty1, tzl1;
		 //       getRPCId(xyzId[nll], tm1, tx1, ty1, tzl1);
		 //       if(tzl!=tzl1 && tzl1==ij && totLayHit[nll][nj]>0) {
		 // 	 triggerinfo1t[nj]+=1;}
		 //     }
		 //   } // 
		 // }   // for(int nj=0;nj<nside;nj++) {

#ifdef isDebug
		 // cout << " \ttrigger1 " << bitset<10>(triggerinfo1[0])
		 //      << " " << bitset<10>(triggerinfo1[1])
		 //      << " " << bitset<10>(triggerinfo1[0]|triggerinfo1[1])
		 //      << endl;
#endif  // #ifdef isDebug
		 // if(triggerinfo1[0]!=triggerinfo_ref1 && triggerinfo1[1]!=triggerinfo_ref1) {continue;}
		 // if(triggerinfo1t[0]<MinLayHit && triggerinfo1t[1]<MinLayHit) {continue;}
#ifdef isDebug
		 cout << " \t Not Continued"<< endl;
#endif  // #ifdef isDebug
		 // cout<<" \t Not Continued "<<iev<<" "<<triggerinfo1t[0]<<" "<<triggerinfo1t[1]<<endl;
		 
		 
		 vector<bool> tmpOccu = occulay;
		 tmpOccu[nl] = false;
		 LinearVectorFit(0,xyzpos1,xyerr,tmpOccu,
				 slope,inter,xychi2,xyext,xyexterr);
		 
		 TVector3 tmpExt1 = (xyext[nl] -
				     getRPCpos(tm,tx,ty,tzl)*(1./strpwidth));
		 
		 laymul_2DmTotFine[tm][tx][ty][tzl][0]->Fill(tmpExt1.X(),tmpExt1.Y());
		 laymul_2DmTotFine[tm][tx][ty][tzl][1]->Fill(tmpExt1.X(),tmpExt1.Y());
#ifdef is2dpos
		 laymul_2DmTot[tm][tx][ty][tzl][0]->Fill(tmpExt1.X(),tmpExt1.Y());
		 laymul_2DmTot[tm][tx][ty][tzl][1]->Fill(tmpExt1.X(),tmpExt1.Y());
#endif	// #ifdef is2dpos
		 for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		   if(TimeGroup[isj].alllayer[ij].module==tm &&
		      TimeGroup[isj].alllayer[ij].xrow==tx &&
		      TimeGroup[isj].alllayer[ij].yrow==ty &&
		      TimeGroup[isj].alllayer[ij].zlay==tzl) {
		     for(int nj=0;nj<nside;nj++) {
		       int tttothit = int(TimeGroup[isj].alllayer[ij].hit[nj].size());
		       for(int iji=0;iji<tttothit;iji++) {
			 int ttstrp = TimeGroup[isj].alllayer[ij].hit[nj][iji].strp;
			 laymul_2Dstrp[tm][tx][ty][tzl][nj][ttstrp]->Fill(tmpExt1.X(),tmpExt1.Y());
#ifdef is2dpos
			 laymul_2Dpos[tm][tx][ty][tzl][nj]->Fill(tmpExt1.X(),tmpExt1.Y(),ttstrp-int(tmpExt1[nj]));
#endif	// #ifdef is2dpos
		       } // }
		     }	// for(int nj=0;nj<nside;nj++) {
		     break;
		   }
		 } // for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		 

		 // cout<<posmulti[nl].X()<<" "<<posmulti[nl].Y()<<endl;
		 if(posmulti[nl].X()>nmxhits || posmulti[nl].Y()>nmxhits) {continue;}
		 
		 // reverting to original position
		 TVector3 tmpExt = (xyext[nl] -
				    getRPCpos(tm,tx,ty,tzl)*(1./strpwidth) +
				    posCorr[nl]);
		 
		 if(tmpExt.X()<0 || tmpExt.X()>nstrip ||
		    tmpExt.Y()<0 || tmpExt.Y()>nstrip) {continue;}
		 int blockMx = int(tmpExt.X()*blockM/nstrip);
		 int blockMy = int(tmpExt.Y()*blockM/nstrip);
		 				     
		 if(posCorr[nl].Mag()==0) {
		   TVector3 tmpPosCorr(blockposoff[tm][tx][ty][tzl][0][blockMx][blockMy],
				       blockposoff[tm][tx][ty][tzl][1][blockMx][blockMy],0.);
		   tmpExt += tmpPosCorr;;
		 }
		 // TVector3 tmpExt = (xyext[nl] -
		 // 		    getRPCpos(tm,tx,ty,tzl)*(1./strpwidth) +
		 // 		    posCorr[nl]);
		 
		 
		 // cout << " " << tm
		 //      << " " << tx
		 //      << " " << ty
		 //      << " l" << tzl
		 //      << " occu " << tmpOccu[nl]
		 //      << " poscorr " << posCorr[nl].X()
		 //      << " " << posCorr[nl].Y()
		 //      << " x " << xyzpos1[nl].X()
		 //      << " y " << xyzpos1[nl].Y()
		 //      << " m " << posmulti[nl].X()
		 //      << " " << posmulti[nl].Y()
		 //      << " ext " << xyext[nl].X()
		 //      << " " << xyext[nl].Y()
		 //      << " ori " << tmpExt.X()
		 //      << " " << tmpExt.Y()
		 //      << " dev " << xyzpos1[nl].X()-xyext[nl].X()
		 //      << " " << xyzpos1[nl].Y()-xyext[nl].Y()
		 //      << endl;
		 
#ifndef isSimData
		 for(int ij=0;ij<int(xyzpos1.size());ij++) {
		   layPosExtErr[tm][tx][ty][tzl]->Fill(xyexterr[nl][0]);}
		 
		 if(posmulti[nl].X()>0 && posmulti[nl].Y()>0) {
		   layPosR[tm][tx][ty][tzl][0][int(posmulti[nl].X())-1]->Fill(xyzpos1[nl].X()-xyext[nl].X()); // X
		   layPosR[tm][tx][ty][tzl][1][int(posmulti[nl].Y())-1]->Fill(xyzpos1[nl].Y()-xyext[nl].Y()); // Y
		 }
#endif	// #ifndef isSimData
		 
		 inefficiency_tot[tm][tx][ty][tzl]->Fill(tmpExt.X(),tmpExt.Y());
		 // inefficiency_tot30[tm][tx][ty][tzl]->Fill(tmpExt.X(),tmpExt.Y());
		 
		 int ttmul[nside] = {0};
		 double ttpos[nside] = {0};
		 bool isSigAbsent[nside] = {1,1};
		 for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		   // cout << " " << ij << endl;
		   if(TimeGroup[isj].alllayer[ij].module==tm &&
		      TimeGroup[isj].alllayer[ij].xrow==tx &&
		      TimeGroup[isj].alllayer[ij].yrow==ty &&
		      TimeGroup[isj].alllayer[ij].zlay==tzl) {
		     for(int nj=0;nj<nside;nj++) {
		       // 		       // if(int(TimeGroup[isj].alllayer[ij].cluster[nj].size())>0) {
		       // 		       int tttothit = int(TimeGroup[isj].alllayer[ij].hit[nj].size());
		       // 		       // cout<<" "<<sideMark[nj]<<" "<<tttothit<<endl;
		       // 		       // if(tttothit<MaxTotHit) {
		       // 		       for(int iji=0;iji<tttothit;iji++) {
		       // 			 int ttstrp = TimeGroup[isj].alllayer[ij].hit[nj][iji].strp;
		       // 			 laymul_2Dstrp[tm][tx][ty][tzl][nj][ttstrp]->Fill(tmpExt1.X(),tmpExt1.Y());
		       // #ifdef is2dpos
		       // 			 laymul_2Dpos[tm][tx][ty][tzl][nj]->Fill(tmpExt1.X(),tmpExt1.Y(),ttstrp);
		       // #endif	// #ifdef is2dpos
		       // 			 // laymul_2Dpos[tm][tx][ty][tzl][nj]->Fill(tmpExt.X(),tmpExt.Y(),ttstrp-int(tmpExt[nj]));
		       // 			 // cout<<" "<<iev<<" tzl "<<tzl<<" "<<sideMark[nj]<<" "<<tmpExt.X()<<" "<<tmpExt.Y()<<" "<<ttstrp-int(tmpExt[nj])<<endl;
		       // 			 // cout<<" "<<iev<<" tzl "<<tzl<<" "<<sideMark[nj]<<" "<<tmpExt1.X()<<" "<<tmpExt1.Y()<<" "<<ttstrp-int(tmpExt1[nj])<<endl;
		       // 		       } // }
		       for(int iji=0;iji<int(TimeGroup[isj].alllayer[ij].cluster[nj].size());iji++) {
			 int tttmul = int(TimeGroup[isj].alllayer[ij].cluster[nj][iji].size());
			 for(int ijj=0;(ijj<tttmul && isSigAbsent[nj]);ijj++) {
			   int ttstrp = TimeGroup[isj].alllayer[ij].cluster[nj][iji][ijj].strp;
			   if(ttstrp==int(tmpExt[nj])) {
			     isSigAbsent[nj] = false;}
			   // laymul_2Dpos[tm][tx][ty][tzl][nj]->Fill(tmpExt.X(),tmpExt.Y(),ttstrp-int(tmpExt[nj]));
			 } // for(int ijj=0;ijj<tttmul;ijj++) {
			 double tttpos = TimeGroup[isj].alllayer[ij].cluster[nj][iji][0].strp + tttmul*0.5;
			 // cout << " something " << tmpExt[nj] << " " << tttpos << endl;
			 if(ttmul[nj]==0 && fabs(tttpos-tmpExt[nj])<maxPosDev_Effi) {
			   ttmul[nj] = tttmul;
			   ttpos[nj] = tttpos;
			   // break;
			 }
		       } // for(int iji=0;iji<int(TimeGroup[isj].alllayer[ij].cluster[nj].size());iji++) {
		     }	// for(int nj=0;nj<nside;nj++) {
		     break;
		   }
		 } // for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		 // cout<<" isSigAbsent "<<tzl<<" "<<isSigAbsent[0]<<" "<<isSigAbsent[1]<<endl;
		 // if(ttmul[0]==0 && ttmul[1]==0) {
		 //   cout << " " << tzl << " " << ttmul[0] << " " << ttpos[0] << " " << ttmul[1] << " " << ttpos[1] << endl;
		 // }
		 // if(ttmul[0]==3 || ttmul[1]==3) {
		 //   inefficiency_tot3[tm][tx][ty][tzl]->Fill(tmpExt.X(),tmpExt.Y());
		 // }
		 if(ttmul[0]>0) {
		   inefficiency_tot1[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
		   // if(ttmul[0]<=nmxhits) {
		   //   laymul_2D[tm][tx][ty][tzl][0][ttmul[0]-1]->Fill(tmpExt.X(),tmpExt.Y());}
		   // if(ttmul[0]==1) {
		   //   laymul_2Dm1[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());}
		 }
		 if(ttmul[1]>0) {
		   inefficiency_tot1[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
		   // if(ttmul[1]<=nmxhits) {
		   //   laymul_2D[tm][tx][ty][tzl][1][ttmul[1]-1]->Fill(tmpExt.X(),tmpExt.Y());}
		   // if(ttmul[1]==1) {
		   //   laymul_2Dm1[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());}
		 }
		 if(ttmul[0]==0 && ttmul[1]==0) {
		   inefficiency_cor[tm][tx][ty][tzl]->Fill(tmpExt.X(),tmpExt.Y());
		   laymul_2DcTotFine[tm][tx][ty][tzl][0]->Fill(tmpExt1.X(),tmpExt1.Y());
		   laymul_2DcTotFine[tm][tx][ty][tzl][1]->Fill(tmpExt1.X(),tmpExt1.Y());
#ifdef is2dpos
		   laymul_2DcTot[tm][tx][ty][tzl][0]->Fill(tmpExt1.X(),tmpExt1.Y());
		   laymul_2DcTot[tm][tx][ty][tzl][1]->Fill(tmpExt1.X(),tmpExt1.Y());
#endif	// #ifdef is2dpos
		 }
		 if(isSigAbsent[0]==1 && isSigAbsent[1]==1) {
		   inefficiency_cor_str[tm][tx][ty][tzl]->Fill(tmpExt.X(),tmpExt.Y());
		 }
		 if(ttmul[0]>0 && ttmul[1]>0 && isSigAbsent[0]==1 && isSigAbsent[1]==1) {
		   inefficiency_cor_strS[tm][tx][ty][tzl]->Fill(tmpExt.X(),tmpExt.Y());
		 }
		 if(ttmul[0]==0 && ttmul[1]>0) {
		   inefficiency_unc[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
		 } // unc x
		 if(ttmul[0]>0 && ttmul[1]==0) {
		   inefficiency_unc[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
		 } // unc y
		 if(ttmul[0]==0) {
		   inefficiency_unc_tot[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
#ifdef is2dpos
		   laymul_2DuTot[tm][tx][ty][tzl][0]->Fill(tmpExt1.X(),tmpExt1.Y());
#endif	// #ifdef is2dpos
		   laymul_2DuTotFine[tm][tx][ty][tzl][0]->Fill(tmpExt1.X(),tmpExt1.Y());
		 } // unc x
		 if(ttmul[1]==0) {
		   inefficiency_unc_tot[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
#ifdef is2dpos
		   laymul_2DuTot[tm][tx][ty][tzl][1]->Fill(tmpExt1.X(),tmpExt1.Y());
#endif	// #ifdef is2dpos
		   laymul_2DuTotFine[tm][tx][ty][tzl][1]->Fill(tmpExt1.X(),tmpExt1.Y());
		 } // unc y
		 // if(isSigAbsent[0]==1 && ttmul[0]>0) {
		 //   inefficiency_unc_strS_tot[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
		 // } // unc x
		 // if(ttmul[1]>0 && isSigAbsent[1]==1) {
		 //   inefficiency_unc_strS_tot[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
		 // } // unc y
		 if(isSigAbsent[0]==1 && ttmul[1]>0) {
		   inefficiency_unc_str[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
		 } // unc x
		 if(ttmul[0]>0 && isSigAbsent[1]==1) {
		   inefficiency_unc_str[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
		 } // unc y
		 if(isSigAbsent[0]==1) {
		   inefficiency_unc_str_tot[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
		 } // unc x
		 if(isSigAbsent[1]==1) {
		   inefficiency_unc_str_tot[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
		 } // unc y
		 if(totLayHit[nl].X()>0) {
		   triggereffi_evt[tm][tx][ty][tzl][0]->Fill(tmpExt.X(),tmpExt.Y());
		 } // X
		 if(totLayHit[nl].Y()>0) {
		   triggereffi_evt[tm][tx][ty][tzl][1]->Fill(tmpExt.X(),tmpExt.Y());
		 } // Y
		 
		 // if(tmpExt.X()<0 || tmpExt.X()>nstrip ||
		 //    tmpExt.Y()<0 || tmpExt.Y()>nstrip) {continue;}
		 
		 // int blockMx = int(tmpExt.X()*blockM/nstrip);
		 // int blockMy = int(tmpExt.Y()*blockM/nstrip);

#ifndef isIter
#ifdef isStrpMulti
		 layblocktotmul[tm][tx][ty][tzl][0]->Fill(int(tmpExt.X()),int(tmpExt.Y()),posmulti[nl].X()); // X
		 layblocktotmul[tm][tx][ty][tzl][1]->Fill(int(tmpExt.X()),int(tmpExt.Y()),posmulti[nl].Y()); // Y
		 laymul[tm][tx][ty][tzl][0]->Fill(tmpExt.X()-(int(tmpExt.X())+0.5),posmulti[nl].X()); // X
		 laymul[tm][tx][ty][tzl][1]->Fill(tmpExt.Y()-(int(tmpExt.Y())+0.5),posmulti[nl].Y()); // Y
		 layblockmul[tm][tx][ty][tzl][0][blockMx][blockMy]->Fill(tmpExt.X()-(int(tmpExt.X())+0.5),posmulti[nl].X()); // X
		 layblockmul[tm][tx][ty][tzl][1][blockMx][blockMy]->Fill(tmpExt.Y()-(int(tmpExt.Y())+0.5),posmulti[nl].Y()); // Y
#endif	// #ifdef isStrpMulti
#endif	// #ifndef isIter
		 
		 if(posmulti[nl].X()>0 && posmulti[nl].Y()>0) {
#ifdef isIter
		   blockPosR[tm][tx][ty][tzl][0][blockMx][blockMy]->Fill(xyzpos1[nl].X()-xyext[nl].X()); // X
		   blockPosR[tm][tx][ty][tzl][1][blockMx][blockMy]->Fill(xyzpos1[nl].Y()-xyext[nl].Y()); // Y
#endif	// #ifdef isIter
#ifndef isIter
#ifdef isStrpMulti
		   // laymul[tm][tx][ty][tzl][0][int(posmulti[nl].X())-1]->Fill(tmpExt.X()-(int(tmpExt.X())+0.5)); // X
		   // laymul[tm][tx][ty][tzl][1][int(posmulti[nl].Y())-1]->Fill(tmpExt.Y()-(int(tmpExt.Y())+0.5)); // Y
		   // layblockmul[tm][tx][ty][tzl][0][blockMx][blockMy][int(posmulti[nl].X())-1]->Fill(tmpExt.X()-(int(tmpExt.X())+0.5)); // X
		   // layblockmul[tm][tx][ty][tzl][1][blockMx][blockMy][int(posmulti[nl].Y())-1]->Fill(tmpExt.Y()-(int(tmpExt.Y())+0.5)); // Y
		   // layblocktotmul[tm][tx][ty][tzl][0][int(tmpExt.X())][int(tmpExt.Y())]->Fill(int(posmulti[nl].X())); // X
		   // layblocktotmul[tm][tx][ty][tzl][1][int(tmpExt.X())][int(tmpExt.Y())]->Fill(int(posmulti[nl].Y())); // Y
		   // layblocktotmul[tm][tx][ty][tzl][0][blockMx][blockMy]->Fill(int(posmulti[nl].X())); // X
		   // layblocktotmul[tm][tx][ty][tzl][1][blockMx][blockMy]->Fill(int(posmulti[nl].Y())); // Y
		   
		   for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		     // cout << " " << ij << endl;
		     if(TimeGroup[isj].alllayer[ij].module==tm &&
			TimeGroup[isj].alllayer[ij].xrow==tx &&
			TimeGroup[isj].alllayer[ij].yrow==ty &&
			TimeGroup[isj].alllayer[ij].zlay==tzl) {
		       
		       for(int jk=0;jk<int(allrawlay.size());jk++) {
		         if(allrawlay[jk].module==tm &&
		       	    allrawlay[jk].xrow==tx &&
		       	    allrawlay[jk].yrow==ty &&
		       	    allrawlay[jk].zlay==tzl) {
			   
			   for(int nj=0;nj<nside;nj++) {
			     // if(int(TimeGroup[isj].alllayer[ij].cluster[nj].size())>0) {
			     for(int iji=0;iji<int(TimeGroup[isj].alllayer[ij].cluster[nj].size());iji++) {
			       // cout << " " << int(TimeGroup[isj].alllayer[ij].cluster[nj].size()) << endl;
			       int ttmul = int(TimeGroup[isj].alllayer[ij].cluster[nj][iji].size());
			       double ttpos = TimeGroup[isj].alllayer[ij].cluster[nj][iji][0].strp + ttmul*0.5;
			       // cout << " something " << iji << " " << tmpExt[nj] << " " << ttpos << endl;
			       
			       bool isSameContinue = 0;
			       for(int ijk=0;ijk<ttmul;ijk++) {
				 int tmpStrp = TimeGroup[isj].alllayer[ij].cluster[nj][iji][ijk].strp;
				 int tmpTDC = tmpStrp%nTDC;
				 if(int(allrawlay[jk].rawTDC[nj][tmpTDC][0].size())==1) {
				   strpCorrel2D[tm][tx][ty][tzl][nj]->Fill(tmpExt[nj]-0.5,tmpStrp);
				   strpCorrel3D[tm][tx][ty][tzl][nj]->Fill(tmpExt[0]-0.5,tmpExt[1]-0.5,tmpStrp-(tmpExt[nj]-0.5));}
				 if(fabs(tmpExt[nj]-0.5-tmpStrp)<(ttmul*0.5+1.)) {continue;}
				 if(int(allrawlay[jk].rawTDC[nj][tmpTDC][0].size())!=1) {
				   isSameContinue = 1;}
				 // if(int(allrawlay[jk].rawTDC[nj][tmpTDC][0].size())==1) {
				 //   strpCorrel2D[tm][tx][ty][tzl][nj]->Fill(tmpExt[nj]-0.5,tmpStrp);
				 // } else {
				 //   // cout<<" "<<iev<<" "<<sideMark[nj]<<tzl<<" "<<tmpStrp<<" tdc "<<tmpTDC<<endl;
				 //   isSameContinue = 1;}
			       } // for(int ijk=0;ijk<ttmul;ijk++) {
			       
			       // if(fabs(xyext[nl][nj]-ttpos)>(ttmul+1.)*0.5) {
			       if(fabs(tmpExt[nj]-ttpos)>(ttmul+posmulti[nl][nj])*0.5+2.) {
				 // cout << " noiseMul " << sideMark[nj] << tzl << " " << xyext[nl][nj] << " " << ttpos << " " << ttmul << endl;
				 // cout << " noiseMul " << sideMark[nj] << tzl << " " << tmpExt[nj] << " " << ttpos << " " << (ttmul+posmulti[nl][nj]) << endl;
				 // if(int(fabs(tmpExt[nj]-ttpos))%8==0) {cout<<iev<<endl;}
				 // noiseStrpMul2D[tm][tx][ty][tzl][nj]->Fill(ttmul,fabs(xyext[nl][nj]-ttpos));
				 if(!isSameContinue) {
				   // cout<<" noiseMul "<<iev<<" "<<sideMark[nj]<<tzl<<" "<<tmpExt[nj]<<" "<<ttpos<<" "<<(ttmul+posmulti[nl][nj])<<endl;
				   noiseStrpMul2D[tm][tx][ty][tzl][nj]->Fill(ttmul,fabs(tmpExt[nj]-ttpos));
				   // if(fabs(tmpExt[nj]-ttpos)>5. && fabs(tmpExt[nj]-ttpos)<6.) {
				   //   noiseStrpMul[tm][tx][ty][tzl][nj][blockMx][blockMy]->Fill(ttmul);}
				 }
			       }
			     } // for(int iji=0;iji<int(TimeGroup[isj].alllayer[ij].cluster[nj].size());iji++) {
			   } // for(int nj=0;nj<nside;nj++) {
			   break;
			 }
		       } // for(int jk=0;jk<int(allrawlay.size());jk++) {
		       break;
		     }
		   } // for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		   
		   
		   // for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {
		   //   // cout << " " << ij << endl;
		   //   if(TimeGroup[isj].alllayer[ij].module==tm &&
		   // 	TimeGroup[isj].alllayer[ij].xrow==tx &&
		   // 	TimeGroup[isj].alllayer[ij].yrow==ty &&
		   // 	TimeGroup[isj].alllayer[ij].zlay==tzl) {
		   //     for(int nj=0;nj<nside;nj++) {
		   // 	 // if(int(TimeGroup[isj].alllayer[ij].cluster[nj].size())>0) {
		   // 	 for(int iji=0;iji<int(TimeGroup[isj].alllayer[ij].cluster[nj].size());iji++) {
		   // 	   // cout << " " << int(TimeGroup[isj].alllayer[ij].cluster[nj].size()) << endl;
		   // 	   int ttmul = int(TimeGroup[isj].alllayer[ij].cluster[nj][iji].size());
		   // 	   double ttpos = TimeGroup[isj].alllayer[ij].cluster[nj][iji][0].strp + ttmul*0.5;
		   // 	   // cout << " something " << iji << " " << tmpExt[nj] << " " << ttpos << endl;
		   // 	   // if(fabs(xyext[nl][nj]-ttpos)>(ttmul+1.)*0.5) {
		   // 	   if(fabs(tmpExt[nj]-ttpos)>(ttmul+posmulti[nl][nj])*0.5+2.) {
		   // 	     // cout << " noiseMul " << sideMark[nj] << tzl << " " << xyext[nl][nj] << " " << ttpos << " " << ttmul << endl;
		   // 	     // cout << " noiseMul " << sideMark[nj] << tzl << " " << tmpExt[nj] << " " << ttpos << " " << (ttmul+posmulti[nl][nj]) << endl;
		   // 	     // if(int(fabs(tmpExt[nj]-ttpos))%8==0) {cout<<iev<<endl;}
		   // 	     noiseStrpMul[tm][tx][ty][tzl][nj][blockMx][blockMy]->Fill(ttmul);
		   // 	     // noiseStrpMul2D[tm][tx][ty][tzl][nj]->Fill(ttmul,fabs(xyext[nl][nj]-ttpos));
		   // 	     noiseStrpMul2D[tm][tx][ty][tzl][nj]->Fill(ttmul,fabs(tmpExt[nj]-ttpos));
		   // 	   }
		   // 	 } // for(int iji=0;iji<int(TimeGroup[isj].alllayer[ij].cluster[nj].size());iji++) {
		   //     } // for(int nj=0;nj<nside;nj++) {
		   //     break;
		   //   }
		   // } // for(int ij=0;ij<int(TimeGroup[isj].alllayer.size());ij++) {


#endif	      // #ifdef isStrpMulti
#endif	      // #ifndef isIter
		 }   // if(posmulti[nl].X()>0 && posmulti[nl].Y()>0) {
	       } // for(int ij=0;ij<int(xyext.size());ij++) {
	
#endif	// #ifdef isCorrection
// #endif	// #ifndef isSimData
	

	       
	/* Final Position fit */
	       LinearVectorFit(0,xyzpos1,xyerr,occulay,
			       slope,inter,xychi2,xyext,xyexterr);
	       thetad = atan(sqrt(pow(slope.X(),2.)+pow(slope.Y(),2.)));
	       phid = atan2(slope.Y(),slope.X());
	       // for(int ij=0;ij<int(xyzpos1.size());ij++) {
	       // 	 cout<<" xyexterr "<<ij<<" "<<xyexterr[ij][0]<<" "<<xyexterr[ij][1]<<endl;}
	       // cout << " " << iev
	       //      << " theta " << thetad*180./TMath::Pi()
	       //      << " phi " << phid*180./TMath::Pi()
	       //      << endl;
	       hxychi2[0]->Fill(xychi2.X()/(TimeGroup[isj].HoughCluster[0].size()-2.));
	       hxychi2[1]->Fill(xychi2.Y()/(TimeGroup[isj].HoughCluster[0].size()-2.));
	       
	
	
	       /* Time Fit */
	
	       xyztime.clear();
	       xyterr.clear();
	
	       tCluster = TimeGroup[isj].HoughCluster[0];
	       tocculay = occulay;
	       
	       for(int jk=0;jk<int(xyext.size());jk++) {
	  
		 TVector3 tpos(0,0,0);
		 TVector2 terr(1.,1.);
	  
		 int tm, tx, ty, tzl;
		 getRPCId(xyzId[jk], tm, tx, ty, tzl);

		 // cumulative distance travelled by particle
		 if(int(xyztime.size())) {
		   // tpos.SetZ(sqrt(pow(xyext[jk].X()-xyext[jk-1].X(),2.) +
		   // 		   pow(xyext[jk].Y()-xyext[jk-1].Y(),2.) +
		   // 		   pow(xyext[jk].Z()-xyext[jk-1].Z(),2.))*strpwidth +
		   // 	      xyztime.back().Z());}
		   ttxx = xyext[jk];
		   ttyy = xyext[jk-1];
		   tpos.SetZ(calPointDist(ttxx,ttyy)*strpwidth +
			     xyztime.back().Z());}
		 
		 for(int ij=0;ij<int(tCluster.size());ij++) {
		   if(tm==tCluster[ij].rawhitInfo.module &&
		      tx==tCluster[ij].rawhitInfo.xrow &&
		      ty==tCluster[ij].rawhitInfo.yrow &&
		      tzl==tCluster[ij].rawhitInfo.zlay) {
		     TVector3 tmpExt = (xyext[jk] -
					getRPCpos(tm,tx,ty,tzl)*(1./strpwidth)+
					posCorr[jk]);
		     if(tmpExt.X()<0 || tmpExt.X()>nstrip ||
			tmpExt.Y()<0 || tmpExt.Y()>nstrip) {continue;}
		     int blkx = int(tmpExt.X()*blockM/nstrip);
		     int blky = int(tmpExt.Y()*blockM/nstrip);
		     
		     // X-time
		     tmpCluster = tCluster[ij].rawhitInfo.cluster[0][0];
		     for(int jki=0;jki<int(tmpCluster.size());jki++) {
		       if(tmpCluster[jki].strp==int(tmpExt.X())
			  && int(tmpCluster[jki].tdc[0].size())) {
			 tpos.SetX(tmpCluster[jki].tdc[0][0]-
				   tCluster[ij].rawhitInfo.tdc_ref[0]-
				   tmpExt.Y()*strpwidth*spdl-
				   blocktimeoff[tm][tx][ty][tzl][0][blkx][blky]-
				   multimeoff[tm][tx][ty][tzl][0][int(posmulti[jk].X())-1]
				   );
		       }
		     }
		     // Y-time
		     tmpCluster = tCluster[ij].rawhitInfo.cluster[1][0];
		     for(int jki=0;jki<int(tmpCluster.size());jki++) {
		       if(tmpCluster[jki].strp==int(tmpExt.Y())
			  && int(tmpCluster[jki].tdc[0].size())) {
			 tpos.SetY(tmpCluster[jki].tdc[0][0]-
				   tCluster[ij].rawhitInfo.tdc_ref[0]-
				   tmpExt.X()*strpwidth*spdl-
				   blocktimeoff[tm][tx][ty][tzl][1][blkx][blky]-
				   multimeoff[tm][tx][ty][tzl][1][int(posmulti[jk].Y())-1]
				   );
		       }
		     }
		     // terr.SetX(timeErr_ref[int(posmulti[jk].X())-1]);
		     // terr.SetY(timeErr_ref[int(posmulti[jk].Y())-1]);
		     terr.SetX(timeerrsq[tm][tx][ty][tzl][0][int(posmulti[jk].X())-1]);
		     terr.SetY(timeerrsq[tm][tx][ty][tzl][1][int(posmulti[jk].Y())-1]);
		     tCluster.erase(tCluster.begin()+ij);
		     break;
		   }
		 } // for(int ij=0;ij<int(tCluster.size());ij++) {
		 xyztime.push_back(tpos);
		 xyterr.push_back(terr);
		 if(xyztime.back().X()==0 || xyztime.back().Y()==0) {
		   tocculay[jk] = false;}
		 // cout << " " << jk
		 //      << " " << xyztime.back().X()
		 //      << " " << xyztime.back().Y()
		 //      << " " << xyztime.back().Z()
		 //      << endl;
	       } // for(int jk=0;jk<int(xyext.size());jk++) {
	
#ifndef isSimData
#ifdef isCorrection
	
	       // time iteration over all points
	       for(int nl=0;nl<int(xyztime.size());nl++) {

#ifdef isIter
		 if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=xyzId[nl]) {continue;}
#endif	// #ifdef isIter
	  
		 int tm, tx, ty, tzl;
		 getRPCId(xyzId[nl], tm, tx, ty, tzl);

		 vector<bool> tmpOccu = tocculay;
		 tmpOccu[nl] = false;
		 LinearVectorFit(1,xyztime,xyterr,tmpOccu,
				 tslope,tinter,xytchi2,xytext,xytexterr);
		 if(posmulti[nl].X()>nmxhits || posmulti[nl].Y()>nmxhits) {continue;}
	  
		 // reverting to original position
		 TVector3 tmpExt = (xyext[nl] -
				    getRPCpos(tm,tx,ty,tzl)*(1./strpwidth) +
				    posCorr[nl]);

		 for(int ij=0;ij<int(xyztime.size());ij++) {
		   layTimeExtErr[tm][tx][ty][tzl]->Fill(xytexterr[nl][0]);}
	  
		 if(posmulti[nl].X()>0 && posmulti[nl].Y()>0) {
		   layTimeR[tm][tx][ty][tzl][0][int(posmulti[nl].X())-1]->Fill(xyztime[nl].X()-xytext[nl].X()); // X
		   layTimeR[tm][tx][ty][tzl][1][int(posmulti[nl].Y())-1]->Fill(xyztime[nl].Y()-xytext[nl].Y()); // Y
		 }
	  
		 if(tmpExt.X()<0 || tmpExt.X()>nstrip ||
		    tmpExt.Y()<0 || tmpExt.Y()>nstrip) {continue;}
		 if(posmulti[nl].X()!=1 && posmulti[nl].Y()!=1) {continue;}
	  
		 int blockMx = int(tmpExt.X()*blockM/nstrip);
		 int blockMy = int(tmpExt.Y()*blockM/nstrip);
#ifdef isIter
		 blockTimeR[tm][tx][ty][tzl][0][blockMx][blockMy]->Fill(xyztime[nl].X()-xytext[nl].X()); // X
		 blockTimeR[tm][tx][ty][tzl][1][blockMx][blockMy]->Fill(xyztime[nl].Y()-xytext[nl].Y()); // Y
#endif	// #ifdef isIter
	  
	       } // for(int nl=0;nl<int(xyztime.size());nl++) {

#endif	// #ifdef isCorrection
#endif	// #ifndef isSimData
	
	
	/* Final Time Fit  */
	       LinearVectorFit(0,xyztime,xyterr,tocculay,
			       tslope,tinter,xytchi2,xytext,xytexterr);
	       // for(int ij=0;ij<int(xyztime.size());ij++) {
	       // cout<<" xytexterr "<<ij<<" "<<xytexterr[ij][0]<<" "<<xytexterr[ij][1]<<endl;}
	       // cout << " " << tslope.X()*cval
	       //      << " " << tslope.Y()*cval << endl;
	       if(tslope.X()!=0) {xybeta[0]->Fill(tslope.X()*cval);}
	       if(tslope.Y()!=0) {xybeta[1]->Fill(tslope.Y()*cval);}
	




#ifndef isCorrection
	
	
	       /* Reconstuction */
		
	       TrackInfo inPoints;
	       
	       inPoints.xyzpos.clear();
	       inPoints.xyzId.clear();
	       inPoints.posCorr.clear();
	       inPoints.posmulti.clear();
	       inPoints.xyerr.clear();
	       inPoints.xyztime.clear();
	       inPoints.xyterr.clear();
	       
	       tCluster = TimeGroup[isj].HoughCluster[0];
	       for(int nl=0;nl<int(tCluster.size());nl++) {
		 
		 // cout << " nl0 " << nl << endl;
	  
		 int tm = tCluster[nl].rawhitInfo.module;
		 int tx = tCluster[nl].rawhitInfo.xrow;
		 int ty = tCluster[nl].rawhitInfo.yrow;
		 int tzl = tCluster[nl].rawhitInfo.zlay;
		 TVector3 tpos = (tCluster[nl].RawHitPos1 -
				  tCluster[nl].PosCorrection);
		 tpos *= strpwidth;
		 UShort_t tid = constructRPCId(tm, tx, ty, tzl);
		 TVector3 tposCorr = tCluster[nl].PosCorrection*strpwidth;
		 TVector3 tmul(tCluster[nl].rawhitInfo.cluster[0][0].size(),
			       tCluster[nl].rawhitInfo.cluster[1][0].size(),0);
		 // TVector2 terr(poserrsq_ref[int(tmul.X())-1],
		 // 	       poserrsq_ref[int(tmul.Y())-1]);
		 TVector2 terr(poserrsq[tm][tx][ty][tzl][0][int(tmul.X())-1],
			       poserrsq[tm][tx][ty][tzl][1][int(tmul.Y())-1]);
		 terr *= strpwidth*strpwidth;
	  
		 inPoints.xyzpos.push_back(tpos);
		 inPoints.xyzId.push_back(tid);
		 inPoints.posCorr.push_back(tposCorr);
		 inPoints.posmulti.push_back(tmul);
		 inPoints.xyerr.push_back(terr);
	  
	  
		 TVector3 ttpos(0,0,0);
		 
		 if(nl>0) {
		   ttxx = inPoints.xyzpos[nl];
		   ttyy = inPoints.xyzpos[nl-1];
		   ttpos.SetZ(calPointDist(ttxx,ttyy) +
			      inPoints.xyztime.back().Z());
		   // ttpos.SetZ(sqrt(pow(inPoints.xyzpos[nl].X()-
		   // 			inPoints.xyzpos[nl-1].X(),2.) +
		   // 		    pow(inPoints.xyzpos[nl].Y()-
		   // 			inPoints.xyzpos[nl-1].Y(),2.) +
		   // 		    pow(inPoints.xyzpos[nl].Z()-
		   // 			inPoints.xyzpos[nl-1].Z(),2.)) +
		   // 	       inPoints.xyztime.back().Z());
		 }
	  
		 TVector3 tmpExt = (tpos -
				    getRPCpos(tm,tx,ty,tzl) +
				    tposCorr);
		 tmpExt *= 1./strpwidth;
		 if(tmpExt.X()<0 || tmpExt.X()>nstrip ||
		    tmpExt.Y()<0 || tmpExt.Y()>nstrip) {continue;}
		 int blkx = int(tmpExt.X()*blockM/nstrip);
		 int blky = int(tmpExt.Y()*blockM/nstrip);
	  
		 // X-time
		 tmpCluster = tCluster[nl].rawhitInfo.cluster[0][0];
		 for(int jki=0;jki<int(tmpCluster.size());jki++) {
		   if(tmpCluster[jki].strp==int(tmpExt.X())
		      && int(tmpCluster[jki].tdc[0].size())) {
		     ttpos.SetX(tmpCluster[jki].tdc[0][0]-
				tCluster[nl].rawhitInfo.tdc_ref[0]-
				tmpExt.Y()*strpwidth*spdl-
				blocktimeoff[tm][tx][ty][tzl][0][blkx][blky]-
				multimeoff[tm][tx][ty][tzl][0][int(tmul.X())-1]);
		   }
		 }
		 // Y-time
		 tmpCluster = tCluster[nl].rawhitInfo.cluster[1][0];
		 for(int jki=0;jki<int(tmpCluster.size());jki++) {
		   if(tmpCluster[jki].strp==int(tmpExt.Y())
		      && int(tmpCluster[jki].tdc[0].size())) {
		     ttpos.SetY(tmpCluster[jki].tdc[0][0]-
				tCluster[nl].rawhitInfo.tdc_ref[0]-
				tmpExt.X()*strpwidth*spdl-
				blocktimeoff[tm][tx][ty][tzl][1][blkx][blky]-
				multimeoff[tm][tx][ty][tzl][1][int(tmul.Y())-1]);
		   }
		 } // [tm][tx][ty][tzl][0]
		 // TVector2 tterr(timeErr_ref[int(tmul.X())-1],
		 // 		timeErr_ref[int(tmul.Y())-1]);
		 TVector2 tterr(timeerrsq[tm][tx][ty][tzl][0][int(tmul.X())-1],
				timeerrsq[tm][tx][ty][tzl][1][int(tmul.Y())-1]);
		 inPoints.xyztime.push_back(ttpos);
		 inPoints.xyterr.push_back(tterr);
		 
		 // cout << " " << nl
		 //      << " " << inPoints.xyzpos.back().X()
		 //      << " " << inPoints.xyzpos.back().Y()
		 //      << " " << inPoints.xyzpos.back().Z()
		 //      // << " " << inPoints.xyzpos.back()[0]
		 //      // << " " << inPoints.xyzpos.back()[1]
		 //      // << " " << inPoints.xyzpos.back()[2]
		 //   // << " " << inPoints.xyzId.back().X()
		 //      << " " << inPoints.posCorr.back().X()
		 //      << " " << inPoints.posCorr.back().Y()
		 //      << " " << inPoints.posmulti.back().X()
		 //      << " " << inPoints.posmulti.back().Y()
		 //      << " " << inPoints.xyerr.back().X()
		 //      << " " << inPoints.xyerr.back().Y()
		 //      << " " << inPoints.xyztime.back().X()
		 //      << " " << inPoints.xyztime.back().Y()
		 //      << " " << inPoints.xyztime.back().Z()
		 //      << " " << inPoints.xyterr.back().X()
		 //      << " " << inPoints.xyterr.back().Y()
		 //      << endl;
	  
	       } // for(int nl=0;nl<int(tCluster.size());nl++) {
	       
	
	       trkLen = inPoints.xyztime.back().Z();
	
	       /*** checking up or down going ***/
	       vector<bool> tOccu; tOccu.clear();
	       LinearVectorFit(0,inPoints.xyztime,inPoints.xyterr,tOccu,
			       tslope,tinter,xytchi2,xytext,xytexterr);
	
	       if(tslope.X()>0 && tslope.Y()>0) { // reversing the track
		 for(int jk=0;jk<int(inPoints.xyzpos.size());jk++) {
		   inPoints.xyzpos.insert(inPoints.xyzpos.begin()+jk,
					  inPoints.xyzpos.back());
		   inPoints.xyzpos.pop_back();
		   inPoints.xyzId.insert(inPoints.xyzId.begin()+jk,
					 inPoints.xyzId.back());
		   inPoints.xyzId.pop_back();
		   inPoints.posCorr.insert(inPoints.posCorr.begin()+jk,
					   inPoints.posCorr.back());
		   inPoints.posCorr.pop_back();
		   inPoints.posmulti.insert(inPoints.posmulti.begin()+jk,
					    inPoints.posmulti.back());
		   inPoints.posmulti.pop_back();
		   inPoints.xyerr.insert(inPoints.xyerr.begin()+jk,
					 inPoints.xyerr.back());
		   inPoints.xyerr.pop_back();
		   inPoints.xyztime.insert(inPoints.xyztime.begin()+jk,
					   inPoints.xyztime.back());
		   inPoints.xyztime.pop_back();
		   inPoints.xyterr.insert(inPoints.xyterr.begin()+jk,
					  inPoints.xyterr.back());
		   inPoints.xyterr.pop_back();

		   // cout << " up-down " << jk
		   //      << " " << inPoints.xyzpos[jk].X()
		   //      << " " << inPoints.xyzpos[jk].Y()
		   //      << " " << inPoints.xyzpos[jk].Z()
		   //   // << " " << inPoints.xyzId[jk].X()
		   //      << " " << inPoints.posCorr[jk].X()
		   //      << " " << inPoints.posCorr[jk].Y()
		   //      << " " << inPoints.posmulti[jk].X()
		   //      << " " << inPoints.posmulti[jk].Y()
		   //      << " " << inPoints.xyerr[jk].X()
		   //      << " " << inPoints.xyerr[jk].Y()
		   //      << " " << inPoints.xyztime[jk].X()
		   //      << " " << inPoints.xyztime[jk].Y()
		   //      << " " << inPoints.xyztime[jk].Z()
		   //      << " " << inPoints.xyterr[jk].X()
		   //      << " " << inPoints.xyterr[jk].Y()
		   //      << endl;
	    
		 }
	       } // if(tslope.X()>0 && tslope.Y()>0) {
	       
	
	       const int ndfi = inPoints.xyzpos.size();
	       nhits_finder = ndfi;
	       // cout<<" nhits_finder "<<nhits_finder<<endl;
	       
	

	       // /* circle fit */
	       // const int ndfi = inPoints.xyzpos.size();
	       // reals DataX[ndfi];
	       // reals DataY[ndfi];
	       // reals ErrY[ndfi];
	       // for(int jk=0;jk<ndfi;jk++) {
	       //   DataX[jk] = inPoints.xyzpos[jk].Z();
	       //   DataY[jk] = inPoints.xyzpos[jk].X();
	       //   ErrY[jk]  = inPoints.xyerr[jk].X();
	       //   // cout << " " << jk << " " << DataX[jk] << " " << DataY[jk] << endl;
	       // }
	       // reals LambdaIni=0.00001;
	       // Data data1(ndfi,DataX,DataY,ErrY);
	
	       // double ttang = atan(-(DataX[0]-DataX[ndfi-1])/(DataY[0]-DataY[ndfi-1]));
	       // // cout << " " << iev << " ang " << ttang
	       // //      << " A " << 0.5*(DataX[0]+DataX[ndfi-1])+4.*cos(ttang)
	       // //      << " " << 0.5*(DataY[0]+DataY[ndfi-1])+4.*sin(ttang)
	       // //      << " B " << 0.5*(DataX[0]+DataX[ndfi-1])-4.*cos(ttang)
	       // //      << " " << 0.5*(DataY[0]+DataY[ndfi-1])-4.*sin(ttang)
	       // //      << endl;
	       // Circle circle,circleIni(0.5*(DataX[0]+DataX[ndfi-1])+4.*cos(ttang),0.5*(DataY[0]+DataY[ndfi-1])+4.*sin(ttang),4.);
	       // int code0 = CircleFitByChernovLesort(data1,circleIni,LambdaIni,circle);
	       // Circle circle1,circleIni1(0.5*(DataX[0]+DataX[ndfi-1])-4.*cos(ttang),0.5*(DataY[0]+DataY[ndfi-1])-4.*sin(ttang),4.);
	       // int code1 = CircleFitByChernovLesort(data1,circleIni1,LambdaIni,circle1);
	       // // cout.precision(7);
	       // // if ((code == 1)||(code==2)) cout << "\n Geometric circle by Chernov-Lesort did not converge. Iterations maxed out.\n";
	       // // if (code == 3) cout << "\n Geometric circle by Chernov-Lesort did not converge. Fitting circle too big.\n";
	       // // if(code == 0) cout << " " << iev << " center ("  << circle.a <<","<< circle.b <<")  radius "
	       // // 		   << circle.r << "  sigma " << circle.s << "  iterations: " << circle.i << endl;
	       // // if(code1 == 0) cout << " " << iev << " center ("  << circle1.a <<","<< circle1.b <<")  radius "
	       // // 		   << circle1.r << "  sigma " << circle1.s << "  iterations: " << circle1.i << endl;
	
	       // double cradii = circle.s>circle1.s? circle1.r : circle.r;
	       // double ca = circle.s>circle1.s? circle1.a : circle.a;
	       // double cb = circle.s>circle1.s? circle1.b : circle.b;
	       // int code = circle.s>circle1.s? code1 : code0;
	       // double cs = circle.s>circle1.s? circle1.s : circle.s;
	       // // cout << " " << iev << " cradii " << cradii << " " << ca << " " << cb << endl;

	       // circleS[0]->Fill(cs);
	
	       // if(code==0
	       //    && xychi2.Y()/(TimeGroup[isj].HoughCluster[0].size()-2.)<2.
	       //    && cs<0.5
	       //    // && cradii>minRadii
	       //    // && cradii<maxRadii
	       //    // && circle.s/(strpwidth*strpwidth*ndfi)<5
	       //    ) {
	  
	       //   // cout << " " << iev << " center ("  << ca <<","<< cb
	       //   //      <<")  radius " << cradii << "  sigma " << cs // << "  iterations: " << circle.i
	       //   //      << endl;
	  
	       //   // circleS[0]->Fill(cs);
	  
	  	  
	       // }
	

	

	       /* Getting pt value */
	
	       xyzpos.clear();
	       xyzerr.clear();
	       for(int jk=0;jk<ndfi;jk++) {
		 xyzpos.push_back(inPoints.xyzpos[ndfi-1-jk]);
		 TVector3 xxx(inPoints.xyerr[ndfi-1-jk].X(),
			      inPoints.xyerr[ndfi-1-jk].Y(),0.);
		 xyzerr.push_back(xxx);
		 // cout << " " << jk
		 //      << " " << xyzpos[jk].X()/strpwidth
		 //      << " " << xyzpos[jk].Y()/strpwidth
		 //      << " " << xyzpos[jk].Z()/strpwidth
		 //      << " " << xyzerr[jk].X()
		 //      << " " << xyzerr[jk].Y()
		 //      << " " << xyzerr[jk].Z()
		 //      << endl;
		 if(jk>1 &&
		    // calPointDist(xyzpos[0],xyzpos[jk])>circleLen) {
		    calPointDist(xyzpos[0],xyzpos[jk])>=calPointDist(xyzpos[0],xyzpos[1])*circlePt*1.1) {
		   break;}
	       }
	       const int ndfi3 = int(xyzpos.size());
               // cout<<" ndfi3 "<<ndfi3<<endl;
	       
	       TVector3 MomIniDir(xyzpos[1].X()-xyzpos[0].X(),
				  xyzpos[1].Y()-xyzpos[0].Y(),
				  xyzpos[1].Z()-xyzpos[0].Z());
	       MomIniDir *= 1./MomIniDir.Mag();
	       // cout << " MomIniDir " << MomIniDir.X()
	       //      << " " << MomIniDir.Y()
	       //      << " " << MomIniDir.Z() << endl;
	       
	       TVector3 MagField = GetMagneticField(xyzpos[0].X(),
						    xyzpos[0].Y(),
						    xyzpos[0].Z()+(airGap+ironThickness)*0.5);
	       // cout << " MagField " << MagField.X()
	       //      << " " << MagField.Y()
	       //      << " " << MagField.Z() << endl;
	
	       TVector3 MomAlongMagField = MomIniDir;
	       MomAlongMagField.RotateUz(MagField.Unit());
	       double pZ = MomAlongMagField.Z();
	       double pT = MomAlongMagField.Perp();
	       // cout << " pZ " << pZ << " pT " << pT << endl;
	       
	       TVector3 forceDirection = MomIniDir.Cross(MagField); // F = Q*PxB
	       double rotAngle = -forceDirection.Phi();
	       forceDirection.RotateZ(rotAngle); // Rotating around Z, making Y component zero
	       // cout << " rotAngle " << rotAngle*180./TMath::Pi() << endl;
	       // cout << " forceDirection " << forceDirection.X()
	       //      << " " << forceDirection.Y()
	       //      << " " << forceDirection.Z() << endl;
	       
	       for(int jk=0;jk<ndfi3;jk++) {
		 xyzpos[jk].RotateZ(rotAngle);
		 xyzerr[jk].RotateZ(rotAngle);
		 // cout << " " << jk
		 //      << " " << xyzpos[jk].X()/strpwidth
		 //      << " " << xyzpos[jk].Y()/strpwidth
		 //      << " " << xyzpos[jk].Z()/strpwidth
		 //      << " " << xyzerr[jk].X()
		 //      << " " << xyzerr[jk].Y()
		 //      << " " << xyzerr[jk].Z()
		 //      << endl;
	       }
	       
	       reals DataX3[ndfi3];
	       reals DataY3[ndfi3];
	       reals ErrY3[ndfi3];
	       for(int jk=0;jk<ndfi3;jk++) {
		 DataX3[jk] = xyzpos[jk].Z();
		 DataY3[jk] = xyzpos[jk].X();
		 ErrY3[jk]  = fabs(xyzerr[jk].X());
		 // cout << " " << jk << " " << DataX[jk] << " " << DataY[jk] << endl;
	       }
	       reals LambdaIni3=0.00001;
	       Data data13(ndfi3,DataX3,DataY3,ErrY3);
	       // Circle circle3,circleIni3(0.5,DataY3[0],2.);
	       // int codet = CircleFitByChernovLesort(data13,circleIni3,LambdaIni3,circle3);
	       // if(codet!=0) {continue;}
	       // cout<<" "<<iev<<" cradii "<<circle3.r<<" "<<circle3.a/strpwidth<<" "<<circle3.b/strpwidth<<" s "<<circle3.s<<endl;
	       // double cradii3 = circle3.r;
	       // double ca3 = circle3.a;
	       // double cb3 = circle3.b;
	       
	       double ttang = atan(-(DataX3[0]-DataX3[ndfi3-1])/(DataY3[0]-DataY3[ndfi3-1]));
	       Circle circle31,circleIni31(0.5*(DataX3[0]+DataX3[ndfi3-1])+4.*cos(ttang),0.5*(DataY3[0]+DataY3[ndfi3-1])+4.*sin(ttang),4.);
	       int code31 = CircleFitByChernovLesort(data13,circleIni31,LambdaIni3,circle31);
	       Circle circle32,circleIni32(0.5*(DataX3[0]+DataX3[ndfi3-1])-4.*cos(ttang),0.5*(DataY3[0]+DataY3[ndfi3-1])-4.*sin(ttang),4.);
	       int code32 = CircleFitByChernovLesort(data13,circleIni32,LambdaIni3,circle32);
	       
	       // if(fabs(circle31.s-circle32.s)>0.0001) {
	       // 	 cout<<" "<<iev<<endl;
	       // 	 cout<<" circle31 "<<circle31.r<<" "<<circle31.a/strpwidth<<" "<<circle31.b/strpwidth<<" s "<<circle31.s<<endl;
	       // 	 cout<<" circle32 "<<circle32.r<<" "<<circle32.a/strpwidth<<" "<<circle32.b/strpwidth<<" s "<<circle32.s<<endl;
	       // }
	       
	       int codet = (circle31.s<circle32.s)?code31:code32;
	       if(codet!=0) {continue;}
	       
	       double cradii3 = (circle31.s<circle32.s)?circle31.r:circle32.r;
	       double ca3 = (circle31.s<circle32.s)?circle31.a:circle32.a;
	       double cb3 = (circle31.s<circle32.s)?circle31.b:circle32.b;
	       double cs3 = (circle31.s<circle32.s)?circle31.s:circle32.s;
	       // cout<<" "<<iev<<" cradii "<<cradii3<<" "<<ca3/strpwidth<<" "<<cb3/strpwidth<<" s "<<cs3<<endl;
	       
	       double pXZ = 0.3*uniformField*cradii3;
	       double parMom = pXZ*sqrt(1. + (pZ/pT)*(pZ/pT));
	       // cout<<" "<<iev<<" pXZ "<<pXZ<<" parMom "<<parMom<<endl;
	       MomIniDir.SetMag(parMom);
	       
	       // TVector3 startPt(xyzpos.back().X()-cb3,
	       // 		 0.,xyzpos.back().Z()-ca3);
	       // TVector3   endPt(xyzpos.front().X()-cb3,
	       // 		 0.,xyzpos.front().Z()-ca3);
	       // double skangle = endPt.Angle(startPt);
	       // cout << " " << iev
	       //      << " skangle " << skangle*180./TMath::Pi()
	       //      << endl;
	       
	       TVector3   endPt(xyzpos.back().X()-cb3,
				xyzpos.back().Z()-ca3,0.);
	       TVector3 startPt(xyzpos.front().X()-cb3,
				xyzpos.front().Z()-ca3,0.);
	       double skangle = endPt.Phi()-startPt.Phi();
	       // cout << " 0 " << iev
	       //      << " skangle " << skangle*180./TMath::Pi()
	       //      << endl;
	       if(skangle<-TMath::Pi()) {
		 skangle += 2.*TMath::Pi();
	       } else if(skangle>TMath::Pi()) {
		 skangle -= 2.*TMath::Pi();}
	       // cout << " 1 " << iev
	       //      << " skangle " << skangle*180./TMath::Pi()
	       //      << endl;
	       int chargeVal = skangle>0.?1.:-1.;
	       // cout<<" chargeVal "<<chargeVal<<endl;
	       
	       csFill = cs3;
	       crFill = cradii3;
	       cmomFill = parMom*chargeVal;
	       
	       // bool isLeft = atan2(cb3-PosAfterTime.front().posXY[0],ca3-PosAfterTime.front().posZ)<0 ? true : false;
	       // int chargeVal = isLeft?1:-1;
	       
	       
	       // cout << "\n\n\n" << iev << endl;
	
#ifndef isLifetime
	       
	       inPoints1 = inPoints;
	       inPoints1.ipos = inPoints1.xyzpos.back();
	       // inPoints1.ipos[0] = inPoints1.xyzpos.back().X();
	       // inPoints1.ipos[1] = inPoints1.xyzpos.back().Y();
	       // inPoints1.ipos[2] = inPoints1.xyzpos.back().Z();
	       // inPoints1.iMom.SetMagThetaPhi(1.,TMath::Pi()*160./180.,0.);
	       inPoints1.iMom = MomIniDir;
	       inPoints1.charge = chargeVal;
	       
	       vector<Double_t> vstart =
		 {inPoints1.ipos[0],
		  inPoints1.ipos[1],
		  inPoints1.iMom.Mag(),
		  inPoints1.iMom.Theta(),
		  inPoints1.iMom.Phi()};
	       
	       vector<Double_t> vstep =
		 {0.001, 0.001,
		  0.01, 0.001, 0.001};
	       
	       vector<Double_t> vlimL =
		 {inPoints1.ipos[0]-posSpread,
		  inPoints1.ipos[1]-posSpread,
		  // (inPoints1.iMom.Mag()-momSpread<minMom?
		  //  minMom:inPoints1.iMom.Mag()-momSpread),
		  inPoints1.iMom.Mag()*(1.-momSpread),
		  inPoints1.iMom.Theta()-thetaSpread,
		  inPoints1.iMom.Phi()-phiSpread};

	       vector<Double_t> vlimU =
		 {inPoints1.ipos[0]+posSpread,
		  inPoints1.ipos[1]+posSpread,
		  // inPoints1.iMom.Mag()+momSpread,
		  inPoints1.iMom.Mag()*(1.+momSpread),
		  inPoints1.iMom.Theta()+thetaSpread,
		  inPoints1.iMom.Phi()+phiSpread};
	       
	       vector<TString> parnames =
		 {"vx","vy",
		  "pmag","ptheta","pphi"};
	       
	       TMinuit *tMinuit = new TMinuit(5);
	       tMinuit->SetPrintLevel(-1);
	       tMinuit->SetFCN(fcn);
	
	       Double_t arglist[10];
	       Int_t ierflg = 0;
	
	       arglist[0] = 1;
	       tMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
	
	       for(int npr=0;npr<5;npr++) {
		 tMinuit->mnparm(npr, parnames[npr], vstart[npr], vstep[npr],
				 vlimL[npr], vlimU[npr], ierflg);
	       }
	
	
	       // Set starting values and step sizes for parameters
	       // tMinuit->mnparm(0, "vx", vstart[0], vstep[0],
	       // 		vlimL[0], vlimU[0], ierflg);
	       // tMinuit->mnparm(1, "vy", vstart[1], vstep[1],
	       // 		vlimL[1], vlimU[1], ierflg);
	       // tMinuit->mnparm(2, "pmag", vstart[2], vstep[2],
	       // 		vlimL[2], vlimU[2], ierflg);
	       // tMinuit->mnparm(3, "ptheta", vstart[3], vstep[3],
	       // 		vlimL[3], vlimU[3], ierflg);
	       // tMinuit->mnparm(4, "pphi", vstart[4], vstep[4],
	       // 		vlimL[4], vlimU[4], ierflg);
	       // tMinuit->mnparm(3, "ptheta", vstart[3], vstep[3],
	       // 		TMath::Pi()*0.5, TMath::Pi(), ierflg);
	       // tMinuit->mnparm(4,   "pphi", vstart[4], vstep[4],
	       // 		-TMath::Pi(), TMath::Pi(), ierflg);
	
	
	       // Now ready for minimization step
	       arglist[0] = 500.;
	       arglist[1] = 0.001;
	       tMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
	
	       // // Print results
	       // Double_t amin,edm,errdef;
	       // Int_t nvpar,nparx,icstat;
	       // tMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	       // //tMinuit->mnprin(3,amin);
	
	       TString tmpname;
	       Double_t parval[5];
	       Double_t parerr[5];
	       Double_t lupval, llowval;
	       Int_t tmpc;
	
	       for(int npr=0;npr<5;npr++) {
		 tMinuit->mnpout(npr,tmpname,parval[npr],parerr[npr],
				 lupval,llowval,tmpc);
		 // cout << " par " << npr
		 //      << " " << parnames[npr]
		 //      << " " << parval[npr]
		 //      << " " << parerr[npr]
		 //      << endl;
	       }
	       // cout << " " << iev
	       //      << " theta " << parval[3]*180./TMath::Pi()
	       //      << " phi " << parval[4]*180./TMath::Pi()
	       //      << endl;
	
	       /* final propagation */
	       inPoints1.ipos[0] = parval[0];
	       inPoints1.ipos[1] = parval[1];
	       inPoints1.iMom.SetMagThetaPhi(parval[2],parval[3],parval[4]);
	       PropagateTrack(inPoints1);
	       
	       chi2n = inPoints1.chi2;
	       // nhits_finder = ndfi;
	       momout = parval[2]*inPoints1.charge;
	       momerr = parerr[2];
	       theout = parval[3];
	       theerr = parerr[3];
	       phiout = parval[4];
	       phierr = parerr[4];
	       // dEdx   = parval[5];
	       // dEdxerr = parerr[5];
	       eloss = inPoints1.energyLoss;
	       
#ifdef isSimData
	       // if(nhits_finder>5) {
	       // 	 cout<<" iev "<<iev
	       // 	     <<" momin "<<momInFill
	       // 	     <<" momout "<<momout<<" err "<<momerr
	       // 	     <<" eloss "<< eloss
	       // 	     <<" nhits "<<nhits_finder
	       // 	     <<" chi2ndf "<<chi2n/(nhits_finder-5.)
	       // 	     <<endl<<endl;
	       // }
#endif

	       if(tMinuit) {delete tMinuit;}
	       
#endif	// #ifndef isLifetime
	       
#ifdef isLifetime
	       ePosInLayFill[0] = -1; ePosInLayFill[1] = -1;
	       if(isj==0) {
		 muDecayLayIDFill = inPoints.xyzId.front();
		 decayPosInLayFill[0] = inPoints.xyzpos.front().X()/strpwidth;
		 decayPosInLayFill[1] = inPoints.xyzpos.front().Y()/strpwidth;
		 UShort_t tid1 = -1;
		 for(int isjj=1;isjj<int(TimeGroup.size());isjj++) {
		   // if(int(TimeGroup[isjj].PosForTrack.size())==1) {
		   for(int nll=0;nll<int(TimeGroup[isjj].PosForTrack.size());nll++) {
		     int tm = TimeGroup[isjj].PosForTrack[nll].rawhitInfo.module;
		     int tx = TimeGroup[isjj].PosForTrack[nll].rawhitInfo.xrow;
		     int ty = TimeGroup[isjj].PosForTrack[nll].rawhitInfo.yrow;
		     int tzl = TimeGroup[isjj].PosForTrack[nll].rawhitInfo.zlay;
		     UShort_t tid = constructRPCId(tm, tx, ty, tzl);
		     if(((inPoints.xyzId[0]<inPoints.xyzId[1]) &&
			 (tid-inPoints.xyzId[0]==0 || tid-inPoints.xyzId[0]==-1)) ||
			((inPoints.xyzId[0]>inPoints.xyzId[1]) &&
			 (tid-inPoints.xyzId[0]==0 || tid-inPoints.xyzId[0]==1))
			) {
		       // if(fabs(tid-inPoints.xyzId.front())==1) {
		       // rejecting later hits in same layer
		       if(tid==tid1//  &&
			  // ePosInLayFill[0]==TimeGroup[isjj].PosForTrack[nll].RawHitPos.X() &&
			  // ePosInLayFill[1]==TimeGroup[isjj].PosForTrack[nll].RawHitPos.Y()
			  ) {continue;}
		       eLayIDFill = tid;
		       ePosInLayFill[0] = TimeGroup[isjj].PosForTrack[nll].RawHitPos.X();
		       ePosInLayFill[1] = TimeGroup[isjj].PosForTrack[nll].RawHitPos.Y();
		       TVector3 tmpExt(ePosInLayFill[0],
				       ePosInLayFill[1],0);
		       int blkx = int(tmpExt.X()*blockM/nstrip);
		       int blky = int(tmpExt.Y()*blockM/nstrip);
		       
		       // X-time
		       tmpCluster = TimeGroup[isjj].PosForTrack[nll].rawhitInfo.cluster[0][0];
		       eMultiFill[0] = int(tmpCluster.size());
		       for(int jki=0;jki<int(tmpCluster.size());jki++) {
			 if(tmpCluster[jki].strp==int(tmpExt.X())
			    && int(tmpCluster[jki].tdc[0].size())) {
			   eTimeFill[0] = (tmpCluster[jki].tdc[0][0]-
					   TimeGroup[isjj].PosForTrack[nll].rawhitInfo.tdc_ref[0]-
					   tmpExt.Y()*strpwidth*spdl-
					   blocktimeoff[tm][tx][ty][tzl][0][blkx][blky]-
					   multimeoff[tm][tx][ty][tzl][0][int(tmpCluster.size())-1]
					   // -inPoints.xyztime[2].X()
					   -xytext[2].X()
					   );
			 }
		       }
		       // Y-time
		       tmpCluster = TimeGroup[isjj].PosForTrack[nll].rawhitInfo.cluster[1][0];
		       eMultiFill[1] = int(tmpCluster.size());
		       for(int jki=0;jki<int(tmpCluster.size());jki++) {
			 if(tmpCluster[jki].strp==int(tmpExt.Y())
			    && int(tmpCluster[jki].tdc[0].size())) {
			   eTimeFill[1] = (tmpCluster[jki].tdc[0][0]-
					   TimeGroup[isjj].PosForTrack[nll].rawhitInfo.tdc_ref[0]-
					   tmpExt.X()*strpwidth*spdl-
					   blocktimeoff[tm][tx][ty][tzl][1][blkx][blky]-
					   multimeoff[tm][tx][ty][tzl][1][int(tmpCluster.size())-1]
					   // -inPoints.xyztime[2].Y()
					   -xytext[2].Y()
					   );
			 }
		       }
		       // cout<<" xytext "<<xytext[2].X()<<" "<<xytext[2].Y()<<endl;
		       // cout<<" inPoints.xyzpos[1] "<<inPoints.xyzpos[1].X()<<" "<<inPoints.xyzpos[1].Y()<<" "<<inPoints.xyzpos[1].Z()<<" inPoints.xyzpos[0] "<<inPoints.xyzpos[0].X()<<" "<<inPoints.xyzpos[0].Y()<<" "<<inPoints.xyzpos[0].Z()<<endl;
		       // TVector3 pt11(inPoints.xyzpos[1].X()/strpwidth-inPoints.xyzpos[0].X()/strpwidth,
		       // 		     inPoints.xyzpos[1].Y()/strpwidth-inPoints.xyzpos[0].Y()/strpwidth,
		       // 		     inPoints.xyzpos[1].Z()/strpwidth-inPoints.xyzpos[0].Z()/strpwidth);
		       // TVector3 pt22(inPoints.xyzpos[0].X()/strpwidth-TimeGroup[isjj].PosForTrack[nll].RawHitPos.X(),
		       // 		     inPoints.xyzpos[0].Y()/strpwidth-TimeGroup[isjj].PosForTrack[nll].RawHitPos.Y(),
		       // 		     inPoints.xyzpos[0].Z()/strpwidth-TimeGroup[isjj].PosForTrack[nll].RawHitPos.Z());
		       // TVector3 pt11 = inPoints.xyzpos[1]*(1./strpwidth)-inPoints.xyzpos[0]*(1./strpwidth);
		       TVector3 dir01 = inPoints.xyzpos[0]-inPoints.xyzpos[1];
		       muExitAngleFill = dir01.Theta();
		       TVector3 apporxDecayPt = (inPoints.xyzpos[1] +
						 dir01*(1./fabs(dir01.Z()))*(fabs(dir01.Z())+(ironThickness+airGap)*0.5+rpcZShift))*(1./strpwidth);
		       TVector3 pt11 = apporxDecayPt-inPoints.xyzpos[0]*(1./strpwidth);
		       TVector3 pt22 = TimeGroup[isjj].PosForTrack[nll].RawHitPos - apporxDecayPt;
		       // cout<<" pt11 "<<pt11.X()<<" "<<pt11.Y()<<" "<<pt11.Z()<<" pt22 "<<pt22.X()<<" "<<pt22.Y()<<" "<<pt22.Z()<<endl;
		       // cout<<" pt11 "<<pt11.X()<<" "<<pt11.Y()<<" "<<pt11.Z()<<" apporxDecayPt "<<apporxDecayPt.X()<<" "<<apporxDecayPt.Y()<<" "<<apporxDecayPt.Z()<<endl;
		       eExitAngleFill = pt22.Angle(pt11);
		       ironInteractionFill = apporxDecayPt.Mag()*(1./fabs(apporxDecayPt.Z()))*(ironThickness*0.5/strpwidth);
		       // cout<<" "<<iev<<" isForward "<<isForwardFill<<" decayLay "<<muDecayLayIDFill<<" pos "<<int(decayPosInLayFill[0])<<" "<<int(decayPosInLayFill[1])<<" eLay "<<eLayIDFill<<" ePos "<<int(ePosInLayFill[0])<<" "<<int(ePosInLayFill[1])<<" angle "<<eExitAngleFill*180./TMath::Pi()<<" iron "<<ironInteractionFill<<" time "<<eTimeFill[0]<<" "<<eTimeFill[1]<<endl;
		     } // for(int nll=0;nll<int(TimeGroup[isjj].PosForTrack.size());nll) {
		       // } // if(fabs(tid-inPoints.xyzID.back())==1) {
		     tid1 = tid; // rejecting later hits in same layer
		   } // if(int(TimeGroup[isjj].PosForTrack.size())==1) {
		 } // for(int isjj=0;isjj<int(TimeGroup.size());isjj++) {
	       }   // if(isj==0) {
	       if(ePosInLayFill[0]==-1) {continue;}
#endif	// #ifdef isLifetime
	       
#ifndef isSimData
	       recosepFill = eventTime.AsDouble() - recotimeFill;
	       recotimeFill = eventTime.AsDouble();
#endif	// #ifndef isSimData
	       
	       momTree->Fill();
	       
#endif	//#ifndef isCorrection
	
	
	     }	// for(int isj=0;isj<int(TimeGroup.size());isj++) {
      
      
#endif	// #ifdef isAnalysis
      
      
      
      
      
      
      
      
      
      
    } // for(Long64_t iev=nentrymn;iev<nentry;iev++) {
    
  } // if(!fileIn->IsZombie()) {

  fileIn->Close();
  
  f1->cd();
  TriggerTree->Write();
#ifdef isSimData
  SimTree->Write();
#endif	// #ifdef isSimData
#ifdef isAnalysis
  outtreeh->Write();
  momTree->Write();
  for(int nj=0;nj<nside;nj++) {
    xybeta[nj]->Write();
  } // for(int nj=0;nj<nside;nj++) {
  for(int nj=0;nj<nside;nj++) {
    hxychi2[nj]->Write();
  } // for(int nj=0;nj<nside;nj++) {
  for(int nj=0;nj<nside;nj++) {
    circleS[nj]->Write();
  } // for(int nj=0;nj<nside;nj++) {
#endif	// #ifdef isAnalysis
  f1->Purge();
  f1->Close();
  
#ifdef isCorrection

  cout << "corr file writting" << endl;
  f2->cd();

#ifndef isSimData

  // for(int nj=0;nj<nside;nj++) {
  //   for(int nl=0;nl<nlayer;nl++) {
  //     for(int ns=0;ns<nstrip;ns++) {
  // 	rawtdcStrp[nj][nl][ns]->Write();
  //     }}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  rawtdcStrpProf[nm][nx][ny][nl]->Write();
	}}}}
  
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  layPosExtErr[nm][nx][ny][nl]->Write();
	}}}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nmxhits;nh++) {
	      layPosR[nm][nx][ny][nl][nj][nh]->Write();
	    }}}}}}
#ifdef isIter
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  for(int nj=0;nj<nside;nj++) {
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		blockPosR[nm][nx][ny][nl][nj][nsx][nsy]->Write();
	      }}
	  }}}}}
#endif	// #ifdef isIter
  
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  layTimeExtErr[nm][nx][ny][nl]->Write();
	}}}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  for(int nj=0;nj<nside;nj++) {
	  for(int nh=0;nh<nmxhits;nh++) {
	  layTimeR[nm][nx][ny][nl][nj][nh]->Write();
	}}}}}}
#ifdef isIter
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
#ifdef isIter
	  if(constructRPCId(itermodule,iterxrow,iteryrow,iterlayer)!=constructRPCId(nm,nx,ny,nl)) {continue;}
#endif	// #ifdef isIter
	  for(int nj=0;nj<nside;nj++) {
	    for(int nsx=0;nsx<blockM;nsx++) {
	      for(int nsy=0;nsy<blockM;nsy++) {
		blockTimeR[nm][nx][ny][nl][nj][nsx][nsy]->Write();
	      }}}}}}}
#endif	// #ifdef isIter

#ifndef isIter
#ifdef isStrpMulti
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
  	for(int nl=0;nl<nlayer;nl++) {
  	  for(int nj=0;nj<nside;nj++) {
	    strpCorrel2D[nm][nx][ny][nl][nj]->Write();
	    strpCorrel3D[nm][nx][ny][nl][nj]->Write();
	  }}}}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
  	for(int nl=0;nl<nlayer;nl++) {
  	  for(int nj=0;nj<nside;nj++) {
	    laymul[nm][nx][ny][nl][nj]->Write();
	  }}}}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
  	for(int nl=0;nl<nlayer;nl++) {
  	  for(int nj=0;nj<nside;nj++) {
  	    for(int nsx=0;nsx<blockM;nsx++) {
  	      for(int nsy=0;nsy<blockM;nsy++) {
		layblockmul[nm][nx][ny][nl][nj][nsx][nsy]->Write();
	      }}}}}}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
  	for(int nl=0;nl<nlayer;nl++) {
  	  for(int nj=0;nj<nside;nj++) {
	    layblocktotmul[nm][nx][ny][nl][nj]->Write();
	  }}}}}
  // for(int nm=0;nm<nmodule;nm++) {
  //   for(int nx=0;nx<nxrow;nx++) {
  //     for(int ny=0;ny<nyrow;ny++) {
  // 	for(int nl=0;nl<nlayer;nl++) {
  // 	  for(int nj=0;nj<nside;nj++) {
  // 	    for(int nsx=0;nsx<nstrip;nsx++) {
  // 	      for(int nsy=0;nsy<nstrip;nsy++) {
  // 		layblocktotmul[nm][nx][ny][nl][nj][nsx][nsy]->Write();
  // 	      }}}}}}}
  // for(int nm=0;nm<nmodule;nm++) {
  //   for(int nx=0;nx<nxrow;nx++) {
  //     for(int ny=0;ny<nyrow;ny++) {
  // 	for(int nl=0;nl<nlayer;nl++) {
  // 	  for(int nj=0;nj<nside;nj++) {
  // 	    for(int nsx=0;nsx<blockM;nsx++) {
  // 	      for(int nsy=0;nsy<blockM;nsy++) {
  // 		layblocktotmul[nm][nx][ny][nl][nj][nsx][nsy]->Write();
  // 	      }}}}}}}
  // for(int nm=0;nm<nmodule;nm++) {
  //   for(int nx=0;nx<nxrow;nx++) {
  //     for(int ny=0;ny<nyrow;ny++) {
  // 	for(int nl=0;nl<nlayer;nl++) {
  // 	  for(int nj=0;nj<nside;nj++) {
  // 	    for(int nsx=0;nsx<blockM;nsx++) {
  // 	      for(int nsy=0;nsy<blockM;nsy++) {
  // 		noiseStrpMul[nm][nx][ny][nl][nj][nsx][nsy]->Write();
  // 	      }}}}}}}
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
  	for(int nl=0;nl<nlayer;nl++) {
  	  for(int nj=0;nj<nside;nj++) {
	    noiseStrpMul2D[nm][nx][ny][nl][nj]->Write();
	  }}}}}
#endif	      // #ifdef isStrpMulti

#endif	// #ifndef isIter
#endif	// #ifndef isSimData
  
#ifndef isIter
  
  h_ndfout->Write();
  for(int nm=0;nm<nmodule;nm++) {
    for(int nx=0;nx<nxrow;nx++) {
      for(int ny=0;ny<nyrow;ny++) {
	for(int nl=0;nl<nlayer;nl++) {
	  inefficiency_tot[nm][nx][ny][nl]->Write();
	  // inefficiency_tot30[nm][nx][ny][nl]->Write();
	  // inefficiency_tot3[nm][nx][ny][nl]->Write();
	  // inefficiency_cor[nm][nx][ny][nl]->Divide(inefficiency_tot[nm][nx][ny][nl]);
	  inefficiency_cor[nm][nx][ny][nl]->Write();
	  inefficiency_cor_str[nm][nx][ny][nl]->Write();
	  inefficiency_cor_strS[nm][nx][ny][nl]->Write();
	  for(int nj=0;nj<nside;nj++) {
	    for(int nh=0;nh<nstrip;nh++) {
	      laymul_2Dstrp[nm][nx][ny][nl][nj][nh]->Write();
	    } // for(int nh=0;nh<nmxhits;nh++) {
	    // for(int nh=0;nh<nmxhits;nh++) {
	    //   laymul_2D[nm][nx][ny][nl][nj][nh]->Write();} // for(int nh=0;nh<nmxhits;nh++) {
	    // laymul_2Dm1[nm][nx][ny][nl][nj]->Write();
	    laymul_2DmTotFine[nm][nx][ny][nl][nj]->Write();
	    laymul_2DuTotFine[nm][nx][ny][nl][nj]->Write();
	    laymul_2DcTotFine[nm][nx][ny][nl][nj]->Write();
#ifdef is2dpos
	    laymul_2DmTot[nm][nx][ny][nl][nj]->Write();
	    laymul_2DuTot[nm][nx][ny][nl][nj]->Write();
	    laymul_2DcTot[nm][nx][ny][nl][nj]->Write();
	    laymul_2Dpos[nm][nx][ny][nl][nj]->Write();
#endif	// #ifdef is2dpos
	    inefficiency_tot1[nm][nx][ny][nl][nj]->Write();
	    // triggereffi_evt[nm][nx][ny][nl][nj]->Divide(inefficiency_tot[nm][nx][ny][nl]);
	    triggereffi_evt[nm][nx][ny][nl][nj]->Write();
	    inefficiency_unc[nm][nx][ny][nl][nj]->Write();
	    inefficiency_unc_tot[nm][nx][ny][nl][nj]->Write();
	    // inefficiency_unc_strS_tot[nm][nx][ny][nl][nj]->Write();
	    inefficiency_unc_str[nm][nx][ny][nl][nj]->Write();
	    inefficiency_unc_str_tot[nm][nx][ny][nl][nj]->Write();
	  }}}}}
#endif	// #ifndef isIter
  
  cout << "corr file purging" << endl;
  
  // f2->Purge();
  f2->Close();
  
  cout << "corr file wrote" << endl;

#endif	// #ifdef isCorrection
  
  cout << "end of program" << endl;
  
  return 0;
  
  cout << "after return" << endl;
}; // main

