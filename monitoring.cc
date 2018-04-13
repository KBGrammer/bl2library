// detsigave - Kyle Grammer
// Returns histograms of adc and gage card pulse height spectra
// and saves them to the specified file.
// Parameters - run_i, run_f, histogram size in mV, stepsize,
//              file name for storage, limit for number of triggers
//              per run.
//
// ./detsigave 4038 4040 0.005 1 vac_25kv_20140925 500
//
//


#include "TauData.h"
#include "taun_vars.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TVectorD.h>
#include <Math/Minimizer.h>
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

using namespace std;

const int nhists = 4;
int histsize = 2.1/(.0025);

const int nmax = 1821;

int countlimit;

double xp[nmax];
double yp[nmax];

int curr_i;

double detsignalave(int run, int series);
double squ(double a);
double inv(double a);

double func(const double * params);
double fcn(const double *params);

//#ifndef __CINT__
int main(int argc, char** argv)
{
  int r0,rf,stepsize, series;
  if(argc < 3) {
    cout << "Not enough parameters!" << endl;
    cout << "Parameters - run_i, run_f, histogram size in mV, stepsize,"
	 << endl 
	 << " file name for storage, limit for number of triggers"
	 << " per run." << endl
	 << "./detsigave 4038 4040 0.005 1 vac_25kv_20140925 500" 
	 << endl << "Last parameter defaults to 1e6 if left blank."
	 << endl;
    return 0;
  }
  r0 = atoi(argv[1]);
  if(argc > 2) {
    rf = atoi(argv[2]);
    series = atoi(argv[3]);
    stepsize = atoi(argv[4]);
    
    // histsize = 2.1/(atof(argv[3]));
    // stepsize = atoi(argv[4]);
    // if(argc == 7) {
    //   countlimit = atoi(argv[6]);
    // }
    // else {
    //   countlimit = 1e6;
    // }
  }
  else {
    rf = r0;
    stepsize = 10;
  }

  
  int counter = 0;

  cout << "#run cur(uA) volt(kV) T4K(K) T77K(K) TTrap(K) P1 P2 L1 L2 A1 A2 A3 A4 A5 A6" << endl;

  for(int j = r0; j <= rf; j+=stepsize) {
    //cout << j << endl;
    int temp =  detsignalave(j, series);
      //cout << temp << endl;
      counter += temp;
  }


  cout << counter << endl;
  
  return(0);
}

double squ(double a)
{
  return a*a;
}

double inv(double a)
{
  return 1/a;
}

double detsignalave(int run, int series)
{
  TauData t(run, series);
  int llength = t.c_tree->GetEntries();
  int glength = t.g_tree->GetEntries();
  int hlength = t.h_tree->GetEntries();
  int mlength = t.m_tree->GetEntries();
  //cout << "len " << glength << ' ' << llength << ' ' << hlength << endl;
  if (t.GoodRun()) {

    int d0 = 0;

    int x = 0;
    int nomisscounter;

    double marray[NMONITOR];
    int aarray[NMONITOR];
    short adcarray[HBINS];
    short tdcarray[HBINS];
    short larray[2];
    float gvarray[TBINS];
    float gbarray[TBINS];
    short garray[TBINS];
    int tsarray[1];
    double lhearray[1];
    double ln2array[1];

    TBranch *gvbranch  = t.g_tree->GetBranch("g_v");
    gvbranch->SetAddress(gvarray);
    TBranch *gbbranch  = t.g_tree->GetBranch("g_b");
    gbbranch->SetAddress(gbarray);
    TBranch *gbranch  = t.g_tree->GetBranch("g");
    gbranch->SetAddress(garray);
    TBranch *lbranch  = t.c_tree->GetBranch("list");
    lbranch->SetAddress(larray);
    TBranch *adchbranch  = t.h_tree->GetBranch("adch");
    adchbranch->SetAddress(adcarray);
    TBranch *tdchbranch  = t.h_tree->GetBranch("tdch");
    tdchbranch->SetAddress(tdcarray);
    TBranch *mbranch  = t.m_tree->GetBranch("m");
    mbranch->SetAddress(marray);
    TBranch *alphabranch  = t.m_tree->GetBranch("alpha");
    alphabranch->SetAddress(aarray);
    TBranch *timebranch  = t.m_tree->GetBranch("ts");
    timebranch->SetAddress(tsarray);
    TBranch *ln2branch  = t.m_tree->GetBranch("ln2");
    ln2branch->SetAddress(ln2array);
    TBranch *lhebranch  = t.m_tree->GetBranch("lhe");
    lhebranch->SetAddress(lhearray);

    int triggercount = 0;

    int cc = 0;
    int gc = 0;

    cout << run << " " << series << ' ';
    t.m_tree->GetEntry(0);
    cout << tsarray[0] << ' ';
    for(int n= 0 ; n < NMONITOR; n++) {
      if(n != 7 && n != 8) {
	cout << marray[n] << ' ';
      }
      else {
	cout << ln2array[0] << ' ' << lhearray[0] << ' ';
	n++;
      }
    }
    for(int n= 0 ; n < NMON_INTS; n++) {
      cout << aarray[n] << ' ';
    }
    cout << endl;


    //delete gvbranch, gbranch, lbranch, adchbranch, tdchbranch, mbranch;
    //  delete t;

    if(((cc*1.5 > gc && cc*0.5 < gc) || (cc + 100 > gc && cc - 100 < gc)) && 
       (gc < countlimit && cc < countlimit)) {
      return 1;
    }
    return 0;
  }
  else {
    return 0;
  }
}

double func(const double * params)
{
  double value=params[0] + params[1] * xp[curr_i] + params[2] * xp[curr_i] * xp[curr_i];
  return value;
}

//______________________________________________________________________________
double fcn(const double *params)
{
   const int nbins = nmax;
   int i;

//calculate chisquare
   double chisq = 0;
   double delta;
   for (i = 0; i < nbins; i++) {
     curr_i = i;
     delta  = (yp[i]-func(params))/1.0;
     chisq += delta*delta;
   }
   return chisq;
}



