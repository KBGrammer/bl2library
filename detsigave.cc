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
const int nhists2d = 4;
int histsize = 2.1/(.0025);

const int nmax = 1821;

int countlimit;

double xp[nmax];
double yp[nmax];

int curr_i;

double detsignalave(int run, TH1D *hasyms[], int series);
double squ(double a);
double inv(double a);

double func(const double * params);
double fcn(const double *params);

//#ifndef __CINT__
int main(int argc, char** argv)
{
  //TauData t(150, 0, 29);

  //  return 0;
  int r0,rf,stepsize, series;
  if(argc < 6) {
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
    histsize = 2.1/(atof(argv[3]));
    stepsize = atoi(argv[4]);
    series = atoi(argv[6]);
    if(argc == 8) {
      countlimit = atoi(argv[7]);
    }
    else {
      countlimit = 1e6;
    }
  }
  else {
    rf = r0;
    stepsize = 10;
  }


  TH1D *hasyms[nhists];
  hasyms[0] = new TH1D("gage_spectrum","Gage Spectrum",histsize*2,0,3.75);
  hasyms[1] = new TH1D("camac_spectrum","Camac Spectrum",histsize*4,0,8.4);
  hasyms[2] = new TH1D("ctime_spectrum","Camac Timing Spectrum",8192*5.0/8.0,0,8192*0.04*5.0/8.0-0.04);
  hasyms[3] = new TH1D("gtime_spectrum","Gage Timing Spectrum",2048,0,2047* 1.0/1.0e1);
  
  int counter = 0;

  for(int j = r0; j <= rf; j+=stepsize) {
      cout << j << endl;
      int temp =  detsignalave(j, hasyms, series);
      cout << temp << endl;
      counter += temp;
  }

  
  TFile rootfile(argv[5], "RECREATE");


  hasyms[0]->SetLineColor(kRed);
  hasyms[1]->SetLineColor(kBlue);
  hasyms[1]->SetLineColor(kBlack);
  

  for(int i = 0; i < nhists; i++) {
    hasyms[i]->Write();
  }
  TVectorD v(1);
  v[0] = counter;
  v.Write("v");
  for(int i = 0; i < nhists; i++) {
    delete hasyms[i];
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

double detsignalave(int run, TH1D *hasyms[], int series)
{
  TH1D *hasyms_local[nhists];
  hasyms_local[0] = new TH1D("","",histsize*2,0,3.75);
  hasyms_local[1] = new TH1D("","",histsize*4,0,8.4);
  hasyms_local[2] = new TH1D("","",8192*5.0/8.0,0,8192*0.04*5.0/8.0-0.04);
  hasyms_local[3] = new TH1D("","",2048,0,2047* 1.0/1.0e1);
  TauData t(run, series);
  int llength = t.c_tree->GetEntries();
  int glength = t.g_tree->GetEntries();
  int hlength = t.h_tree->GetEntries();
  int mlength = t.m_tree->GetEntries();
  cout << "len " << glength << ' ' << llength << ' ' << hlength << endl;
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


    //t.h_tree->Draw("adch:hb","");

  
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

    int triggercount = 0;

    int cc = 0;
    int gc = 0;

    cout << run << " header ";
    t.m_tree->GetEntry(0);
    for(int n= 0 ; n < NMONITOR; n++) {
      cout << marray[n] << ' ';
    }
    for(int n= 0 ; n < NMON_INTS; n++) {
      cout << aarray[n] << ' ';
    }
    cout << tsarray[0];
    cout << endl;
    int adcsum = 0;
    t.h_tree->GetEntry(0);
    for(int n= 0 ; n < HBINS; n++) {
      //t.h_tree->GetEntry(n);
      short val = adcarray[n];
      adcsum += val;
      for(int i = 0; i < val; i++) {
	hasyms_local[1]->Fill(n*8.0/HBINS-0.03);
	//cout << n << ' ' << val << endl;
      }
      val = tdcarray[n];
      for(int i = 0; i < val; i++) {
	hasyms_local[2]->Fill(n*0.04);
	//cout << n << ' ' << val << endl;
      }
    }
    cout << adcsum << endl;

    cc = adcsum;

    // times:
    // det               = 10000
    // door down/ramp on = 10023
    // mirror down       = 10128
    // close/det off     = 10160

 
    for(int n = 0; n < glength; n++) {
      t.g_tree->GetEntry(n);
      double amigos[100];
      double amigoscount[100];
      int index_val = 0;
      for(int i = 0; i < 100; i++) {
	amigos[i] = 0;
	amigoscount[i] = 0;
      }
      int delta = 20;
      int bgcount = 0;
      double sum = 0;
      double max = 0;
      double b_prev = 0;
      bool triggered = false;
      for(int i = 0; i < TBINS; i++) {
	//float val = (-16.0 - gvarray[i])*4.0 / (-32768.0 * 2.0)+1.75;
	float val = gvarray[i]+1.75;
	float val2 = gbarray[i];

	if(val2 == 2 && b_prev < 2 && !triggered) {
	  hasyms_local[3]->Fill((i-256) * 1.0/1.0e1);
	  triggered = true;
	}
	if(val2 < 2 && triggered) {
	  //triggered = false;
	}
	b_prev = val2;
	//short val2 = garray[i];

	//cout << 1000 << ' ' << i << ' ' << val << endl;

	if(val > max) {
	  max = val;
	  index_val = i;
	}

	if(i > 0 && i < 150) {
	  sum += val;
	  int index = bgcount;
	  bgcount ++;
	}
	int tb_door = 22/204.8;
	int tb_mirr = 98/204.8;
	int tb_rest = 128/204.8;
	if(i/delta < 100) {
	  amigos[i/delta] += val;
	  amigoscount[i/delta] ++;
	}

      }
      hasyms_local[0]->Fill(max - (sum/bgcount));
      cout << series << ' ' << run << ' ' << n;
      for(int i = 0; i < 100; i++) {
	double mult = 0;
	if(index_val/delta == i) {
	  mult = 1.0;
	}
	amigos[i] /= amigoscount[i];
	cout << ' ' << amigos[i]*mult;
      }
      cout << endl;
      gc++;
    }


    cout << "camac counts " << cc << endl
	 << "gage counts  " << gc << endl;
  

    for(int i = 0; i < nhists; i++) {

      if((((cc*1.5 > gc && cc*0.5 < gc) || (cc + 100 > gc && cc - 100 < gc))) && 
	 (gc < countlimit && cc < countlimit)) {
	hasyms[i]->Add(hasyms_local[i]);
      }
      else {
	cout <<"bad" << endl;
      }
    }
    for(int i = 0; i < nhists; i++) {
      delete hasyms_local[i];
    }

    //delete gvbranch, gbranch, lbranch, adchbranch, tdchbranch, mbranch;
    //  delete t;

    if(((cc*1.5 > gc && cc*0.5 < gc) || (cc + 100 > gc && cc - 100 < gc)) && 
       (gc < countlimit && cc < countlimit)) {
      return 1;
    }
    return 0;
  }
  else {
    for(int i = 0; i < nhists; i++) {
      delete hasyms_local[i];
    }

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



