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
#include <TNtuple.h>
#include <TFile.h>
#include <TVectorD.h>
#include <Math/Minimizer.h>
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"

using namespace std;

const int nhists = 19;
const int nhists2d = 4;
int histsize = 2.1/(.0025);

const int nmax = 1821;
const int mode = 2;

int countlimit;

double xp[nmax];
double yp[nmax];

int curr_i;

double detsignalave(int run, TH1D *hasyms[], int series, TH2D *hasyms2d[]);//, TNtuple *ntuple);
double squ(double a);
double inv(double a);


void get_derivative(float values[], double deriv[], int smooth);
void get_derivative(double values[], double deriv[], int smooth);
void get_smoothed(double deriv[], double tderiv[], int smooth);
void get_smoothed(float deriv[], double tderiv[], int smooth);

double func(const double * params);
double fcn(const double *params);

double checkhighvolts(int run, int series);

//void TopDownMergeSort(A[], B[], n);
//void TopDownSplitMerge(A[], iBegin, iEnd, B[]);
//void TopDownMerge(A[], iBegin, iMiddle, iEnd, B[]);
//void CopyArray(B[], iBegin, iEnd, A[]);
//void InvertArray(B[], A[]);


//#ifndef __CINT__
int main(int argc, char** argv)
{
  //TauData t(150, 0, 29);

  //  return 0;
  int r0,rf,stepsize, series;
  /*
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
    }*/


  TH1D *hasyms[nhists];
  hasyms[0] = new TH1D("gage_spectrum","Gage Spectrum",1000,4,75);
  hasyms[1] = new TH1D("gage_spectrum_cut","Gage Spectrum Cut at High Voltage",1000,4,75);
  hasyms[2] = new TH1D("gtime_spectrum_cut","Gage Timing Spectrum Cut at High Voltage",2048,-25.6,2047*1.0/1.0e1-25.6);
  hasyms[3] = new TH1D("gtime_spectrum","Gage Timing Spectrum",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms[4] = new TH1D("gage_spectrum_cut_time","Gage Spectrum Cut around time peak",1000,4,75);
  hasyms[5] = new TH1D("bg_hist", "Gage Trace Pedestal" ,2000, -1, 1);
  hasyms[6] = new TH1D("gtime_high","Gage Timing E above 75keV",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms[7] = new TH1D("gspec_high","Gage Spectrum E above 75keV",1000,75, 3.8*60/1.6);
  hasyms[8] = new TH1D("gspec_low","Gage Spectrum E below 4keV",1000,-0.25*60/1.6, 4);
  hasyms[9] = new TH1D("gtime_low","Gage Timing E below 4keV",2048,-25.6,2047* 1.0/1.0e1-25.6);

  hasyms[10] = new TH1D("gspec_rank1","Gage Spectrum 1st Highest",1000,4,3.8*60/1.6);
  hasyms[11] = new TH1D("gtime_rank1","Gage Timing 1st Highest",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms[12] = new TH1D("gspec_rank2","Gage Spectrum 2nd Highest",1000,4,3.8*60/1.6);
  hasyms[13] = new TH1D("gtime_rank2","Gage Timing 2nd Highest",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms[14] = new TH1D("gspec_rank3","Gage Spectrum 3rd Highest",1000,4,3.8*60/1.6);
  hasyms[15] = new TH1D("gtime_rank3","Gage Timing 3rd Highest",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms[16] = new TH1D("gspec_volts","Gage Spectrum Volts",1000, 0, 4);
  hasyms[17] = new TH1D("gtime_all","Gage Timing All",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms[18] = new TH1D("ni_time","NI Timing All",200, 0, 200);
  
  TH2D *hasyms2d[nhists2d];
  hasyms2d[0] = new TH2D("gage_spectrum_series","Gage Spectrum",1000,4,75, 100, 103, 202);
  hasyms2d[1] = new TH2D("gtime_spectrum_series","Gage Timing Spectrum",2048,0,2047* 1.0/1.0e1, 100, 103, 202);
  hasyms2d[2] = new TH2D("gage_spectrum_series_cut","Gage Spectrum Cut at High Voltage",1000,4,75, 100, 103, 202);
  hasyms2d[3] = new TH2D("gtime_spectrum_series_cut","Gage Timing Spectrum Cut at High Voltage",2048,0,2047* 1.0/1.0e1, 100, 103, 202);
  int counter = 0;

  //TNtuple *ntuple = new TNtuple("taudata", "Lifetime Data", "series:run:Entry:bits:bg:gv:En:t:kV:pa:pp:trig", 12*32);



  std::ifstream fin;
  stepsize = 1;
  fin.open(argv[1]);
  int rcount = 0;
  do {
    fin >> series >> r0 >> rf;
    if (series != -1) {
      for(int j = r0; j <= rf; j+=stepsize) {
	if(j % 25 == 0) {
	  cout << "Currently on " << series << ' ' << j << endl;
	}
	int temp =  detsignalave(j, hasyms, series, hasyms2d);//, ntuple);
	//double temp = checkhighvolts(j, series);
	//cout << temp << endl;
	counter += temp;
	rcount ++;
      }
    }
  }while(series != -1);
  //cout << rcount << endl;
  fin.close();
  
  TFile rootfile(argv[2], "RECREATE");


  hasyms[0]->SetLineColor(kRed);
  hasyms[1]->SetLineColor(kBlue);
  hasyms[1]->SetLineColor(kBlack);
  

  for(int i = 0; i < nhists; i++) {
    hasyms[i]->Write();
  }
  for(int i = 0; i < nhists2d; i++) {
    hasyms2d[i]->Write();
  }
  TVectorD v(1);
  v[0] = counter;
  v.Write("v");

  //ntuple->Write();
  for(int i = 0; i < nhists; i++) {
    delete hasyms[i];
  }
  for(int i = 0; i < nhists2d; i++) {
    delete hasyms2d[i];
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

double detsignalave(int run, TH1D *hasyms[], int series, TH2D *hasyms2d[])//, TNtuple *ntuple)
{
  TH1D *hasyms_local[nhists];
  hasyms_local[0] = new TH1D("","",1000,4,75);
  hasyms_local[1] = new TH1D("","",1000,4,75);
  hasyms_local[2] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[3] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[4] = new TH1D("", "" ,1000,4,75);
  hasyms_local[5] = new TH1D("", "" ,2000, -1, 1);
  hasyms_local[6] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[7] = new TH1D("","",1000,75, 3.8*60/1.6);
  hasyms_local[8] = new TH1D("","",1000,-0.25*60/1.6, 4);
  hasyms_local[9] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[10] = new TH1D("","",1000,4,3.8*60/1.6);
  hasyms_local[11] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[12] = new TH1D("","",1000,4,3.8*60/1.6);
  hasyms_local[13] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[14] = new TH1D("","",1000,4,3.8*60/1.6);
  hasyms_local[15] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[16] = new TH1D("","",1000, 0, 4);
  hasyms_local[17] = new TH1D("","",2048,-25.6,2047* 1.0/1.0e1-25.6);
  hasyms_local[18] = new TH1D("","",200, 0, 200);

  TH2D *hasyms2d_local[nhists2d];
  hasyms2d_local[0] = new TH2D("","",1000,4,75, 100, 103, 202);
  hasyms2d_local[1] = new TH2D("","",2048,0, 2047*1.0/1.0e1, 100, 103, 202);
  hasyms2d_local[2] = new TH2D("","",1000,4,75, 100, 103, 202);
  hasyms2d_local[3] = new TH2D("","",2048,0, 2047*1.0/1.0e1, 100, 103, 202);



  TauData t(1, run, series);
  int llength = t.c_tree->GetEntries();
  int glength = t.g_tree->GetEntries();
  int hlength = t.h_tree->GetEntries();
  int mlength = t.m_tree->GetEntries();
  int nlength = t.n_tree->GetEntries();
  //cout << "len " << glength << ' ' << llength << ' ' << hlength << endl;
  int total_hit_count = 0;
  int total_trig_count = 0;
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

    int narray[200];


    //t.h_tree->Draw("adch:hb","");

  
    TBranch *gvbranch  = t.g_tree->GetBranch("g_v");
    gvbranch->SetAddress(gvarray);
    TBranch *gbbranch  = t.g_tree->GetBranch("g_b");
    gbbranch->SetAddress(gbarray);
    TBranch *gbranch  = t.g_tree->GetBranch("g");
    gbranch->SetAddress(garray);
    //TBranch *lbranch  = t.c_tree->GetBranch("list");
    //lbranch->SetAddress(larray);
    //TBranch *adchbranch  = t.h_tree->GetBranch("adch");
    // adchbranch->SetAddress(adcarray);
    //TBranch *tdchbranch  = t.h_tree->GetBranch("tdch");
    //tdchbranch->SetAddress(tdcarray);
    TBranch *mbranch  = t.g_tree->GetBranch("m");
    mbranch->SetAddress(marray);
    TBranch *alphabranch  = t.g_tree->GetBranch("alpha");
    alphabranch->SetAddress(aarray);
    TBranch *timebranch  = t.g_tree->GetBranch("ts");
    timebranch->SetAddress(tsarray);

    TBranch * nbranch = t.n_tree->GetBranch("pcounter");
    nbranch->SetAddress(narray);

    int triggercount = 0;

    int cc = 0;
    int gc = 0;

     //  cout << run << " header ";
    t.g_tree->GetEntry(0);
    // for(int n= 0 ; n < NMONITOR; n++) {
    //   cout << marray[n] << ' ';
    // }
    cerr << series << ' ' << run << ' ' << marray[0] 
	 << ' ' << marray[1] << ' ' << marray[5] 
	 << ' ' << marray[6];

    double high_volts = marray[1];

    // for(int n= 0 ; n < NMON_INTS; n++) {
    //   cout << aarray[n] << ' ';
    // }
    ///cout << tsarray[0];
    //cout << endl;
    int adcsum = 0;
    // t.h_tree->GetEntry(0);
    // for(int n= 0 ; n < HBINS; n++) {
    //   //t.h_tree->GetEntry(n);
    //   short val = adcarray[n];
    //   adcsum += val;
    //   for(int i = 0; i < val; i++) {
    // 	hasyms_local[1]->Fill(n*8.0/HBINS-0.03);
    // 	//cout << n << ' ' << val << endl;
    //   }
    //   val = tdcarray[n];
    //   for(int i = 0; i < val; i++) {
    // 	hasyms_local[2]->Fill(n*0.04);
    // 	//cout << n << ' ' << val << endl;
    //   }
    // }
    //cout << adcsum << endl;

    cc = adcsum;

    // times:
    // det               = 10000
    // door down/ramp on = 10023
    // mirror down       = 10128
    // close/det off     = 10160

    double gv_derivative[2048];
    double gv_derivative1[2048];
    double gv_derivative2[2048];
    double gv_derivative3[2048];
    double gv_derivative4[2048];
    double gv_smoothed[2048];
    int bgcount = 0;
    double sum = 0;
    double sum2 = 0;
    for(int n = 0; n < nlength; n++) {
      t.n_tree->GetEntry(n);
      for(int i = 0; i < 200; i++) {
	if(i == 0) {
	  if(narray[i] > 0) {
	    //hasyms_local[18]->AddAt(narray[i], i);
	    hasyms_local[18]->AddBinContent(i, narray[i]);
	    //for(int j = 0; j < narray[i]; j++) {
	    //  hasyms_local[18]->Fill(i);
	    //}
	  }
	}
	else {
	  if(narray[i] - narray[i-1] > 0) {
	    //hasyms_local[18]->AddAt(narray[i], i);
	    hasyms_local[18]->AddBinContent(i, narray[i] - narray[i-1]);
	    //for(int j = 0; j < narray[i] - narray[i-1]; j++) {
	    //  hasyms_local[18]->Fill(i);
	    //}	    
	  }
	}
      }

    }
    for(int n = 0; n < glength; n++) {
      t.g_tree->GetEntry(n);
      double tsum = 0;
      double tsum2 = 0;
      int tbgcount = 0;
      bool tbadness = false;
      for(int i = 0; i < 2047; i++) {
	float val = gvarray[i]+1.75;
	if((i > 0 && i < 150) || (i > 2047 -150)) {
	  tsum += val;
	  tsum2 += val * val;
	  int index = bgcount;
	  tbgcount ++;
	  if(val > 0.5) {
	    tbadness = true;
	  }
	}
      }
      hasyms_local[5]->Fill(tsum/tbgcount);
      if(!tbadness) {
	bgcount += tbgcount;
	sum += tsum;
	sum2 += tsum2;
      }
    }
    // cerr << run << ' ' << series << ' ' << sqrt(sum2/bgcount - (sum/bgcount)*(sum/bgcount)) << ' ' << bgcount << ' '
    // 	 << sum << ' ' << sum2 << endl;
 
    for(int n = 0; n < glength; n++) {
      int hit_count = 0;
      double max_list[2048];
      double max_time_list[2048];

      for(int i = 0; i < 2048; i++) {
	max_list[i] = 0;
	max_time_list[i] = 0;
      }

      // wave_event *root;
      // root = new wave_event;
      // root->pPrev = NULL;
      // root->pNext = NULL;
      // root->peak_height = 0;
      // root->peak_time = 0;
      // wave_event *conductor;
      // conductor = root;

      t.g_tree->GetEntry(n);
      double amigos[100];
      double amigoscount[100];
      int index_val = 0;
      for(int i = 0; i < 100; i++) {
	amigos[i] = 0;
	amigoscount[i] = 0;
      }
      int delta = 20;
      double max = 0;
      double min = 0;
      double b_prev = 0;
      bool triggered = false;
      bool badness = false;
      int smooth = 3;
      get_smoothed(gvarray, gv_smoothed, 1);
      get_derivative(gv_smoothed, gv_derivative, 0);
      //get_smoothed(gv_derivative, gv_derivative1, 1);
      //get_smoothed(gv_derivative1, gv_derivative2, 1);
      //get_smoothed(gv_derivative2, gv_derivative3, 1);
      //get_smoothed(gv_derivative3, gv_derivative4, 1);
      for(int i = 0; i < 5; i++) {
	get_smoothed(gv_derivative, gv_derivative4, 1);
	for(int i = 0; i < 2048; i++) {
	  gv_derivative[i] = gv_derivative4[i];
	}
      }
      int big_count = 0;
      bool new_max = false;
      for(int i = 0; i < 2047; i++) {
	//float val = (-16.0 - gvarray[i])*4.0 / (-32768.0 * 2.0)+1.75;
	float val = gvarray[i]+1.75;
	float val2 = gbarray[i];

	//if(val2 == 2 && b_prev < 2 && !triggered) {
	//hasyms_local[3]->Fill((i-256) * 1.0/1.0e1);
	//triggered = true;
	  //}
	if(val2 < 2 && triggered) {
	  //triggered = false;
	}
	b_prev = val2;
	//short val2 = garray[i];

	//cout << 1000 << ' ' << i << ' ' << val << endl;

	if(i > 4 && i < 2044 && gv_derivative4[i-1] > 0 
	   && gv_derivative4[i+1] < 0 
	   && max_time_list[hit_count-1] != i-1) {
	  bool good_peak = true;
	  int shape_time = 10;
	  for(int qq = i - shape_time; qq < i; qq++) {
	    good_peak = good_peak && gv_derivative4[qq] > 0;
	  }
	  for(int qq = i+1; qq <= i+shape_time; qq++) {
	    good_peak = good_peak && gv_derivative4[qq] < 0;
	  }
	  if(good_peak) {
	    max_list[hit_count] = val;
	    max_time_list[hit_count] = i;
	    hit_count++;
	    //conductor = root;
	    // wave_event * curr_event;
	    // curr_event = new wave_event;
	    // curr_event->pPrev = NULL;
	    // curr_event->pNext = NULL;
	    // curr_event->peak_height = val;
	    // curr_event->peak_time = i;
	    // if(hit_count == 1){
	    //   root->peak_height = val;
	    //   root->peak_time = i;
	    //   root->pNext = new wave_event;
	    // }
	    // else {
	    //   for(int i = 0; i < hit_count; i++) {
	    // 	if(curr_event->peak_height > conductor->peak_height) {
	    // 	  curr_event->pNext = conductor;
	    // 	  curr_event->pPrev = conductor->pPrev;
	    // 	  conductor->pPrev = curr_event;
	    // 	  i = hit_count;
	    // 	}
	    // 	conductor = conductor->pNext;
	    //   }
	  }
	    // if(run == 5 && series == 103 && n == 5) {
	    //   cerr << i << ' ' << val << endl;
	    // }
	}
	  // else {
	  //   if(run == 5 && series == 103 && n == 5) {
	  //     cerr << i << ' ' << val << " bad" <<  endl;
	  //   }
	  // }
	//	if(min < val 

	int tb_door = 22/204.8;
	int tb_mirr = 98/204.8;
	int tb_rest = 128/204.8;
	if(i/delta < 100) {
	  amigos[i/delta] += val;
	  amigoscount[i/delta] ++;
	}
	// if(gv_derivative[i] > 0.02 || gv_derivative[i] < 0.02) {
	//   badness = true;
	// }
	if(val > 0.5) { // badness value for signals that are too high above pedestal value
	  big_count++;
	}
      }
      if(mode == 2) {
	double tmax_list[2047];
	int tmax_time_list[2047];
	int tsort_count = 0;
	// if(run == 5 && series == 103 && n == 5) {
	//   cout << "pre-sort dump" << endl;
	//   for(int qq = 0; qq < hit_count; qq++) {
	//     //cout << max_list[qq] << ' ' << max_time_list[qq] << endl;  
	//     tmax_list[qq] = -1;
	//     tmax_time_list[qq] = -1;
	//   }
	// }
	for(int qq = 0; qq < hit_count; qq++) {
	  double tmax = -1.0;
	  double tmaxt = 0;
	  for(int pp = 0; pp < hit_count; pp++) {
	    if(max_list[pp] > tmax && (qq == 0 || max_list[pp] < tmax_list[qq-1])) {
		 tmax = max_list[pp];
		 tmaxt = max_time_list[pp];
	    }
	  }
	  tmax_list[qq] = tmax;
	  tmax_time_list[qq] = tmaxt;
	}
	total_trig_count++;
	for(int qq = 0; qq < hit_count; qq++) {
	  total_hit_count++;
	  double mv = (max_list[qq] - (sum/bgcount))*60.0/1.6;
	  double mt = (max_time_list[qq]-256) * 1.0/1.0e1;
	  double tmv = (tmax_list[qq] - (sum/bgcount))*60.0/1.6;
	  double tmt = (tmax_time_list[qq]-256) * 1.0/1.0e1;
	  if(max_list[qq] - (sum/bgcount) > 0.08) {
	    hasyms_local[16]->Fill(max_list[qq] - (sum/bgcount));
	    hasyms_local[17]->Fill(mt);
	  }
	  //  TNtuple *ntuple = new TNtuple("taudata", "Lifetime Data", "series:run:Entry:bits:bg:gv:En:t:kV:pa:pp:trig","");
	  // ntuple->Fill(float(series), float(run), float(n), 0, (sum/bgcount), 
	  // 	       float(max_list[qq]), float(mv), 
	  // 	       float(max_time_list[qq]), 0, 0, float(qq));
	  if (tmv > 4 && tsort_count < 3) {
	    hasyms_local[10+tsort_count*2]->Fill(tmv); // main hist
	    hasyms_local[11+tsort_count*2]->Fill(tmt); // main hist
	    tsort_count++;
	    // cout << tsort_count << ' ' << tmax_list[qq] << ' ' << tmax_time_list[qq] << 
	    //   ' ' << tmv << ' ' << tmt << endl;  
	  }
	  // if(run == 5 && series == 103 && n == 5) {
	  //   cout << tmax_list[qq] << ' ' << tmax_time_list[qq] << endl;  
	  // }

	  if((mv > 4 && mv < 75)) {
	    hasyms_local[0]->Fill(mv); // main hist
	    hasyms_local[3]->Fill(mt); // main hist
	    hasyms2d_local[0]->Fill(mv, series); // main hist
	    hasyms2d_local[1]->Fill(mt, series); // main hist
	    
	    if(mv > 10 && mv < high_volts ) {
	      hasyms_local[1]->Fill(mv); // cut hist	
	      hasyms_local[2]->Fill(mt); // cut hist
	      hasyms2d_local[2]->Fill(mv, series); // cut hist	
	      hasyms2d_local[3]->Fill(mt, series); // cut hist
	      if(mt > 31.73 - 2*1.427 && mt < 31.73 + 2*1.427) {
		hasyms_local[4]->Fill(mv); // cut hist	
	      }
	    }
	  }
	  else if (mv >= 75){
	    hasyms_local[6]->Fill(mt); // main hist
	    hasyms_local[7]->Fill(mv); // main hist

	  }
	  else if(mv < 4) {
	    hasyms_local[9]->Fill(mt); // main hist
	    hasyms_local[8]->Fill(mv); // main hist

	  }
	  
	}
      }
      // if(n < 100) {
      // 	int tcount = 0;
      // 	for(int i = 0; i < 2048; i++) {
      // 	  cout << run << ' ' << series << ' ' << n << ' ' << i << ' ' << gvarray[i]+1.75 << ' ';
      // 	  if(i == max_time_list[tcount]) {
      // 	    cout << max_list[tcount] << ' ' << gv_derivative[i] <<
      // 	      ' ' << gv_derivative1[i] << ' ' << gv_derivative2[i] << ' ' <<
      // 	      ' ' << gv_derivative3[i] << ' ' << gv_derivative4[i] <<
      // 	      ' ' << gv_smoothed[i] <<  endl;
      // 	    tcount++;
      // 	  }
      // 	  else {
      // 	    cout << 0 << ' ' << gv_derivative[i] <<
      // 	      ' ' << gv_derivative1[i] << ' ' << gv_derivative2[i] << ' ' <<
      // 	      ' ' << gv_derivative3[i] << ' ' << gv_derivative4[i] <<
      // 	      ' ' << gv_smoothed[i] <<  endl;
      // 	  }

      // 	}
      // }
      // if(run == 5 && series == 103 && n == 5) {
      // 	for(int i = 0; i < 2048; i++) {
      // 	  cerr << max_time_list[i] << ' ' << max_list[i] << endl;
      // 	  if(max_time_list[i] == 0) {
      // 	    i = 2048;
      // 	  }
      // 	}
      // }
      // if(big_count > 700) {
      // 	badness = true;
      // }
      if(mode == 1) {
	if(index_val < 5 || index_val > 2042) {
	  badness = true;
	}
	if(!badness) {
	  hasyms_local[0]->Fill((max - (sum/bgcount))*60.0/1.6); // main hist

	  hasyms_local[3]->Fill((index_val-256) * 1.0/1.0e1); // main hist

	  hasyms2d_local[0]->Fill((max - (sum/bgcount))*60.0/1.6, series); // main hist

	  hasyms2d_local[1]->Fill((index_val-256) * 1.0/1.0e1, series); // main hist


	  if((max - (sum/bgcount))*60.0/1.6 > 10 && (max - (sum/bgcount))*60.0/1.6 < high_volts ) {
	    hasyms_local[1]->Fill((max - (sum/bgcount))*60.0/1.6); // cut hist	
	    hasyms_local[2]->Fill((index_val-256) * 1.0/1.0e1); // cut hist
	    hasyms2d_local[2]->Fill((max - (sum/bgcount))*60.0/1.6, series); // cut hist	
	    hasyms2d_local[3]->Fill((index_val-256) * 1.0/1.0e1, series); // cut hist
	    if((index_val-256) * 1.0/1.0e1 > 31.73 - 2*1.427 && (index_val-256) * 1.0/1.0e1 < 31.73 + 2*1.427) {
	      hasyms_local[4]->Fill((max - (sum/bgcount))*60.0/1.6); // cut hist	
	    }
	  }
	}
      }
      //	  triggered = true;
      //cout << series << ' ' << run << ' ' << n << ' ' << max;
      for(int i = 0; i < 100; i++) {
	double mult = 0;
	if(index_val/delta == i) {
	  mult = 1.0;
	}
	amigos[i] /= amigoscount[i];
	//cout << ' ' << amigos[i]*mult;
      }
      // cout << endl;
      gc++;
    }


    //cout << "camac counts " << cc << endl
    //	 << "gage counts  " << gc << endl;
  

    for(int i = 0; i < nhists; i++) {

      //if((((cc*1.5 > gc && cc*0.5 < gc) || (cc + 100 > gc && cc - 100 < gc))) && 
      //	 (gc < countlimit && cc < countlimit)) {
	hasyms[i]->Add(hasyms_local[i]);
	// }
	// else {
	//	cout <<"bad" << endl;
	//}
    }
    for(int i = 0; i < nhists; i++) {
      delete hasyms_local[i];
    }
    for(int i = 0; i < nhists2d; i++) {

      //if((((cc*1.5 > gc && cc*0.5 < gc) || (cc + 100 > gc && cc - 100 < gc))) && 
      //	 (gc < countlimit && cc < countlimit)) {
	hasyms2d[i]->Add(hasyms2d_local[i]);
	// }
	// else {
	//	cout <<"bad" << endl;
	//}
    }
    for(int i = 0; i < nhists2d; i++) {
      delete hasyms2d_local[i];
    }
      cerr << ' ' << total_hit_count << ' ' << total_trig_count << endl;

    //delete gvbranch, gbranch, lbranch, adchbranch, tdchbranch, mbranch;
    //  delete t;

    //if(((cc*1.5 > gc && cc*0.5 < gc) || (cc + 100 > gc && cc - 100 < gc)) && 
    //   (gc < countlimit && cc < countlimit)) {
      return 1;
      //}

    return 0;
  }
  else {
    for(int i = 0; i < nhists; i++) {
      delete hasyms_local[i];
    }
    for(int i = 0; i < nhists2d; i++) {
      delete hasyms2d_local[i];
    }

    return 0;

  }
}

double func(const double * params)
{
  double value=params[0] + params[1] * xp[curr_i] + params[2] * xp[curr_i] * xp[curr_i];
  return value;
}

void get_derivative(double values[], double deriv[], int smooth)
{
  double tderiv[2048];
  for(int i = 0; i < 2047; i++) {
    if (i > 3 && i < 2044) {
      deriv[i] = (values[i+4]-32.0/3.0*values[i+3]+56.0*values[i+2]-224.0*values[i-1]+224*values[i+1]-56.0*values[i-2]+32.0/3.0*values[i-3]-values[i-4])/280.0;
      tderiv[i] = deriv[i];
    }
    else {
      deriv[i] = 0;
      tderiv[i] = deriv[i];
    }
  }
  if(smooth == 0) {
    return;
  }
  return;
}

void get_derivative(float values[], double deriv[], int smooth)
{
  double temp[3048];
  for(int i = 0; i < 2048; i++) {
    temp[i] = values[i];
  }
  get_derivative(temp, deriv, smooth);
  return;
}

void get_smoothed(double deriv[], double tderiv[], int smooth)
{
  if(smooth == 1) {
    for(int i = 0; i < 2047; i++) {
      if (i > 3 && i < 2044) {
	tderiv[i] = (deriv[i-1]+deriv[i]+deriv[i+1])/3.0;
      }
      else if (i < 4){
	tderiv[i] = deriv[4];
      }
      else {
	tderiv[i] = deriv[2043];
      }
    }
  }
  if(smooth == 2) {
    for(int i = 0; i < 2047; i++) {
      if (i > 3 && i < 2044) {
	tderiv[i] = (deriv[i-2]+2*deriv[i-1]+3*deriv[i]+2*deriv[i+1]+deriv[i+2])/9.0;
      }
      else if (i < 4){
	tderiv[i] = deriv[4];
      }
      else {
	tderiv[i] = deriv[2043];
      }
    }
  }
  if(smooth == 3) {
    for(int i = 0; i < 2047; i++) {
      if (i > 3 && i < 2044) {
	tderiv[i] = (deriv[i-3]+3*deriv[i-2]+6*deriv[i-1]+7*deriv[i]+6*deriv[i+1]+3*deriv[i+2]+deriv[i+3])/17.0;
      }
      else if (i < 4){
	tderiv[i] = deriv[4];
      }
      else {
	tderiv[i] = deriv[2043];
      }
    }
  }
  return;
}

void get_smoothed(float deriv[], double tderiv[], int smooth)
{
  double temp[3048];
  for(int i = 0; i < 2048; i++) {
    temp[i] = deriv[i];
  }
  get_smoothed(temp, tderiv, smooth);
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
/*
// Array A[] has the items to sort; array B[] is a work array.
void TopDownMergeSort(A[], B[], n)
{
    TopDownSplitMerge(A, 0, n, B);
}

// iBegin is inclusive; iEnd is exclusive (A[iEnd] is not in the set).
void TopDownSplitMerge(A[], iBegin, iEnd, B[])
{
    if(iEnd - iBegin < 2)                       // if run size == 1
        return;                                 //   consider it sorted
    // recursively split runs into two halves until run size == 1,
    // then merge them and return back up the call chain
    iMiddle = (iEnd + iBegin) / 2;              // iMiddle = mid point
    TopDownSplitMerge(A, iBegin,  iMiddle, B);  // split / merge left  half
    TopDownSplitMerge(A, iMiddle,    iEnd, B);  // split / merge right half
    TopDownMerge(A, iBegin, iMiddle, iEnd, B);  // merge the two half runs
    CopyArray(B, iBegin, iEnd, A);              // copy the merged runs back to A
}

//  Left half is A[iBegin :iMiddle-1].
// Right half is A[iMiddle:iEnd-1   ].
void TopDownMerge(A[], iBegin, iMiddle, iEnd, B[])
{
    i = iBegin, j = iMiddle;
    
    // While there are elements in the left or right runs...
    for (k = iBegin; k < iEnd; k++) {
        // If left run head exists and is <= existing right run head.
        if (i < iMiddle && (j >= iEnd || A[i] <= A[j])) {
            B[k] = A[i];
            i = i + 1;
        } else {
            B[k] = A[j];
            j = j + 1;    
        }
    } 
}

void CopyArray(B[], iBegin, iEnd, A[])
{
    for(k = iBegin; k < iEnd; k++)
        A[k] = B[k];
}


void InvertArray(A[], int nmax)
{
  double t[2048];
  for(int k = 0; k < 2048; k++) {
    t[k] = A[k];
  }
  for(int k = 0; k < 2048; k++) {
    A[2047-k] = t[k];
  }
}
*/


double checkhighvolts(int run, int series)
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


  
    TBranch *mbranch  = t.m_tree->GetBranch("m");
    mbranch->SetAddress(marray);
    TBranch *alphabranch  = t.m_tree->GetBranch("alpha");
    alphabranch->SetAddress(aarray);
    TBranch *timebranch  = t.m_tree->GetBranch("ts");
    timebranch->SetAddress(tsarray);

    int triggercount = 0;

    int cc = 0;
    int gc = 0;

     //  cout << run << " header ";
    t.m_tree->GetEntry(0);
    //for(int n= 0 ; n < NMONITOR; n++) {
    //  cout << marray[n] << ' ';
    // }
    
    double high_volts = marray[1];
    cout << series << ' ' << run << ' ' << high_volts << endl;
    return high_volts;
  }
}
