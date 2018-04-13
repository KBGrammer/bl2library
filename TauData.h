// $Id:$   -*-c++-*-

#ifndef TAUDATA_H
#define TAUDATA_H

#include "TTree.h"
#include "taun_vars.h"

class TauData {

public:
  const char* fDataFileName[7]; //!
  int fNruns; //!
  int *fRunList; //!
  int *fSerList; //!
  //int *fRunIndex;
  //int fMaxRunEvents;


  char* MCAFileName;


  TTree *g_tree;
  TTree *n_tree;
  TTree *c_tree;
  TTree *h_tree;
  TTree *m_tree;

  TauData(int run=0, int series = 0);
  TauData(int nruns, int run, int series);
  TauData(int nruns, int run, int series, const char * mcafile);
  TauData(const char * mcafile);
  //TauData(const char* path);
  //TauData(int nruns, int *runlist, const char* path=0, int maxrp=0);
  ~TauData();
  void Clear();
  void ReLoad(int run=0, int series = 0);
  void ReLoad(int nruns, int run, int series);
  void ReLoad(int nruns, int run, int series, const char * mcafile);
  void ReLoad(const char * mcafile);
  void Init(const char* path=0);
  //void Init3(const char* path=0);
  char* NewString(const char* init);

  bool GoodRun();

  void GetFileNames(const char * filename[], const char * path, int r, int s);

  //ClassDef(TauData,1);
  
private:
  
  //short *gadc;
  void get_derivative(float values[], float deriv[], int smooth);
  void get_smoothed(float smooth[], float tsmooth[], int smoothness);

  float convert_to_volts(short g_signal);

  bool debug;
  int run_mode;
  int calib;

};

//void libTAUn();

#endif
