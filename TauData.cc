// $Id: TauData.cc,v 1.1 2006/01/28 08:34:21\ Exp $   -*-c++-*-

// reads raw TAUn data files in TTree format

#include "TauData.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "utility.h"
#include "TTree.h"
#include "TSystem.h"
#include "taun_vars.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

//void libTAUn() {};
//using namespace std;

//ClassImp(TauData);

//______________________________________________________________________
TauData::TauData(int run, int series)
  : debug(false), fNruns(1)
{
  int nruns = fNruns; 
  run_mode = 1;
  calib = 0;
  debug = false;
  fRunList = new int [nruns];
  //fSerList = new int(series);
  fSerList = new int[1];
  fSerList[0] = series;
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
  /*
  debug = false;
  fRunList = new int(run);
  fSerList = new int(series);
  Init();
  */
}

//______________________________________________________________________
TauData::TauData(int nruns, int run, int series)
  : debug(false), fNruns(nruns)
{
  debug = false;
  run_mode = 1;
  calib = 0;
  fRunList = new int [nruns];
  fSerList = new int[1];
  fSerList[0] = series;
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
}

//______________________________________________________________________
TauData::TauData(int nruns, int run, int series, const char * mcafile)
  : debug(false), fNruns(nruns)
{
  debug = false;
  run_mode = 15;
  calib = 1;
  MCAFileName = new char;
  MCAFileName = NewString(Form("/home/lifetime/data/%s", mcafile));
  fRunList = new int [nruns];
  fSerList = new int(series);
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
}

//______________________________________________________________________
TauData::TauData(const char * mcafile)
  : debug(false)
{
  fNruns = 1;
  int nruns = 1;
  int series = 0;
  int run = 0;
  debug = false;
  run_mode = 15;
  calib = 1;
  MCAFileName = new char;
  MCAFileName = NewString(Form("/home/lifetime/data/%s", mcafile));
  fRunList = new int [nruns];
  fSerList = new int(series);
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
}



//______________________________________________________________________
TauData::~TauData()
{
  delete fRunList;
  if(debug) std::cout<<"<< Checkpoint #1 >>"<<std::endl;
  for (int i=0; i<7; ++i) {
    if(debug) std::cout<<"<< Checkpoint #2 >>"<<std::endl;
    delete[] fDataFileName[i];
  }
  delete MCAFileName;
  delete g_tree;
  //delete c_tree;
  //delete h_tree;
  delete m_tree;
  delete n_tree;
  if(debug) std::cout<<"<< Checkpoint #3 >>"<<std::endl;
  //delete fRunIndex;
  if(debug) std::cout<<"<< Checkpoint #4 >>"<<std::endl;
}

//______________________________________________________________________
void TauData::ReLoad(int run, int series)
{
  int nruns = fNruns; 
  run_mode = 1;
  calib = 0;
  debug = false;
  fRunList = new int [nruns];
  //fSerList = new int(series);
  fSerList = new int[1];
  fSerList[0] = series;
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
  /*
  debug = false;
  fRunList = new int(run);
  fSerList = new int(series);
  Init();
  */
}

//______________________________________________________________________
void TauData::ReLoad(int nruns, int run, int series)
{
  fNruns = nruns;
  debug = false;
  run_mode = 1;
  calib = 0;
  fRunList = new int [nruns];
  fSerList = new int[1];
  fSerList[0] = series;
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
}

//______________________________________________________________________
void TauData::ReLoad(int nruns, int run, int series, const char * mcafile)
{
  fNruns = nruns;
  debug = false;
  run_mode = 15;
  calib = 1;
  MCAFileName = new char;
  MCAFileName = NewString(Form("/home/lifetime/data/%s", mcafile));
  fRunList = new int [nruns];
  fSerList = new int(series);
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
}

//______________________________________________________________________
void TauData::ReLoad(const char * mcafile)
{
  fNruns = 1;
  int nruns = 1;
  int series = 0;
  int run = 0;
  debug = false;
  run_mode = 15;
  calib = 1;
  MCAFileName = new char;
  MCAFileName = NewString(Form("/home/lifetime/data/%s", mcafile));
  fRunList = new int [nruns];
  fSerList = new int(series);
  for (int r=0; r < nruns; r++) {
    fRunList[r] = run+r;
  }

  Init();
}

//______________________________________________________________________
void TauData::Clear()
{
  delete fRunList;
  if(debug) std::cout<<"<< Checkpoint #1 >>"<<std::endl;
  for (int i=0; i<7; ++i) {
    if(debug) std::cout<<"<< Checkpoint #2 >>"<<std::endl;
    delete[] fDataFileName[i];
  }
  delete MCAFileName;
  /*  TBranch *b = g_tree->GetBranch("g");
  g_tree->GetListOfBranches()->Remove(b);
  TLeaf *l = g_tree->GetLeaf("g");
  g_tree->GetListOfLeaves()->Remove(l);
  b = g_tree->GetBranch("g_v");
  g_tree->GetListOfBranches()->Remove(b);
  l = g_tree->GetLeaf("g_v");
  g_tree->GetListOfLeaves()->Remove(l);
  b = g_tree->GetBranch("dg_v");
  g_tree->GetListOfBranches()->Remove(b);
  l = g_tree->GetLeaf("dg_v");
  g_tree->GetListOfLeaves()->Remove(l);
  b = g_tree->GetBranch("g_vsmooth");
  g_tree->GetListOfBranches()->Remove(b);
  l = g_tree->GetLeaf("g_vsmooth");
  g_tree->GetListOfLeaves()->Remove(l);
  b = g_tree->GetBranch("dg_vsmooth");
  g_tree->GetListOfBranches()->Remove(b);
  l = g_tree->GetLeaf("dg_vsmooth");
  g_tree->GetListOfLeaves()->Remove(l);
  b = g_tree->GetBranch("g_b");
  g_tree->GetListOfBranches()->Remove(b);
  l = g_tree->GetLeaf("g_b");
  g_tree->GetListOfLeaves()->Remove(l);
  b = g_tree->GetBranch("t");
  g_tree->GetListOfBranches()->Remove(b);
  l = g_tree->GetLeaf("t");
  g_tree->GetListOfLeaves()->Remove(l);*/
  delete g_tree;
  //delete c_tree;
  //delete h_tree;
  delete m_tree;
  delete n_tree;
  if(debug) std::cout<<"<< Checkpoint #3 >>"<<std::endl;
  //delete fRunIndex;
  if(debug) std::cout<<"<< Checkpoint #4 >>"<<std::endl;
}

void TauData::GetFileNames(const char * filename[], const char * path, int r, int s)
{
  int run_no = r;
  int series = s;
  for (int i=0;i<7;++i) {
    if(debug)  std::cout << i << ' ';
    if (!path) {
      if(i+1 == 1 && 1 == 0) {
	filename[i] = NewString(Form("s_%04d_r_%06d_%d.dat", 
				     series, run_no, 4+1));

      }
      else {
	filename[i] = NewString(Form("s_%04d_r_%06d_%d.dat", 
				     series, run_no, i+1));
      }
      if(series == 0) {
	filename[i] = NewString(Form("run_%d_%d.dat", 
				     run_no, i+1));
      }
    }
    else {
      if(i+1 == 1 && 1 == 0) {
	filename[i] = NewString(Form("%s/s_%04d_r_%06d_%d.dat", 
				     path, series, run_no, 4+1));

      }
      else {
	filename[i] = NewString(Form("%s/s_%04d_r_%06d_%d.dat", 
				     path, series, run_no, i+1));
      }
      if(series == 0) {
	filename[i] = NewString(Form("%s/run_%d_%d.dat", 
				     path, run_no, i+1));
      }
    }

    if(debug)  std::cout << i << ' ' << filename[i] << std::endl;
  }
}

char* TauData::NewString(const char* init)
{
  char* str=new char[strlen(init)+1];
  strcpy(str,init);
  return str;
}

float TauData::convert_to_volts(short g_signal)
{
  const float toVolts = 4.0 / (-32768.0 * 2.0);
  return (-16.0 - g_signal) * toVolts + 0.0;
}

bool TauData::GoodRun()
{
  return (g_tree->GetEntries() > 0)
    && (m_tree->GetEntries() > 0);

  return (g_tree->GetEntries() > 0)
    && (h_tree->GetEntries() > 0)
    && (c_tree->GetEntries() > 0)
    && (m_tree->GetEntries() > 0);
}

void TauData::Init(const char* path)
{
  if(debug)  std::cout << "Begin Init " << fNruns<< std::endl;
  const char* ipath=path;
  for (int i=0;i<5;++i) {
    fDataFileName[i]=new char;
  }

  if(debug)  std::cout << "Declare gage tree" << std::endl;
  g_tree=new TTree("g_tree","gage card signal");
  n_tree=new TTree("n_tree","ni daq signal");
  short tBuffer [TBINS];
  short gBuffer [TBINS];
  short gBuffer2 [TBINS];
  float gvBuffer [TBINS];
  float dgvBuffer [TBINS];
  float smoothgvBuffer [TBINS];
  float dsmoothgvBuffer [TBINS];
  float gbBuffer [TBINS];
  double tsBuffer [2];
  double pBuffer [5];
  short gentry[1];

  int alphaBuffer[22];
  int pcounterBuffer[200];
  int npcounterBuffer[200];
  int ntbBuffer[200];
  int npBuffer[17];

  for(int i = 0; i < TBINS; i++) {
    tBuffer[i] = i;
  }
  g_tree->Branch("g", gBuffer, Form("g[%i]/S", TBINS));
  g_tree->Branch("g_v", gvBuffer, Form("g_v[%i]/F", TBINS));
  g_tree->Branch("dg_v", dgvBuffer, Form("dg_v[%i]/F", TBINS));
  g_tree->Branch("g_vsmooth", smoothgvBuffer, Form("g_vsmooth[%i]/F", TBINS));
  g_tree->Branch("dg_vsmooth", dsmoothgvBuffer, Form("dg_vsmooth[%i]/F", TBINS));
  g_tree->Branch("g_b", gbBuffer, Form("g_b[%i]/F", TBINS));

  g_tree->Branch("alpha", alphaBuffer, Form("alpha[%i]/I", 22));
  g_tree->Branch("pcounter", pcounterBuffer, Form("pcounter[%i]/I", 200));

  n_tree->Branch("alpha", alphaBuffer, Form("alpha[%i]/I", 22));
  n_tree->Branch("pcounter", npcounterBuffer, Form("pcounter[%i]/I", 200));
  n_tree->Branch("ntb", ntbBuffer, Form("ntb[%i]/I", 200));
  n_tree->Branch("p", npBuffer, Form("p[%i]/I", 17));

  g_tree->Branch("t", tBuffer, Form("t[%i]/S", TBINS));
  g_tree->Branch("ntb", ntbBuffer, Form("ntb[%i]/I", 200));
  g_tree->Branch("ts", tsBuffer, Form("ts[%i]/D", 2));
  g_tree->Branch("p", pBuffer, Form("p[%i]/D", 5));
  g_tree->Branch("entry", gentry, "entry/S"); //not needed
  if(debug)  std::cout << "Read file" << std::endl;

  m_tree=new TTree("m_tree","monitoring data");
  double mBuffer[NMONITOR];
  int mBuffer2[NMON_INTS];
  unsigned int mBuffer3[1];
  double mBuffer4[1];
  double mBuffer5[1];
  m_tree->Branch("m", mBuffer, Form("m[%i]/D", NMONITOR));
  m_tree->Branch("alphac", mBuffer2, Form("alphac[%i]/I", NMON_INTS));
  m_tree->Branch("ts", mBuffer3, "ts/I");
  m_tree->Branch("lhe", mBuffer4, "lhe/D");
  m_tree->Branch("ln2", mBuffer5, "ln2/D");
  g_tree->Branch("m", mBuffer, Form("m[%i]/D", NMONITOR));
  g_tree->Branch("alphac", mBuffer2, Form("alphac[%i]/I", NMON_INTS));
  g_tree->Branch("tsc", mBuffer3, "tsc/I");
  g_tree->Branch("lhe", mBuffer4, "lhe/D");
  g_tree->Branch("ln2", mBuffer5, "ln2/D");

  std::ifstream fin_camac;


  //std::ofstream fout2;
  //fout.open("test.dat");

  for (int irun=0; irun<fNruns; irun++) {
    std::ifstream fin_gage;

    //std::ofstream fout;
    std::ifstream fin_alpha;

    //std::ofstream fout6;
    std::ifstream fin_proton;

    //std::ofstream fout7;

    std::ifstream fin_monitor;
    int run_no = fRunList[irun];
    if(run_no %10 == 0) {
      std::cout << run_no << std::endl;
    }
    int series = fSerList[0];
    if(debug)  std::cout << run_no << std::endl;
    if (!path) {
      ipath="/home/lifetime/data";
    }
    if(debug)  std::cout << ipath << std::endl;
    GetFileNames(fDataFileName, ipath, run_no, series);
    fin_monitor.open(fDataFileName[3], std::ifstream::binary);
    //  short tBuffer [TBINS];
    //float gvBuffer [TBINS];

    if(debug)  std::cout << "monitoring read" << std::endl;
    if(debug)  std::cout << "Declare tree" << std::endl;

    if (fin_monitor != NULL) {
      fin_monitor.seekg(0, fin_monitor.end);
      long fileSize = fin_monitor.tellg();
      fin_monitor.seekg (0, fin_monitor.beg);
  
      if(debug)  std::cout << fileSize << std::endl;
      char * cBuf = new char [NMONITOR*8];

      // Read the file in to the buffer
      //for(int i = 0; i < fileSize / (4); i++) {
      //fin.read(cBuf, NMONITOR*8);

      if(fileSize == 100) {
	fin_monitor.read(reinterpret_cast<char*> (&mBuffer3[0]), sizeof(mBuffer3[0]));
      }
      else {
	fin_monitor.read(reinterpret_cast<char*> (&mBuffer3[0]), sizeof(mBuffer3[0]));
	//mBuffer3[0] = 0;
      }
      for(int i = 0; i < NMONITOR; i++) {
	//fin >> mBuffer[i];
	fin_monitor.read(reinterpret_cast<char*> (&mBuffer[i]), sizeof(mBuffer[i]));
	if(debug)  std::cout << mBuffer[i] << std::endl;
      }
      for(int i = 0; i < NMON_INTS; i++) {
	fin_monitor.read(reinterpret_cast<char*> (&mBuffer2[i]), sizeof(mBuffer2[i]));
	if(debug)  std::cout << mBuffer2[i] << std::endl;
      }

      mBuffer5[0] = (1-(mBuffer[7]-0.01)/0.53)*100.0;
      mBuffer4[0] = (1-(mBuffer[8]-0.0165)/0.53)*100.0;
    
      // current, voltage, t4, t77, ttrap, povc, pbore
    
      m_tree->Fill();
      //delete cBuf;
    }
  
    //fin_monitor.close();

    //int run_no = fRunList[irun];
    //if(debug) std::cout << run_no << std::endl;
    //int series = fSerList[0];
    if(debug)  std::cout << run_no << std::endl;
    if (!path) {
      ipath="/home/lifetime/data";
    }
    if(debug)  std::cout << ipath << std::endl;
    GetFileNames(fDataFileName, ipath, run_no, series);
    if(run_mode == 5) {
      fin_gage.open(fDataFileName[4], std::ifstream::binary);
      std::cout << "mode 5" << std::endl;
    }
    else if(calib == 1) {
      std::cout << "MCA mode" << std::endl;
      fin_gage.open(MCAFileName, std::ifstream::binary);
      //std::cout << "MCA mode" << std::endl;
    }
    else {
      fin_gage.open(fDataFileName[0], std::ifstream::binary);
    }
    fin_alpha.open(fDataFileName[5], std::ifstream::binary); //alpha
    fin_proton.open(fDataFileName[6], std::ifstream::binary); //proton

    //fout6.open("temp.dat");

    char * abuffer = new char [10000*4*22];
    //abuffer = new char [10000*4*22];
    int * abufferi = new int ;//[10000*22];
    char * pbuffer = new char [10000*4*200];
    //pbuffer = new char [10000*4*200];
    int * pbufferi = new int ;//[10000*200];
    bool good6 = false;
    bool good7 = false;
    //gv_tree->Branch("t", tBuffer, "t[2048]/S");
    

    for(int i = 0; i < 200; i++) {
      ntbBuffer[i] = i;
    }


    fin_alpha.seekg(0, fin_alpha.end);
    long tfileSize = fin_alpha.tellg();
    fin_alpha.seekg (0, fin_alpha.beg);
    if(fin_alpha != NULL) {
      if(debug) std::cout << "file 6 " << tfileSize << std::endl;
      fin_alpha.read(abuffer, tfileSize);
      abufferi = reinterpret_cast<int*> (abuffer);
      good6 = true;
      //for(int i = 0; i < 10000; i++) {
      //std::cout << abuffer[i] << std::endl;
      //int delta = 0;
      //if(i > 0) delta = abufferi[i] - abufferi[i-1];
      //fout6 << abufferi[i] << ' ' << delta << std::endl;
      //}
    }
    int alphasize = tfileSize;
    //fin_alpha.close();
    fin_proton.seekg(0, fin_proton.end);
    tfileSize = fin_proton.tellg();
    fin_proton.seekg (0, fin_proton.beg);
    if(fin_proton != NULL) {
      if(debug) std::cout << "file 7 " << tfileSize << std::endl;
      fin_proton.read(pbuffer, tfileSize);
      pbufferi = reinterpret_cast<int*> (pbuffer);
      good7 = true;
    }
    for(int ii = 0; ii < 10000; ii++) {
      //std::cout << ii << std::endl;
      if(good6) {
      	for(int j = 0; j < 22; j++) {
      	  if( j < 16 || alphasize > 16*10000*4) {
      	    alphaBuffer[j] = abufferi[ii + j * 10000];
      	  }
      	  else {
      	    alphaBuffer[j] = 0;
      	  }
      	  //fout6 << ii << ' ' << j << ' ' << alphaBuffer[j] << std::endl;
      	}
      }
      if(good7) {
      	for(int j = 0; j < 200; j++) {
      	  pcounterBuffer[j] = pbufferi[j + ii * 200];
      	  if(j > 0) {
      	    if(pcounterBuffer[j] - pcounterBuffer[j-1] != 0) {
      	      npBuffer[0] = j;
      	    }
      	  }
      	  else {
      	    if(pcounterBuffer[j] == 1) {
      	      npBuffer[0] = j;
      	    }
      	  }
      	}
      }
      n_tree->Fill();
    }
    //fin_proton.close();
    //fout6.close();
    gentry[0] = 0;
    const float toVolts = 4.0 / (-32768.0*2.0);
    if (fin_gage != NULL) {
      fin_gage.seekg(0, fin_gage.end);
      long fileSize = fin_gage.tellg();
      fin_gage.seekg (0, fin_gage.beg);
  
      if(debug)  std::cout << fileSize << std::endl;
      int temptbins = 2048;
      if(calib == 1) {
	temptbins = 512;
      }
      char * cBuffer = new char [TBINS*2];
      char * dBuffer1 = new char [8];
      char * dBuffer2 = new char [8];
      short * sbuf = new short ;//[TBINS];
      double * dbuf1 = new double [1];
      double * dbuf2 = new double [1];
      double * bgbuf = new double [2];
      if(debug) std::cout << fileSize << std::endl;
      // Read the file in to the buffer
      for(int i = 0; i < fileSize / ((temptbins*2+8)*2); i++) {
	//std::cout << i << std::endl;
	if(good6) {
	  for(int j = 0; j < 22; j++) {
	    alphaBuffer[j] = abufferi[i + j * 10000];
	    //fout6 << i << ' ' << j << ' ' << alphaBuffer[j] << std::endl;
	  }
	}
	if(good7) {
	  for(int j = 0; j < 200; j++) {
	    pcounterBuffer[j] = pbufferi[j + i * 200];
	  }
	}
	for(int j = 0; j < 2; j++) {
	  bgbuf[j] = 0;
	}
	for(int k = 0; k < 2; k++) {
	  if(k == 0) {
	    fin_gage.read(dBuffer1, 8);
	    dbuf1 = reinterpret_cast<double*> (dBuffer1);
	  }
	  else {
	    fin_gage.read(dBuffer2, 8);
	    dbuf2 = reinterpret_cast<double*> (dBuffer2);
	  }
	  fin_gage.read(cBuffer, temptbins*2);
	  sbuf = reinterpret_cast<short*> (cBuffer);
	  if(k == 0) {
	    double max = -4;
	    double tmax = 0;
	    for(int j = 0; j < TBINS; j++) {
	      float temp;
	      if(j < temptbins) {
		gBuffer[j] = sbuf[j];
		temp = gBuffer[j];
	      }
	      else {
		gBuffer[j] = 0;
		temp = 0;
	      }
	      if(j < 256 && calib == 0) {
		bgbuf[0] += gBuffer[j];
		bgbuf[1] += (-16.0 - temp) * toVolts + 0.0;
	      }
	      //gvBuffer[j] = convert_to_volts(gBuffer[j]);
	      //gvBuffer[j] = (-16.0 - sbuf[j]) * toVolts + 0.0;
	      gvBuffer[j] = (-16.0 - temp) * toVolts + 0.0;
	      //std::cout << j << ' ' << gvBuffer[j] << std::endl;
	      if(gvBuffer[j] > max) {
		max = gvBuffer[j];
		tmax = j;
	      }
	      //fout << temp << std::endl;
	    }
	    pBuffer[2] = max; // max volts
	    pBuffer[3] = tmax; // max arrival time
	    get_smoothed(gvBuffer, smoothgvBuffer, 15);
	    get_derivative(smoothgvBuffer, dgvBuffer, 0);
	    get_smoothed(dgvBuffer, dsmoothgvBuffer, 15);
	  }
	  else {
	    for(int j = 0; j < TBINS; j++) {
	      float temp;
	      if(j < temptbins) {
		gBuffer2[j] = sbuf[j];
		temp = gBuffer2[j];
	      }
	      else {
		gBuffer2[j] = 0;
		temp = 0;
	      }
	      //gvBuffer[j] = convert_to_volts(gBuffer[j]);
	      //gvBuffer[j] = (-16.0 - sbuf[j]) * toVolts + 0.0;
	      gbBuffer[j] = (-16.0 - temp) * toVolts * 0.2/4.0 + 0.0;
	      //if(i == 0 && j < 512) std::cout << temp << std::endl;
	    }
	  }
	}
	double tlast = tsBuffer[0]+10000.0;
	tsBuffer[0] = *dbuf1;
	tsBuffer[1] = *dbuf2;
	int tempbin = 200*int(tsBuffer[1]);
	//uncomment below to make delta T file
	//fout << tsBuffer[0] << ' ' << tsBuffer[1] << std::endl;
	//  fout << (tsBuffer[0]+10000.0) << ' ' 
	//    << (tsBuffer[0]+10000.0) / 10200.5 << ' ' 
	//   << tempbin << ' ' << tsBuffer[1] 
	//   << ' ' << tlast << ' '
	//   << (tsBuffer[0]+10000.0) - tlast << std::endl;
	
	for(int i = 0; i < 200; i++) {
	  //npcounterBuffer[i] = pbufferi[i + tempbin];
	}
	//if(debug) std::cout << tsBuffer[0] / 10200.0 << ' ' << npcounterBuffer[0] << std::endl;
	//npcounterBuffer[0] = int(tsBuffer[0] / 10200.0);

	//if(debug) std::cout << pBuffer[3] << std::endl;
	pBuffer[0] = bgbuf[0] / 256.0;
	pBuffer[1] = bgbuf[1] / 256.0;
	pBuffer[4] = run_no;
	gentry[0] = fileSize / ((TBINS*2+8)*2);
	g_tree->Fill();
	//gv_tree->Fill();
      }
      //delete[]fileBuf;
      //fin_gage.close();
      //if(debug) std::cout << "gage " << irun << std::endl;
      delete [] cBuffer;
      delete [] dBuffer1;
      delete [] dBuffer2;
    }
    if(fin_gage) fin_gage.close();
    if(fin_alpha) fin_alpha.close();
    if(fin_proton) fin_proton.close();
    if(fin_monitor) fin_monitor.close();
    //delete [] pbufferi;
    delete [] pbuffer;
    //delete [] abufferi;
    delete [] abuffer;
    //delete [] cBuffer;

  }
  //fout2.close();
  return; 
  
  c_tree=new TTree("c_tree","camac list");
  short lBuffer[2];
  short entry[1];
  c_tree->Branch("list", lBuffer, "list[2]/S");
  entry[0] = 0;
  for (int irun=0; irun<fNruns; ++irun) {
    int run_no = fRunList[irun];
    if(debug) std::cout << run_no << std::endl;
    int series = fSerList[0];
    if(debug)  std::cout << run_no << std::endl;
    if (!path) {
      ipath="/home/lifetime/data";
    }
    if(debug)  std::cout << ipath << std::endl;
    GetFileNames(fDataFileName, ipath, run_no, series);
    fin_camac.open(fDataFileName[1], std::ifstream::binary);

    if(debug)  std::cout << "Camac list read" << std::endl;
    if(debug)  std::cout << "Declare camac tree" << std::endl;
    short entry[1];
    entry[0] = 0;

    if (fin_camac != NULL) {
      fin_camac.seekg(0, fin_camac.end);
      long fileSize = fin_camac.tellg();
      fin_camac.seekg (0, fin_camac.beg);
  
      if(debug)  std::cout << fileSize << std::endl;
      char * cBuf = new char [4];

      short * lbuf = new short [2];
      //c_tree->Branch("entry", entry, "entry/S"); not needed
      //g_tree->Branch("g_v", gvBuffer, "g_v[TBINS]/F");
      //g_tree->Branch("t", tBuffer, "t[TBINS]/S");

      // Read the file in to the buffer
      for(int i = 0; i < fileSize / (4); i++) {
	fin_camac.read(cBuf, 4);

	lbuf = reinterpret_cast<short*> (cBuf);
	if(lbuf[0] < 8192) {
	  lBuffer[0] = lbuf[0];
	  lBuffer[1] = lbuf[1];
	}
	else {
	  lBuffer[0] = lbuf[1];
	  lBuffer[1] = lbuf[0];
	}
      
	c_tree->Fill();

	//if(debug) std::cout << entry[0] << ' ' << lBuffer[0] << ' ' << lBuffer[1] << std::endl;

	//delete cBuffer;
	//delete sbuf;
      }
      //delete lbuf;
      //delete cBuf;
    }
    fin_camac.close();
  }

  h_tree=new TTree("h_tree","camac histogram");

  short adcBuffer[HBINS*8192/8192];
  short tdcBuffer[HBINS*8192/8192];
  short hb[HBINS*8192/8192];
  h_tree->Branch("adch", adcBuffer, 
		 Form("adch[%i]/S", HBINS*8192/8192));
  h_tree->Branch("tdch", tdcBuffer, 
		 Form("tdch[%i]/S", HBINS*8192/8192));
  for (int irun=0; irun<fNruns; ++irun) {
    int run_no = fRunList[irun];
    if(debug) std::cout << run_no << std::endl;
    int series = fSerList[0];
    if(debug)  std::cout << run_no << std::endl;
    if (!path) {
      ipath="/home/lifetime/data";
    }
    if(debug)  std::cout << ipath << std::endl;
    GetFileNames(fDataFileName, ipath, run_no, series);
    fin_camac.open(fDataFileName[2], std::ifstream::binary);
    if(debug)  std::cout << "Camac hist read" << std::endl;
    if(debug)  std::cout << "Declare hist tree" << std::endl;

    if (fin_camac != NULL) {
      fin_camac.seekg(0, fin_camac.end);
      long fileSize = fin_camac.tellg();
      fin_camac.seekg (0, fin_camac.beg);
  
      if(debug)  std::cout << fileSize << std::endl;
      char * hBuf = new char [HBINS*4];

      //h_tree->Branch("hist", adcBuffer, Form("hist[%i]/S", HBINS));
      //h_tree->Branch("hb", hb, Form("hb[%i]/S", HBINS*8192/8192));
      //h_tree->Branch("hb", hb, "hb[8192]/S");

      // Read the file in to the buffer
      //for(int i = 0; i < fileSize / (4); i++) {
      fin_camac.read(hBuf, HBINS*4);
      short * tbuf = new short [HBINS*2];
      tbuf = reinterpret_cast<short*> (hBuf);
      for(int qq = 0; qq < HBINS*8192/8192; qq++) {
	adcBuffer[qq] = tbuf[qq];
	tdcBuffer[qq] = tbuf[qq+8192];
	//hb[qq] = qq;
	//if(tdcBuffer[qq] != 0)
	//  if(debug) std::cout << qq << ' ' << tdcBuffer[qq] << std::endl;
	//h_tree->Fill();
      }
      h_tree->Fill();
      //delete hBuf;
    }
    fin_camac.close();
  }

  // 8192 entry version
  //  fin.open(fDataFileName[2], std::ifstream::binary);
  // if(debug)  std::cout << "Camac hist read" << std::endl;
  // if(debug)  std::cout << "Declare hist tree" << std::endl;
  // h_tree=new TTree("h_tree","camac histogram");
  // if (fin) {
  //   short adcBuffer[1];
  //   short tdcBuffer[1];
  //   short hb[HBINS*8191/8192];
  //   fin.seekg(0, fin.end);
  //   long fileSize = fin.tellg();
  //   fin.seekg (0, fin.beg);
  
  //   if(debug)  std::cout << fileSize << std::endl;
  //   char * cBuf = new char [HBINS*4];

  //   h_tree->Branch("adch", adcBuffer, Form("adch[%i]/S", 1));
  //   h_tree->Branch("tdch", tdcBuffer, Form("tdch[%i]/S", 1));
  //   //h_tree->Branch("hb", hb, Form("hb[%i]/S", HBINS*8191/8192));
  //   //h_tree->Branch("hb", hb, "hb[8192]/S");

  //   // Read the file in to the buffer
  //   //for(int i = 0; i < fileSize / (4); i++) {
  //     fin.read(cBuf, HBINS*4);
  //     short * tbuf = new short [HBINS*2];
  //     tbuf = reinterpret_cast<short*> (cBuf);
  //     for(int qq = 0; qq < HBINS; qq++) {
  // 	adcBuffer[0] = tbuf[qq];
  // 	tdcBuffer[0] = tbuf[qq+8192];
  // 	//hb[qq] = qq;
  // 	//if(debug) std::cout << qq << ' ' << hb[qq] << std::endl;
  // 	h_tree->Fill();
  //     }
  //     //h_tree->Fill();
  // }
  // fin.close();

}

void TauData::get_derivative(float values[], float deriv[], int smooth)
{
  double tderiv[2048];
  for(int i = 0; i < 2048; i++) {
    if (i > 6 && i < 2041) {
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

void TauData::get_smoothed(float smooth[], float tsmooth[], int smoothness)
{
  int nweights = 0;
  float weights[11];
  if(smoothness == 1) { //3 point triangle
    weights[0] = 1/3; 
    weights[1] = 1/3; 
    weights[2] = 1/3;
    nweights = 3;
  }
  if(smoothness == 2) { //5 point triangle
    weights[0] = 1/9; 
    weights[1] = 2/9; 
    weights[2] = 1/3;
    weights[3] = 2/9; 
    weights[4] = 1/9; 
    nweights = 5;
  }
  if(smoothness == 3) { //7 point triangle
    weights[0] = 1/16; 
    weights[1] = 2/16; 
    weights[2] = 3/16;
    weights[3] = 4/16; 
    weights[4] = 3/16;
    weights[5] = 2/16; 
    weights[6] = 1/16; 
    nweights = 7;
  }
  if(smoothness == 11) { //3 point gaussian
    weights[0] = 0.1842; 
    weights[1] = 0.6316; 
    weights[2] = 0.1842;
    nweights = 3;
  }
  if(smoothness == 12) { //5 point gaussian
    weights[0] = 0.0647; 
    weights[1] = 0.2447; 
    weights[2] = 0.3813;
    weights[3] = 0.2447; 
    weights[4] = 0.0647; 
    nweights = 5;
  }
  if(smoothness == 13) { //7 point gaussian
    weights[0] = 0.0356; 
    weights[1] = 0.1104; 
    weights[2] = 0.2176;
    weights[3] = 0.2729; 
    weights[4] = 0.2176; 
    weights[5] = 0.1104; 
    weights[6] = 0.0356; 
    nweights = 7;
  }
  if(smoothness == 14) { //9 point gaussian
    weights[0] = 0.0238; 
    weights[1] = 0.062; 
    weights[2] = 0.1228;
    weights[3] = 0.1852; 
    weights[4] = 0.2124; 
    weights[5] = 0.1852; 
    weights[6] = 0.1228; 
    weights[7] = 0.062;
    weights[8] = 0.0238;  
    nweights = 9;
  }
  if(smoothness == 15) { //11 point gaussian
    weights[0] = 0.017587;
    weights[1] = 0.040124;
    weights[2] = 0.076209;
    weights[3] = 0.120507;
    weights[4] = 0.15864;
    weights[5] = 0.173865;
    weights[6] = 0.15864;
    weights[7] = 0.120507;
    weights[8] = 0.076209;
    weights[9] = 0.040124;
    weights[10] = 0.017587;
    nweights = 11;
  }
  for(int i = 0; i < 2048; i++) {
    tsmooth[i] = 0;
    float weightsum = 0;
    for(int j = -1*(nweights-1)/2; j < (nweights-1)/2; j++) {
      if(i + j > 0 && i + j < 2048) {
	tsmooth[i] += smooth[i+j]*weights[j+(nweights-1)/2];
	weightsum += weights[j+(nweights-1)/2];
      }
    }
    tsmooth[i] /= weightsum;
  }

  return;
}



/*

derivative expressions

t3.g_tree->Draw("(g_v[t+1]-g_v[t-1])/2.0:t","p[3]<5 && t > 0 && t < 2047");

t3.g_tree->Draw("(-g_v[t+2]+8.0*g_v[t+1]-8.0*g_v[t-1]+g_v[t-2])/12.0:t","p[3]<5 && t > 1 && t <2046");

t3.g_tree->Draw("(g_v[t+3]-9.0*g_v[t+2]+45*g_v[t+1]-45*g_v[t-1]+9.0*g_v[t-2]-g_v[t-3])/60.0:t","p[3]<5 && t > 2 && t < 2045");

t3.g_tree->Draw("(-g_v[t+4]+32.0/3.0*g_v[t+3]-56.0*g_v[t+2]+224.0*g_v[t+1]-224*g_v[t-1]+56.0*g_v[t-2]-32.0/3.0*g_v[t-3]+g_v[t-4])/280.0:t","p[3]<5 && t > 3 && t < 2044");

2nd order

t3.g_tree->Draw("-1/560*g_v[t-4]+8/315*g_v[t-3]-1/5*g_v[t-2]+8/5*g_v[t-1]-205/72*g_v[t]+8/5*g_v[t+1]-1/5*g_v[t+2]+8/315*g_v[t+3]-1/560*g_v[t+4]:t","p[3]<5 && t > 3 && t < 2044");

cuts on deriv

t3.g_tree->Draw("g_v:t","t > 3 && t < 2044 (-g_v[t+4]+32.0/3.0*g_v[t+3]-56.0*g_v[t+2]+224.0*g_v[t-1]-224*g_v[t+1]+56.0*g_v[t-2]-32.0/3.0*g_v[t-3]+g_v[t-4])/280.0 < 0.05")

.L plot_histograms.cxx

plot_histograms("series_103_to_105_300mmpips")
plot_histograms("series_106_to_112_300mmpips")
plot_histograms("series_113_pressuretest_300mmpips")
plot_histograms("series_115_to_117_600mm")
plot_histograms("series_119_to_129_300mmpips")
plot_histograms("series_130_to_138_300mmpips")
plot_histograms("series_103_to_138_300mmpips")


plot_histograms("series_103_to_105_300mmpips_2")
plot_histograms("series_106_to_112_300mmpips_2")
plot_histograms("series_113_pressuretest_300mmpips_2")
plot_histograms("series_115_to_117_600mm_2")
plot_histograms("series_119_to_129_300mmpips_2")
plot_histograms("series_130_to_138_300mmpips_2")
plot_histograms("series_103_to_138_300mmpips_2")
plot_histograms("series_106_to_138_300mmpips_2")
plot_histograms("series_139_to_141_300mmpips_olddet_2")


*/
