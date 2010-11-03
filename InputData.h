#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

using namespace std;

class InputData
{

public:

  InputData(const char*);
  ~InputData();
  bool Read();
  void Print();
  void PrintparSyst();
  void PrintSystIntLim();
  void PrintSystTyp();
  void PrintparCV();
  int GetNDataSet();

  void SetNobs(vector<double>);
  vector< map<string,vector<double> > >* GetData();
  map<string,vector<double> >* GetGammaPar();
  map<string,vector<double> >* GetMeasPar();
  map<string,vector<double> >* GetMeasLim();
  map<string,vector<int> >* GetMeasTyp();
  vector< map<string,double> >* GetparCV();
  vector< map<string, map<string,vector<double> > > >* GetparSyst();
  map<string,vector<double> >* GetSystIntLim();
  map<string,vector<int> >* GetSystTyp();
  vector<int>* GetLikelihood();
private:
  void readkey(string, string, vector<string>*);  
  void readkey(string, string, vector<double>*);  
  void readkey(string, string, vector<int>*);  
  bool enddataset(string, string);  
  void StoreMeasPar();
  void StoreData();
  void StoreCentralValues();
  void StoreSystPar();
  void StoreSystType();
  void AddGammaPar();
  bool keyisfound(map<string,vector<double> >, string);

  ifstream ifs;
  vector<string> meas, data;
  vector<double> measlimd, measpard, datavd;
  vector < vector<double> > datavd_vec;
  vector<int> meastypei, likelihood;  
  vector<double> parvd, systvd, systlimd;
  vector< vector<double> > parvd_vec, systvd_vec, systlimd_vec;
  vector<int> systtypei;
  vector< vector<int> > systtypei_vec;
  vector<string> par, parv, syst, systv, systlim, systtype;
  vector< vector<string> > parv_vec, syst_vec, systv_vec, systlim_vec, systtype_vec;
  int Nds;
  vector< map<string,double> > parcv_vec;
  vector< map<string, map<string,vector<double> > > > parsyst_vec;
  vector< map<string, vector<double> > > datav_vec;
  map<string,vector<double> > systintlim;
  map<string,vector<int> > systtyp;
  map<string,vector<double> > measlim, measpar;
  map<string,vector<int> > meastyp;
  map<string,vector<double> > gammapar;
};
