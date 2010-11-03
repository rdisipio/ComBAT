#include "InputData.h"

InputData::InputData(const char* file)
{
  // default constructor
  ifs.open(file);
}
InputData::~InputData()
{
  // default destructor
}

bool InputData::Read()
{
  Nds=0;
  string line;
  if(!ifs){
    cout << "Error in opening input card file !!! " << endl;
    return false;
  }
  while(!ifs.eof())
  {
    std::getline(ifs,line);
    readkey(line,"<Likelihood>", &likelihood); 
    readkey(line,"<Data>", &data); 
    readkey(line,"<Data-value>", &datavd); 
    readkey(line,"<Measure-lim>", &measlimd); 
    readkey(line,"<Measure-prior>", &meastypei); 
    readkey(line,"<Measure-par>", &measpard); 
    readkey(line,"<Measure>", &meas); 
    readkey(line,"<Parameters>", &par); 
    readkey(line,"<Parameters-values>", &parvd); 
    readkey(line,"<Syst>", &syst); 
    readkey(line,"<Syst-Prior>", &systtypei); 
    readkey(line,"<Syst-lim>", &systlimd); 
    for(int i=0; i<par.size(); i++)
    {
      readkey(line,"<Syst-"+par.at(i)+">", &systvd);
    }

    
    if(enddataset(line,"<EndDataSet>")){
      Nds++;
      datavd_vec.push_back(datavd);
      parvd_vec.push_back(parvd);
      syst_vec.push_back(syst);
      systvd_vec.push_back(systvd);
      systlimd_vec.push_back(systlimd);
      systtypei_vec.push_back(systtypei);
      if(systvd.size() != 2*par.size()*syst.size())
      {
        cout << " Error: inconsistency between Npar and Nsyst !!! " << endl; 
        return false;
      }
      datavd.clear();
      parvd.clear();
      syst.clear();
      systvd.clear();
      systlimd.clear();
      systtypei.clear();
    }
  }
  // Prepare final maps to be used as input for BAT
  StoreMeasPar();
  StoreData();
  StoreCentralValues();
  StoreSystType();
  StoreSystPar();
  AddGammaPar();
  map<string, vector<double> >::iterator itg;
  //  cout << "gammapar size " << gammapar.size() << endl; 
  for(itg=gammapar.begin(); itg != gammapar.end(); itg++)
    cout << " test test test test " << itg->first << " " << itg->second.size() << " " << itg->second.at(0) << " " << itg->second.at(1) << endl;
  ifs.close();
  return true;
}

void InputData::StoreData()
{
  map<string, vector<double> > datatv;
  for (int k=0; k<Nds; k++) {
    datatv.insert(make_pair(data.at(k),datavd_vec.at(k)));
    datav_vec.push_back(datatv);
    datatv.clear();
  }
}

void InputData::StoreMeasPar()
{
  measlim.insert(make_pair(meas.at(0),measlimd));
  meastyp.insert(make_pair(meas.at(0),meastypei));
  measpar.insert(make_pair(meas.at(0),measpard));
}



void InputData::StoreCentralValues()
{
  map<string,double> parcv;
  for (int k=0; k<Nds; k++) {
    for(int i=0; i<par.size(); ++i)
      parcv.insert(make_pair(par.at(i),parvd_vec.at(k).at(i)));
    parcv_vec.push_back(parcv);
    parcv.clear();
  }
}

vector< map<string,double> >* InputData::GetparCV()
{
  return &parcv_vec;
}

int InputData::GetNDataSet()
{
  return Nds;
}

void InputData::StoreSystPar()
{
  map<string, map<string,vector<double> > > parsyst;
  map<string,vector<double> > systsig;
  vector<double> v, vlim;
  vector<int> vtyp;
  map<string,vector<int> >::iterator itt;
  int id=0;
  for (int k=0; k<Nds; k++) {
    int kk1=0;
    for(int i=0; i<par.size(); ++i){      
      for(int j=0; j<syst_vec.at(k).size(); j++){
	itt=systtyp.find(syst_vec.at(k).at(j));
        if(itt->second.size() > 0) {
         v.push_back(systvd_vec.at(k).at(kk1++));
         v.push_back(systvd_vec.at(k).at(kk1++));
         systsig.insert(make_pair(syst_vec.at(k).at(j),v));
         v.clear();
	}
	else{kk1++; kk1++;}
      }
      parsyst.insert(make_pair(par.at(i),systsig));
      systsig.clear();
    }
    parsyst_vec.push_back(parsyst);
    parsyst.clear();
  }
}

void InputData::StoreSystType()
{
  vector<double> vlim1, vlim2;
  vector<int> vtyp;
  int id=0,kk2;
  map<string,vector<double> >::iterator itg;
  for (int k=0; k<Nds; k++) {
    kk2=0;
    for(int j=0; j<syst_vec.at(k).size(); j++){
      if( systtypei_vec.at(k).at(j) != -1) {
        if( !keyisfound(systintlim,syst_vec.at(k).at(j))){
          vlim1.push_back(systlimd_vec.at(k).at(kk2++));
	  vlim1.push_back(systlimd_vec.at(k).at(kk2++));
          systintlim.insert(make_pair(syst_vec.at(k).at(j),vlim1));
	  //          cout << " DEBUG " << syst_vec.at(k).at(j) << " " << vlim1.at(0) << " " << vlim1.at(1) << endl;
          if(systtyp.find(syst_vec.at(k).at(j)) == systtyp.end()){
            vtyp.push_back(++id);
	    vtyp.push_back(systtypei_vec.at(k).at(j));
            systtyp.insert(make_pair(syst_vec.at(k).at(j),vtyp));
            vtyp.clear();
	  }
          vlim1.clear();
	}
        else {kk2++; kk2++;}
      }
      else{kk2++; kk2++;}
    }
  }
}

bool InputData::keyisfound(map<string,vector<double> > mp, string s)
{
  if(mp.find(s) != mp.end()) return true;
  return false;
}
void InputData::PrintparCV()
{
  map<string,double>::iterator it;
  for(int i=0; i<parcv_vec.size(); i++)
  {
    cout << " ==========> Next Dataset <=========== " << endl;
    for (it=parcv_vec.at(i).begin(); it!=parcv_vec.at(i).end(); it++)
    {
      cout << it->first << " " << it->second << endl;
    } 
  }
}

void InputData::PrintSystIntLim()
{
  map<string, vector<double> >::iterator it1;
  cout << " ==========> Systematics Integration limits <=========== " << endl;
  for(it1=systintlim.begin(); it1!=systintlim.end(); it1++)
  {
    if(it1->second.size()==4)
    cout << it1->first << "...... " << it1->second.at(0) << " " << it1->second.at(1) << 
      "            " << it1->second.at(2) << " " << it1->second.at(3) << endl;
    else cout << it1->first << "...... " << it1->second.at(0) << " " << it1->second.at(1) << endl;
  }
}

void InputData::PrintSystTyp()
{
  map<string, vector<int> >::iterator it1;
  cout << " ==========> Systematics id and type <=========== " << endl;
  for(it1=systtyp.begin(); it1!=systtyp.end(); it1++)
  {
    cout << it1->first << "...... " << it1->second.at(0) << "  " << it1->second.at(1) << endl;
  }
}

void InputData::PrintparSyst()
{
  map<string, vector<double> >::iterator it1;
  map<string, map<string, vector<double> > >::iterator it2;
  for(int i=0; i<parsyst_vec.size(); i++){
    cout << " ==========> Next Dataset <=========== " << endl;
    for(it2=parsyst_vec.at(i).begin(); it2!=parsyst_vec.at(i).end(); it2++){
      cout << " Parameter: " << it2->first << endl;
      for(it1=it2->second.begin(); it1!=it2->second.end(); it1++){
	cout << "------------------> " << it1->first << ": " << it1->second.at(0) << " " << it1->second.at(1) << endl;
      }
    }
  }
}

void InputData::AddGammaPar()
{
  map<string, vector<int> >::iterator it0;
  map<string, vector<double> >::iterator it1;
  map<string, map<string, vector<double> > >::iterator it2;
  vector<double> gammap;
  for(it0=systtyp.begin(); it0!=systtyp.end(); it0++){
    if(it0->second.at(1) == 2){
      for(int i=0; i<parsyst_vec.size(); i++){
        for(it2=parsyst_vec.at(i).begin(); it2!=parsyst_vec.at(i).end(); it2++){
          it1=it2->second.find(it0->first);
	  if(it1!=it2->second.end()) {
	    if(it1->second.at(0) != 0){
	      //  	      cout << " -+-+-+->      " << it1->first << " " << it1->second.at(0) << " " << it1->second.at(1) << endl; 
              gammap.push_back(it1->second.at(0));
              gammap.push_back(it1->second.at(1));
              gammapar.insert(make_pair(it0->first,gammap));
              gammap.clear();
	    }
	  }
        }
      }
    }    
  }
}

vector< map<string, vector<double> > >* InputData::GetData()
{
  return &datav_vec;
}

void InputData::SetNobs(vector<double> obs)
{
  map<string, vector<double> >::iterator it;
  for(int is=0; is<obs.size(); is++){
    it=datav_vec.at(is).begin();
    it->second.at(0)=obs.at(is);
  }
}

map<string, vector<double> >* InputData::GetMeasPar()
{
  return &measpar;
}

map<string, vector<double> >* InputData::GetGammaPar()
{
  return &gammapar;
}

map<string, vector<int> >* InputData::GetMeasTyp()
{
  return &meastyp;
}

map<string, vector<double> >* InputData::GetMeasLim()
{
  return &measlim;
}

vector< map<string, map<string, vector<double> > > >* InputData::GetparSyst()
{
  return &parsyst_vec;
}

map<string, vector<double> >* InputData::GetSystIntLim()
{
  return &systintlim;
}

map<string, vector<int> >* InputData::GetSystTyp()
{
  return &systtyp;
}

vector<int>* InputData::GetLikelihood()
{
  return &likelihood;
}

void InputData::Print()
{
  for(int i=0; i<par.size(); ++i)
  {
    cout << "<Parameters> " << par.at(i) << endl;
  }  
  for (int k=0; k<Nds; k++) {
    for(int i=0; i<parv_vec.at(k).size(); ++i)
    {
      cout << "<Parameters> " << k << " " << parv_vec.at(k).at(i) << endl;
    }  
    for(int i=0; i<syst_vec.at(k).size(); ++i)
    {
      cout << "<Syst> " << syst_vec.at(k).at(i) << endl;
    }  

    for(int i=0; i<systtypei_vec.at(k).size(); ++i)
    {
      cout << "<Syst-type> " << k << " " << systtypei_vec.at(k).at(i) << endl;
    }  

    for(int i=0; i<systvd_vec.at(k).size(); ++i)
    {
      cout << "<Systv> " << k << " " << systvd_vec.at(k).at(i) << endl;
    }  

    for(int i=0; i<systlimd_vec.at(k).size(); ++i)
    {
      cout << "<Syst-lim> " << k << " " << systlimd_vec.at(k).at(i) << endl;
    }  
  }
}

void InputData::readkey(string p, string key, vector<string>* str){
  stringstream ss;
  string a;
  string b;
  ss << p;
  ss >> a; 
  if( a == key ) {
    while (!ss.eof()) {
      b.erase();
      ss >> b;
      if (b.find(a,0) == string::npos && b.find(a,b.size()-1) == string::npos) 
      {
	str->push_back(b);
      }
    }
  }
}

void InputData::readkey(string p, string key, vector<double>* dbl){
  stringstream ss;
  string a;
  string b;
  ss << p;
  ss >> a; 
  if( a == key ) {
    while (!ss.eof()) {
      b.erase();
      ss >> b;
      if (b.find(a,0) == string::npos && b.find(a,b.size()-1) == string::npos) 
      {
	dbl->push_back(std::atof(b.c_str()));
      }
    }
  }
}

void InputData::readkey(string p, string key, vector<int>* integ){
  stringstream ss;
  string a;
  string b;
  ss << p;
  ss >> a; 
  if( a == key ) {
    while (!ss.eof()) {
      b.erase();
      ss >> b;
      if (b.find(a,0) == string::npos && b.find(a,b.size()-1) == string::npos) 
      {
	integ->push_back(std::atoi(b.c_str()));
      }
    }
  }
}

bool InputData::enddataset(string p, string key){
  stringstream ss;
  string a;
  string b;
  ss << p;
  ss >> a; 
  if( a == key ) return true;
  else return false;  
}
