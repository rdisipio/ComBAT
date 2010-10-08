#ifndef __COMMONS__
#define __COMMONS__

#include <iostream>
#include <map>

typedef const unsigned int nparam_t;
typedef const double       param_t;
typedef const double       syst_t;

using namespace std;


//typedef map< const string, vector< double > > SigmaArray;


///////////////////////////////////////////



namespace CH {
  nparam_t ELE = 0;
  nparam_t MU  = 1;
};

namespace DATAPOINT {
  nparam_t ELE_0 = 0;
  nparam_t MU_0  = 1;
};

namespace VAL {
  enum QUANTITY {
    Nobs, 
    NsigMC, 
    NbkgW_datadriven,
    NbkgQCD_datadriven,
    NbkgST,
    NbkgZ,
    NbkgOthers,
    W_c_frac,
    W_b_frac,
    eff_btag,
    eff_mistag_c,
    eff_mistag_lq,
    eff_sig,
    eff_bkg,
    iLumi
  };
};

/*
namespace VAL {
  nparam_t Nobs    = 0;
  nparam_t NsigMC  = 1;
  nparam_t NbkgW_datadriven = 2;
  nparam_t NbkgQCD_datadriven = 3;
  nparam_t NbkgST;
  nparam_t NbkgZ;
  nparam_t NbkgOthers;
  nparam_t W_c_frac = 4;
  nparam_t W_b_frac = 5;
  nparam_t eff_btag = 6;
  nparam_t eff_mistag_c = 7;
  nparam_t eff_mistag_lq = 8;
  nparam_t eff_sig = 9;
  nparam_t eff_bkg = 10;
  nparam_t iLumi   = 11;
};
*/

namespace PARAM {
  nparam_t xsec        = 0;

  nparam_t Nbkg_ele    = 1;
  nparam_t Nbkg_mu     = 2;
  nparam_t eff_sig_ele = 3;  
  nparam_t eff_sig_mu  = 4;

  nparam_t N_PARAMS    = 5;
};


namespace SYST {
  nparam_t JES         = 3;
  nparam_t Wnorm       = 4;
  nparam_t iLumi       = 5;
  nparam_t ISRFSR      = 6;
  nparam_t QCD         = 6;
  nparam_t bTagging    = 7;
  nparam_t HFW         = 8;
  
  nparam_t N_SYSTS     = 7;
};


#endif /* __COMMONS__ */
