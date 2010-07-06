#ifndef __COMMONS__
#define __COMMONS__

#include <iostream>
#include <map>

typedef const unsigned int nparam_t;
typedef const double       param_t;
typedef const double       syst_t;

using namespace std;


typedef map< const string, vector< double > > SigmaArray;


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
  nparam_t Nobs    = 0;
  nparam_t NsigMC  = 1;
  nparam_t Nbkg    = 2;
  nparam_t eff_sig = 3;
  nparam_t eff_bkg = 4;
  nparam_t iLumi   = 5;
};


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
