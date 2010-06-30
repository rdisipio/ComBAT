#ifndef __COMMONS__
#define __COMMONS__

typedef const unsigned int nparam_t;
typedef const double       param_t;
typedef const double       syst_t;

namespace CH {
  nparam_t ELE = 0;
  nparam_t MU  = 1;
};

namespace DATAPOINT {
  nparam_t ELE_0 = 0;
  nparam_t MU_0  = 1;
  nparam_t ELE_JES_UP = 1;
  nparam_t ELE_JES_DOWN = 2;
  nparam_t ELE_WNORM_UP = 3;
  nparam_t ELE_WNORM_DOWN = 4;
  
 // nparam_t MU_0 = 5;
  nparam_t MU_JES_UP = 6;
  nparam_t MU_JES_DOWN = 7;
  nparam_t MU_WNORM_UP = 8;
  nparam_t MU_WNORM_DOWN = 9;
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

  nparam_t eff_sig_ele = 1;  
  nparam_t eff_sig_mu  = 2;
};


namespace SYST {
  nparam_t JES         = 3;
  nparam_t Wnorm       = 4;
  nparam_t iLumi       = 5;
};

typedef double** syst_matrix_t;

#endif /* __COMMONS__ */
