#ifndef __CALIBRATE_H__
#define __CALIBRATE_H__

#include "globals.h"
#include "solver.h"

typedef struct
{
  // ......................................................................
  // final demand
  // ......................................................................

  //combine final demand from different sectors... cons is CES, inv is Cobb-Douglas. Only inv uses construction
  // Cons = [eps1*up^rho + eps2*dn^rho + eps3*svcs^rho]^(1/rho)
  // Inv = G * [up^eps1 * dn*^eps2 * svcs^eps3 * cons^eps4]
  double rho;
  double eps[NC][NF][NS];
  double G[NC];
  
  // combine final demand from different countries into sector-specific bundles
  // f_s = H * [sum_{j=1}^{NC} (theta_j * f_j)^sig]^(1/sig)
  double sig[NC][NS-1];
  double theta[NC][NS-1][NC];
  double H[NC][NS-1];

  // ......................................................................
  // gross output parameters
  // ......................................................................

  // combine intermediates from different countries into sector-specific bundles
  // M_s = C * [sum_{j=1}^{NC} (mu_j * m_j)^zeta]^(1/zeta)
  double zeta[NC][NS-1];
  double mu[NC][NS-1][NC];
  double M[NC][NS-1];

  // combine value added and intermediate bundles from different sectors... Leontief
  // Gross output = min[VA/lam_va, M_up/lam_up, M_dn/lam_dn, M_svcs/lam_svcs]
  double lam_va[NC][NS];
  double lam[NC][NS][NS-1];
  double B[NC][NS];

  // value added... Cobb-Douglas
  // VA = A * [k^alpha * (gam * ell)^(1-alpha)]
  double alpha[NC][NS];
  double A[NC][NS];

  // ......................................................................
  // households
  // ......................................................................
  
  // capital formation
  double delta; // depreciation rate
  double tauk[NC];  // capital tax rate
  double rss; // steady-state real interest rate
  
  // household preferences
  double beta[NC]; // discount factors
  double psi; // intertemporal elasticity
  double phi[NC]; // consumption share

  // endowments
  double lbar[NC];
  double kk0[NC];
  double b0[NC];

  // ......................................................................
  // time series parameters
  // ......................................................................
  // sector-level productivities
  double a_ts[NT+1][NC][NS];

  // import tariffs: destination-sector-source
  double tau_m_ts[NT+1][NC][NS][NC];
  double tau_f_ts[NT+1][NC][NS][NC];

  // ......................................................................
  // base-period equilibrium values
  // ......................................................................
  double iomat[NS*NC+2][NS*NC + NF*NC + 1];
  double r0[NC];
  double ii0[NC];
  double ll0[NC];
  double y0[NC][NS];
  double va0[NC][NS];
  double k0[NC][NS];
  double l0[NC][NS];
  double md0[NC][NS][NS-1];
  double md0_[NC][NS][NS-1];
  double ex0[NC][NC];
  double im0[NC][NC];
  double nx0[NC][NC];
  double c0[NC][NS-1];
  double i0[NC][NS];
  double m0[NC][NS-1];
  double m02[NC][NS-1][NC];
  double q0[NC][NS-1];
  double q02[NC][NS-1][NC];
  double im02[NC][NS-1][NC];
  double ex02[NC][NS-1][NC];
  double nx02[NC][NS-1][NC];
  double lshare0[NC][NS];

  // ......................................................................
  // adjustment cost parameters
  // ......................................................................
  double etaM;
  double etaF;
  double etaK;
  double etaL;
 
  // ......................................................................
  // firm dynamics parameters
  // ......................................................................
  //double eta;
  //double sig_z[NC][NS];
  //double kappa0[NC][NS][NC-1];
  //double kappa1[NC][NS][NC-1];
  //uint Ji[NC][NC-1];
 
}params;

params ppp0[NTH];
params ppp1[NTH];

//exporter_vars ev[NC][NS][NC-1];

uint copy_params(params * dest, const params * src);
uint set_nontargeted_params(params * p);
uint load_iomat(params * p);
void load_ts_params(params * p);
void set_tariffs(params * p, uint scenario);
uint store_base_period_values(params * p);
uint calibrate_prod_params(params * p);
uint calibrate_fin_params(params * p);
uint calibrate_hh_params(params * p);
//uint calibrate_firm_params();
uint stack_calvars();
uint calibrate();
uint write_params();

int calfunc_f(const gsl_vector * x, void * data, gsl_vector * f);
int calfunc_df(const gsl_vector * x, void * data, gsl_matrix * J);
int calfunc_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);

// Gross output = min[VA/lam_va, M_up/lam_up, M_dn/lam_dn, M_svcs/lam_svcs]
static inline double prod_go(double va, const double md[NS], double lam_va, const double lam[NS-1])
{
  if(cobb_douglas_flag==0)
    {
      double mm = HUGE_VAL;
      uint i;
      for(i=0; i<NS-1; i++)
	{
	  mm = fmin(mm,md[i]/lam[i]);
	}
      return fmin( va/lam_va, mm );
    }
  else
    {
      double tmp = pow(va,lam_va);
      uint i;
      for(i=0; i<NS-1; i++) // note we do not include construction, which is last sector by design
	{
	  tmp = tmp * pow(md[i],lam[i]);
	}

      return tmp;
    }
}

static inline double prod_va(double k, double l, double A, double alpha)
{
  return A * pow(k,alpha) * pow(l,(1.0 - alpha));
}

// Inv = G * [up^eps1 * dn*^eps2 * svcs^eps3 * cons^eps4]
static inline double prod_inv(const double x[NS], const double eps[NS], double G)
{
  double tmp=1.0;
  for(int s=0; s<NS; s++)
    tmp = tmp * pow(x[s],eps[s]);

  return G*tmp;
  //return G * pow(x[0],eps[0]) * pow(x[1],eps[1]) * pow(x[2],eps[2]) * pow(x[3],eps[3]);
}

// M_s = C * [sum_{j=1}^{NC} (mu_j * m_j)^zeta]^(1/zeta)
static inline double prod_m(const double m2[NC], double M, const double mu[NC], double zeta)
{
  double tmp = 00;
  for(int c=0; c<NC; c++)
    tmp = tmp + mu[c]*pow(m2[c],zeta);

  return M*pow(tmp,1.0/zeta);

  /*
  return M * pow( mu[0]*pow(m2[0],zeta) + 
		  mu[1]*pow(m2[1],zeta) + 
		  mu[2]*pow(m2[2],zeta), 1.0/zeta );
  */
}

// f_s = H * [sum_{j=1}^{NC} (theta_j * f_j)^sig]^(1/sig)
static inline double prod_q(const double q2[NC], double H, const double theta[NC], double sig)
{
  double tmp = 0.0;
  for(int c=0; c<NC; c++)
    tmp = tmp + theta[c]*pow(q2[c],sig);

  return H*pow(tmp,1.0/sig);

  /*
  return H * pow( theta[0]*pow(q2[0],sig) + 
		  theta[1]*pow(q2[1],sig) + 
		  theta[2]*pow(q2[2],sig), 1.0/sig );
  */
}

// d u(c_up, c_dn, c_svcs, leis) / d c_s
static inline double muc(const double c[NS-1], double l, double lbar, const double eps[NS], double rho, double phi, double psi, uint s)
{
  
  double leisure;
  if(lbar-l > 0.0001)
    {
      leisure = lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(lbar-l));
    }
  
  return phi * eps[s] * pow(c[s],rho-1.0) * 
    pow(DOT_PROD_EX(c,eps,NS-1,rho),psi*phi/rho-1.0) * 
    pow(leisure,(1.0-phi)*psi);
}

static inline double mul(const double c[NS-1], double l, double lbar, const double eps[NS], double rho, double phi, double psi)
{
  double leisure;
  if(lbar-l > 0.0001)
    {
      leisure = lbar-l;
    }
  else
    {
      leisure = 0.0001 / log(0.0001-(lbar-l));
    }

  return (1.0-phi) * 
    pow(DOT_PROD_EX(c,eps,NS-1,rho),psi*phi/rho) * 
    pow(leisure,(1.0-phi)*psi - 1.0);
}

static inline double phiK(double x, double delta, double etaK)
{
  return (pow(delta,1.0-etaK) * pow(x,etaK) - (1.0-etaK)*(delta))/etaK;
}

static inline double dphiK(double x, double delta, double etaK)
{
  return pow(delta,1.0-etaK) * pow(x,etaK-1.0);
}

#endif
