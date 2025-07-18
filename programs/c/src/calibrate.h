#ifndef __CALIBRATE_H__
#define __CALIBRATE_H__

#include "globals.h"
#include "exporters.h"
#include "solver.h"

typedef struct
{
  // ......................................................................
  // final demand
  // ......................................................................

  //combine final demand from different sectors... cons is CES, inv is Cobb-Douglas
  // Cons = [eps * goods^rho + (1-eps) * services^rho]^(1/rho)
  // Inv = G * [goods^eps * services^(1-eps)]
  double rho;
  double eps[NC][NF][NS];
  double G[NC];

  // combine final demand from different countries into sector-specific bundles
  // f_s = H * [sum_{j=1}^{NC} (theta_j * f_j)^sig]^(1/sig)
  double sig[NC][NS];
  double theta[NC][NS][NC];
  double H[NC][NS];
  double H2[NC][NS][NC];

  // ......................................................................
  // gross output parameters
  // ......................................................................

  // combine intermediates from different countries into sector-specific bundles
  // M_s = C * [sum_{j=1}^{NC} (mu_j * m_j)^zeta]^(1/zeta)
  double zeta[NC][NS];
  double mu[NC][NS][NC];
  double M[NC][NS];
  double M2[NC][NS][NC];

  // combine value added and intermediate bundles from different sectors... Leontief
  // Gross output = min[VA/lam_va, M_goods/lam_goods, M_svcs/lam_svcs]
  // CD alt       = B * VA^lam_va * M_goods^lam_goods * M_csvs^lam_svcs
  double lam_va[NC][NS];
  double lam[NC][NS][NS];
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
  //double tariff_data[NC][NS][NC];
  double tau_m_ts[NT+1][NC][NS][NC];
  double tau_f_ts[NT+1][NC][NS][NC];

  // import NTB: destination-sector-source
  double ntb_m_ts[NT+1][NC][NS][NC];
  double ntb_f_ts[NT+1][NC][NS][NC];
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
  double md0[NC][NS][NS];
  double md0_[NC][NS][NS];
  double ex0[NC][NC];
  double im0[NC][NC];
  double nx0[NC][NC];
  double c0[NC][NS];
  double i0[NC][NS];
  double m0[NC][NS];
  double m02[NC][NS][NC];
  double q0[NC][NS];
  double q02[NC][NS][NC];
  double im02[NC][NS][NC];
  double ex02[NC][NS][NC];
  double nx02[NC][NS][NC];
  double lshare0[NC][NS];
  double n0[NC][NS][NC-1];

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
  double eta;
  double sig_z[NC][NS];
  double kappa0[NC][NS][NC-1];
  double kappa1[NC][NS][NC-1];
  uint Ji[NC][NC-1];
 
}params;

params ppp0[NTH];
params ppp1[NTH];

exporter_vars ev[NC][NS][NC-1];

uint copy_params(params * dest, const params * src);
uint set_nontargeted_params(params * p);
uint load_iomat(params * p);
void load_ts_params(params * p);
//uint load_tariffs(params * p);
void set_tariffs(params * p, uint scenario);
uint store_base_period_values(params * p);
uint calibrate_prod_params(params * p);
uint calibrate_fin_params(params * p);
uint calibrate_hh_params(params * p);
uint calibrate_firm_params();
uint stack_calvars();
uint calibrate();
uint write_params();

int calfunc_f(const gsl_vector * x, void * data, gsl_vector * f);
int calfunc_df(const gsl_vector * x, void * data, gsl_matrix * J);
int calfunc_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);

static inline double prod_go(double va, const double md[NS], double lam_va, const double lam[NS])
{
  if(cobb_douglas_flag==0)
    {
      double mm = HUGE_VAL;
      uint i;
      for(i=0; i<NS; i++)
	{
	  mm = fmin(mm,md[i]/lam[i]);
	}
      return fmin( va/lam_va, mm );
    }
  else
    {
      double tmp = pow(va,lam_va);
      uint i;
      for(i=0; i<NS; i++)
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

static inline double prod_inv(const double x[NS], const double eps[NS], double G)
{
  return G * pow(x[1],eps[1]) * 
    pow(x[2],eps[2]) * 
    pow(x[3],eps[3]) *
    pow(x[4],eps[4]);
}

static inline double prod_m(const double m2[NC], double M, const double mu[NC], double zeta)
{
  return M * pow( mu[0]*pow(m2[0],zeta) + 
		  mu[1]*pow(m2[1],zeta) + 
		  mu[2]*pow(m2[2],zeta) +
		  mu[3]*pow(m2[3],zeta) +
		  mu[4]*pow(m2[4],zeta) + 
		  mu[5]*pow(m2[5],zeta), 1.0/zeta );
}

static inline double prod_q(const double q2[NC], double H, const double theta[NC], double sig)
{
  return H * pow( theta[0]*pow(q2[0],sig) + 
		  theta[1]*pow(q2[1],sig) + 
		  theta[2]*pow(q2[2],sig) +
		  theta[3]*pow(q2[3],sig) +
		  theta[4]*pow(q2[4],sig) + 
		  theta[5]*pow(q2[5],sig), 1.0/sig );
}

static inline double muc(const double c[NS], double l, double lbar, const double eps[NS], double rho, double phi, double psi, uint s)
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
    pow(DOT_PROD_EX(c,eps,NS,rho),psi*phi/rho-1.0) * 
    pow(leisure,(1.0-phi)*psi);
}

static inline double mul(const double c[NS], double l, double lbar, const double eps[NS], double rho, double phi, double psi)
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
    pow(DOT_PROD_EX(c,eps,NS,rho),psi*phi/rho) * 
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
