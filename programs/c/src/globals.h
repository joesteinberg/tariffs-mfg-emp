#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_cdf.h>
#include <lapacke.h>

// Regions. NC = 4 (baseline): USA, China, low-tariff partners, high-tariff partners (non-China
// partners split at the import-weighted median 2025 tariff increase). NC = 3 (make model_nc3, or
// -DNC=3; the original version, outputs _nc3): USA, China, rest of world. In both, the last region
// holds the residual bonds (market clearing) and "non-China partners" means regions 2..NC-1.
#ifndef NC
#define NC 4
#endif
#define USA 0
#define CHN 1
#if NC==3
#define ROW 2
#elif NC==4
#define LTP 2 // low-tariff partners
#define HTP 3 // high-tariff partners
#else
#error "NC must be 3 or 4"
#endif
static const char * const CNAME[NC] = {"usa","chn",
#if NC==3
				       "row"
#else
				       "ltp","htp"
#endif
};


#define NS 6 // 4 sectors (0 = upstream-hi, 1 = upstream-hi, 2 = downstreadm-hi, 3 - downstream-lo, 4 = services, 5 = construction)
#define UPH 0
#define UPL 1
#define DNH 2
#define DNL 3
#define SVC 4
#define CNS 5

#define NF 2 // 2 final uses (0=Cons, 1=Inv)
#define FIN 0
#define INT 1

#ifndef NT
#define NT 100 // time periods until convergence (override at compile time: make model_nt200)
#endif

#define TINY 1.0e-12
#define TINYSQ 1.0e-6

#ifdef _OPENMP
#define NTH  128
#include <omp.h>
#else
#define NTH 1
#endif


typedef unsigned int uint;
FILE * logfile;
uint eval_eqm_once_flag;
uint eval_bgp_once_flag;
uint read_seed;
uint write_seed;
uint target_sector_flag;
uint target_country_flag;
// -c: 0 China only; 1 all non-China partners; 2 all partners; with NC = 4 also 3 low-tariff
// partners only and 4 high-tariff partners only. Partner j (1..NC-1) is targeted:
static inline uint target_partner(uint j)
{
  switch(target_country_flag)
    {
    case 0: return j==CHN;
    case 1: return j!=CHN && j!=USA;
    case 2: return j!=USA;
#if NC==4
    case 3: return j==LTP;
    case 4: return j==HTP;
#endif
    default: return 0;
    }
}
uint retaliation_flag;
uint adjustment_flag;
uint duration_flag;
uint tariff_policy_flag;
uint m_adj_cost;
uint f_adj_cost;
uint k_adj_cost;
uint l_adj_cost;
uint fixl;
double psib; // portfolio adjustment cost psi_B (-b <value>; appendix); 0 in the baseline
char psib_suf[32]; // output-name suffix when psi_B > 0: "_psib<value>"
uint cobb_douglas_flag;
uint cobb_douglas_flag2;
uint par;
uint ghh_prefs;
uint eval_calfn_once;
double tariff;
double phiz; // sourcing adjustment cost phi_Z (etaM = etaF); default 2.9, set with -z
char phiz_suf[32]; // output-name suffix when -z is given: "_phiZ_lo", "_phiZ_hi", or "_pz<value>"
// Alternative phi_Z values from calibrate_phiz.py: short-run Armington elasticity (US, impact,
// symmetric 25% all-goods tariff) of 1.0 (lo, weak frictions) and 0.2 (hi, strong frictions).
// The baseline 2.9 gives 0.473 (4 regions). Keep in sync with PHIZ in programs/python/model_io.py.
#define PHIZ_LO 0.6645
#define PHIZ_HI 11.4

// Long-run Armington elasticity vectors (-e NAME), by goods sector (Oil = UPH, Steel = UPL,
// Toys = DNH, Cars = DNL) for intermediate (zeta) and final (sig) sourcing; services stay at 1.05.
// Output suffix "_sig_<name>". SIG_PHIZ is phi_Z re-targeted (calibrate_phiz.py --sig NAME) so the
// impact Armington elasticity in the symmetric experiment equals the baseline's 0.4733; 0 means not
// yet calibrated, in which case -z must be given.
//   cp:    Caliendo-Parro (2015) theta + 1, trade-weighted by cluster (icio.py, 4-region weights)
//   krs:   Kehoe-Ruhl-Steinberg (2018): 3 for intermediates, 2 for final goods, every sector
//   half:  baseline with half the dispersion: sbar*(s/sbar)^0.5, rescaled to the baseline
//          trade-weighted mean sbar = 3.25 (weights as in the `ae` statistic)
//   u2, umean, u5: the same elasticity in every goods sector: 2, 3.25 (= sbar), 5
#define NSIGV 7
static const char * const SIG_NAMES[NSIGV] = {"base","cp","krs","half","u2","umean","u5"};
static const double SIG_INT[NSIGV][4] = {
  {8.0, 2.0, 4.0, 1.05},
  {20.13, 3.87, 12.82, 3.14},
  {3.0, 3.0, 3.0, 3.0},
  {5.50, 2.75, 3.89, 1.99},
  {2.0, 2.0, 2.0, 2.0},
  {3.25, 3.25, 3.25, 3.25},
  {5.0, 5.0, 5.0, 5.0}};
static const double SIG_FIN[NSIGV][4] = {
  {8.0, 2.0, 4.0, 1.05},
  {20.13, 3.87, 12.82, 3.14},
  {2.0, 2.0, 2.0, 2.0},
  {5.50, 2.75, 3.89, 1.99},
  {2.0, 2.0, 2.0, 2.0},
  {3.25, 3.25, 3.25, 3.25},
  {5.0, 5.0, 5.0, 5.0}};
static const double SIG_PHIZ[NSIGV] = {2.9, 4.568, 2.95, 3.423, 2.522, 3.648, 4.578};
uint sig_vec; // 0 = baseline (8/2/4/1.05)
uint phiz_given; // 1 if -z was given (overrides the vector's re-targeted phi_Z)
char sig_suf[32];
double tariff_scale; // scales all tariffs in set_tariffs (1 = full policy); used for continuation
int TSHOCK;

static inline void set_all_v(double * v, uint n, double x)
{
  int i;
  for(i=0; i<n; i++)
    {
      v[i] = x;
    }
}
#define SET_ALL_V(v,n,x) set_all_v( (double *)(v), (n), (x) )

static inline double sum(const double * v, uint n)
{
  uint i;
  double mysum = 0.0;
  for(i=0; i<n; i++)
    {
      mysum = mysum + v[i];
    }
  return mysum;
}
#define SUM(v,n) (sum ( (double *)(v), (n) ))

static inline double dot_prod(const double * v1, const double * v2, uint n)
{
  uint i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + v1[i]*v2[i];
    }
  return sum;
}
#define DOT_PROD(v1,v2,n) (dot_prod( (double *)(v1), (double *)(v2), (n) ))

static inline double dot_prod_ex(const double * v1, const double * v2, uint n, double ex)
{
  uint i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + pow(v1[i],ex)*v2[i];
    }
  return sum;  
}
#define DOT_PROD_EX(v1,v2,n,ex) (dot_prod_ex( (double *)(v1), (double *)(v2), (n), (ex) ))

static inline double dot_prod_ex2(const double * v1, const double * v2, uint n, double ex1, double ex2)
{
  uint i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + pow(v1[i],ex1)*pow(v2[i],ex2);
    }
  return sum;  
}
#define DOT_PROD_EX2(v1,v2,n,ex1,ex2) (dot_prod_ex2( (double *)(v1), (double *)(v2), (n), (ex1), (ex2) ))

static inline double vec_mul_ex(const double * v, const double * ex, uint n)
{
  uint i;
  double prod  = 0.0;
  for(i=0; i<n; i++)
    {
      prod = prod * pow(v[i],ex[i]);
    }
  return prod;
}
#define VEC_MUL_EX(v,ex,n) (vec_mul_ex( (double *)(v), (double *)(ex), (n) ))

static inline void copy_subvector(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = src[i];
    }
}
#define COPY_SUBVECTOR(d,s,n) ( copy_subvector( (double *)(d), (const double *)(s), (n) ) )

static inline void copy_subvector_log(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = log(src[i]);
    }
}
#define COPY_SUBVECTOR_LOG(d,s,n) ( copy_subvector_log( (double *)(d), (const double *)(s), (n) ) )

static inline void copy_subvector_log_plus_one(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = log(src[i]+1.0);
    }
}
#define COPY_SUBVECTOR_LOG_PLUS_ONE(d,s,n) ( copy_subvector_log_plus_one( (double *)(d), (const double *)(s), (n) ) )

static inline void copy_subvector_exp(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = exp(src[i]);
    }
}
#define COPY_SUBVECTOR_EXP(d,s,n) ( copy_subvector_exp( (double *)(d), (const double *)(s), (n) ) )

static inline void copy_subvector_exp_minus_one(double * dest, const double * src, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      dest[i] = exp(src[i]) - 1.0;
    }
}
#define COPY_SUBVECTOR_EXP_MINUS_ONE(d,s,n) ( copy_subvector_exp_minus_one( (double *)(d), (const double *)(s), (n) ) )


static inline void write_array(double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"wb");
  uint i;    
  for(i=0; i<n; i++)
    {
      fprintf(file,"%.15f\n",v[i]);
    }    
  fclose(file);
}

static inline void write_vec_bin(const double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"wb");
  fwrite(v,n*sizeof(double),1,file);
  fclose(file);
}

static inline void write_vec_txt(const double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"wb");
  uint i;
  for(i=0; i<n; i++)
    {
      fprintf(file,"%0.15f\n",v[i]);
    }
  fclose(file);
}

static inline uint read_vec_bin(double * v, uint n, char * fname)
{
  FILE * file = fopen(fname,"rb");

  if(file)
    {
      uint bytes_read = fread(v,n*sizeof(double),1,file);
      fclose(file);
      if(bytes_read == 0)
        {   
	  return 1;
        }

      uint i;
      for(i=0; i<n; i++)
	{
	  if(!isfinite(v[i]) || isnan(v[i]))
	    {
	      fprintf(logfile,"Inf or NaN dected in binary data!\n");
	      return 1;
	    }
	}
      return 0;
    }
  else
    {
      fprintf(logfile,"File %s doesn't exist!\n",fname);
      return 1;
    }
}

static inline void linspace(double lo, double hi, uint n, double * v)
{
  double d=(hi-lo)/(n-1.0);  
  v[0]=lo;
  int i;
  for(i=1;i<n;i++)
    {
      v[i] = v[i-1]+d;
    }
  return;
}
#define LINSPACE(l,h,n,v) ( linspace((l),(h),(n),(double *)(v)) )

static inline void linspace_2d(const double * lo, const double * hi, uint n1, uint n2, double * v)
{
  uint j;
  for(j=0; j<n2; j++)
    {
      double d=(hi[j]-lo[j])/(n1-1.0);
      v[0*n2+j]=lo[j];
      int i;
      for(i=1;i<n1;i++)
	{
	  v[i*n2+j] = v[(i-1)*n2+j]+d;
	}
    }
  return;
}
#define LINSPACE_2D(l,h,n,m,v) ( linspace_2d( ( (double *)l),((double *)h),(n),(m),(double *)(v) ) )

static inline void logspace(double lo, double hi, uint n, double * v)
{
  uint i;
  linspace(log(lo),log(hi),n,v);
  for(i=0; i<n; i++)
    {
      v[i] = exp(v[i]);
    }
  return;
}
#define LOGSPACE(l,h,n,v) ( logspace((l),(h),(n),(double *)(v)) )

static inline void logspace_2d(const double * lo, const double * hi, uint n1, uint n2, double * v)
{
  uint i,j;

  double llo[n2], lhi[n2];
  for(j=0; j<n2; j++)
    {
      llo[j] = log(lo[j]);
      lhi[j] = log(hi[j]);
    }
  linspace_2d(llo,lhi,n1,n2,v);
  
  for(j=0; j<n1; j++)
    {
      for(i=0; i<n2; i++)
	{
	  v[i*n2+j] = exp(v[i*n2+j]);
	}
    }
  return;
}
#define LOGSPACE_2D(l,h,n,m,v) ( linspace_2d( ( (double *)l),((double *)h),(n),(m),(double *)(v) ) )

static inline char* concat(const char *s1, const char *s2)
{
  char *result = malloc(strlen(s1)+strlen(s2)+1);
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}

static inline double fabs_max(const double * x, uint n)
{
  uint i;
  double retval=0.0;
  for(i=0; i<n; i++)
    {
      if(fabs(x[i])>retval)
	{
	  retval=fabs(x[i]);
	}
    }
  return retval;
}

static inline uint fabs_imax(const double * x, uint n)
{
  uint i;
  double max=0.0;
  uint retval=0;
  for(i=0; i<n; i++)
    {
      if(fabs(x[i])>max)
	{
	  max=fabs(x[i]);
	  retval=i;
	}
    }
  return retval;
}

static inline double enorm(const gsl_vector * f) 
{
  double e2 = 0 ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    double fi= gsl_vector_get(f, i);
    e2 += fi * fi ;
  }
  return sqrt(e2);
}

static inline double vsum(const gsl_vector * f) 
{
  double sum = 0 ;
  size_t i, n = f->size ;
  for (i = 0; i < n ; i++) {
    double fi= gsl_vector_get(f, i);
    sum += fabs(fi) ;
  }
  return sum;
}

static inline void fprintf_vec(FILE * file, const double * v, uint n)
{
  uint i;
  for(i=0; i<n; i++)
    {
      fprintf(file," %0.4f",v[i]);
    }
  fprintf(file,"\n");
}

static inline void fprintf_mat_1(FILE * file, const double mat[][NS-1], uint n)
{
  uint i,j;
  for(i=0; i<n; i++)
    {
      for(j=0; j<NS-1; j++)
	{
	  fprintf(file," %0.4f",mat[i][j]);
	}
      fprintf(file,"\n");
      
    }
  fprintf(file,"\n");
}

static inline void fprintf_mat_2(FILE * file, const double mat[][NS], uint n)
{
  uint i,j;
  for(i=0; i<n; i++)
    {
      for(j=0; j<NS; j++)
	{
	  fprintf(file," %0.4f",mat[i][j]);
	}
      fprintf(file,"\n");
      
    }
  fprintf(file,"\n");
}


static inline void fprintf_3d_1(FILE * file, const double x[][NF][NS], uint n)
{
  uint i,j,k;
  for(i=0; i<n; i++)
    {
      for(j=0; j<NF; j++)
	{
	  for(k=0; k<NS; k++)
	    {
	      fprintf(file," %0.4f",x[i][j][k]);
	    }
	  if(j==(NF-1))
	    {
	      //fprintf(file,"\n");
	    }
	  else
	    {
	      fprintf(file,"\t");
	    }
	}
      fprintf(file,"\n");
    }
  fprintf(file,"\n");
}

static inline void fprintf_3d_2(FILE * file, const double x[][NS-1][NC], uint n)
{
  uint i,j,k;
  for(i=0; i<n; i++)
    {
      for(j=0; j<NS-1; j++)
	{
	  for(k=0; k<NC; k++)
	    {
	      fprintf(file," %0.4f",x[i][j][k]);
	    }
	  if(j==(NS-2))
	    {
	      //fprintf(file,"\n");
	    }
	  else
	    {
	      fprintf(file,"\t");
	    }
	}
      fprintf(file,"\n");
    }
  fprintf(file,"\n");
}

static inline void fprintf_3d_3(FILE * file, const double x[][NS][NS-1], uint n)
{
  uint i,j,k;
  for(i=0; i<n; i++)
    {
      for(j=0; j<NS; j++)
	{
	  for(k=0; k<NS-1; k++)
	    {
	      fprintf(file," %0.4f",x[i][j][k]);
	    }
	  if(j==(NS-1))
	    {
	      //fprintf(file,"\n");
	    }
	  else
	    {
	      fprintf(file,"\t");
	    }
	}
      fprintf(file,"\n");
    }
  fprintf(file,"\n");
}

/*
static inline void fprintf_3d_4(FILE * file, const double x[][NS][NC-1], uint n)
{
  uint i,j,k;
  for(i=0; i<n; i++)
    {
      for(j=0; j<NS; j++)
	{
	  for(k=0; k<(NC-1); k++)
	    {
	      fprintf(file," %0.4f",x[i][j][k]);
	    }
	  if(j==(NS-1))
	    {
	      //fprintf(file,"\n");
	    }
	  else
	    {
	      fprintf(file,"\t");
	    }
	}
      fprintf(file,"\n");
    }
  fprintf(file,"\n");
}
*/

#endif
