#ifndef __CALIBRATE_C__
#define __CALIBRATE_C__

#include "calibrate.h"

double export_participation_rate_target[NC][NS][NC-1];
double exit_rate_target[NC][NS][NC-1];
double exit_rate_target2[NC][NS][NC-1];
double ratio_target[NC][NS];
double tau_k_tmp[NC];
double sig_z_tmp[NC][NS];

uint homotopy_times = 15;

uint copy_params(params * dest, const params * src)
{
  dest->rho = src->rho;
  memcpy((double *)(dest->eps),(const double *)(src->eps),sizeof(double)*NC*NF*NS);
  memcpy((double *)(dest->G),(const double *)(src->G),sizeof(double)*NC);

  memcpy((double *)(dest->sig),(const double *)(src->sig),sizeof(double)*NC*NS);
  memcpy((double *)(dest->theta),(const double *)(src->theta),sizeof(double)*NC*NS*NC);
  memcpy((double *)(dest->H),(const double *)(src->H),sizeof(double)*NC*NS);

  memcpy((double *)(dest->zeta),(const double *)(src->zeta),sizeof(double)*NC*NS);
  memcpy((double *)(dest->mu),(const double *)(src->mu),sizeof(double)*NC*NS*NC);
  memcpy((double *)(dest->M),(const double *)(src->M),sizeof(double)*NC*NS);
  
  memcpy((double *)(dest->lam_va),(const double *)(src->lam_va),sizeof(double)*NC*NS);
  memcpy((double *)(dest->lam),(const double *)(src->lam),sizeof(double)*NC*NS*NS);
  memcpy((double *)(dest->B),(const double *)(src->B),sizeof(double)*NC*NS);
  memcpy((double *)(dest->alpha),(const double *)(src->alpha),sizeof(double)*NC*NS);
  memcpy((double *)(dest->A),(const double *)(src->A),sizeof(double)*NC*NS);

  dest->delta = src->delta;
  memcpy((double *)(dest->tauk),(const double *)(src->tauk),sizeof(double)*NC);
  dest->rss = src->rss;
  memcpy((double *)(dest->beta),(const double *)(src->beta),sizeof(double)*NC);
  dest->psi = src->psi;
  memcpy((double *)(dest->phi),(const double *)(src->phi),sizeof(double)*NC);
  memcpy((double *)(dest->lbar),(const double *)(src->lbar),sizeof(double)*NC);
  memcpy((double *)(dest->kk0),(const double *)(src->kk0),sizeof(double)*NC);
  memcpy((double *)(dest->b0),(const double *)(src->b0),sizeof(double)*NC);

  memcpy((double *)(dest->a_ts),(const double *)(src->a_ts),sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(dest->tau_m_ts),(const double *)(src->tau_m_ts),sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(dest->tau_f_ts),(const double *)(src->tau_f_ts),sizeof(double)*(NT+1)*NC*NS*NC);
  
  memcpy((double *)(dest->iomat),(const double *)(src->iomat),sizeof(double)*(NS*NC+2)*(NS*NC+NF*NC+1));
  memcpy((double *)(dest->r0),(const double *)(src->r0),sizeof(double)*NC);
  memcpy((double *)(dest->ii0),(const double *)(src->ii0),sizeof(double)*NC);
  memcpy((double *)(dest->ll0),(const double *)(src->ll0),sizeof(double)*NC);
  memcpy((double *)(dest->y0),(const double *)(src->y0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->va0),(const double *)(src->va0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->k0),(const double *)(src->k0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->l0),(const double *)(src->l0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->md0),(const double *)(src->md0),sizeof(double)*NC*NS*NS);
  memcpy((double *)(dest->ex0),(const double *)(src->ex0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->im0),(const double *)(src->im0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->nx0),(const double *)(src->nx0),sizeof(double)*NC*NC);
  memcpy((double *)(dest->c0),(const double *)(src->c0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->i0),(const double *)(src->i0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->m0),(const double *)(src->m0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->m02),(const double *)(src->m02),sizeof(double)*NC*NS*NC);
  memcpy((double *)(dest->q0),(const double *)(src->q0),sizeof(double)*NC*NS);
  memcpy((double *)(dest->q02),(const double *)(src->q02),sizeof(double)*NC*NS*NC);
  memcpy((double *)(dest->im02),(const double *)(src->im02),sizeof(double)*NC*NS*NC);
  memcpy((double *)(dest->ex02),(const double *)(src->ex02),sizeof(double)*NC*NS*NC);

  dest->etaM = src->etaM;
  dest->etaF = src->etaF;
  dest->etaK = src->etaK;
  dest->etaL = src->etaL;
  
  return 0;
}

uint set_nontargeted_params(params * p)
{
  // parameters common across countries
  p->rss = 0.02;
  p->delta = 0.06;

  if(k_adj_cost==0)
    {
      p->etaK = 0.0001;
    }
  else
    {
      p->etaK = 6.6;
    }
  if(l_adj_cost==0)
    {
      p->etaL = 6.5;
    }
  else
    {
      p->etaL = 1.0;
    }
  if(f_adj_cost==0)
    {
      p->etaF = 0.0;
    }
  else
    {
      p->etaF = 1.0;
    }
  if(m_adj_cost==0)
    {
      p->etaM = 0.0;
    }
  else
    {
      p->etaM = 3.5;
    }
  
  if(cobb_douglas_flag2)
    {
      p->rho = 1.0-1.0/0.9;
    }
  else
    {
      p->rho = 1.0-1.0/0.65;
    }
  p->psi = -1.0;

  SET_ALL_V(p->r0,NC,p->rss);
  SET_ALL_V(p->tauk,NC,0.25);
  SET_ALL_V(p->alpha,NC*NS,0.34);
  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);

  // need to adjust these!
  int i;
  for(i=0; i<NC; i++)
    {
      p->zeta[i][0] = 1.0-1.0/4.0;
      p->sig[i][0] = 1.0-1.0/4.0;
      
      p->zeta[i][1] = 1.0-1.0/4.0;
      p->sig[i][1] = 1.0-1.0/4.0;
      
      p->zeta[i][2] = 1.0-1.0/4.0;
      p->sig[i][2] = 1.0-1.0/4.0;
      
      p->zeta[i][3] = 1.0-1.0/4.0;
      p->sig[i][3] = 1.0-1.0/4.0;

      p->zeta[i][3] = 1.0-1.0/0.0;
      p->sig[i][3] = 1.0-1.0/0.0;

    }     

  return 0;
}

uint load_iomat(params * p)
{
  uint i, j, got;
  double tmp;
  FILE * file;

  file = fopen("../python/output/iomat.txt","rb");

  if(file)
    {
      for(i=0; i<(NS*NC+2); i++)
	{
	  for(j=0; j<(NS*NC+NF*NC+1); j++)
	    {
	      got = fscanf(file,"%lf",&tmp);
	      if(got != 1)
		{
		  fprintf(logfile,"Error reading IO matrix!\n");
		  fclose(file);
		  return 1;
		}
	      p->iomat[i][j] = tmp;
	    }
	}
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,"Error loading IO matrix!\n");
      return 1;
    }
}

void load_ts_params(params * p)
{
  uint i, s, t;
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->a_ts[0][i][s] = 1.0;
	  for(t=0; t<NT; t++)
	    {
	      p->a_ts[t+1][i][s] = 1.0;
	    }
	}
    }
}


void set_tariffs(params * p, double tau, uint scenario)
{
  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);

  uint s, t;
  
  if(scenario>=1)
    {
      for(t=TSHOCK; t<(NT+1); t++)
	{
	  for(s=0; s<NS-2; s++) // don't include services or construction
	    {
	      double x = tau;

	      // temporary tariffs turn off after 4 years
	      if(duration_flag==1 && t>=4)
		x=0.0;

	      // turn off tariffs for nontargeted sectors
	      if(target_sector_flag==0 && s==DNS)
		x=0.0;

	      if(target_sector_flag==1 && s==UPS)
		x=0.0;

	      // tariffs on China
	      if(target_country_flag==0 || target_country_flag==2) 
		{
		  p->tau_m_ts[t][USA][s][CHN] = x;
		  p->tau_f_ts[t][USA][s][CHN] = x;

		  if(retaliation_flag)
		    {
		      p->tau_m_ts[t][CHN][s][USA] = x;
		      p->tau_f_ts[t][CHN][s][USA] = x;
		    }
		}

	      // tariffs on ROW
	      if(target_country_flag==1 || target_country_flag==2)
		{
		  p->tau_m_ts[t][USA][s][ROW] = x;
		  p->tau_f_ts[t][USA][s][ROW] = x;

		  if(retaliation_flag)
		    {
		      p->tau_m_ts[t][ROW][s][USA] = x;
		      p->tau_f_ts[t][ROW][s][USA] = x;
		    }
		}
	    }
	}
    }
}

uint store_base_period_values(params * p)
{
  double mkt_clear_tol = 1.0e-7;
  uint varow = NC*NS;
  uint gorow = NC*NS+1;


  SET_ALL_V(p->y0,NC*NS,0.0);
  SET_ALL_V(p->va0,NC*NS,0.0);
  SET_ALL_V(p->k0,NC*NS,0.0);
  SET_ALL_V(p->l0,NC*NS,0.0);
  SET_ALL_V(p->md0,NC*NS*NS,0.0);
  SET_ALL_V(p->m0,NC*NS,0.0);
  SET_ALL_V(p->m02,NC*NS*NC,0.0);
  SET_ALL_V(p->q0,NC*NS,0.0);
  SET_ALL_V(p->q02,NC*NS*NC,0.0);
  SET_ALL_V(p->ex0,NC*NC,0.0);
  SET_ALL_V(p->im0,NC*NC,0.0);
  SET_ALL_V(p->nx0,NC*NC,0.0);
  SET_ALL_V(p->c0,NC*NS,0.0);
  SET_ALL_V(p->i0,NC*NS,0.0);
  SET_ALL_V(p->ii0,NC,0.0);

  uint i, s, j, r;

  for(i=0; i<NC; i++)
    {
      uint ccol = NC*NS+i;
      uint icol = NC*NS+NC+i;

      //double rdky = p->alpha[i][0]*(p->va0[i][0]+p->va0[i][1]);
      //double dky = p->delta*p->kk0[i];
      //p->tauk[i] = 1.0 - ( (dky+p->r0[i]*p->kk0[i])/rdky );
      //double rky = rdky - dky;
      //p->r0[i] = ((1.0-p->tauk[i])*rdky-dky)/p->kk0[i];

      for(s=0; s<NS; s++)
	{  
	  // first get value added and factors
	  uint scol = i*NS + s;
	  p->y0[i][s] = p->iomat[gorow][scol];
	  p->va0[i][s] = p->iomat[varow][scol];

	  p->l0[i][s] = (1.0 - p->alpha[i][s]) * p->va0[i][s];
	  //p->k0[i][s] = p->alpha[i][s] * p->va0[i][s] / ((p->r0[i] + p->delta) / (1.0 - p->tauk[i]));

	  // now get demand for products from different source countries and sectors 
	  for(j=0; j<NC; j++)
	    {
	      p->c0[i][s] = p->c0[i][s] + p->iomat[j*NS+s][ccol];
	      p->i0[i][s] = p->i0[i][s] + p->iomat[j*NS+s][icol];
	      p->q02[i][s][j] = p->iomat[j*NS+s][ccol] + p->iomat[j*NS+s][icol];

	      for(r=0; r<NS; r++)
		{
		  uint rcol = i*NS + r;
		  p->m02[i][s][j] = p->m02[i][s][j] + p->iomat[j*NS+s][rcol];
		  p->md0[i][r][s] = p->md0[i][r][s] + p->iomat[j*NS+s][rcol];
		  p->md0_[i][r][s] = p->md0_[i][r][s] + p->iomat[j*NS+s][rcol];
		}
	    }
	  p->q0[i][s] = sum(p->q02[i][s],NC);
	  p->m0[i][s] = sum(p->m02[i][s],NC);

	}

      p->ll0[i] = sum(p->l0[i],NS);
      p->ii0[i] = sum(p->i0[i],NS);

      p->tauk[i] = 1.0 - ((p->r0[i]+p->delta)*(p->ii0[i]/p->delta))/(p->alpha[i][0]*SUM(p->va0[i],NS));
      for(s=0; s<NS; s++)
	{
	  p->k0[i][s] = p->alpha[i][s] * p->va0[i][s] / ((p->r0[i] + p->delta) / (1.0 - p->tauk[i]));
	  if(p->k0[i][s]<0.0)
	    {
	      fprintf(logfile,"negative capital for country/sector %d/%d\n",i,s);
	      printf("%0.4f %0.4f %0.4f %0.4f\n",p->alpha[i][s],p->va0[i][s],p->r0[i] + p->delta,1.0 - p->tauk[i]);
	      return 1;
	    }
	}
      p->kk0[i] = sum(p->k0[i],NS);

    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(j != i)
		{
		  p->im02[i][s][j] = p->q02[i][s][j] + p->m02[i][s][j];
		  p->ex02[i][s][j] = p->q02[j][s][i] + p->m02[j][s][i];
		  p->nx02[i][s][j] = p->ex02[i][s][j] - p->im02[i][s][j];
		}
	      else
		{
		  p->im02[i][s][j] = 0.0;
		  p->ex02[i][s][j] = 0.0;
		  p->nx02[i][s][j] = 0.0;
		}
	    }
	  p->im0[i][j] = p->im02[i][0][j] + 
	    p->im02[i][1][j] + 
	    p->im02[i][2][j] + 
	    p->im02[i][3][j];
	  
	  p->ex0[i][j] = p->ex02[i][0][j] + 
	    p->ex02[i][1][j] + 
	    p->ex02[i][2][j] + 
	    p->ex02[i][3][j];

	  p->nx0[i][j] = p->nx02[i][0][j] + 
	    p->nx02[i][1][j] + 
	    p->nx02[i][2][j] + 
	    p->nx02[i][3][j];
	}
    }

  double tmp=0.0;
  for(i=0; i<NC; i++)
    {
      tmp = (sum(p->va0[i],NS) - (sum(p->q0[i],NS) + sum(p->ex0[i],NC) - sum(p->im0[i],NC)))/sum(p->va0[i],NS);
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,"GDP != C+I+NX for country %d, error = %f\n",i,tmp);
	  return 1;
	}

      for(s=0; s<NS; s++)
	{
	  tmp = p->y0[i][s];
	  for(j=0; j<NC; j++)
	    {
	      tmp = tmp - p->q02[j][s][i] - p->m02[j][s][i];
	    }
	  tmp = tmp/p->y0[i][s];
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,"supply != demand for country/sector %d/%d, error = %f\n",i,s,tmp);
	      return 1;
	    }
	}

      for(s=0; s<NS; s++)
	{
	  tmp = p->y0[i][s] - (p->va0[i][s] + SUM(p->md0[i][s],NS));
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,"go != va + m for country/sector %d/%d, error = %f\n",i,s,tmp);
	      return 1;
	    }

	  tmp = p->m0[i][s] - SUM(p->m02[i][s],NC);
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,"m != sum(m2) for country/sector %d/%d, error = %f\n",i,s,tmp);
	      return 1;
	    }
	}
    }

  // from IMF BoP dataset
  double usa_gdp_usdmm = 27000000;
  double usa_nfa_usdmm = -19853153.3747712;
  double chn_nfa_usdmm = 2908203;

  p->b0[0] = 100*usa_nfa_usdmm/usa_gdp_usdmm;
  p->b0[1] = 100*chn_nfa_usdmm/usa_gdp_usdmm;
  p->b0[2] = -sum(p->b0,2);

  printf("%0.4f %0.4f %0.4f\n",p->b0[0],p->b0[1],p->b0[2]);

  SET_ALL_V(ratio_target,NC*NS,0.6);

  return 0;

}

uint calibrate_prod_params(params * p)
{
  uint i,s,r;
  double tmp;
  
  SET_ALL_V(p->lam,NC*NS*NS,0.0);
  SET_ALL_V(p->A,NC*NS,0.0);
  SET_ALL_V(p->B,NC*NS,1.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  p->A[i][s] = p->va0[i][s] / ( pow(p->k0[i][s],p->alpha[i][s])*pow(p->l0[i][s],1.0-p->alpha[i][s]) );
	  p->lam_va[i][s] = p->va0[i][s] / p->y0[i][s];
	  if(p->lam_va[i][s]<0.0)
	    {
	      printf("Negative lam_va! i/s = %d/%d\n",i,s);
	      return 1;
	    }
	  
	  for(r=0; r<NS-1; r++)
	    {
	      p->lam[i][s][r] = p->md0[i][s][r] / p->y0[i][s];
	    }

	  if(cobb_douglas_flag)
	    {
	      p->B[i][s] = p->y0[i][s]/prod_go(p->va0[i][s],p->md0[i][s],p->lam_va[i][s],p->lam[i][s]);
	    }
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  tmp = p->B[i][s]*prod_go(p->va0[i][s],p->md0[i][s],p->lam_va[i][s],p->lam[i][s])
	    - p->y0[i][s];
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"prod_go != y0 for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }

	  tmp = prod_va(p->k0[i][s],p->l0[i][s],p->A[i][s],p->alpha[i][s]) - p->va0[i][s];
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"prod_va != va0 for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }

	  tmp = p->y0[i][s] - p->va0[i][s] - sum(p->md0[i][s],NS);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"nonzero profits for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }

	  if(cobb_douglas_flag==0)
	    {
	      tmp = (1.0-sum(p->lam[i][s],NS)) * (1.0-p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
		pow(p->k0[i][s],p->alpha[i][s]) * pow(p->l0[i][s],-(p->alpha[i][s])) - 1.0;
	    }
	  else
	    {
	      double tmp = pow(p->A[i][s]*pow(p->k0[i][s],p->alpha[i][s]),p->lam_va[i][s])*
		pow(p->l0[i][s],p->lam_va[i][s]*(1.0-p->alpha[i][s])-1.0);
	      uint r;
	      for(r=0; r<NS-1; r++)
		{
		  tmp = tmp * pow(p->md0[i][s][r],p->lam[i][s][r]);
		}
	      tmp =  (1.0-p->alpha[i][s]) * p->B[i][s] * p->lam_va[i][s] * tmp - 1.0;
	      
	    }
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"labor FOC for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }

	  if(cobb_douglas_flag==0)
	    {
	      tmp = (1.0-sum(p->lam[i][s],NS-1)) * (p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
		pow(p->k0[i][s],p->alpha[i][s]-1.0) * pow(p->l0[i][s],1.0-p->alpha[i][s]) - 
		(p->r0[i] + p->delta)/(1.0-p->tauk[i]);
	    }
	  else
	    {
	      double tmp = pow(p->A[i][s]*pow(p->l0[i][s],1.0-p->alpha[i][s]),p->lam_va[i][s])*
		pow(p->k0[i][s],p->lam_va[i][s]*p->alpha[i][s]-1.0);

	      uint r;
	      for(r=0; r<NS-1; r++)
		{
		  tmp = tmp * pow(p->md0[i][s][r],p->lam[i][s][r]);
		}
	      tmp =  p->alpha[i][s] * p->A[i][s] * p->B[i][s] * p->lam_va[i][s] * tmp
		-(p->r0[i] + p->delta)/(1.0-p->tauk[i]);
	    }
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"capital FOC for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }
	}
    }

  return 0;
}

uint calibrate_fin_params(params * p)
{
  uint i,s,j,jj,cnt,idx;
  double tmp;
  double tmp1[NC];
  double tmp2[NS];

  SET_ALL_V(p->mu,NC*NS*NC,0.0);
  SET_ALL_V(p->M,NC*NS,0.0);
  SET_ALL_V(p->G,NC,0.0);
  SET_ALL_V(p->H,NC*NS,0.0);
  SET_ALL_V(p->eps,NC*NF*NS,0.0);
  SET_ALL_V(p->theta,NC*NS*NC,0.0);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS-1; s++) // we don't do this for construction!
	{
	  idx=USA;
	  for(j=0; j<NC; j++)
	    {
	      tmp1[j] = pow(p->m02[i][s][j]/p->m02[i][s][idx],1.0 - p->zeta[i][s]);
	    }
	  p->mu[i][s][idx] = 1.0/sum(tmp1,NC);
	  cnt=0;
	  for(j=0; j<NC; j++)
	    {
	      cnt=cnt+1;
	      if(j != idx)
		{
		  if(cnt<NC)
		    {
		      p->mu[i][s][j] = p->mu[i][s][idx]*tmp1[j];
		    }
		  else
		    {
		      p->mu[i][s][j] = 1.0 - sum(p->mu[i][s],NC-1);
		    }
		}
	    
	    }
	  tmp = pow( DOT_PROD_EX(p->m02[i][s], p->mu[i][s], NC, p->zeta[i][s]), 1.0/p->zeta[i][s] );
	  p->M[i][s] = p->m0[i][s]/tmp;
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS-1; s++)  // we don't do this for construction!
	{
	  idx=USA;
	  for(j=0; j<NC; j++)
	    {
	      tmp1[j] = pow(p->q02[i][s][j]/p->q02[i][s][idx],1.0 - p->sig[i][s]);
	    }
	  p->theta[i][s][idx] = 1.0/sum(tmp1,NC);
	  cnt=0;
	  for(j=0; j<NC; j++)
	    {
	      cnt=cnt+1;
	      if(j != idx)
		{
		  if(cnt<NC)
		    {
		      p->theta[i][s][j] = p->theta[i][s][idx]*tmp1[j];
		    }
		  else
		    {
		      p->theta[i][s][j] = 1.0 - sum(p->theta[i][s],NC-1);
		    }
		}
	    
	    }
	  tmp = pow( DOT_PROD_EX(p->q02[i][s],p->theta[i][s],NC,p->sig[i][s]), 1.0/p->sig[i][s] );
	  p->H[i][s] = p->q0[i][s]/tmp;
	}

      p->theta[i][CNS][i]=1.0;
      p->H[i][CNS] = 1.0;      
    }

  idx=SVC;
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS-1; s++)  // we don't do this for construction!
	{
	  tmp2[s] = pow(p->c0[i][s]/p->c0[i][idx],1.0 - p->rho);
	}
      p->eps[i][0][idx] = 1.0/sum(tmp2,NS-1);
      cnt=0;
      for(s=0; s<NS-1; s++) // we don't do this for construction!
	{
	  cnt=cnt+1;
	  if(s != idx)
	    {
	      if(cnt<NS-1)  // we don't do this for construction!
		{
		  p->eps[i][0][s] = p->eps[i][0][idx]*tmp2[s];
		}
	      else
		{
		  p->eps[i][0][s] = 1.0 - sum(p->eps[i][0],NS-1);
		}
	    }
	}
    }

  for(i=0; i<NC; i++) 
    {
      for(s=0; s<NS; s++) // this is the only one where we do it for construction
	{
	  p->eps[i][1][s] = p->i0[i][s] / p->ii0[i];
	}
      p->G[i] = p->ii0[i]/ ( pow(p->i0[i][0],p->eps[i][1][0]) * 
			     pow(p->i0[i][1],p->eps[i][1][1]) *
			     pow(p->i0[i][2],p->eps[i][1][2]) * 
			     pow(p->i0[i][3],p->eps[i][1][3]));
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS-1; s++)
	{
	  tmp = p->m0[i][s] - prod_m(p->m02[i][s], p->M[i][s], p->mu[i][s], p->zeta[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"Intermediate Armington production function for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }

	  tmp = p->q0[i][s] - prod_q(p->q02[i][s],p->H[i][s],p->theta[i][s],p->sig[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"Final Armington production function for country/sector %d/%d, error = %f",i,s,tmp);
	      return 1;
	    }

	  for(j=0; j<NC; j++)
	    {
	      
	      tmp = 1.0 - p->mu[i][s][j] * pow(p->M[i][s],p->zeta[i][s]) * 
		pow(p->m0[i][s]/p->m02[i][s][j],1.0-p->zeta[i][s]);
	      
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,"Intermediate Armington FOC for country/sectorcountry %d/%d/%d, error = %f", i,s,j,tmp);
		  return 1;
		}
	      
	      tmp = 1.0 - p->theta[i][s][j] * pow(p->H[i][s],p->sig[i][s]) * 
		pow(p->q0[i][s]/p->q02[i][s][j],1.0-p->sig[i][s]);
	      
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,"Final Armington FOC for country/sector/country %d/%d/%d, error = %f", i,s,j,tmp);
		  return 1;
		}
	    }

	  for(j=0; j<NC; j++)
	    {
	      for(jj=0; jj<NC; jj++)
		{
		  tmp = 1.0 - (p->mu[i][s][j] / p->mu[i][s][jj]) * 
		    pow(p->m02[i][s][jj]/p->m02[i][s][j],1.0-p->zeta[i][s]);
		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,"Intermediate Armington FOC v2 for country/sector/country/country %d/%d/%d/%d, error = %f",
			      i,s,j,j,tmp);
		      return 1;
		    }

		  tmp = 1.0 - (p->theta[i][s][j] / p->theta[i][s][jj]) * 
		    pow(p->q02[i][s][jj]/p->q02[i][s][j],1.0-p->sig[i][s]);
		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,"Final Armington FOC v2 for country/sector/country/country %d/%d/%d/%d, error = %f",
			      i,s,j,jj,tmp);
		      return 1;
		    }
		}
	    }

	  for(s=0; s<(NS-2); s++)
	    {
	      tmp = muc(p->c0[i],p->ll0[i],p->lbar[i],p->eps[i][0],p->rho, p->phi[i], p->psi, s) /
	      muc(p->c0[i],p->ll0[i],p->lbar[i],p->eps[i][0],p->rho, p->phi[i], p->psi, NS-2) - 1.0;
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,"HH intratemp FOC 1 for country %d, sector %d, error = %f",i,s,tmp);
		  return 1;
		}
	    }

	  tmp = p->ii0[i] - prod_inv(p->i0[i],p->eps[i][1],p->G[i]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,"ii0 != prod_inv for country %d, error = %f",i,tmp);
	      return 1;
	    }

	  for(s=0; s<NS; s++)
	    {
	      tmp = 1.0 - p->eps[i][1][s]*(p->ii0[i]/p->i0[i][s]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,"Investment FOC for country/sector %d/%d, error = %f",i,s,tmp);
		  return 1;
		}
	    }
	}
    }

  return 0;
  
}

uint calibrate_hh_params(params * p)
{
  uint i, s;
  double tmp1[NS];
  
  for(i=0; i<NC; i++)
    {
      //tmp = pow(p->c0[i][0]/p->c0[i][1],1.0 - p->rho);
      //p->eps[i][0][0] = tmp/(1.0+tmp);
      //p->eps[i][0][1] = 1.0 - p->eps[i][0][0];

      p->lbar[i] = 3.0 * p->ll0[i];
      p->phi[i] = 1.0;

      for(s=0; s<NS-1; s++) // we don't do this for construction!
	{
	  tmp1[s] = p->c0[i][s];
	}
      // note: beta = (1.0+rss)/gbgp^(phi*psi-1) > 1
      // but, as long as beta*gbgp^(phi*psi) < 1, we can calculate welfare just fine
      p->beta[i] = muc(p->c0[i],p->ll0[i],p->lbar[i],p->eps[i][0],p->rho, p->phi[i], p->psi, 0) /
	muc(tmp1,p->ll0[i],p->lbar[i],p->eps[i][0],p->rho, p->phi[i], p->psi, 0) / (1.0 + p->rss);
    }

  return 0;
}

uint calibrate()
{
  params * p = &(ppp0[0]);
  
  if(set_nontargeted_params(p))
    {
      return 1;
    }

  if(load_iomat(p))
    {
      return 1;
    }

  load_ts_params(p);


  if(store_base_period_values(p))
    {
      return 1;
    }


  if(calibrate_prod_params(p))
    {
      return 1;
    }


  if(calibrate_fin_params(p))
    {
      return 1;
    }


  if(calibrate_hh_params(p))
    {
      return 1;
    }


  write_params();
  
  uint it;
  for(it=0; it<NTH; it++)
    {
      if(it>0 && copy_params(&(ppp0[it]),&(ppp0[0])))
	{
	  fprintf(logfile, "\nFailed to copy ppp0!\n");
	  return 1;
	}
      if(copy_params(&(ppp1[it]),&(ppp0[0])))
	{
	  fprintf(logfile, "\nFailed to copy ppp1!\n");
	  return 1;
	}
    }

  return 0;
}

uint write_params()
{
  const params * p = &(ppp0[0]);
  
  FILE * file = fopen("output/params.txt","wb");
  if(file)
    {
      fprintf(file,"Scalar parameters:\n");
      fprintf(file,"rho: %0.4f\n",p->rho);
      fprintf(file,"delta: %0.4f\n",p->delta);
      fprintf(file,"rss: %0.4f\n",p->rss);
      fprintf(file,"psi: %0.4f\n",p->psi);

      fprintf(file,"\nVECTOR PARAMETERS (1 x NC):\n");

      fprintf(file,"G:");
      fprintf_vec(file,p->G,NC);

      fprintf(file,"beta:");
      fprintf_vec(file,p->beta,NC);

      fprintf(file,"phi:");
      fprintf_vec(file,p->phi,NC);

      fprintf(file,"lbar:");
      fprintf_vec(file,p->lbar,NC);

      fprintf(file,"kk0:");
      fprintf_vec(file,p->kk0,NC);

      fprintf(file,"tauk:");
      fprintf_vec(file,p->tauk,NC);

      fprintf(file,"b0:");
      fprintf_vec(file,p->b0,NC);

      fprintf(file,"\nMATRIX PARAMETERS (NC x NS):\n\n");

      fprintf(file,"sig:\n");
      fprintf_mat(file,p->sig,NC);

      fprintf(file,"H:\n");
      fprintf_mat(file,p->H,NC);

      fprintf(file,"zeta:\n");
      fprintf_mat(file,p->zeta,NC);

      fprintf(file,"M:\n");
      fprintf_mat(file,p->M,NC);

      fprintf(file,"lam_va:\n");
      fprintf_mat(file,p->lam_va,NC);

      fprintf(file,"alpha:\n");
      fprintf_mat(file,p->alpha,NC);

      fprintf(file,"A:\n");
      fprintf_mat(file,p->A,NC);
	
      fprintf(file,"\n3D PARAMETERS:\n\n");

      fprintf(file,"eps (NC x 2 x NS):\n");
      fprintf_3d_1(file,p->eps,NC);

      fprintf(file,"theta (NC x NS x NC):\n");
      fprintf_3d_2(file,p->theta,NC);

      fprintf(file,"mu (NC x NS x NC):\n");
      fprintf_3d_2(file,p->mu,NC);

      fprintf(file,"lam (NC x NS x NS):\n");
      fprintf_3d_3(file,p->lam,NC);

      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,"Error opening file to write parameters!\n");
      return 1;
    }
}

#endif
