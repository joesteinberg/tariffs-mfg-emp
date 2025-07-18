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
  memcpy((double *)(dest->H2),(const double *)(src->H2),sizeof(double)*NC*NS*NC);

  memcpy((double *)(dest->zeta),(const double *)(src->zeta),sizeof(double)*NC*NS);
  memcpy((double *)(dest->mu),(const double *)(src->mu),sizeof(double)*NC*NS*NC);
  memcpy((double *)(dest->M),(const double *)(src->M),sizeof(double)*NC*NS);
  memcpy((double *)(dest->M2),(const double *)(src->M2),sizeof(double)*NC*NS*NC);
  
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
  memcpy((double *)(dest->ntb_m_ts),(const double *)(src->ntb_m_ts),sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(dest->ntb_f_ts),(const double *)(src->ntb_f_ts),sizeof(double)*(NT+1)*NC*NS*NC);
  
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
  //memcpy((double *)(dest->tariff_data),(const double *)(src->tariff_data),sizeof(double)*NC*NS*NC);

  dest->etaM = src->etaM;
  dest->etaF = src->etaF;
  dest->etaK = src->etaK;
  dest->etaL = src->etaL;
  
  dest->eta = src->eta;
  memcpy((double *)(dest->sig_z),(const double *)(src->sig_z),sizeof(double)*NC*NS);
  memcpy((double *)(dest->kappa0),(const double *)(src->kappa0),sizeof(double)*NC*NS*(NC-1));
  memcpy((double *)(dest->kappa1),(const double *)(src->kappa1),sizeof(double)*NC*NS*(NC-1));
  memcpy((double *)(dest->n0),(const double *)(src->n0),sizeof(double)*NC*NS*(NC-1));
  memcpy((uint *)(dest->Ji),(const uint *)(src->Ji),sizeof(uint)*NC*(NC-1));

  return 0;
}

uint set_nontargeted_params(params * p)
{
  // parameters common across countries
  p->rss = 0.02;
  p->delta = 0.06;
  p->eta = 5.0;

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

  if(no_k_flag)
    {
      SET_ALL_V(p->alpha,NC*NS,0.0);
    }
  else
    {
      SET_ALL_V(p->alpha,NC*NS,0.34);
    }
  
  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);

  if(sym_te_flag)
    {
      SET_ALL_V(p->zeta,NC*NS,1.0-1.0/5.0);
      SET_ALL_V(p->sig,NC*NS,1.0-1.0/5.0);
    }
  else if(ltp_flag)
    {
      int i;
      for(i=0; i<NC; i++)
	{
	  p->zeta[i][0] = 1.0-1.0/14.78;
	  p->sig[i][0] = 1.0-1.0/14.78;

	  p->zeta[i][1] = 1.0-1.0/33.75;
	  p->sig[i][1] = 1.0-1.0/33.75;

	  p->zeta[i][2] = 1.0-1.0/0.94;
	  p->sig[i][2] = 1.0-1.0/0.94;

	  p->zeta[i][3] = 1.0-1.0/7.5;
	  p->sig[i][3] = 1.0-1.0/7.5;

	  p->zeta[i][4] = 1.0-1.0/5.00;
	  p->sig[i][4] = 1.0-1.0/5.00;
	}
    }
  else
    {
      int i;
      for(i=0; i<NC; i++)
	{
	  p->zeta[i][0] = 1.0-1.0/8.11;
	  p->sig[i][0] = 1.0-1.0/8.11;

	  p->zeta[i][1] = 1.0-1.0/31.82;
	  p->sig[i][1] = 1.0-1.0/31.82;

	  p->zeta[i][2] = 1.0-1.0/0.88;
	  p->sig[i][2] = 1.0-1.0/0.88;

	  p->zeta[i][3] = 1.0-1.0/5.17;
	  p->sig[i][3] = 1.0-1.0/5.17;

	  p->zeta[i][4] = 1.0-1.0/5.00;
	  p->sig[i][4] = 1.0-1.0/5.00;
	}     
    }

  uint i;
  for(i=0; i<NC; i++)
    {
      if(i==0)
	{
	  p->Ji[i][0]=1;
	  p->Ji[i][1]=2;
	  p->Ji[i][2]=3;
	  p->Ji[i][3]=4;
	  p->Ji[i][4]=5;
	}
      else if(i==1)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=2;
	  p->Ji[i][2]=3;
	  p->Ji[i][3]=4;
	  p->Ji[i][4]=5;
	}
      else if(i==2)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=1;
	  p->Ji[i][2]=3;
	  p->Ji[i][3]=4;
	  p->Ji[i][4]=5;
	}
      else if(i==3)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=1;
	  p->Ji[i][2]=2;
	  p->Ji[i][3]=4;
	  p->Ji[i][4]=5;
	}
      else if(i==4)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=1;
	  p->Ji[i][2]=2;
	  p->Ji[i][3]=3;
	  p->Ji[i][4]=5;
	}
      else if(i==5)
	{
	  p->Ji[i][0]=0;
	  p->Ji[i][1]=1;
	  p->Ji[i][2]=2;
	  p->Ji[i][3]=3;
	  p->Ji[i][4]=4;
	}
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
		  fprintf(logfile,KRED "Error reading IO matrix!\n" RESET);
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
      fprintf(logfile,KRED "Error loading IO matrix!\n" RESET);
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

/*
uint load_tariffs(params * p)
{
  SET_ALL_V(p->tariff_data,NC*NS*NC,0.0);

  uint i, j, s, got;
  double tmp;
  FILE * file;
  if(old_mfn_flag==1)
    {
      file = fopen("../python/output/tariffs_old.txt","rb");
    }
  else if(us_ht_flag==1 || us_ht_flag2==1)
    {
      file = fopen("../python/output/tariffs_alt1.txt","rb");
    }
  else
    {
      file = fopen("../python/output/tariffs.txt","rb");
    }

  if(file)
    {
      for(i=0; i<NC-1; i++)
	{
	  for(j=0; j<NC; j++)
	    {
	      if(j!=i)
		{
		  for(s=0; s<NS-1; s++)
		    {
		      int tmpi, tmps, tmpj;
		      got = fscanf(file,"%d %d %d %lf",&tmpi,&tmps,&tmpj,&tmp);
		      if(got != 4)
			{
			  fprintf(logfile,KRED "Error reading tariffs file!\n" RESET);
			  fclose(file);
			  return 1;
			}
		      else
			{
			  p->tariff_data[i][s][j] = tmp/100.0;
			}
		    }
		}
	    }	    
	}
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error loading tariffs file!\n" RESET);
      return 1;
    }
}
*/

void set_tariffs(params * p, uint scenario)
{
  SET_ALL_V(p->tau_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->tau_f_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_m_ts,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(p->ntb_f_ts,(NT+1)*NC*NS*NC,0.0);

  uint s, t;
  
  if(scenario>=1)
    {
      for(t=TNAFTA; t<(NT+1); t++)
	{
	  for(s=0; s<NS-1; s++) // don't include services
	    {
	      // USA to Canada
	      if(scenario2%2==0 || t<4)
		{
		  double x=0.0;
		  
		  if(t<=2)
		    {
		      x = 0.25*0.35;
		  
		      if(s==CARS)
			x = 0.25*0.3;
		      else if(s==OIL)
			x=0.1;
		    }
		  else
		    {
		      x = 0.12*0.35;
		      
		      if(s==CARS)
			x = 0.25*0.3;
		      else if(s==OIL)
			x=0.0;
		    }
		    
		      
		  p->tau_m_ts[t][USA][s][CAN] = x;
		  p->tau_f_ts[t][USA][s][CAN]= x;
		}
	      
	      // USA to Mexico
	      if(scenario2%2==0 || t<4)
		{
		  double x=0.0;
		  
		  if(t<=2)
		    {
		      x = 0.25*0.35;
		  
		      if(s==CARS)
			x = 0.25*0.3;
		      else if(s==OIL)
			x=0.1;
		    }
		  else
		    {
		      x = 0.1*0.35;
		      
		      if(s==CARS)
			x = 0.25*0.3;
		      else if(s==OIL)
			x=0.1;
		    }
		    	  
		  p->tau_m_ts[t][USA][s][MEX] = x;
		  p->tau_f_ts[t][USA][s][MEX]= x;
		}

	      // USA to China
	      if(scenario2%2==0 || t<4)
		{
		  p->tau_m_ts[t][USA][s][CHN] = 0.54;
		  p->tau_f_ts[t][USA][s][CHN]= 0.54;
		}

	      // USA to EU
	      if(scenario2%2==0 || t<4)
		{
		  //double x = s==CARS ? 0.25 : 0.1;
		  double x = 0.5;
		  p->tau_m_ts[t][USA][s][EU] = x;
		  p->tau_f_ts[t][USA][s][EU]= x;
		}

	      // USA to RW
	      if(scenario2%2==0 || t<4)
		{
		  double x = s==CARS ? 0.25 : 0.1;
		  p->tau_m_ts[t][USA][s][ROW] = x;
		  p->tau_f_ts[t][USA][s][ROW]= x;
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
	      if(no_k_flag)
		{
		  p->c0[i][s] += p->i0[i][s];
		  p->i0[i][s] = 0.0;
		}
	      p->q02[i][s][j] = p->iomat[j*NS+s][ccol] + p->iomat[j*NS+s][icol];

	      for(r=0; r<NS; r++)
		{
		  uint rcol = i*NS + r;
		  if(noio_flag)
		    {
		      p->m02[i][s][j] = 0.0;
		      p->md0[i][r][s] = 0.0;
		      p->md0_[i][r][s] = 0.0;
		    }
		  else
		    {
		      p->m02[i][s][j] = p->m02[i][s][j] + p->iomat[j*NS+s][rcol];
		      p->md0[i][r][s] = p->md0[i][r][s] + p->iomat[j*NS+s][rcol];
		      p->md0_[i][r][s] = p->md0_[i][r][s] + p->iomat[j*NS+s][rcol];
		    }
		}
	    }
	  p->q0[i][s] = sum(p->q02[i][s],NC);
	  p->m0[i][s] = sum(p->m02[i][s],NC);

	}

      p->ll0[i] = sum(p->l0[i],NS);
      p->ii0[i] = sum(p->i0[i],NS);

      if(!no_k_flag)
	{
	  p->tauk[i] = 1.0 - ((p->r0[i]+p->delta)*(p->ii0[i]/p->delta))/(p->alpha[i][0]*SUM(p->va0[i],NS));
	  for(s=0; s<NS; s++)
	    {
	      p->k0[i][s] = p->alpha[i][s] * p->va0[i][s] / ((p->r0[i] + p->delta) / (1.0 - p->tauk[i]));
	      if(p->k0[i][s]<0.0)
		{
		  fprintf(logfile,KRED "negative capital for country/sector %d/%d\n" RESET,i,s);
		  printf("%0.4f %0.4f %0.4f %0.4f\n",p->alpha[i][s],p->va0[i][s],p->r0[i] + p->delta,1.0 - p->tauk[i]);
		  return 1;
		}
	    }
	  p->kk0[i] = sum(p->k0[i],NS);
	}
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
	    p->im02[i][3][j] + 
	    p->im02[i][4][j];
	  
	  p->ex0[i][j] = p->ex02[i][0][j] + 
	    p->ex02[i][1][j] + 
	    p->ex02[i][2][j] + 
	    p->ex02[i][3][j] + 
	    p->ex02[i][4][j];

	  p->nx0[i][j] = p->nx02[i][0][j] + 
	    p->nx02[i][1][j] + 
	    p->nx02[i][2][j] + 
	    p->nx02[i][3][j] + 
	    p->nx02[i][4][j];
	}
    }

  double tmp=0.0;
  for(i=0; i<NC; i++)
    {
      tmp = (sum(p->va0[i],NS) - (sum(p->q0[i],NS) + sum(p->ex0[i],NC) - sum(p->im0[i],NC)))/sum(p->va0[i],NS);
      if(fabs(tmp)>mkt_clear_tol)
	{
	  fprintf(logfile,KRED "GDP != C+I+NX for country %d, error = %f\n" RESET,i,tmp);
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
	      fprintf(logfile,KRED "supply != demand for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }
	}

      for(s=0; s<NS; s++)
	{
	  tmp = p->y0[i][s] - (p->va0[i][s] + SUM(p->md0[i][s],NS));
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "go != va + m for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }

	  tmp = p->m0[i][s] - SUM(p->m02[i][s],NC);
	  if(fabs(tmp)>mkt_clear_tol)
	    {
	      fprintf(logfile,KRED "m != sum(m2) for country/sector %d/%d, error = %f\n" RESET,i,s,tmp);
	      return 1;
	    }
	}
    }

  
  // from WIOD

  if(fix_tb_flag2)
    {
      double xx = (1.0+p->rss)/p->rss;
      p->b0[0] = -SUM(p->nx0[0],NC) * xx;
      p->b0[1] = -SUM(p->nx0[1],NC) * xx;
      p->b0[2] = -SUM(p->nx0[2],NC) * xx;
      p->b0[3] = -SUM(p->nx0[3],NC) * xx;
      p->b0[4] = -SUM(p->nx0[4],NC) * xx;
    }
  else
    {
      // from IMF BoP dataset
      double usa_gdp_usdmm = 27000000;
      double usa_nfa_usdmm = -19853153.3747712;
      double can_nfa_usdmm = 1073984.7;
      double mex_nfa_usdmm = -758198;
      double chn_nfa_usdmm = 2908203;
      double eu_nfa_usdmm = 434695;

      p->b0[0] = 100*usa_nfa_usdmm/usa_gdp_usdmm;
      p->b0[1] = 100*can_nfa_usdmm/usa_gdp_usdmm;
      p->b0[2] = 100*mex_nfa_usdmm/usa_gdp_usdmm;
      p->b0[3] = 100*chn_nfa_usdmm/usa_gdp_usdmm;
      p->b0[4] = 100*eu_nfa_usdmm/usa_gdp_usdmm;

    }
  p->b0[5] = -sum(p->b0,5);

  printf("%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",p->b0[0],p->b0[1],p->b0[2],p->b0[3],p->b0[4],p->b0[5]);

  SET_ALL_V(ratio_target,NC*NS,0.6);
    
  if(nokappa)
    {
      SET_ALL_V(p->n0,NC*NS*(NC-1),1.0);
    }
  else
    {
      uint i,s,ii;
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {     
	      for(ii=0; ii<(NC-1); ii++)
		{
		  if(s==1)
		    {
		      export_participation_rate_target[i][s][ii]=1.0;
		      exit_rate_target[i][s][ii] = 0.0;
		    }
		  else
		    {
		      export_participation_rate_target[i][s][ii]=0.25;
		      exit_rate_target[i][s][ii] = 0.45;
		    }
		  p->n0[i][s][ii]=export_participation_rate_target[i][s][ii];
		}
	    }
	}
    }

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
	  if(!no_k_flag)
	    {
	      p->A[i][s] = p->va0[i][s] / ( pow(p->k0[i][s],p->alpha[i][s])*pow(p->l0[i][s],1.0-p->alpha[i][s]) );
	    }
	  else
	    {
	      p->A[i][s] = p->va0[i][s]/p->l0[i][s];
	    }
	  p->lam_va[i][s] = p->va0[i][s] / p->y0[i][s];
	  if(p->lam_va[i][s]<0.0)
	    {
	      printf("Negative lam_va! i/s = %d/%d\n",i,s);
	      return 1;
	    }
	  
	  for(r=0; r<NS; r++)
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
	  if(noio_flag==0)
	    {
	      tmp = p->B[i][s]*prod_go(p->va0[i][s],p->md0[i][s],p->lam_va[i][s],p->lam[i][s])
		- p->y0[i][s];
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "prod_go != y0 for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }

	  if(!no_k_flag)
	    {
	      tmp = prod_va(p->k0[i][s],p->l0[i][s],p->A[i][s],p->alpha[i][s]) - p->va0[i][s];
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "prod_va != va0 for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }

	  tmp = p->y0[i][s] - p->va0[i][s] - sum(p->md0[i][s],NS);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "nonzero profits for country/sector %d/%d, error = %f" RESET,i,s,tmp);
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
	      for(r=0; r<NS; r++)
		{
		  tmp = tmp * pow(p->md0[i][s][r],p->lam[i][s][r]);
		}
	      tmp =  (1.0-p->alpha[i][s]) * p->B[i][s] * p->lam_va[i][s] * tmp - 1.0;
	      
	    }
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "labor FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  if(!no_k_flag)
	    {
	      if(cobb_douglas_flag==0)
		{
		  tmp = (1.0-sum(p->lam[i][s],NS)) * (p->alpha[i][s]) * p->A[i][s]/p->lam_va[i][s] *
		    pow(p->k0[i][s],p->alpha[i][s]-1.0) * pow(p->l0[i][s],1.0-p->alpha[i][s]) - 
		    (p->r0[i] + p->delta)/(1.0-p->tauk[i]);
		}
	      else
		{
		  double tmp = pow(p->A[i][s]*pow(p->l0[i][s],1.0-p->alpha[i][s]),p->lam_va[i][s])*
		    pow(p->k0[i][s],p->lam_va[i][s]*p->alpha[i][s]-1.0);

		  uint r;
		  for(r=0; r<NS; r++)
		    {
		      tmp = tmp * pow(p->md0[i][s][r],p->lam[i][s][r]);
		    }
		  tmp =  p->alpha[i][s] * p->A[i][s] * p->B[i][s] * p->lam_va[i][s] * tmp
		    -(p->r0[i] + p->delta)/(1.0-p->tauk[i]);
		}
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "capital FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
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

  if(noio_flag==0)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      idx=3;
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
	      tmp = pow( DOT_PROD_EX(p->m02[i][s],p->mu[i][s],NC,p->zeta[i][s]), 1.0/p->zeta[i][s] );
	      p->M[i][s] = p->m0[i][s]/tmp;
	    }
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  idx=3;
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
    }

  idx=3;
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  tmp2[s] = pow(p->c0[i][s]/p->c0[i][idx],1.0 - p->rho);
	}
      p->eps[i][0][idx] = 1.0/sum(tmp2,NS);
      cnt=0;
      for(s=0; s<NS; s++)
	{
	  cnt=cnt+1;
	  if(s != idx)
	    {
	      if(cnt<NS)
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

  if(!no_k_flag)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      p->eps[i][1][s] = p->i0[i][s] / p->ii0[i];
	    }
	  p->G[i] = p->ii0[i]/ ( pow(p->i0[i][1],p->eps[i][1][1]) * 
				 pow(p->i0[i][2],p->eps[i][1][2]) *
				 pow(p->i0[i][3],p->eps[i][1][3]) * 
				 pow(p->i0[i][4],p->eps[i][1][4]));
	}
    }

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  if(noio_flag==0)
	    {
	      tmp = p->m0[i][s] - prod_m(p->m02[i][s],p->M[i][s],p->mu[i][s],p->zeta[i][s]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Intermediate Armington production function for country/sector %d/%d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }

	  tmp = p->q0[i][s] - prod_q(p->q02[i][s],p->H[i][s],p->theta[i][s],p->sig[i][s]);
	  if(fabs(tmp)>TINY)
	    {
	      fprintf(logfile,KRED "Final Armington production function for country/sector %d/%d, error = %f" RESET,i,s,tmp);
	      return 1;
	    }

	  for(j=0; j<NC; j++)
	    {
	      if(noio_flag==0)
		{
		  tmp = 1.0 - p->mu[i][s][j] * pow(p->M[i][s],p->zeta[i][s]) * 
		    pow(p->m0[i][s]/p->m02[i][s][j],1.0-p->zeta[i][s]);

		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,KRED "Intermediate Armington FOC for country/sectorcountry %d/%d/%d, error = %f" RESET, i,s,j,tmp);
		      return 1;
		    }
		}

	      tmp = 1.0 - p->theta[i][s][j] * pow(p->H[i][s],p->sig[i][s]) * 
		pow(p->q0[i][s]/p->q02[i][s][j],1.0-p->sig[i][s]);

	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "Final Armington FOC for country/sector/country %d/%d/%d, error = %f" RESET, i,s,j,tmp);
		  return 1;
		}

	    }

	  for(j=0; j<NC; j++)
	    {
	      for(jj=0; jj<NC; jj++)
		{
		  if(noio_flag==0)
		    {
		      tmp = 1.0 - (p->mu[i][s][j] / p->mu[i][s][jj]) * 
			pow(p->m02[i][s][jj]/p->m02[i][s][j],1.0-p->zeta[i][s]);
		      if(fabs(tmp)>TINY)
			{
			  fprintf(logfile,KRED "Intermediate Armington FOC v2 for country/sector/country/country %d/%d/%d/%d, error = %f" RESET,
				  i,s,j,j,tmp);
			  return 1;
			}
		    }

		  tmp = 1.0 - (p->theta[i][s][j] / p->theta[i][s][jj]) * 
		    pow(p->q02[i][s][jj]/p->q02[i][s][j],1.0-p->sig[i][s]);
		  if(fabs(tmp)>TINY)
		    {
		      fprintf(logfile,KRED "Final Armington FOC v2 for country/sector/country/country %d/%d/%d/%d, error = %f" RESET,
			      i,s,j,jj,tmp);
		      return 1;
		    }
		}
	    }

	  for(s=0; s<(NS-1); s++)
	    {
	      tmp = muc(p->c0[i],p->ll0[i],p->lbar[i],p->eps[i][0],p->rho, p->phi[i], p->psi, s) /
	      muc(p->c0[i],p->ll0[i],p->lbar[i],p->eps[i][0],p->rho, p->phi[i], p->psi, NS-1) - 1.0;
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "HH intratemp FOC 1 for country %d, sector %d, error = %f" RESET,i,s,tmp);
		  return 1;
		}
	    }

	  if(!no_k_flag)
	    {
	      tmp = p->ii0[i] - prod_inv(p->i0[i],p->eps[i][1],p->G[i]);
	      if(fabs(tmp)>TINY)
		{
		  fprintf(logfile,KRED "ii0 != prod_inv for country %d, error = %f" RESET,i,tmp);
		  return 1;
		}

	      for(s=0; s<NS; s++)
		{
		  if(s!=0) // agriculture share of investment is zero
		    {
		      tmp = 1.0 - p->eps[i][1][s]*(p->ii0[i]/p->i0[i][s]);
		      if(fabs(tmp)>TINY)
			{
			  fprintf(logfile,KRED "Investment FOC for country/sector %d/%d, error = %f" RESET,i,s,tmp);
			  return 1;
			}
		    }
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

      for(s=0; s<NS; s++)
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

uint calibrate_firm_params()
{
  params * p = &(ppp0[0]);
  
  uint i,s,r,j,ii;
  double MC[NC][NS];
  double Df[NC][NS][NC];
  double Dbf[NC][NS][NC];
  double Dtf[NC][NS][NC];
  double Dm[NC][NS][NC];
  double Dbm[NC][NS][NC];
  double Dtm[NC][NS][NC];
  double Zi[NC][NS];

  // initial guess for dispersion (George's China paper)
  if(nokappa)
    {
      SET_ALL_V(p->sig_z,NC*NS,0.6);
      //SET_ALL_V(p->sig_z,NC*NS,5.0);
    }
  else
    {
      SET_ALL_V(p->sig_z,NC*NS,0.4);
      //SET_ALL_V(p->sig_z,NC*NS,5.0);
    }

  // modify value added share
  /*
  if(!noio_flag)
    {
      if(cobb_douglas_flag)
	{
	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  double tmp = p->lam_va[i][s];
		  p->lam_va[i][s] = 1.0 - (p->eta/(p->eta-1.0))*SUM(p->md0_[i][s],NS)/p->y0[i][s];
		  double tmp2 = (1.0-p->lam_va[i][s])/(1.0-tmp);
		  for(r=0; r<NS; r++)
		    {
		      p->lam[i][s][r] = p->lam[i][s][r]*tmp2;
		    }
		}
	    }
	}
      else
	{
	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  double tmp = p->lam_va[i][s];
		  double sq = (p->eta/(p->eta-1.0))*SUM(p->md0_[i][s],NS)/p->y0[i][s];
		  double br = 0.0;
		  if(no_k_flag)
		    {
		      br = 1.0;
		    }
		  else
		    {
		      br = pow((p->r0[i]+p->delta)/(1.0-p->tauk[i])/p->alpha[i][s],p->alpha[i][s])*
			pow(1.0/(1.0-p->alpha[i][s]),1.0-p->alpha[i][s]);
		    }
		  p->lam_va[i][s] = (1.0-sq)/(1.0-sq+br*sq);

		  if(p->lam_va[i][s]<0.0)
		    {
		      printf("Negative lam_va! i/s = %d/%d\n",i,s);
		      printf("%0.4f %0.4f %0.4f\n",sq,(p->eta/(p->eta-1.0)),SUM(p->md0_[i][s],NS)/p->y0[i][s]);
		      return 1;
		    }

		  double tmp2 = (1.0-p->lam_va[i][s])/(1.0-tmp);
		  for(r=0; r<NS; r++)
		    {
		      p->lam[i][s][r] = p->lam[i][s][r]*tmp2;
		    }
		}
	    }
	}
    }
  else
    {
      SET_ALL_V(p->lam_va,NC*NS,1.0);
      SET_ALL_V(p->lam,NC*NS*NS,0.0);
    }
  */


  
  // initialize some other stuff
  SET_ALL_V(p->kappa0,NC*NS*(NC-1),0.0);
  SET_ALL_V(p->kappa1,NC*NS*(NC-1),0.0);

  // first calculate the domestic demand stuff
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  Zi[i][s] = exp(p->sig_z[i][s]*p->sig_z[i][s]*(p->eta-1.0)*(p->eta-1.0)/2.0);
	  //Zi[i][s] = p->sig_z[i][s]/(p->sig_z[i][s] + 1.0 - p->eta);

	  if(no_k_flag)
	    {
	      MC[i][s] = p->lam_va[i][s] + SUM(p->lam[i][s],NS);
	    }
	  else if(cobb_douglas_flag==0)
	    {
	      MC[i][s] = (p->lam_va[i][s]*(
					   pow((p->r0[i]+p->delta)/(1.0-p->tauk[i])/p->alpha[i][s],p->alpha[i][s]) *
					   pow(1.0/(1.0-p->alpha[i][s]),1.0-p->alpha[i][s])
					   )
			  ) + SUM(p->lam[i][s],NS);
		
	    }
	  else
	    {
	      MC[i][s] = (pow( (p->r0[i]+p->delta)/(1.0-p->tauk[i])/(p->alpha[i][s]*p->lam_va[i][s]),
			       p->lam_va[i][s]*p->alpha[i][s]) *
			  pow(1.0/(1.0-p->alpha[i][s])/p->lam_va[i][s],(1.0-p->alpha[i][s])*p->lam_va[i][s]));
	      if(!noio_flag)
		{
		  for(r=0; r<NS; r++)
		    {
		      MC[i][s] = MC[i][s] * pow(1.0/p->lam[i][s][r],p->lam[i][s][r]);
		    }
		}
	      //pow(1.0/(p->A[i][s]*p->lam_va[i][s]),p->lam_va[i][s]));
	    }


	  p->H2[i][s][i] = pow(p->eta/(p->eta-1.0),p->eta-1.0) *
	    pow(MC[i][s],p->eta-1.0) * (1.0/Zi[i][s]);
	  p->H2[i][s][i] = pow(p->H2[i][s][i],1.0/(p->eta-1.0));

	  Df[i][s][i] = pow(p->H2[i][s][i], p->eta-1.0) * p->q02[i][s][i];
	  Dbf[i][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Df[i][s][i] * pow(MC[i][s],1.0-p->eta);
	  Dtf[i][s][i] = Dbf[i][s][i]/p->eta;

	  double PhomeF = (1.0/p->H2[i][s][i]) * (p->eta/(p->eta-1.0))*MC[i][s] * pow(Zi[i][s],1.0/(1.0-p->eta));
	  if(fabs(PhomeF-1.0)>TINYSQ)
	    {
	      fprintf(logfile,"Error! PF[%d][%d] = %0.2f\n\n",i,i,PhomeF);
	      return 1;
	    }

	  double YFhome = Dbf[i][s][i] * Zi[i][s];
	  if(fabs(YFhome-p->q02[i][s][i])>TINYSQ)
	    {
	      fprintf(logfile,"Error! PF*YF2[%d][%d] = %0.2f, PF*YF20[%d][%d] = %0.2f\n\n",
		      i,i,YFhome,i,i,p->q02[i][s][i]);
	      return 1;
	    }

  
	  YFhome = p->H2[i][s][i] * pow( MC[i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
	    Df[i][s][i] * pow(Zi[i][s],p->eta/(p->eta-1.0));
	  
	  if(fabs(YFhome-p->q02[i][s][i])>TINYSQ)
	    {
	      fprintf(logfile,"Error! YF2[%d][%d] = %0.2f, YF20[%d][%d] = %0.2f\n\n",
		      i,i,YFhome,i,i,p->q02[i][s][i]);
	      return 1;
	    }

	  if(!noio_flag)
	    {
	      p->M2[i][s][i] = pow(p->eta/(p->eta-1.0),p->eta-1.0) *
		pow(MC[i][s],p->eta-1.0) * (1.0/Zi[i][s]);
	      p->M2[i][s][i] = pow(p->M2[i][s][i],1.0/(p->eta-1.0));

	      Dm[i][s][i] = pow(p->M2[i][s][i], p->eta-1.0) * p->m02[i][s][i];
	      Dbm[i][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Dm[i][s][i] * pow(MC[i][s],1.0-p->eta);
	      Dtm[i][s][i] = Dbm[i][s][i]/p->eta;

	      double PhomeM = (1.0/p->M2[i][s][i]) * (p->eta/(p->eta-1.0))*MC[i][s] * pow(Zi[i][s],1.0/(1.0-p->eta));

	      if(fabs(PhomeM-1.0)>TINYSQ)
		{
		  fprintf(logfile,"Error! PM[%d][%d] = %0.2f\n\n",i,i,PhomeM);
		  return 1;
		}

	      double YMhome = Dbm[i][s][i] * Zi[i][s];
	      if(fabs(YMhome-p->m02[i][s][i])>TINYSQ)
		{
		  fprintf(logfile,"Error! PM*YM2[%d][%d] = %0.2f, PM*YM20[%d][%d] = %0.2f\n\n",
			  i,i,YMhome,i,i,p->m02[i][s][i]);
		  return 1;
		}

	      YMhome = p->M2[i][s][i] * pow( MC[i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		Dm[i][s][i] * pow(Zi[i][s],p->eta/(p->eta-1.0));
	      
	      if(fabs(YMhome-p->m02[i][s][i])>TINYSQ)
		{
		  fprintf(logfile,"Error! YM2[%d][%d] = %0.2f, YM20[%d][%d] = %0.2f\n\n",
			  i,i,YMhome,i,i,p->m02[i][s][i]);
		  return 1;
		}
	    }
	  else
	    {
	      Dm[i][s][i] = 0.0;
	      Dbm[i][s][i] = 0.0;
	      Dtm[i][s][i] = 0.0;
	    }
	}
    }

    

  // calculate export demand stuff and fixed costs in static model where kappa0 and kappa1 are assumed to
  // be equal
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      j=p->Ji[i][ii];

	      if(nokappa || s==1)
		{
		  ev[i][s][ii].zp = -HUGE_VAL;
		  ev[i][s][ii].zm = -HUGE_VAL;
		  ev[i][s][ii].Zp = Zi[i][s];
		  ev[i][s][ii].Zm = Zi[i][s];
		  ev[i][s][ii].Z = Zi[i][s];
		  ev[i][s][ii].n = 1.0;
		}
	      else
		{
		  ev[i][s][ii].zp = p->sig_z[i][s]*gsl_cdf_ugaussian_Pinv(1.0-p->n0[i][s][ii]);
		  //ev[i][s][ii].zp = pow(1.0-p->n0[i][s][ii],-1.0/p->sig_z[i][s]);
		  ev[i][s][ii].zm = ev[i][s][ii].zp;

		  ev[i][s][ii].Zp = Zi[i][s] *
		    gsl_cdf_ugaussian_P( (p->sig_z[i][s]*p->sig_z[i][s]*(p->eta-1.0) - ev[i][s][ii].zp)/p->sig_z[i][s] );
		  //ev[i][s][ii].Zp = Zi[i][s] * pow(ev[i][s][ii].zp,p->eta-1.0-p->sig_z[i][s]);

		  ev[i][s][ii].Zm = Zi[i][s] *
		    gsl_cdf_ugaussian_P( (p->sig_z[i][s]*p->sig_z[i][s]*(p->eta-1.0) - ev[i][s][ii].zm)/p->sig_z[i][s]);
		  //ev[i][s][ii].Zm = Zi[i][s] * pow(ev[i][s][ii].zm,p->eta-1.0-p->sig_z[i][s]);

		  ev[i][s][ii].Z = p->n0[i][s][ii] * ev[i][s][ii].Zm + (1.0-p->n0[i][s][ii]) * ev[i][s][ii].Zp;
		}

	      p->H2[j][s][i] = pow(p->eta/(p->eta-1.0),p->eta-1.0) * pow(MC[i][s],p->eta-1.0) * (1.0/ev[i][s][ii].Z);
	      p->H2[j][s][i] = pow(p->H2[j][s][i],1.0/(p->eta-1.0));

	      Df[j][s][i] = pow(p->H2[j][s][i], p->eta-1.0) * p->q02[j][s][i];
	      Dbf[j][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Df[j][s][i] * pow(MC[i][s],1.0-p->eta);
	      Dtf[j][s][i] = Dbf[j][s][i]/p->eta;


	      double YFex = Dbf[j][s][i] * ev[i][s][ii].Z;
	      if(fabs(YFex-p->q02[j][s][i])>TINYSQ)
		{
		  fprintf(logfile,"Error! PF*YF2[%d][%d] = %0.2f, PF*YF20[%d][%d] = %0.2f\n\n",j,i,YFex,j,i,p->q02[j][s][i]);
		  return 1;
		}

	      YFex = p->H2[j][s][i] * pow( MC[i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		Df[j][s][i] * pow(ev[i][s][ii].Z,p->eta/(p->eta-1.0));
	      
	      if(fabs(YFex-p->q02[j][s][i])>TINYSQ)
		{
		  fprintf(logfile,"Error! YF2[%d][%d] = %0.2f, YF20[%d][%d] = %0.2f\n\n",j,i,YFex,j,i,p->q02[j][s][i]);
		  return 1;
		}

	      if(!noio_flag)
		{
		  p->M2[j][s][i] = pow(p->eta/(p->eta-1.0),p->eta-1.0) * pow(MC[i][s],p->eta-1.0) * (1.0/ev[i][s][ii].Z);
		  p->M2[j][s][i] = pow(p->M2[j][s][i],1.0/(p->eta-1.0));
		  
		  Dm[j][s][i] = pow(p->M2[j][s][i], p->eta-1.0) * p->m02[j][s][i];
		  Dbm[j][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Dm[j][s][i] * pow(MC[i][s],1.0-p->eta);
		  Dtm[j][s][i] = Dbm[j][s][i]/p->eta;

		  double YMex = Dbm[j][s][i] * ev[i][s][ii].Z;
		  if(fabs(YMex-p->m02[j][s][i])>TINYSQ)
		    {
		      fprintf(logfile,"Error! PM*YM2[%d][%d] = %0.2f, PM*YM20[%d][%d] = %0.2f\n\n",j,i,YMex,j,i,p->m02[j][s][i]);
		      return 1;
		    }
 
		  YMex = p->M2[j][s][i] * pow( MC[i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		    Dm[j][s][i] * pow(ev[i][s][ii].Z,p->eta/(p->eta-1.0));
		  
		  if(fabs(YMex-p->m02[j][s][i])>TINYSQ)
		    {
		      fprintf(logfile,"Error! YM2[%d][%d] = %0.2f, YM20[%d][%d] = %0.2f\n\n",j,i,YMex,j,i,p->m02[j][s][i]);
		      return 1;
		    }
		}
	      else
		{
		  Dm[j][s][i] = 0.0;
		  Dbm[j][s][i] = 0.0;
		  Dtm[j][s][i] = 0.0;
		}     

	      if(nokappa || s==1)
		//if(nokappa)
		{
		  p->kappa0[i][s][ii]=0.0;
		  p->kappa1[i][s][ii]=0.0;
		}
	      else
		{
		  p->kappa0[i][s][ii] = (Dtf[j][s][i] + Dtm[j][s][i]) * exp((p->eta-1.0)*ev[i][s][ii].zp);
		  //p->kappa0[i][s][ii] = (Dtf[j][s][i] + Dtm[j][s][i]) * pow(ev[i][s][ii].zp,p->eta-1.0);
		  p->kappa1[i][s][ii] = p->kappa0[i][s][ii];

		  p->kappa0[i][s][ii] = p->kappa0[i][s][ii]*1.25;
		  p->kappa1[i][s][ii] = p->kappa1[i][s][ii]*1.1;
		}
	    }
	}
    }

    

  solver_n = NC*NS + NC*NS*(NC-1)*2 + NC*(NS)*(NC-1)*2 + NC;
  //solver_n = NC*NS + NC*NS*(NC-1)*2 + NC*(NS)*(NC-1)*2;

  alloc_solver_mem();
  gsl_vector * solver_x_copy = gsl_vector_alloc(solver_n);
  uint status = 0;


  if(stack_calvars())
    {
      status=1;
    }
  else if(eval_calfn_once)
    {
      status = calfunc_f(solver_x,NULL,f0[0]);
    }
  else
    {
      solver_verbose=1;
      par=0;
      gsl_multiroot_function_fdf f = {&calfunc_f,&calfunc_df,&calfunc_fdf,solver_n,NULL};

      if(eqkappa || nokappa)
	{
	  fprintf(logfile,KBLU "Calibrating firm parameters...\n" RESET);
	  status = find_root_deriv_mkl(&f);
	  if(status)
	    {
	      fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
	      status=1;
	    }
	}
      else
	{	  
	  eqkappa=1;
	  
	  fprintf(logfile,KBLU "\nOne-shot calibration for kappa0=kappa1...\n" RESET);

	  
	  status = find_root_deriv_mkl(&f);

	  
	  if(status)
	    {
	      fprintf(logfile,KRED "\nError solving calibration function!\n" RESET);
	      status=1;
	    }
	  else
	    {
	      //write_params();
	      COPY_SUBVECTOR(tau_k_tmp,p->tauk,NC);
	      COPY_SUBVECTOR(sig_z_tmp,p->sig_z,NC*NS);
	      eqkappa=0;
	      for(i=0; i<NC; i++)
		{
		  for(s=0; s<NS; s++)
		    {
		      for(ii=0; ii<(NC-1); ii++)
			{
			  exit_rate_target2[i][s][ii] = 1.0-export_participation_rate_target[i][s][ii];
			}
		    }
		}

	      fprintf(logfile,KBLU "\nHomotopy process for kappa0!=kappa1 beginning...\n\n" RESET);

	      double grid[NC][NS][NC-1][homotopy_times];
	      for(i=0; i<NC; i++)
		{
		  for(s=0; s<NS; s++)
		    {
		      for(ii=0; ii<(NC-1); ii++)
			{
			  linspace(exit_rate_target[i][s][ii],
				   exit_rate_target2[i][s][ii],
				   homotopy_times,
				   grid[i][s][ii]);
			}
		    }
		}
	      int h;
	      for(h=homotopy_times-1;h>=0;h--)
		{
		  for(i=0; i<NC; i++)
		    {
		      for(s=0; s<NS; s++)
			{
			  for(ii=0; ii<(NC-1); ii++)
			    {
			      if(h==(homotopy_times-2))
				{
				  p->kappa1[i][s][ii] = p->kappa1[i][s][ii]*0.95;
				}
			      exit_rate_target2[i][s][ii] = grid[i][s][ii][h];
			    }
			}
		    }

		  if(h==(homotopy_times-2))
		    {
		      stack_calvars();
		    }
		  else if(h<(homotopy_times-2))
		    {
		      gsl_vector_memcpy(solver_x,solver_x_copy);
		    }
		  fprintf(logfile,"\tExit rate target: %0.3f",exit_rate_target2[0][0][0]);
		  status = find_root_deriv_mkl(&f);
		  if(!status)
		    {
		      gsl_vector_memcpy(solver_x_copy,solver_x);
		    }
		  //fprintf(logfile,"%0.16f, %0.16f\n",p->kk0[0],SUM(p->k0[0],NS));
		  fprintf(logfile,"\n");
		  
		}
	    }
	}

    }
  solver_verbose=1;

  free_solver_mem();
  gsl_vector_free(solver_x_copy);

  if(!no_k_flag)
    {
      for(i=0; i<NC; i++)
	{
	  double tmp = SUM(p->k0[i],NS)/p->kk0[i];
	  for(s=0; s<NS; s++)
	    {
	      p->i0[i][s] = tmp * p->i0[i][s];
	      p->c0[i][s] = p->q0[i][s] - p->i0[i][s];

	      if(p->c0[i][s]<0.0)
		{
		  fprintf(logfile,KRED "Negative consumption!\n" RESET);
		  return 1;
		}
	    }
	  p->kk0[i] = SUM(p->k0[i],NS);
	}
    }
  
  return status;
}

uint stack_calvars()
{
  params * p = &(ppp0[0]);
  uint i,s,ii;
  uint nx=0;

  COPY_SUBVECTOR_LOG(solver_x->data+nx,p->sig_z,NC*NS);
  nx=nx+NC*NS;

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      uint j = p->Ji[i][ii];
	      
	      solver_x->data[nx] = log(p->H2[j][s][i]);
	      nx=nx+1;

	      if(!noio_flag)
		{
		  solver_x->data[nx] = log(p->M2[j][s][i]);
		}
	      else
		{
		  solver_x->data[nx] = p->M2[j][s][i];
		}
	      nx=nx+1;

	      if(nokappa || s==1)
		{
		  solver_x->data[nx] = p->kappa0[i][s][ii];
		  nx = nx+1;
		  
		  solver_x->data[nx] = p->kappa1[i][s][ii];
		  nx = nx+1;
		}
	      else
	      {
		solver_x->data[nx] = log(p->kappa0[i][s][ii]);
		nx = nx+1;
		  
		solver_x->data[nx] = log(p->kappa1[i][s][ii]);
		nx = nx+1;
	      }
	    }
	}

      solver_x->data[nx] = p->tauk[i];
      nx=nx+1;
    }

  if(nx!=solver_n)
    {
      printf("Wrong number of parameters got stacked!\n");
      return 1;
    }
  
  return 0;
}

uint eval_exporter_moments(const double * calpars, double * myf)
{
  params * p = &(ppp0[0]);
  uint i, s, ii, j, r;
  
  uint nx = 0;

  COPY_SUBVECTOR_EXP(p->sig_z,calpars+nx,NC*NS);
  nx=nx+NC*NS;

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      j = p->Ji[i][ii];
	      
	      p->H2[j][s][i] = exp(calpars[nx]);
	      nx=nx+1;

	      if(!noio_flag)
		{
		  p->M2[j][s][i] = exp(calpars[nx]);
		}
	      else
		{
		  p->M2[j][s][i] = calpars[nx];
		}
	      nx=nx+1;

	      if(nokappa || s==1)
		{
		  p->kappa0[i][s][ii] = calpars[nx];
		  nx=nx+1;

		  p->kappa1[i][s][ii] = calpars[nx];
		  nx=nx+1;
		}
	      else
		{
		  p->kappa0[i][s][ii] = exp(calpars[nx]);
		  nx=nx+1;

		  p->kappa1[i][s][ii] = exp(calpars[nx]);
		  nx=nx+1;
		}
	    }
	}

      p->tauk[i] = calpars[nx];
      nx=nx+1;

    }
  
  if(nx!=solver_n)
    {
      printf("Wrong number of parameters got stacked!\n");
      return 1;
    }

  // modify value added share
  /*
  if((eqkappa || nokappa) && !noio_flag)
    {
      if(cobb_douglas_flag)
	{
	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  double tmp = p->lam_va[i][s];
		  p->lam_va[i][s] = 1.0 - (p->eta/(p->eta-1.0))*SUM(p->md0_[i][s],NS)/p->y0[i][s];
		  double tmp2 = (1.0-p->lam_va[i][s])/(1.0-tmp);
		  for(r=0; r<NS; r++)
		    {
		      p->lam[i][s][r] = p->lam[i][s][r]*tmp2;
		    }
		}
	    }
	}
      else
	{
	  for(i=0; i<NC; i++)
	    {
	      for(s=0; s<NS; s++)
		{
		  double tmp = p->lam_va[i][s];
		  double sq = (p->eta/(p->eta-1.0))*SUM(p->md0_[i][s],NS)/p->y0[i][s];
		  double br = 0.0;
		  if(!no_k_flag)
		    {
		      br = pow((p->r0[i]+p->delta)/(1.0-p->tauk[i])/p->alpha[i][s],p->alpha[i][s])*
			pow(1.0/(1.0-p->alpha[i][s]),1.0-p->alpha[i][s]);
		    }
		  else
		    {
		      br = 1.0;
		    }
	    
		  p->lam_va[i][s] = (1.0-sq)/(1.0-sq+br*sq);

		  double tmp2 = (1.0-p->lam_va[i][s])/(1.0-tmp);
		  for(r=0; r<NS; r++)
		    {
		      p->lam[i][s][r] = p->lam[i][s][r]*tmp2;
		    }
		}
	    }
	}
    }
  */


  double MC[NC][NS];
  double Df[NC][NS][NC];
  double Dbf[NC][NS][NC];
  double Dhf[NC][NS][NC];
  double Dtf[NC][NS][NC];
  double Dm[NC][NS][NC];
  double Dbm[NC][NS][NC];
  double Dhm[NC][NS][NC];
  double Dtm[NC][NS][NC];
  
  // recalculate export demand scale factors
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  if(no_k_flag)
	    {
	      MC[i][s] = p->lam_va[i][s] + SUM(p->lam[i][s],NS);
	    }
	  else if(cobb_douglas_flag==0)
	    {
	      MC[i][s] = (p->lam_va[i][s]*(
					   pow((p->r0[i]+p->delta)/(1.0-p->tauk[i])/p->alpha[i][s],p->alpha[i][s]) *
					   pow(1.0/(1.0-p->alpha[i][s]),1.0-p->alpha[i][s])
					   )
			  ) + SUM(p->lam[i][s],NS);
	      
	    }
	  else
	    {
	      MC[i][s] = (pow( (p->r0[i]+p->delta)/(1.0-p->tauk[i])/(p->alpha[i][s]*p->lam_va[i][s]),
			       p->lam_va[i][s]*p->alpha[i][s]) *
			  pow(1.0/(1.0-p->alpha[i][s])/p->lam_va[i][s],(1.0-p->alpha[i][s])*p->lam_va[i][s]));
	      if(!noio_flag)
		{
		  for(r=0; r<NS; r++)
		    {
		      MC[i][s] = MC[i][s] * pow(1.0/p->lam[i][s][r],p->lam[i][s][r]);
		    }
		}
	    }
	
	  for(j=0; j<NC; j++)
	    {
   
	      Df[j][s][i] = pow(p->H2[j][s][i], p->eta-1.0) * p->q02[j][s][i];
	      Dbf[j][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Df[j][s][i] * pow(MC[i][s],1.0-p->eta);
	      Dhf[j][s][i] = pow(p->eta/(p->eta-1.0),-p->eta) * Df[j][s][i] * pow(MC[i][s],1.0-p->eta);
	      Dtf[j][s][i] = Dbf[j][s][i]/p->eta;

	      if(!noio_flag)
		{
		  Dm[j][s][i] = pow(p->M2[j][s][i], p->eta-1.0) * p->m02[j][s][i];
		  Dbm[j][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Dm[j][s][i] * pow(MC[i][s],1.0-p->eta);
		  Dhm[j][s][i] = pow(p->eta/(p->eta-1.0),-p->eta) * Dm[j][s][i] * pow(MC[i][s],1.0-p->eta);
		  Dtm[j][s][i] = Dbm[j][s][i]/p->eta;
		}
	      else
		{
		  Dm[j][s][i] = 0.0;
		  Dbm[j][s][i] = 0.0;
		  Dhm[j][s][i] = 0.0;
		  Dtm[j][s][i] = 0.0;
		}
	    }
	}
    }


    for(i=0; i<NC; i++)
      {
	for(s=0; s<NS; s++)
	  {
	    double Zi=0.0;
	    Zi = exp(p->sig_z[i][s]*p->sig_z[i][s]*(p->eta-1.0)*(p->eta-1.0)/2.0);
	    //Zi = p->sig_z[i][s]/(p->sig_z[i][s]+1.0-p->eta);

	    p->H2[i][s][i] = pow(p->eta/(p->eta-1.0),p->eta-1.0) *
	      pow(MC[i][s],p->eta-1.0) * (1.0/Zi);
	    p->H2[i][s][i] = pow(p->H2[i][s][i],1.0/(p->eta-1.0));

	    Df[i][s][i] = pow(p->H2[i][s][i], p->eta-1.0) * p->q02[i][s][i];
	    Dbf[i][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Df[i][s][i] * pow(MC[i][s],1.0-p->eta);
	    Dhf[i][s][i] = pow(p->eta/(p->eta-1.0),-p->eta) * Df[i][s][i] * pow(MC[i][s],1.0-p->eta);
	    Dtf[i][s][i] = Dbf[i][s][i]/p->eta;

	    double PhomeF = (1.0/p->H2[i][s][i]) * (p->eta/(p->eta-1.0))*MC[i][s] * pow(Zi,1.0/(1.0-p->eta));

	    if(fabs(PhomeF-1.0)>TINYSQ)
	      {
		fprintf(logfile,"Error! PF[%d][%d] = %0.2f\n\n",i,i,PhomeF);
		return 1;
	      }
	  
	    double YFhome = Dbf[i][s][i] * Zi;
	    if(fabs(YFhome-p->q02[i][s][i])>TINYSQ)
	      {
		fprintf(logfile,"Error! PF*YF2[%d][%d] = %0.2f, PF*YF20[%d][%d] = %0.2f\n\n",
			i,i,YFhome,i,i,p->q02[i][s][i]);
		return 1;
	      }

	    YFhome = p->H2[i][s][i] * pow( MC[i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
	      Df[i][s][i] * pow(Zi,p->eta/(p->eta-1.0));
	  
	    if(fabs(YFhome-p->q02[i][s][i])>TINYSQ)
	      {
		fprintf(logfile,"Error! YF2[%d][%d] = %0.2f, YF20[%d][%d] = %0.2f\n\n",
			i,i,YFhome,i,i,p->q02[i][s][i]);
		return 1;
	      }

	    if(!noio_flag)
	      {
		p->M2[i][s][i] = pow(p->eta/(p->eta-1.0),p->eta-1.0) *
		  pow(MC[i][s],p->eta-1.0) * (1.0/Zi);
		p->M2[i][s][i] = pow(p->M2[i][s][i],1.0/(p->eta-1.0));

		Dm[i][s][i] = pow(p->M2[i][s][i], p->eta-1.0) * p->m02[i][s][i];
		Dbm[i][s][i] = pow(p->eta/(p->eta-1.0),1.0-p->eta) * Dm[i][s][i] * pow(MC[i][s],1.0-p->eta);
		Dhm[i][s][i] = pow(p->eta/(p->eta-1.0),-p->eta) * Dm[i][s][i] * pow(MC[i][s],1.0-p->eta);
		Dtm[i][s][i] = Dbm[i][s][i]/p->eta;

		double PhomeM = (1.0/p->M2[i][s][i]) * (p->eta/(p->eta-1.0))*MC[i][s] * pow(Zi,1.0/(1.0-p->eta));

		if(fabs(PhomeM-1.0)>TINYSQ)
		  {
		    fprintf(logfile,"Error! PM[%d][%d] = %0.2f\n\n",i,i,PhomeM);
		    return 1;
		  }

		double YMhome = Dbm[i][s][i] * Zi;
		if(fabs(YMhome-p->m02[i][s][i])>TINYSQ)
		  {
		    fprintf(logfile,"Error! PM*YM2[%d][%d] = %0.2f, PM*YM20[%d][%d] = %0.2f\n\n",
			    i,i,YMhome,i,i,p->m02[i][s][i]);
		    return 1;
		  }

		YMhome = p->M2[i][s][i] * pow( MC[i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		  Dm[i][s][i] * pow(Zi,p->eta/(p->eta-1.0));
	      
		if(fabs(YMhome-p->m02[i][s][i])>TINYSQ)
		  {
		    fprintf(logfile,"Error! YM2[%d][%d] = %0.2f, YM20[%d][%d] = %0.2f\n\n",
			    i,i,YMhome,i,i,p->m02[i][s][i]);
		    return 1;
		  }
	      }
	  }
      }
  

    // steady state exporter dynamics
    for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      j=p->Ji[i][ii];

	      if(ev_steady_state(p->sig_z[i][s],
				 p->eta,
				 p->kappa0[i][s][ii],
				 p->kappa1[i][s][ii],
				 1.0,
				 p->beta[i],
				 Dtf[j][s][i] + Dtm[j][s][i],
				 p->n0[i][s][ii],
				 &(ev[i][s][ii])))
		{
		  fprintf(logfile, KRED "Failed to solve for steady state exporter moments for country/sector/dest %d/%d/%d!\n" RESET,i,s,j);
		  return 1;
		}
	    }
	}
    }

  // factor demand
  double exrate[NC][NS][NC-1];
  double exitrate[NC][NS][NC-1];
  double exportdiffF[NC][NS][NC-1];
  double exportdiffM[NC][NS][NC-1];
  double ratio[NC][NS];
  double Kd[NC][NS];
  double Ld[NC][NS];
  double Md[NC][NS][NS];

  
  for(i=0; i<NC; i++)
    {      
      for(s=0; s<NS; s++)
	{
	  if(no_k_flag)
	    {
	      Ld[i][s] = (p->lam_va[i][s]/MC[i][s]) * (Dhf[i][s][i]+Dhm[i][s][i])*ev[i][s][0].Zi;

	      for(r=0; r<NS; r++)
		{
		  Md[i][s][r] = (p->lam[i][s][r]/MC[i][s]) * (Dhf[i][s][i]+Dhm[i][s][i]) * ev[i][s][0].Zi;
		}
	    }
	  
	  else if(cobb_douglas_flag==0)
	    {

	      Kd[i][s] = (p->lam_va[i][s]/MC[i][s]) *
		pow( (p->r0[i]+p->delta)/(1.0-p->tauk[i])*(1.0-p->alpha[i][s])/1.0/p->alpha[i][s],p->alpha[i][s]-1.0)*(Dhf[i][s][i]+Dhm[i][s][i])*ev[i][s][0].Zi;

	      if(Kd[i][s]<0.0)
		{
		  printf("Error in Kd[i][s]! %d %d\n",i,s);
		  printf("%0.4f\n",p->lam_va[i][s]);
		  printf("%0.4f\n",MC[i][s]);
		  printf("%0.4f\n",(p->r0[i]+p->delta));
		  printf("%0.4f\n",(1.0-p->tauk[i]));
		  printf("%0.4f\n",(Dhf[i][s][i]+Dhm[i][s][i])*ev[i][s][0].Zi);
		  exit(1);
		}

	      Ld[i][s] = (p->lam_va[i][s]/MC[i][s]) *
		pow( (p->r0[i]+p->delta)/(1.0-p->tauk[i])*(1.0-p->alpha[i][s])/1.0/p->alpha[i][s],p->alpha[i][s])*(Dhf[i][s][i]+Dhm[i][s][i])*ev[i][s][0].Zi;

	      for(r=0; r<NS; r++)
		{
		  Md[i][s][r] = (p->lam[i][s][r]/MC[i][s]) * (Dhf[i][s][i]+Dhm[i][s][i]) * ev[i][s][0].Zi;
		}
	    }
	  else
	    {
	      Kd[i][s] = (p->lam_va[i][s]*p->alpha[i][s]/(p->r0[i]+p->delta)*(1.0-p->tauk[i])) * (Dhf[i][s][i]+Dhm[i][s][i]) * ev[i][s][0].Zi;
	      Ld[i][s] = (p->lam_va[i][s]*(1.0-p->alpha[i][s])/1.0) * (Dhf[i][s][i]+Dhm[i][s][i]) * ev[i][s][0].Zi;
	      for(r=0; r<NS; r++)
		{
		  Md[i][s][r] = (p->lam[i][s][r]/1.0) * (Dhf[i][s][i]+Dhm[i][s][i]) * ev[i][s][0].Zi;
		}
	    }

	  double lf = 0.0;
	  double weight = 0.0;
	  ratio[i][s]=0.0;
	  
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      //double ynon = 0.0;
	      //double yex = 0.0;

	      j=p->Ji[i][ii];

	      exrate[i][s][ii] = ev[i][s][ii].n;
	      exitrate[i][s][ii] = ev[i][s][ii].Fzm;


	      if(no_k_flag)
		{
		  Ld[i][s] += (p->lam_va[i][s]/MC[i][s]) * (Dhf[j][s][i]+Dhm[j][s][i]) * ev[i][s][ii].Z;
		  
		  for(r=0; r<NS; r++)
		    {
		      Md[i][s][r] += (p->lam[i][s][r]/MC[i][s]) * (Dhf[j][s][i]+Dhm[j][s][i]) * ev[i][s][ii].Z;
		    }

		}
	      else if(cobb_douglas_flag==0)
		{
		  Kd[i][s] += (p->lam_va[i][s]/MC[i][s]) *
		    pow( (p->r0[i]+p->delta)/(1.0-p->tauk[i])*(1.0-p->alpha[i][s])/1.0/p->alpha[i][s],p->alpha[i][s]-1.0)*(Dhf[j][s][i]+Dhm[j][s][i]) * ev[i][s][ii].Z;
		  
		  if(Kd[i][s]<0.0)
		{
		  printf("Error in Kd[i][s]! %d %d %d\n",i,s,j);
		  exit(1);
		}

		  
		  Ld[i][s] += (p->lam_va[i][s]/MC[i][s]) *
		    pow( (p->r0[i]+p->delta)/(1.0-p->tauk[i])*(1.0-p->alpha[i][s])/1.0/p->alpha[i][s],p->alpha[i][s]) * (Dhf[j][s][i] + Dhm[j][s][i]) * ev[i][s][ii].Z;

		  for(r=0; r<NS; r++)
		    {
		      Md[i][s][r] += (p->lam[i][s][r]/MC[i][s]) * (Dhf[j][s][i]+Dhm[j][s][i]) * ev[i][s][ii].Z;
		    }
		}
	      else
		{
		  Kd[i][s] += (p->lam_va[i][s]*p->alpha[i][s]/(p->r0[i]+p->delta)*(1.0-p->tauk[i])) * (Dhf[j][s][i] +Dhm[j][s][i]) * ev[i][s][ii].Z;
		  Ld[i][s] += (p->lam_va[i][s]*(1.0-p->alpha[i][s])/1.0) * (Dhf[j][s][i]+Dhm[j][s][i]) * ev[i][s][ii].Z;
		  for(r=0; r<NS; r++)
		    {
		      Md[i][s][r] += (p->lam[i][s][r]/1.0) * (Dhf[j][s][i]+Dhm[j][s][i]) * ev[i][s][ii].Z;
		    }
		}

	      if(!nokappa)
		{
		  lf += (
			 p->kappa0[i][s][ii]*(1.0-ev[i][s][ii].n)*(1.0-ev[i][s][ii].Fzp) +
			 p->kappa1[i][s][ii]*ev[i][s][ii].n*(1.0-ev[i][s][ii].Fzm)
			 );
		}

	      double PYFj = Dbf[j][s][i] * ev[i][s][ii].Z;
	      exportdiffF[i][s][ii] = (PYFj - p->q02[j][s][i])/p->q02[j][s][i];

	      if(!noio_flag)
		{
		  double PYMj = Dbm[j][s][i] * ev[i][s][ii].Z;
		  exportdiffM[i][s][ii] = (PYMj - p->m02[j][s][i])/p->m02[j][s][i];
		}

	      //double Znon1 = ev[i][s][ii].Zi *
	      //gsl_cdf_ugaussian_P((ev[i][s][ii].zm - p->sig_z[i][s]*p->sig_z[i][s]*(p->eta-1.0)) / p->sig_z[i][s]);

	      //double Znon2 = ev[i][s][ii].Zi *
	      //gsl_cdf_ugaussian_P((ev[i][s][ii].zp - p->sig_z[i][s]*p->sig_z[i][s]*(p->eta-1.0)) / p->sig_z[i][s]);

	      //ynon = (Dbf[i][s][i] + Dbm[i][s][i]) * (ev[i][s][ii].n*Znon1 + (1.0-ev[i][s][ii].n)*Znon2) / (1.0-ev[i][s][ii].n);

	      //yex = (Dbf[i][s][i] + Dbm[i][s][i] + Dbf[j][s][i] + Dbm[j][s][i]) * ev[i][s][ii].Z / ev[i][s][ii].n;

	      //ratio[i][s] += ev[i][s][ii].n * yex/ynon;
	      //weight += ev[i][s][ii].n;
	      //ratio[i][s] += 0.5 * yex/ynon;
	      //weight+=0.5;
	      
	      // compute top 5 share
	      double log_top5_cutoff = gsl_cdf_ugaussian_Pinv(1.0-0.05*ev[i][s][ii].n)*p->sig_z[i][s];
	      double tmp1_5 = (exp(p->sig_z[i][s]*p->sig_z[i][s]*
				   (p->eta-1.0)*(p->eta-1.0) / 2.0) *
			       gsl_cdf_ugaussian_P( (p->sig_z[i][s]*
						     p->sig_z[i][s]*
						     (p->eta-1.0) -
						     log_top5_cutoff) /
						     p->sig_z[i][s] ));

	      //double top5_cutoff = pow(1.0-(1.0-0.05*p->n0[i][s][ii]),-1.0/p->sig_z[i][s]);
	      //double tmp1_5 = (p->sig_z[i][s]/(p->sig_z[i][s]+1.0-p->eta)) * pow(top5_cutoff,p->eta-p->sig_z[i][s]-1.0);
	      
	      double Z_5 = tmp1_5;
	      double top5_exports = Z_5/ev[i][s][ii].Z;

	      ratio[i][s] += top5_exports*
	      (p->q02[j][s][i]+p->m02[j][s][i]);
	      
	      weight+=(p->q02[j][s][i]+p->m02[j][s][i]);	      
	    }
	  
	  ratio[i][s] = ratio[i][s]/weight;
	  
	  // recompute factors, investment, and consumption as necessary
	  p->k0[i][s] = Kd[i][s];
	  p->l0[i][s] = Ld[i][s];
	  for(r=0; r<NS; r++)
	    {
	      p->md0[i][s][r] = Md[i][s][r];
	    }
	  Ld[i][s] = Ld[i][s] + lf;
	}    

      p->lbar[i] = 3.0*SUM(Ld[i],NS);
    }

  uint cnt=0;
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  //if(s!=1)
	  {
	    for(ii=0; ii<(NC-1); ii++)
	      {
		if(s==1 || nokappa)
		  {
		    myf[cnt] = p->kappa0[i][s][ii];
		  }
		else
		  {
		    myf[cnt] = (exrate[i][s][ii] -
				export_participation_rate_target[i][s][ii]);
		  }
		cnt = cnt+1;
	      }

	    for(ii=0; ii<(NC-1); ii++)
	      {
		if(s==1 || nokappa)
		  {
		    myf[cnt] = p->kappa1[i][s][ii];
		  }
		else if(eqkappa)
		  {
		    myf[cnt] = p->kappa0[i][s][ii]-p->kappa1[i][s][ii];
		  }
		else
		  {
		    myf[cnt] = (exitrate[i][s][ii] - exit_rate_target2[i][s][ii]);
		  }
		cnt=cnt+1;
	      }
	  }

	  for(ii=0; ii<(NC-1); ii++)
	    {
	      myf[cnt] = exportdiffF[i][s][ii];
	      cnt = cnt+1;

	      if(!noio_flag)
		{
		  myf[cnt] = exportdiffM[i][s][ii];
		}
	      else
		{
		  myf[cnt] = p->M2[i][s][p->Ji[i][ii]];
		}
	      cnt = cnt+1;
	    }

	  if(eqkappa || nokappa)
	    {
	      myf[cnt] = (ratio[i][s] - ratio_target[i][s]);
	      cnt = cnt+1;
	    }
	  else
	    {
	      myf[cnt] = p->sig_z[i][s] - sig_z_tmp[i][s];
	      cnt=cnt+1;
	    }
	}
      
      if(no_k_flag)
	{
	  myf[cnt] = p->tauk[i];
	  cnt = cnt+1;
	}
      else if(eqkappa || nokappa)
	{
	  myf[cnt] = SUM(p->k0[i],NS) - p->kk0[i];
	  cnt = cnt+1;
	}
      else
	{
	  myf[cnt] = p->tauk[i] - tau_k_tmp[i];
	  cnt=cnt+1;
	}
    }

  if(cnt != solver_n)
    {
      printf("Wrong number of calibration equations!\n");
      return 1;
    }
      
  for(i=0; i<cnt; i++)
    {
      if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	{
	  return 1;
	}
    }
      

  return 0;
}

int calfunc_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  if(eval_exporter_moments(x->data,f->data))
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

int calfunc_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&calfunc_f, x, J, 1))
    {
      return 1;
    }
  else
    {
      return 0;
    }
}

int calfunc_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(calfunc_f(x,NULL,f))
    {
      return 1;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&calfunc_f, x, J, 0))
	{
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
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

  /*
  if(load_tariffs(p))
    {
      return 1;
    }
  */

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


  if(calibrate_firm_params())
    {
      return 1;
    }


  write_params();
  
  uint it;
  for(it=0; it<NTH; it++)
    {
      if(it>0 && copy_params(&(ppp0[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp0!\n" RESET);
	  return 1;
	}
      if(copy_params(&(ppp1[it]),&(ppp0[0])))
	{
	  fprintf(logfile, KRED "\nFailed to copy ppp1!\n" RESET);
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

      fprintf(file,"sig_z:\n");
      fprintf_mat(file,p->sig_z,NC);
	
      fprintf(file,"\n3D PARAMETERS:\n\n");

      fprintf(file,"eps (NC x 2 x NS):\n");
      fprintf_3d_1(file,p->eps,NC);

      fprintf(file,"theta (NC x NS x NC):\n");
      fprintf_3d_2(file,p->theta,NC);

      fprintf(file,"mu (NC x NS x NC):\n");
      fprintf_3d_2(file,p->mu,NC);

      fprintf(file,"lam (NC x NS x NS):\n");
      fprintf_3d_3(file,p->lam,NC);

      fprintf(file,"kappa0 (NC x NS x NC-1):\n");
      fprintf_3d_4(file,p->kappa0,NC);

      fprintf(file,"kappa1 (NC x NS x NC-1):\n");
      fprintf_3d_4(file,p->kappa1,NC);

      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write parameters!\n" RESET);
      return 1;
    }
}

#endif
