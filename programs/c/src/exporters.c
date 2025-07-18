#ifndef __EXPORTERS_C__

#define __EXPORTERS_C__

#include "exporters.h"

void copy_exporter_vars(exporter_vars * dest, const exporter_vars * src)
{
  dest->zp = src->zp;
  dest->zm = src->zm;
  dest->Zi = src->Zi;
  dest->Zp = src->Zp;
  dest->Zm = src->Zm;
  dest->Z = src->Z;
  dest->Fzp = src->Fzp;
  dest->Fzm = src->Fzm;
  dest->dV = src->dV;
  dest->n = src->n;
}

void calc_z_Z_Fz(double sigma,
		 double theta,
		 double kappa0,
		 double kappa1,
		 double W,
		 double Dtj,
		 double dVp,
		 exporter_vars * ev)
{
  // now calculate conditional productivities
  if(nokappa || kappa0<1.e-8 || W*kappa0<dVp)
    {
      ev->zp=-HUGE_VAL;
      //ev->zp = 1.0;
    }
  else
    {
      ev->zp = (1.0/(theta-1.0)) * (log(W*kappa0 - dVp) - log(Dtj));
      //ev->zp = fmax(pow( (W*kappa0-dVp)/Dtj , 1.0/(theta-1.0) ),1.0);
    }

  if(nokappa || kappa1<1.0e-8 || W*kappa1 < dVp)
    {
      ev->zm = -HUGE_VAL;
      //ev->zm  = 1.0;
    }
  else
    {
      ev->zm = (1.0/(theta-1.0)) * (log(W*kappa1 - dVp) - log(Dtj));
      //ev->zm = fmax(pow( (W*kappa1-dVp)/Dtj , 1.0/(theta-1.0) ),1.0);
    }

  // now calculate conditional productivities
  ev->Zi = exp(sigma*sigma*(theta-1.0)*(theta-1.0)/2.0);
  //ev->Zi = sigma/(sigma+1.0-theta);
  
  if(nokappa || kappa0<1.0e-8 || gsl_isinf(ev->zp)==-1)
    //if(nokappa || kappa0<1.0e-8 || ev->zp<1.0)
    {
      ev->Zp = ev->Zi;
    }
  else if(gsl_isinf(ev->zp)==1)
    {
      ev->Zp = 0.0;
    }
  else
    {
      ev->Zp = ev->Zi * gsl_cdf_ugaussian_P( ((theta-1.0)*sigma*sigma - ev->zp)/sigma );
      //ev->Zp = ev->Zi * pow(ev->zp,theta-sigma-1.0);
    }

  if(nokappa || kappa1<1.0e-8 || gsl_isinf(ev->zm)==-1)
    //if(nokappa || kappa1<1.0e-8 || ev->zm<1.0)
    {
      ev->Zm = ev->Zi;
    }
  else if(gsl_isinf(ev->zm)==1)
    {
      ev->Zm = 0.0;
    }
  else
    {
      ev->Zm = ev->Zi * gsl_cdf_ugaussian_P( ((theta-1.0)*sigma*sigma - ev->zm)/sigma );
      //ev->Zm = ev->Zi * pow(ev->zm,theta-sigma-1.0);
    }

  // now calculate probabilities
  if(nokappa)
    {
      ev->Fzp = 0.0;
      ev->Fzm = 0.0;
    }
  else
    {
      ev->Fzp = gsl_cdf_ugaussian_P(ev->zp/sigma);
      ev->Fzm = gsl_cdf_ugaussian_P(ev->zm/sigma);
      //ev->Fzp = 1.0-pow(ev->zp,-sigma);
      //ev->Fzm = 1.0-pow(ev->zm,-sigma);
    }

  return;
}

void init_n(double n, exporter_vars * ev)
{
  ev->n = n;
}

void update_n(double nm,
	      exporter_vars * ev)
{
  if(nokappa)
    {
      ev->Z = ev->Zi;
      ev->n = 1.0;
    }
  else
    {
      ev->Z = nm * ev->Zm + (1.0-nm) * ev->Zp;
      //ev->Z = ev->Zi;
      ev->n = nm * (1.0-ev->Fzm) + (1.0-nm) * (1.0-ev->Fzp);
    }
}

void init_dV(double sigma,
	     double theta,
	     double kappa1,
	     double W,
	     double Lambda,
	     double Dtj,
	     exporter_vars * ev)
{
  if(eqkappa)
    {
      ev->dV = 0.0;
    }
  else
    {
      //double Zi = exp(sigma*sigma*(theta-1.0)*(theta-1.0)/2.0);
      //double Zi = sigma/(sigma+1.0-theta);
      //double V0 = Dti*Zi/(1.0-Lambda);
      //double V1 = (Dti*Zi + Dtj*Zi)/(1.0-Lambda);
      //ev->dV = Dtj*Zi/(1.0-Lambda);
      ev->dV = 0.0;
    }

  return;
}

void update_dV(double kappa0,
	       double kappa1,
	       double W,
	       double Lambda,
	       double Dtj,
	       double dVm,
	       exporter_vars * ev)
{	
  if(nokappa || eqkappa)
    {
      ev->dV=0.0;
    }
  else
    {
      ev->dV = Lambda * (Dtj*(ev->Zm - ev->Zp) 
	+ (ev->Fzp - ev->Fzm) * dVm
			 - (1.0-ev->Fzm)*W*kappa1 + (1.0-ev->Fzp)*W*kappa0);
    }
  
}

uint ev_steady_state(double sigma,
		     double theta,
		     double kappa0,
		     double kappa1,
		     double W,
		     double Lambda,
		     double Dtj,
		     double n0,
		     exporter_vars * ev)
{
  if(!nokappa)
    {
      init_n(n0, ev);
      init_dV(sigma,theta,kappa1,W,Lambda,Dtj,ev);

      uint t=0;
      double ddV = +HUGE_VAL;
      double dn = +HUGE_VAL;
      double dVm=0.0;
      double nm=0.0;
      while(t<10000 && (ddV > 1.0e-10 || dn > 1.0e-10))
	{
	  t++;
	  dVm = ev->dV;
	  nm = ev->n;
	  calc_z_Z_Fz(sigma,theta,kappa0,kappa1,W,Dtj,dVm,ev);
	  update_n(ev->n,ev);
	  update_dV(kappa0,kappa1,W,Lambda,Dtj,ev->dV,ev);
	  ddV = fabs(ev->dV-dVm);
	  dn = fabs(ev->n-nm);
	}

      if(t==10000)
	{
	  fprintf(logfile,KRED "Steady state value function failed to converge! ddV = %0.3g\n" KNRM,ev->dV - dVm);
	  return 1;
	}
    }
  else
    {
      calc_z_Z_Fz(sigma,theta,kappa0,kappa1,W,Dtj,0.0,ev);
      update_n(ev->n,ev);
    }
  
  return 0;
}

#endif
