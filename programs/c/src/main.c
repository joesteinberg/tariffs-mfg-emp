#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals.h"
#include "calibrate.h"
#include "eqm.h"

double tau=0.0;

uint parse_args(int argc, char **argv)
{
  int opt = 0;

  fprintf(logfile,"Model scenario:");
  
  while((opt = getopt(argc, argv, "t:s:c:r:d:")) != -1)
    {
      int o = atoi(optarg);
      
      switch(opt){

      case 't':
	fprintf(logfile,"Tariff rate: %d pct",o);
	tau=o/100.0;
      
      case 's':
	fprintf(logfile,"Targeted sector: ");
	if(o==0)
	  {
	    target_sector_flag=0;
	    fprintf(logfile,"Upstream goods only\n");
	  }
	else if(o==1)
	  {
	    target_sector_flag=1;
	    fprintf(logfile,"Downstream goods only\n");
	  }
	else if(o==2)
	  {
	    target_sector_flag=2;
	    fprintf(logfile,"Upstream + downstream goods\n");
	  }
	else
	  {
	    fprintf(logfile,"Invalid command-line option!\n");
	    return 1;
	  }

	break;
      
      case 'c':
	fprintf(logfile,"Targeted country: ");
	if(o==0)
	  {
	    target_country_flag=0;
	    fprintf(logfile,"China only\n");
	  }
	else if(o==1)
	  {
	    target_country_flag=1;
	    fprintf(logfile,"Rest of world\n");
	  }
	else if(o==2)
	  {
	    target_country_flag=2;
	    fprintf(logfile,"China + rest of world\n");
	  }
	else
	  {
	    fprintf(logfile,"Invalid command-line option!\n");
	    return 1;
	  }

	break;

      case 'r':
	fprintf(logfile,"Retaliation: ");
	if(o==0)
	  {
	    retaliation_flag=0;
	    fprintf(logfile,"No\n");
	  }
	else if(o==1)
	  {
	    retaliation_flag=1;
	    fprintf(logfile,"Yes\n");
	  }
	else
	  {
	    fprintf(logfile,"Invalid command-line option!\n");
	    return 1;
	  }

	break;

      case 'd':
	fprintf(logfile,"Duration: ");
	if(o==0)
	  {
	    duration_flag=0;
	    fprintf(logfile,"Permanent\n");
	  }
	else if(o==1)
	  {
	    duration_flag=1;
	    fprintf(logfile,"Transitory\n");
	  }
	else
	  {
	    fprintf(logfile,"Invalid command-line option!\n");
	    return 1;
	  }

	break;
      

      case 'a':
	fprintf(logfile,"Adjustment costs: ");
	if(o==0)
	  {
	    adjustment_flag=0;
	    fprintf(logfile,"None\n");
	  }
	else if(o==1)
	  {
	    adjustment_flag=1;
	    l_adj_cost=1;
	    fprintf(logfile,"Labor only\n");
	  }
	else if(o==2)
	  {
	    adjustment_flag=2;
	    k_adj_cost=1;
	    fprintf(logfile,"Capital only\n");
	  }
	else if(o==3)
	  {
	    adjustment_flag=3;
	    f_adj_cost=1;
	    m_adj_cost=1;
	    fprintf(logfile,"Supply chains only\n");
	  }
	else if(o==4)
	  {
	    adjustment_flag=4;
	    f_adj_cost=1;
	    m_adj_cost=1;
	    k_adj_cost=1;
	    l_adj_cost=1;
	    fprintf(logfile,"Labor + capital + supply chains\n");
	  }
	else
	  {
	    fprintf(logfile,"Invalid command-line option!\n");
	    return 1;
	  }
	break;

      default:
	fprintf(logfile,"Invalid command-line option!\n");
	return 1;      
      }
    }

  return 0;
}

int quant_exercise()
{
  // -----------------------------------------------------------------------------------------------------------
  // set up variable and parameter structures

  uint it;
  if(calibrate())
    {
      fprintf(logfile, "\nProgram failed!\n");
      return 1;
    }
  
  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      init_vars(&(eee1[it]));
    }

  // -------------------------------------------------------------------------------------------------------
  // free trade

  fprintf(logfile, "----------------------------------------------------------------------\n");

  scenario = 0;
  int fixl_ = fixl;
  fixl=1;
  set_neqm();

  fprintf(logfile,"\nSolving for free-trade benchmark...\n");
  if(solve_eqm())
    {
      fprintf(logfile, "\nProgram failed!\n");
      return 1;
    }
  calc_welfare(&(eee0[0]), &(ppp0[0]));

  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_usa",0);
  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_chn",1);
  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_row",2);
  
  fixl = fixl_;

  // ------------------------------------------------------------------------
  // preference params
  eqm * e = &(eee0[0]);
  int tn=0;
  for(tn=0; tn<NTH; tn++)
    {
      params * p = &(ppp0[tn]);
      int t = 0;
      int i=0;
      for(i=0; i<NC; i++)
	{
	  if(fixl==1)
	    {
	      p->phi[i]=1;
	    }
	  else if(ghh_prefs==1)
	    {
	      p->phi[i] = e->w_t[t][i] / (1.0/3.0);
	      
	      double tmp = ( pow(e->w_t[t][i] / p->phi[i], 1.0)
			    - e->ll_t[t][i]/p->lbar[i] );
	      
	      if(fabs(tmp)>TINYSQ)
		{
		  fprintf(logfile,"HH intratemp FOC 2 for country %d, error = %0.16f\n",i,tmp);
		  exit(1);
		}
	      
	    }
	  else
	    {
	      double Ctmp = (p->eps[i][0][0]*pow(e->c_t[t][i][0], p->rho) +
			     p->eps[i][0][1]*pow(e->c_t[t][i][1], p->rho) + 
			     p->eps[i][0][2]*pow(e->c_t[t][i][2], p->rho));
	  
	      double Ltmp = (p->lbar[i] - e->ll_t[t][i])*
		p->eps[i][0][0] * pow(e->c_t[t][i][0], p->rho-1.0);
	      
	      double tmp = e->p_t[t][i][0]*Ctmp/Ltmp/e->w_t[t][i];
	      p->phi[i] = tmp/(1.0+tmp);
	  
	      tmp = muc(e->c_t[t][i],e->ll_t[t][i],p->lbar[i],p->eps[i][0],p->rho,p->phi[i],p->psi,2)/e->p_t[t][i][2];

	      tmp = tmp - mul(e->c_t[t][i],e->ll_t[t][i],p->lbar[i],p->eps[i][0],p->rho,p->phi[i],p->psi) / e->w_t[t][i];
	      
	      if(fabs(tmp)>TINYSQ)
		{
		  fprintf(logfile,"HH intratemp FOC 2 for country %d, error = %0.16f\n",i,tmp);
		  exit(1);
		}
	    }
	}
      
    }

  // -------------------------------------------------------------------------------------------------------
  // Tariffs
  
  scenario = 1;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs(&(ppp1[it]),tau,scenario);
    }

  fprintf(logfile,"\nSolving for equilibrium with repeal...\n");
  if(solve_eqm())
    {
      fprintf(logfile, "\nProgram failed!\n");
      return 1;
    }


  calc_welfare(&(eee1[0]), &(ppp1[0]));

  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_usa",0);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_chn",1);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_row",2);

  // -------------------------------------------------------------------------------------------------------
  // 1-period average SR trade elasticity
  double tmp = eee1[0].te_t[TSHOCK][0];
  fprintf(logfile,"\nAverage SR trade elasticity: %0.4f\n",tmp);
  
  return 0;
}

int main(int argc, char * argv[])
{
  par = 0;
  solver_verbose=1;
  cobb_douglas_flag=0;
  cobb_douglas_flag2=0;
  f_adj_cost=0;
  m_adj_cost=0;
  k_adj_cost=0;
  l_adj_cost=0;
  fixl=0;
  ghh_prefs=1;
  eval_eqm_once_flag=0;
  eval_bgp_once_flag=0;
  read_seed=0;
  write_seed=1;
  logfile = stdout;

  fprintf(logfile,  "\nTariffs, Manufacturing Employment, and Supply-Chain Adjustment Frictions");
  fprintf(logfile,  "\nJoseph Steinberg, University of Toronto");
  fprintf(logfile, "\n");

  fprintf(logfile, "\n----------------------------------------------------------------------\n");
  fprintf(logfile, "\nSetting up environment...\n\n");

  if(parse_args(argc,argv))
    {
      fprintf(logfile, "\nProgram failed!\n");
      return 1;
    };

  // -----------------------------------------------------------------------
  // set up parallel environment
#ifdef _OPENMP
  omp_set_num_threads(NTH);
  uint nt = omp_get_max_threads();
  fprintf(logfile, "\n\tParallel processing using %d OMP threads\n",nt);
  
#pragma omp parallel num_threads(nt)
  {
    int it = omp_get_thread_num();
    fprintf(logfile,"\t\tHello from thread %d out of %d\n",
	    it, nt);
  }
  fprintf(logfile,"\n");
#endif

  quant_exercise();

  // ------------------------------------------------------------------------

  fprintf(logfile, "\n----------------------------------------------------------------------\n");
  fprintf(logfile,  "\nProgram complete!\n\n");

  return 0;
}

#endif
