#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals.h"
#include "calibrate.h"
#include "eqm.h"

uint parse_args(int argc, char **argv)
{
  int opt = 0;
  uint cnt=0;

  fprintf(logfile,KBLU "Model scenario:" RESET);
  
  while((opt = getopt(argc, argv, "abcdef")) != -1){
    switch(opt){
    case 'a':
      scenario2=1;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Tariffs end after Trump's term" RESET);
      
      break;
    case 'b':
      scenario2 = 2;
      sens=1;
      f_adj_cost=0;
      m_adj_cost=0;
      cnt=cnt+1;
      fprintf(logfile,KBLU " No import adj costs" RESET);
      break;
    case 'c':
      scenario2=3;
      f_adj_cost=0;
      m_adj_cost=0;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Transitory + no adj costs" RESET);
      break;
      /*
    case 'd':
      scenario2=4;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Retaliation by China, Canada + Mexico" RESET);
      break;
    case 'e':
      scenario2=5;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Add tariffs on EU" RESET);
      break;
    case 'f':
      scenario2=6;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Retaliation by China, Canada, Mexico, and EU" RESET);
      break;
      */
    }
  }

  if(cnt>1)
    {
      fprintf(logfile,KRED "\nOnly one sensitivity analysis allowed at a time!\n" RESET);
      return 1;
    }
  else
    {
      if(scenario2==0)
	{
	  fprintf(logfile,KBLU " Baseline analysis" RESET);
	  write_seed = 1;
	}
	
      fprintf(logfile,"\n");
  
      return 0;
    }
}

int quant_exercise()
{
  // -----------------------------------------------------------------------------------------------------------
  // set up variable and parameter structures

  //params * p = &(ppp0[0]);
  
  uint it;
  if(calibrate())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  
  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      init_vars(&(eee1[it]));
    }

  // -------------------------------------------------------------------------------------------------------
  // no repeal

  fprintf(logfile, KNRM "----------------------------------------------------------------------\n" RESET);

  scenario = 0;
  int fixl_ = fixl;
  fixl=1;
  set_neqm();

  fprintf(logfile,KBLU "\nSolving for no-repeal equilibrium...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee0[0]), &(ppp0[0]));

  if(scenario2==0)
    {
      write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_usa",0);
      write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_can",1);
      write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_mex",2);
      write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_chn",3);
      write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_eu",4);
      write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_row",5);
    }

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
		  fprintf(logfile,KRED "HH intratemp FOC 2 for country %d, error = %0.16f\n" RESET,i,tmp);
		  exit(1);
		}
	      
	    }
	  else
	    {
	      double Ctmp = (p->eps[i][0][0]*pow(e->c_t[t][i][0], p->rho) +
			     p->eps[i][0][1]*pow(e->c_t[t][i][1], p->rho) + 
			     p->eps[i][0][2]*pow(e->c_t[t][i][2], p->rho) + 
			     p->eps[i][0][3]*pow(e->c_t[t][i][3], p->rho) + 
			     p->eps[i][0][4]*pow(e->c_t[t][i][4], p->rho));
	  
	      double Ltmp = (p->lbar[i] - e->ll_t[t][i])*
		p->eps[i][0][0] * pow(e->c_t[t][i][0], p->rho-1.0);
	      
	      double tmp = e->p_t[t][i][0]*Ctmp/Ltmp/e->w_t[t][i];
	      p->phi[i] = tmp/(1.0+tmp);
	  
	      tmp = muc(e->c_t[t][i],e->ll_t[t][i],p->lbar[i],p->eps[i][0],p->rho,p->phi[i],p->psi,2)/e->p_t[t][i][2];

	      tmp = tmp - mul(e->c_t[t][i],e->ll_t[t][i],p->lbar[i],p->eps[i][0],p->rho,p->phi[i],p->psi) / e->w_t[t][i];
	      
	      if(fabs(tmp)>TINYSQ)
		{
		  fprintf(logfile,KRED "HH intratemp FOC 2 for country %d, error = %0.16f\n" RESET,i,tmp);
		  exit(1);
		}
	    }
	}
      
    }

  // -------------------------------------------------------------------------------------------------------
  // With repeal
  scenario = 1;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs(&(ppp1[it]),scenario);
    }

  fprintf(logfile,KBLU "\nSolving for equilibrium with repeal...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }


  calc_welfare(&(eee1[0]), &(ppp1[0]));

  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_usa",0);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_can",1);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_mex",2);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_chn",3);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_eu",4);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_row",5);

  // -------------------------------------------------------------------------------------------------------
  // 1-period average SR trade elasticity
  //double tmp = (eee1[0].te_t[TNAFTA][0] + eee1[0].te_t[TNAFTA][1] + eee1[0].te_t[TNAFTA][2])/3.0;
  //fprintf(logfile,KBLU "\nAverage SR trade elasticity: %0.4f\n" RESET,tmp);
  
  return 0;
}

int main(int argc, char * argv[])
{
  noio_flag=0;
  sens=0;
  dom_con_flag=0;
  par = 0;
  iceberg=1;
  solver_verbose=1;
  cobb_douglas_flag=0;
  cobb_douglas_flag2=0;
  f_adj_cost=1;
  m_adj_cost=1;
  k_adj_cost=1;
  l_adj_cost=1;
  fixl=1;
  ghh_prefs=1;
  scenario2=0;
  iceberg_flag=0;
  eval_eqm_once_flag=0;
  eval_bgp_once_flag=0;
  read_seed=0;
  write_seed=0;
  camta_flag=0;
  ucta_flag=0;
  sym_te_flag=0;
  us_ht_flag=0;
  us_ht_flag2=0;
  ltp_flag=0;
  old_mfn_flag=0;
  fix_tb_flag=0;
  fix_tb_flag2=0;
  old_iomat_flag=0;
  eqkappa=0;
  nokappa=0;
  fix_k_flag=0;
  no_k_flag=0;
  eval_calfn_once=0;
  logfile = stdout;

  fprintf(logfile, KGRN "\nThe Macroeconomic Effects of NAFTA Repeal" RESET);
  fprintf(logfile, KGRN "\nJoseph Steinberg, University of Toronto" RESET);
  fprintf(logfile, KNRM "\n" RESET);

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KBLU "\nSetting up environment...\n\n" RESET);

  if(parse_args(argc,argv))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    };

  // -----------------------------------------------------------------------
  // set up parallel environment
#ifdef _OPENMP
  omp_set_num_threads(NTH);
  //mkl_set_num_threads(NTH);
  uint nt = omp_get_max_threads();
  fprintf(logfile, KBLU "\n\tParallel processing using %d OMP threads\n" RESET,nt);
  
#pragma omp parallel num_threads(nt)
  {
    int it = omp_get_thread_num();
    fprintf(logfile,KBLU "\t\tHello from thread %d out of %d\n" RESET,
	    it, nt);
  }
  fprintf(logfile,"\n");
#endif

  quant_exercise();

  // ------------------------------------------------------------------------

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KGRN "\nProgram complete!\n\n" RESET);

  return 0;
}

#endif
