#ifndef __EQM_C__
#define __EQM_C__

#include "eqm.h"

const uint nbgp = NC+4*NC*NS + 4*NC*NS*NC;
const uint nbgp_no_k = NC+3*NC*NS + 4*NC*NS*NC;
const uint nbgp_noio = NC+3*NC*NS + 2*NC*NS*NC;

void set_neqm()
{
  // key eqm vars used in solver: wage (NC), sectoral investment (NC*NS, but not in last period),
  // sectoral labor (NC*NS), gross output prices (NC*NS), sectoral consumption (NC*NS), 
  // bonds (NC-1, but not in first period), bond price (1, but not in last period)

  // deterministic no-Brexit counterfactual
  if(scenario==0)
    {
      if(noio_flag)
	{
	  neqm = (NT+1)*NC + NT*NC*NS + 4*(NT+1)*NC*NS + (NT+1)*NC*NC*NS*2 + NT*(NC-1) + NT;
	}
      else if(no_k_flag)
	{
	  neqm = (NT+1)*NC + 4*(NT+1)*NC*NS + (NT+1)*NC*NC*NS*4 + NT*(NC-1) + NT;
	}
      else
	{
	  neqm = (NT+1)*NC + NT*NC*NS + 5*(NT+1)*NC*NS + (NT+1)*NC*NC*NS*4 + NT*(NC-1) + NT;
	}
    }
  else if(scenario>=1)
    {
      uint nn = NT+1-TNAFTA;
      if(noio_flag)
	{
	  neqm = (nn)*NC + (nn-1)*NC*NS + 4*(nn)*NC*NS + 2*nn*NC*NC*NS + (nn-1)*(NC-1) + (nn-1);
	}
      else if(no_k_flag)
	{
	  neqm = (nn)*NC + 4*(nn)*NC*NS + 4*nn*NC*NC*NS + (nn-1)*(NC-1) + (nn-1);
	}
      else
	{
	  neqm = (nn)*NC + (nn-1)*NC*NS + 5*(nn)*NC*NS + 4*nn*NC*NC*NS + (nn-1)*(NC-1) + (nn-1);
	}
    }
}

void init_vars(eqm * e)
{
  SET_ALL_V(e->b_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cc_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ii_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ll_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->kk_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->cpi_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->pi_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->w_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->rk_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->ngdp_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->rgdp_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->iy_t,(NT+1)*NC,0.0);
  SET_ALL_V(e->lp_agg_t,(NT+1)*NC,0.0);

  SET_ALL_V(e->ex_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->im_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->nx_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->exf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->imf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->nxf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->exm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->imm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->nxm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rer_t,(NT+1)*NC*NC,0.0);

  SET_ALL_V(e->rex_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rim_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rexf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rimf_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rexm_t,(NT+1)*NC*NC,0.0);
  SET_ALL_V(e->rimm_t,(NT+1)*NC*NC,0.0);

  SET_ALL_V(e->y_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->va_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->md_t,(NT+1)*NC*NS*NS,0.0);
  SET_ALL_V(e->k_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->l_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->is_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->lp_t,(NT+1)*NC*NS,0.0);

  SET_ALL_V(e->exs_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->ims_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->nxs_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->rexs_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->rims_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->pm_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->m_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->m2_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->pm2_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->p_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->q_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->q2_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->pq2_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->c_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->i_t,(NT+1)*NC*NS,0.0);

  SET_ALL_V(e->Df_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->Dbf_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->Dhf_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->Dtf_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->Dm_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->Dbm_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->Dhm_t,(NT+1)*NC*NS*NC,0.0);
  SET_ALL_V(e->Dtm_t,(NT+1)*NC*NS*NC,0.0);

  SET_ALL_V(e->exrate_t,(NT+1)*NC*NS*(NC),0.0);
  SET_ALL_V(e->exitrate_t,(NT+1)*NC*NS*(NC),0.0);

  uint t;
  for(t=0; t<NT+1; t++)
    {
      uint i;
      for(i=0; i<NC; i++)
	{
	  uint s;
	  for(s=0; s<NS; s++)
	    {
	      uint ii;
	      for(ii=0; ii<2; ii++)
		{
		  e->ev_t[t][i][s][ii].zp=0.0;
		  e->ev_t[t][i][s][ii].zm=0.0;
		  e->ev_t[t][i][s][ii].Zi=0.0;
		  e->ev_t[t][i][s][ii].Zp=0.0;
		  e->ev_t[t][i][s][ii].Zm=0.0;
		  e->ev_t[t][i][s][ii].Z=0.0;
		  e->ev_t[t][i][s][ii].Fzp=0.0;
		  e->ev_t[t][i][s][ii].Fzm=0.0;
		  e->ev_t[t][i][s][ii].dV=0.0;
		  e->ev_t[t][i][s][ii].n=0.0;
		}
	    }
	}
    }

}

void copy_vars(eqm * e1, const eqm * e0)
{
  memcpy((double *)(e1->pb_t),(const double *)(e0->pb_t),sizeof(double)*(NT+1));

  memcpy((double *)(e1->b_t),(const double *)(e0->b_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cc_t),(const double *)(e0->cc_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ii_t),(const double *)(e0->ii_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ll_t),(const double *)(e0->ll_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->kk_t),(const double *)(e0->kk_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->cpi_t),
	 (const double *)(e0->cpi_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->pi_t),(const double *)(e0->pi_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->w_t),(const double *)(e0->w_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->rk_t),(const double *)(e0->rk_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->ngdp_t),
	 (const double *)(e0->ngdp_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->rgdp_t),
	 (const double *)(e0->rgdp_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->iy_t),(const double *)(e0->iy_t),sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->lp_agg_t),
	 (const double *)(e0->lp_agg_t),
	 sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->ex_t),
	 (const double *)(e0->ex_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->im_t),
	 (const double *)(e0->im_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->nx_t),
	 (const double *)(e0->nx_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->exf_t),
	 (const double *)(e0->exf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->imf_t),
	 (const double *)(e0->imf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->nxf_t),
	 (const double *)(e0->nxf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->exm_t),
	 (const double *)(e0->exm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->imm_t),
	 (const double *)(e0->imm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->nxm_t),
	 (const double *)(e0->nxm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rer_t),
	 (const double *)(e0->rer_t),
	 sizeof(double)*(NT+1)*NC*NC);

  memcpy((double *)(e1->rex_t),
	 (const double *)(e0->rex_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rim_t),
	 (const double *)(e0->rim_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rexf_t),
	 (const double *)(e0->rexf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rimf_t),
	 (const double *)(e0->rimf_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rexm_t),
	 (const double *)(e0->rexm_t),
	 sizeof(double)*(NT+1)*NC*NC);
  memcpy((double *)(e1->rimm_t),
	 (const double *)(e0->rimm_t),
	 sizeof(double)*(NT+1)*NC*NC);


  memcpy((double *)(e1->y_t),
	 (const double *)(e0->y_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->py_t),
	 (const double *)(e0->py_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->va_t),
	 (const double *)(e0->va_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->k_t),
	 (const double *)(e0->k_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->l_t),
	 (const double *)(e0->l_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->is_t),
	 (const double *)(e0->is_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->rks_t),
	 (const double *)(e0->rks_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->ws_t),
	 (const double *)(e0->ws_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->md_t),
	 (const double *)(e0->md_t),
	 sizeof(double)*(NT+1)*NC*NS*NS);
  memcpy((double *)(e1->lp_t),
	 (const double *)(e0->lp_t),
	 sizeof(double)*(NT+1)*NC*NS);

  memcpy((double *)(e1->exs_t),
	 (const double *)(e0->exs_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->ims_t),
	 (const double *)(e0->ims_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->nxs_t),
	 (const double *)(e0->nxs_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->rexs_t),
	 (const double *)(e0->rexs_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->rims_t),
	 (const double *)(e0->rims_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->pm_t),
	 (const double *)(e0->pm_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->m_t),
	 (const double *)(e0->m_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->m2_t),
	 (const double *)(e0->m2_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->pm2_t),
	 (const double *)(e0->pm2_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->p_t),
	 (const double *)(e0->p_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->q_t),
	 (const double *)(e0->q_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->q2_t),
	 (const double *)(e0->q2_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->pq2_t),
	 (const double *)(e0->pq2_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->c_t),
	 (const double *)(e0->c_t),
	 sizeof(double)*(NT+1)*NC*NS);
  memcpy((double *)(e1->i_t),
	 (const double *)(e0->i_t),
	 sizeof(double)*(NT+1)*NC*NS);

  memcpy((double *)(e1->welfare_t),
	 (const double *)(e0->welfare_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare2_t),
	 (const double *)(e0->welfare2_t),
	 sizeof(double)*(NT+1)*NC);
  memcpy((double *)(e1->welfare_cost_t),
	 (const double *)(e0->welfare_cost_t),
	 sizeof(double)*(NT+1)*NC);

  memcpy((double *)(e1->Df_t),
	 (const double *)(e0->Df_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->Dbf_t),
	 (const double *)(e0->Dbf_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->Dhf_t),
	 (const double *)(e0->Dhf_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->Dtf_t),
	 (const double *)(e0->Dtf_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->Dm_t),
	 (const double *)(e0->Dm_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->Dbm_t),
	 (const double *)(e0->Dbm_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->Dhm_t),
	 (const double *)(e0->Dhm_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);
  memcpy((double *)(e1->Dtm_t),
	 (const double *)(e0->Dtm_t),
	 sizeof(double)*(NT+1)*NC*NS*NC);

  memcpy((double *)(e1->exrate_t),
	 (const double *)(e0->exrate_t),
	 sizeof(double)*(NT+1)*NC*NS*(NC-1));
  memcpy((double *)(e1->exitrate_t),
	 (const double *)(e0->exitrate_t),
	 sizeof(double)*(NT+1)*NC*NS*(NC-1));


  uint t, i, s,ii;
  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      for(ii=0; ii<2; ii++)
		{
		  copy_exporter_vars(&(e1->ev_t[t][i][s][ii]),
				     (const exporter_vars *)(&(e0->ev_t[t][i][s][ii])));
		}
	    }
	}
    }
}

uint stack_bgp_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t = NT;
  
  COPY_SUBVECTOR_LOG(myx+nx,e->w_t[t],NC);
  nx=nx+NC;

  if(!no_k_flag)
    {
      COPY_SUBVECTOR_LOG(myx+nx,e->k_t[t],NC*NS);
      nx=nx+NC*NS;
    }
  
  COPY_SUBVECTOR_LOG(myx+nx,e->p_t[t],NC*NS);
  nx=nx+NC*NS;
  
  COPY_SUBVECTOR_LOG(myx+nx,e->pq2_t[t],NC*NS*NC);
  nx=nx+NC*NS*NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->q2_t[t],NC*NS*NC);
  nx=nx+NC*NS*NC;

  if(!noio_flag)
    {
      COPY_SUBVECTOR_LOG(myx+nx,e->pm_t[t],NC*NS);
      nx=nx+NC*NS;

      COPY_SUBVECTOR_LOG(myx+nx,e->pm2_t[t],NC*NS*NC);
      nx=nx+NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,e->m2_t[t],NC*NS*NC);
      nx=nx+NC*NS*NC;
    }
  
  COPY_SUBVECTOR_LOG_PLUS_ONE(myx+nx,e->c_t[t],NC*NS);
  nx=nx+NC*NS;

  if( (noio_flag && nx != nbgp_noio) || (no_k_flag && nx != nbgp_no_k) || (!noio_flag && !no_k_flag && nx != nbgp))
    {
      fprintf(logfile,KRED "Error stacking bgp vars! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }

    return 0;
}

uint unstack_bgp_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t = NT;

  copy_subvector_exp( (double *)(e->w_t[t]), myx+nx, NC);
  nx=nx+NC;

  if(!no_k_flag)
    {
      COPY_SUBVECTOR_EXP(e->k_t[t],myx+nx,NC*NS);
      nx=nx+NC*NS;
    }
  
  COPY_SUBVECTOR_EXP(e->p_t[t],myx+nx,NC*NS);
  nx= nx + NC*NS;

  COPY_SUBVECTOR_EXP(e->pq2_t[t],myx+nx,NC*NS*NC);
  nx= nx + NC*NS*NC;

  COPY_SUBVECTOR_EXP(e->q2_t[t],myx+nx,NC*NS*NC);
  nx= nx + NC*NS*NC;
  
  if(!noio_flag)
    {
      COPY_SUBVECTOR_EXP(e->pm_t[t],myx+nx,NC*NS);
      nx= nx + NC*NS;

      COPY_SUBVECTOR_EXP(e->pm2_t[t],myx+nx,NC*NS*NC);
      nx= nx + NC*NS*NC;
      
      COPY_SUBVECTOR_EXP(e->m2_t[t],myx+nx,NC*NS*NC);
      nx= nx + NC*NS*NC;
    }

  COPY_SUBVECTOR_EXP_MINUS_ONE(e->c_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  if( (noio_flag && nx != nbgp_noio) || (no_k_flag && nx != nbgp_no_k) || (!noio_flag && !no_k_flag && nx != nbgp))
    {
      fprintf(logfile,KRED "Error unstacking bgp vars! nx = %d, nbgp = %d\n" RESET,nx,nbgp);
      return 1;
    }

    return 0;
}

uint stack_eqm_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;
  
  if(scenario>=1)
    {
      nn = NT+1-TNAFTA;
      t0=TNAFTA;
    }

  if(!no_k_flag)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->rks_t[t0]),(nn)*NC*NS);
      nx = nx+(nn)*NC*NS;
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->ws_t[t0]),(nn)*NC*NS);
  nx = nx+(nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  if(!no_k_flag)
    {
      COPY_SUBVECTOR(myx+nx,&(e->is_t[t0]),(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }
  
  COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pq2_t[t0]),(nn)*NC*NS*NC);
  nx = nx + (nn)*NC*NS*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn)*NC*NS*NC);
  nx = nx + (nn)*NC*NS*NC;

  if(!noio_flag)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm2_t[t0]),(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;
    }
  
  COPY_SUBVECTOR_LOG_PLUS_ONE(myx+nx,&(e->c_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  *(myx+nx) = e->b_t[t][i];
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->pb_t[t0]),(nn-1));
  nx = nx + (nn-1);

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error stacking eqm vars! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

    return 0;
}

uint unstack_eqm_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;

  if(scenario>=1)
    {
      nn = NT+1-TNAFTA;
      t0=TNAFTA;
    }

  if(!no_k_flag)
    {
      COPY_SUBVECTOR_EXP(&(e->rks_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx+(nn)*NC*NS;
    }

  COPY_SUBVECTOR_EXP(&(e->ws_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx+(nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;

  if(!no_k_flag)
    {
      COPY_SUBVECTOR(&(e->is_t[t0]),myx+nx,(nn-1)*NC*NS);
      nx = nx + (nn-1)*NC*NS;
    }

  COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->pq2_t[t0]),myx+nx,(nn)*NC*NS*NC);
  nx = nx + (nn)*NC*NS*NC;

  COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn)*NC*NS*NC);
  nx = nx + (nn)*NC*NS*NC;

  if(!noio_flag)
    {
      COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn)*NC*NS);
      nx = nx + (nn)*NC*NS;

      COPY_SUBVECTOR_EXP(&(e->pm2_t[t0]),myx+nx,(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;

      COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn)*NC*NS*NC);
      nx = nx + (nn)*NC*NS*NC;
    }
  
  COPY_SUBVECTOR_EXP_MINUS_ONE(&(e->c_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  uint i = 0;
  uint t = 0;
  for(t=t0+1; t<(NT+1); t++)
    {
      for(i=0; i<(NC-1); i++)
	{
	  e->b_t[t][i] = *(myx+nx);
	  nx = nx+1;
	}
    }

  COPY_SUBVECTOR_EXP(&(e->pb_t[t0]),myx+nx,(nn-1));
  nx = nx + (nn-1);

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error unstacking eqm vars! nx = %d, neqm = %d\n" RESET ,nx,neqm);
      return 1;
    }

  return 0;
}

uint set_initial_bgp_guess()
{
  uint i, t, s, j;
  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);
  t=NT;
  
  for(i=0; i<NC; i++)
    {
      e->w_t[t][i] = 1.0;
      for(s=0; s<NS; s++)
	{
	  e->c_t[t][i][s] = p->c0[i][s];
	  e->l_t[t][i][s] = p->l0[i][s];
	  e->k_t[t][i][s] = p->k0[i][s];
	  e->p_t[t][i][s] = 1.0;
	  e->pm_t[t][i][s] = 1.0;
	  for(j=0; j<NC; j++)
	    {
	      e->pm2_t[t][i][s][j] = 1.0;
	      e->pq2_t[t][i][s][j] = 1.0;
	      e->m2_t[t][i][s][j] = p->m02[i][s][j];
	      e->q2_t[t][i][s][j] = p->q02[i][s][j];

	    }
	}
    }

  if(stack_bgp_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for balanced growth path!\n" RESET);
      return 1;
    }
  else
    {
      return 0;
    }
}

uint set_initial_eqm_guess()
{
  uint i,s,t;
  double bb[NC];

  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);

  bb[0] = p->b0[0];
  bb[1] = p->b0[1];
  bb[2] = p->b0[2];
  bb[3] = p->b0[3];
  bb[4] = p->b0[4];
  bb[5] = p->b0[5];

  free_solver_mem();

  if(noio_flag)
    {
      solver_n = nbgp_noio;
    }
  else
    {
      solver_n = nbgp;
    }
  
  alloc_solver_mem();
  if(solve_bgp(bb))
    {
      fprintf(logfile, KRED "Error solving for balanced growth path!\n");
      return 1;
    }
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

  for(t=0; t<NT+1; t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS; s++)
	    {
	      e->rks_t[t][i][s] = e->rks_t[NT][i][s];
	      e->ws_t[t][i][s] = e->ws_t[NT][i][s];
	    }
	}
    }

  // first construct bond guess... a little awkward to logspace this because we have to deal with
  // absolute values
  double tmpb0[NC];// = {fabs(p->b0[0]),fabs(p->b0[1]),fabs(p->b0[2]),fabs(p->b0[3])};
  double tmpb1[NC];// = {fabs(bb[0]),fabs(bb[1]),fabs(bb[2]),fabs(bb[3])};

  for(int i=0; i<NC; i++)
    {
      tmpb0[i] = fabs(p->b0[i]);
      tmpb1[i] = fabs(bb[i]);
    }
  
  LOGSPACE_2D(tmpb0,tmpb1,NT+1,NC,e->b_t);
  for(i=0; i<(NC-1); i++)
    {
      if(fabs(p->b0[i])<1.0e-6)
	{
	  for(t=0; t<(NT+1); t++)
	    {
	      e->b_t[t][i] = 0.0;
	    }
	}
      else
	{
	  if(p->b0[i] < -TINY)
	    {
	      for(t=0; t<(NT+1); t++)
		{
		  e->b_t[t][i] = -e->b_t[t][i];
		}
	    }
	}
    }
  set_all_v(e->pb_t,NT+1,e->pb_t[NT]);

  // now construct guesses for prices real variables
  double tmpp[NC];
  SET_ALL_V(tmpp,NC,1.0);
  
  double tmpp2[NC][NS];
  SET_ALL_V(tmpp2,NC*NS,1.0);

  double tmpp3[NC][NS][NC];
  SET_ALL_V(tmpp3,NC*NS*NC,1.0);
  
  LOGSPACE_2D(p->k0,e->k_t[NT],NT+1,NC*NS,e->k_t);
  LOGSPACE_2D(p->l0,e->l_t[NT],NT+1,NC*NS,e->l_t);
  LOGSPACE_2D(p->c0,e->c_t[NT],NT+1,NC*NS,e->c_t);
  LOGSPACE_2D(p->q0,e->q_t[NT],NT+1,NC*NS,e->q_t);
  LOGSPACE_2D(p->m0,e->m_t[NT],NT+1,NC*NS,e->m_t);
  LINSPACE_2D(tmpp,e->w_t[NT],NT+1,NC,e->w_t);
  LINSPACE_2D(tmpp,e->rk_t[NT],NT+1,NC,e->rk_t);
  LINSPACE_2D(tmpp2,e->pm2_t[NT],NT+1,NC*NS*NC,e->pm2_t);
  LINSPACE_2D(tmpp2,e->pq2_t[NT],NT+1,NC*NS*NC,e->pq2_t);
  LINSPACE_2D(p->m02,e->m2_t[NT],NT+1,NC*NS*NC,e->m2_t);
  LINSPACE_2D(p->q02,e->q2_t[NT],NT+1,NC*NS*NC,e->q2_t);
  LINSPACE_2D(tmpp2,e->p_t[NT],NT+1,NC*NS,e->p_t);
  LINSPACE_2D(tmpp2,e->pm_t[NT],NT+1,NC*NS,e->pm_t);

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(t=0; t<NT; t++)
	    {
	      e->is_t[t][i][s] = e->k_t[t+1][i][s] - (1.0-p->delta)*e->k_t[t][i][s];
	    }
	}
    }

  if(stack_eqm_vars(solver_x->data,e))
    {
      fprintf(logfile,KRED "Failed to create guess for balanced growth path!\n" RESET);
      return 1;
    }
  else
    {      
      return 0;
    }
}

uint write_eqm_vars(const eqm * e, const params * p, char * fname, uint i)
{
  char * fname2 = concat("output/",fname);

  char * fname3;

  if(scenario2==1)
    {
      fname3 = concat(fname2,"_s1.csv");
    }
  else if(scenario2==2)
    {
      fname3 = concat(fname2,"_s2.csv");
    }
  else if(scenario2==3)
    {
      fname3 = concat(fname2,"_s3.csv");
    }
    else
    {
      fname3 = concat(fname2,".csv");
    }

  /*
  if(scenario2==4)
    {
      fname3 = concat(fname2,"_s4.csv");
    }
  if(scenario2==5)
    {
      fname3 = concat(fname2,"_s5.csv");
    }
  if(scenario2==6)
    {
      fname3 = concat(fname2,"_s6.csv");
    }
  */


  FILE * file = fopen(fname3,"w");
  free(fname2);

  if(file)
    {
      uint s,j,t;
      fprintf(file,"period,rgdp,ngdp,ii,ll,c,W,ae,te");
      for(s=0;s<NS;s++)
	{
	  fprintf(file,",y%d,va%d,i%d,c%d,l%d,py%d,ae%d,te%d",s,s,s,s,s,s,s,s);
	}
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      fprintf(file,",rer%d,ex%d,im%d,nx%d,rex%d,rim%d",
		      j,j,j,j,j,j);
	      for(s=0;s<NS;s++)
		{
		  fprintf(file,
			  ",taus%d-%d,exs%d-%d,ims%d-%d,nxs%d-%d,rexs%d-%d,rims%d-%d,exsm%d-%d,exsf%d-%d,imsm%d-%d,imsf%d-%d,aes%d-%d,tes%d-%d,expart%d-%d",s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j);
		}
	    }
	}
      fprintf(file,"\n");

      for(t=0;t<(NT+1);t++)
	{
	  fprintf(file,"%d,",t);
	  fprintf(file,"%0.16f,",e->rgdp_t[t][i]);
	  fprintf(file,"%0.16f,",e->ngdp_t[t][i]);
	  fprintf(file,"%0.16f,",e->ii_t[t][i]);
	  fprintf(file,"%0.16f,",e->ll_t[t][i]);
	  fprintf(file,"%0.16f,",e->cc_t[t][i]);
	  fprintf(file,"%0.16f,",e->welfare_t[t][i]);
	  fprintf(file,"%0.16f,",e->ae_t[t][i]);
	  fprintf(file,"%0.16f",e->te_t[t][i]);
	  for(s=0; s<NS; s++)
	    {
	      fprintf(file,",%0.16f,",e->y_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->va_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->is_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->c_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->l_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->py_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->aes_t[t][i][s]);
	      fprintf(file,"%0.16f",e->tes_t[t][i][s]);
	    }
	  for(j=0; j<NC; j++)
	    {
	      if(j!=i)
		{
		  fprintf(file,",%0.16f,",e->rer_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->ex_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->im_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->nx_t[t][i][j]);
		  fprintf(file,"%0.16f,",e->rex_t[t][i][j]);
		  fprintf(file,"%0.16f",e->rim_t[t][i][j]);
		  for(s=0; s<NS; s++)
		    {
		      fprintf(file,",%0.16f,",p->tau_m_ts[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->exs_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->ims_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->nxs_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->q2_t[t][j][s][i]+e->m2_t[t][j][s][i]);
		      fprintf(file,"%0.16f,",e->q2_t[t][i][s][j]+e->m2_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->m2_t[t][j][s][i]*e->py_t[t][i][s]);
		      fprintf(file,"%0.16f,",e->q2_t[t][j][s][i]*e->py_t[t][i][s]);
		      fprintf(file,"%0.16f,",e->m2_t[t][i][s][j]*e->py_t[t][j][s]);
		      fprintf(file,"%0.16f,",e->q2_t[t][i][s][j]*e->py_t[t][j][s]);
		      fprintf(file,"%0.16f,",e->aes2_t[t][i][s][j]);
		      fprintf(file,"%0.16f,",e->tes2_t[t][i][s][j]);
		      fprintf(file,"%0.16f",e->exrate_t[t][i][s][j]);
		    }
		}
	    }
	  fprintf(file,"\n");
	}
      
      fclose(file);
      return 0;
    }
  else
    {
      fprintf(logfile,KRED "Error opening file to write equilibrium vars!\n" RESET);
      return 1;
    }
}

uint set_vars1(eqm * e, const params * p, uint t, uint bgp)
{
  uint i,s,j,r;

  SET_ALL_V(e->ngdp_t[t],NC,0.0);
  SET_ALL_V(e->rgdp_t[t],NC,0.0);
  SET_ALL_V(e->ex_t[t],NC*NC,0.0);
  SET_ALL_V(e->im_t[t],NC*NC,0.0);
  SET_ALL_V(e->nx_t[t],NC*NC,0.0);
  SET_ALL_V(e->exf_t[t],NC*NC,0.0);
  SET_ALL_V(e->imf_t[t],NC*NC,0.0);
  SET_ALL_V(e->nxf_t[t],NC*NC,0.0);
  SET_ALL_V(e->exm_t[t],NC*NC,0.0);
  SET_ALL_V(e->imm_t[t],NC*NC,0.0);
  SET_ALL_V(e->nxm_t[t],NC*NC,0.0);
  SET_ALL_V(e->exs_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->ims_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->nxs_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->rex_t[t],NC*NC,0.0);
  SET_ALL_V(e->rim_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexf_t[t],NC*NC,0.0);
  SET_ALL_V(e->rimf_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexm_t[t],NC*NC,0.0);
  SET_ALL_V(e->rimm_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexs_t[t],NC*NS*NC,0.0);
  SET_ALL_V(e->rims_t[t],NC*NS*NC,0.0);

  // bond market clearing
  e->b_t[t][5] = -(e->b_t[t][0]+e->b_t[t][1]+e->b_t[t][2]+e->b_t[t][3]+e->b_t[t][4]);

  // compute sector-level aggregates
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  e->q_t[t][i][s] = prod_q(e->q2_t[t][i][s], p->H[i][s], p->theta[i][s], p->sig[i][s]);

	  if(!noio_flag)
	    {
	      e->m_t[t][i][s] = prod_m(e->m2_t[t][i][s], p->M[i][s], p->mu[i][s], p->zeta[i][s]);
	    }
	}
    }
  
  // investment and labor stuff
  for(i=0; i<NC; i++)
    {
      if(!no_k_flag)
	{
	  if(t<NT)
	    {
	      for(s=0; s<NS; s++)
		{
		  e->k_t[t+1][i][s] = (1.0-p->delta) * e->k_t[t][i][s] + e->is_t[t][i][s];
		}      
	    }
	  else
	    {
	      for(s=0; s<NS; s++)
		{
		  e->is_t[t][i][s] = p->delta * e->k_t[t][i][s];
		}
	    }
      
	  e->kk_t[t][i] = e->k_t[t][i][0] + 
	    e->k_t[t][i][1] + 
	    e->k_t[t][i][2] + 
	    e->k_t[t][i][3] +
	    e->k_t[t][i][4];

	  e->ii_t[t][i] = e->is_t[t][i][0] + 
	    e->is_t[t][i][1] + 
	    e->is_t[t][i][2] + 
	    e->is_t[t][i][3] + 
	    e->is_t[t][i][4];
    
	  e->pi_t[t][i] = 1.0/p->G[i];
	  for(s=1;s<NS; s++)
	    {
	      e->pi_t[t][i] = e->pi_t[t][i] * pow(e->p_t[t][i][s]/p->eps[i][1][s],p->eps[i][1][s]);
	    }
	}
      
      e->cpi_t[t][i] = 0.0;
      for(s=0; s<NS; s++)
	{
	  e->cpi_t[t][i] = e->cpi_t[t][i] + e->p_t[t][i][s]*p->c0[i][s];
	}
      e->cpi_t[t][i] = e->cpi_t[t][i]/SUM(p->c0[i],NS);

      if(t>0)
	{
	  double mutp = muc(
		e->c_t[t][i],
		e->ll_t[t][i],
		p->lbar[i],
		p->eps[i][0],
		p->rho,
		p->phi[i],
		p->psi,
		2);
	  double mut = muc(
		  e->c_t[t-1][i],
		  e->ll_t[t-1][i],
		  p->lbar[i],
		  p->eps[i][0],
		  p->rho,
		  p->phi[i],
		  p->psi,
		  2);

	  e->Q_t[t-1][i] = p->beta[i] * (mutp / e->p_t[t][i][2]) / (mut / e->p_t[t-1][i][2]);
	}

      if(t==NT)
	{
	  e->Q_t[t][i] = p->beta[i];
	}

      if(t == 0)
	{
	  e->rk_t[t][i] = p->r0[i] + p->delta;
	}
      else if(t==NT)
	{
	  if(i==0)
	    {
	      e->pb_t[t] = e->cpi_t[t][i]/(1.0+p->rss);
	    }

	  if(bgp)
	    {
	      e->rk_t[t][i] = e->pi_t[t][i]*e->cpi_t[t][0]/e->pb_t[t] - (1.0-p->delta)*e->pi_t[t][i];
	    }
	  else
	    {
	      e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->Q_t[t-1][i] - (1.0-p->delta)*e->pi_t[t][i];
	    }  
	}
      else
	{
	  e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->Q_t[t-1][i] - (1.0-p->delta)*e->pi_t[t][i];
	}

      e->rk_t[t][i] = e->rk_t[t][i]/(1.0-p->tauk[i]);

      if(bgp && t==NT)
	{
	  for(s=0; s<NS; s++)
	    {
	      e->rks_t[t][i][s] = e->rk_t[t][i];
	      e->ws_t[t][i][s] = e->w_t[t][i];
	    }
	}
    }

    // compute trade flows
  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      e->rer_t[t][i][j] = e->cpi_t[t][j]/e->cpi_t[t][i];
	      
	      for(s=0; s<NS; s++)
		{
		  double m = 0.0;
		  if(!noio_flag)
		    {
		      m = e->pm2_t[t][j][s][i]*e->m2_t[t][j][s][i];
		    }
		  double f = e->pq2_t[t][j][s][i]*e->q2_t[t][j][s][i];
	
		  e->exs_t[t][i][s][j] = m+f;
		  e->exm_t[t][i][j] = e->exm_t[t][i][j] + m;
		  e->exf_t[t][i][j] = e->exf_t[t][i][j] + f;
		  e->ex_t[t][i][j] = e->ex_t[t][i][j] + m+f;

		  m = 0.0;
		  if(!noio_flag)
		    {
		      m = e->m2_t[t][j][s][i];
		    }
		  f = e->q2_t[t][j][s][i];
		  e->rexs_t[t][i][s][j] = m+f;
		  e->rexm_t[t][i][j] = e->rexm_t[t][i][j] + m;
		  e->rexf_t[t][i][j] = e->rexf_t[t][i][j] + f;
		  e->rex_t[t][i][j] = e->rex_t[t][i][j] + m+f;

		  m = 0.0;
		  if(!noio_flag)
		    {
		      m = e->pm2_t[t][i][s][j]*e->m2_t[t][i][s][j];
		    }
		  f = e->pq2_t[t][i][s][j]*e->q2_t[t][i][s][j];
		  e->ims_t[t][i][s][j] = m+f;
		  e->imm_t[t][i][j] = e->imm_t[t][i][j] + m;
		  e->imf_t[t][i][j] = e->imf_t[t][i][j] + f;
		  e->im_t[t][i][j] = e->im_t[t][i][j] + m+f;

		  m = 0.0;
		  if(!noio_flag)
		    {
		      m = e->m2_t[t][i][s][j];
		    }
		  f = e->q2_t[t][i][s][j];
		  e->rims_t[t][i][s][j] = m+f;
		  e->rimm_t[t][i][j] = e->rimm_t[t][i][j] + m;
		  e->rimf_t[t][i][j] = e->rimf_t[t][i][j] + f;
		  e->rim_t[t][i][j] = e->rim_t[t][i][j] + m+f;

		  e->nxs_t[t][i][s][j] = e->exs_t[t][i][s][j] - e->ims_t[t][i][s][j];
		}
	      e->nxm_t[t][i][j] = e->exm_t[t][i][j] - e->imm_t[t][i][j];
	      e->nxf_t[t][i][j] = e->exf_t[t][i][j] - e->imf_t[t][i][j];
	      e->nx_t[t][i][j] = e->ex_t[t][i][j] - e->im_t[t][i][j];
	    }
	  else
	    {
	      e->rer_t[t][i][j] = 1.0;
	    }
	}
    }

    for(i=0; i<NC; i++)
    {
      e->ae_t[t][i] = 0.0;
      e->te_t[t][i] = 0.0;
      double w = 0.0;
      double wt = 0.0;

      for(s=0; s<NS; s++)
	{
	  e->aes_t[t][i][s] = 0.0;
	  e->tes_t[t][i][s] = 0.0;
	  double w2 = 0.0;
	  double w2t = 0.0;

	  for(j=0; j<NC; j++)
	    {
	      if(j!=i)
		{
		  double yt = log(e->m2_t[t][i][s][j] / e->m2_t[t][i][s][i]);
		  double xt = log((e->pm2_t[t][i][s][i]/((1.0+p->tau_m_ts[t][i][s][j])*e->pm2_t[t][i][s][j])));
		  double a0 = (1.0/(1.0-p->zeta[i][s]))*log(p->mu[i][s][j]/p->mu[i][s][i]);
		  e->aes2_m_t[t][i][s][j] = (yt-a0)/xt;

		  yt = log(e->q2_t[t][i][s][j] / e->q2_t[t][i][s][i]);
		  xt = log((e->pq2_t[t][i][s][i]/((1.0+p->tau_f_ts[t][i][s][j])*e->pq2_t[t][i][s][j])));
		  a0 = (1.0/(1.0-p->sig[i][s]))*log(p->theta[i][s][j]/p->theta[i][s][i]);
		  e->aes2_f_t[t][i][s][j] = (yt-a0)/xt;

		  e->aes2_t[t][i][s][j] = ( (e->aes2_m_t[t][i][s][j]*p->m02[i][s][j] + 
					     e->aes2_f_t[t][i][s][j]*p->q02[i][s][j]) /
					    (p->m02[i][s][j] + p->q02[i][s][j]) );

		  double tmp3 = p->m02[i][s][j] + p->q02[i][s][j];
		  e->aes_t[t][i][s] = e->aes_t[t][i][s] + e->aes2_t[t][i][s][j]*tmp3;
		  w2 = w2 + tmp3;

		  if(p->tau_m_ts[t][i][s][j]>0.0)
		    {
		      e->tes2_t[t][i][s][j] = -(log(e->ims_t[t][i][s][j]/eee0[0].ims_t[t][i][s][j])/
						log(1.0+p->tau_m_ts[t][i][s][j]));
		      e->tes_t[t][i][s] = e->tes_t[t][i][s] + e->tes2_t[t][i][s][j]*tmp3;
		      w2t = w2t + tmp3;
		    }
		}
	    }
	  
	  e->aes_t[t][i][s] = e->aes_t[t][i][s]/w2;
	  e->ae_t[t][i] = e->ae_t[t][i] + e->aes_t[t][i][s]*w2;
	  w = w + w2;

	  if(w2t>1.e-6)
	    {
	      e->tes_t[t][i][s] = e->tes_t[t][i][s]/w2t;
	      e->te_t[t][i] = e->te_t[t][i] + e->tes_t[t][i][s]*w2t;
	      wt = wt + w2t;
	    }

	}
      
      e->ae_t[t][i] = e->ae_t[t][i]/w;
      e->te_t[t][i] = e->te_t[t][i]/wt;
    }

  // compute firm-level stuff
  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  if(no_k_flag)
	    {
	      e->MC_t[t][i][s] = p->lam_va[i][s]*e->ws_t[t][i][s] + DOT_PROD(p->lam[i][s],e->pm_t[t][i],NS);
	    }
	  else if(cobb_douglas_flag==0)
	    {
	      e->MC_t[t][i][s] = (p->lam_va[i][s]*(pow(e->rks_t[t][i][s]/p->alpha[i][s],p->alpha[i][s]) *
						   pow(e->ws_t[t][i][s]/(1.0-p->alpha[i][s]),1.0-p->alpha[i][s]))
				  + DOT_PROD(p->lam[i][s],e->pm_t[t][i],NS));
	      
	    }
	  else
	    {
	      e->MC_t[t][i][s] = (pow(e->rks_t[t][i][s]/(p->alpha[i][s]*p->lam_va[i][s]),p->lam_va[i][s]*p->alpha[i][s])*
				  pow(e->ws_t[t][i][s]/((1.0-p->alpha[i][s])*p->lam_va[i][s]),p->lam_va[i][s]*(1.0-p->alpha[i][s])));
				  ///pow(1.0/(p->lam_va[i][s]*p->A[i][s]),p->lam_va[i][s]));

	      if(!noio_flag)
		{
		  for(r=0; r<NS; r++)
		    {
		      e->MC_t[t][i][s] = e->MC_t[t][i][s] * pow(e->pm_t[t][i][r]/p->lam[i][s][r],p->lam[i][s][r]);
		    }
		}
	    }

	  for(j=0; j<NC; j++)
	    {
	      e->Df_t[t][j][s][i] = (pow(p->H2[j][s][i], p->eta-1.0) *
				     pow(e->pq2_t[t][j][s][i],p->eta) *
				     e->q2_t[t][j][s][i]);
	      
	      e->Dbf_t[t][j][s][i] = (pow(p->eta/(p->eta-1.0),1.0-p->eta) *
				      pow(1.0+p->ntb_f_ts[t][j][s][i],1.0-p->eta) *
				     e->Df_t[t][j][s][i] *
				     pow(e->MC_t[t][i][s],1.0-p->eta));
	      
	      e->Dhf_t[t][j][s][i] = (pow(p->eta/(p->eta-1.0),-p->eta) *
				      pow(1.0+p->ntb_f_ts[t][j][s][i],1.0-p->eta) *
				      e->Df_t[t][j][s][i] *
				      pow(e->MC_t[t][i][s],1.0-p->eta));

	      e->Dtf_t[t][j][s][i] = e->Dbf_t[t][j][s][i] / p->eta;

	      if(!noio_flag)
		{
		  e->Dm_t[t][j][s][i] = (pow(p->M2[j][s][i], p->eta-1.0) *
					 pow(e->pm2_t[t][j][s][i],p->eta) *
					 e->m2_t[t][j][s][i]);
	      
		  e->Dbm_t[t][j][s][i] = (pow(p->eta/(p->eta-1.0),1.0-p->eta) *
					  pow(1.0+p->ntb_m_ts[t][j][s][i],1.0-p->eta) *
					  e->Dm_t[t][j][s][i] *
					  pow(e->MC_t[t][i][s],1.0-p->eta));
	      
		  e->Dhm_t[t][j][s][i] = (pow(p->eta/(p->eta-1.0),-p->eta) *
					  pow(1.0+p->ntb_m_ts[t][j][s][i],1.0-p->eta) *
					  e->Dm_t[t][j][s][i] *
					  pow(e->MC_t[t][i][s],1.0-p->eta));

		  e->Dtm_t[t][j][s][i] = e->Dbm_t[t][j][s][i] / p->eta;
		}
	      else
		{
		  e->Dm_t[t][j][s][i] = 0.0;
		  e->Dbm_t[t][j][s][i] = 0.0;
		  e->Dhm_t[t][j][s][i] = 0.0;
		  e->Dtm_t[t][j][s][i] = 0.0;
		}
	    }

	}
    }

  return 0;
  
}

uint set_vars2(eqm * e, const params * p, uint t, uint bgp)
{
  uint i,ii,s,j;

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  for(ii=0; ii<(NC-1); ii++)
	    {
	      j = p->Ji[i][ii];

	      if(t==NT)
		{
		  if(ev_steady_state(p->sig_z[i][s],
				     p->eta,
				     p->kappa0[i][s][ii],
				     p->kappa1[i][s][ii],
				     e->w_t[t][i],
				     e->Q_t[t][i],
				     e->Dtf_t[t][j][s][i] + e->Dtm_t[t][j][s][i],
				     p->n0[i][s][ii],
				     &(e->ev_t[t][i][s][ii])))
		    {
		      fprintf(logfile,
			      KRED "Failed to solve for steady state exporter moments for country/sector/dest %d/%d/%d!\n" RESET,
			      i,s,j);
		      return 1;
		    }
		}
	      else
		{
		  calc_z_Z_Fz(p->sig_z[i][s],
			      p->eta,
			      p->kappa0[i][s][ii],
			      p->kappa1[i][s][ii],
			      e->w_t[t][i],
			      e->Dtf_t[t][j][s][i] + e->Dtm_t[t][j][s][i],
			      e->ev_t[t+1][i][s][ii].dV,
			      &(e->ev_t[t][i][s][ii]));

		  double lam = 0.0;
		  if(t==0)
		    {
		      lam = p->beta[i];
		    }
		  else
		    {
		      lam = e->Q_t[t-1][i];
		    }

		  update_dV(p->kappa0[i][s][ii],
			    p->kappa1[i][s][ii],
			    e->w_t[t][i],
			    lam,
			    e->Dtf_t[t][j][s][i] + e->Dtm_t[t][j][s][i],
			    e->ev_t[t+1][i][s][ii].dV,
			    &(e->ev_t[t][i][s][ii]));
		}	      
	    }
	}
    }

  return 0;
}

uint set_vars3(eqm * e, const params * p, uint t, uint bgp)
{
  uint i, j, s, r, ii;

  
  for(i=0; i<NC; i++)
    {
      e->rgdp_t[t][i] = 0.0;
      e->ngdp_t[t][i] = 0.0;
      
      for(s=0; s<NS; s++)
	{
	  if(no_k_flag)
	    {
	      e->l_t[t][i][s] = (p->lam_va[i][s]/e->MC_t[t][i][s]) *
		(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;
	      
	      for(r=0; r<NS; r++)
		{
		  e->md_t[t][i][s][r] = (p->lam[i][s][r]/e->MC_t[t][i][s]) *
		    (e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;
		}
	    }
	  else if(cobb_douglas_flag==0)
	    {
	      e->kd_t[t][i][s] = (p->lam_va[i][s]/e->MC_t[t][i][s]) *
		pow(e->rks_t[t][i][s]*(1.0-p->alpha[i][s])/e->ws_t[t][i][s]/p->alpha[i][s],p->alpha[i][s]-1.0) *
		(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;

	      e->l_t[t][i][s] = (p->lam_va[i][s]/e->MC_t[t][i][s]) *
		pow(e->rks_t[t][i][s]*(1.0-p->alpha[i][s])/e->ws_t[t][i][s]/p->alpha[i][s],p->alpha[i][s]) *
		(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;

	      if(!noio_flag)
		{
		  for(r=0; r<NS; r++)
		    {
		      e->md_t[t][i][s][r] = (p->lam[i][s][r]/e->MC_t[t][i][s]) *
			(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;
		    }
		}
	    }
	  else
	    {
	      e->kd_t[t][i][s] = (p->lam_va[i][s]*p->alpha[i][s]/e->rks_t[t][i][s]) *
		(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;

	      e->l_t[t][i][s] = (p->lam_va[i][s]*(1.0-p->alpha[i][s])/e->ws_t[t][i][s]) *
		(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;

	      if(!noio_flag)
		{
		  for(r=0; r<NS; r++)
		    {
		      e->md_t[t][i][s][r] = (p->lam[i][s][r]/e->pm_t[t][i][r]) *
			(e->Dhf_t[t][i][s][i] + e->Dhm_t[t][i][s][i]) * e->ev_t[t][i][s][0].Zi;
		    }
		}
	    }

	  e->q2s_t[t][i][s][i] = p->H2[i][s][i] * pow(e->MC_t[t][i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
	    e->Df_t[t][i][s][i] * pow(e->ev_t[t][i][s][0].Zi,p->eta/(p->eta-1.0));   

	  if(!noio_flag)
	    {
	      e->m2s_t[t][i][s][i] = p->M2[i][s][i] * pow(e->MC_t[t][i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		e->Dm_t[t][i][s][i] * pow(e->ev_t[t][i][s][0].Zi,p->eta/(p->eta-1.0));
	    }

	  e->y_t[t][i][s] = e->q2s_t[t][i][s][i] + e->m2s_t[t][i][s][i];
	  e->va_t[t][i][s] = e->q2s_t[t][i][s][i] + e->m2s_t[t][i][s][i];
	  e->rgdp_t[t][i] += e->q2s_t[t][i][s][i] + e->m2s_t[t][i][s][i];
	  e->ngdp_t[t][i] += e->pq2_t[t][i][s][i]*e->q2s_t[t][i][s][i] + e->pm2_t[t][i][s][i]*e->m2s_t[t][i][s][i];
	  
	  e->lf_t[t][i][s] = 0.0;

	  e->exrate_t[t][i][s][i] = 1.0;
	  e->exitrate_t[t][i][s][i] = 0.0;


	  for(ii=0; ii<(NC-1); ii++)
	    {
	      j=p->Ji[i][ii];

	      double n0;

	      if(bgp)
		{
		  n0 = e->ev_t[t][i][s][ii].n;
		}
	      else if(t==0)
		{
		  n0 = p->n0[i][s][ii];
		  update_n(n0,&(e->ev_t[t][i][s][ii]));
		}
	      else
		{
		  n0 = e->ev_t[t-1][i][s][ii].n;
		  update_n(n0,&(e->ev_t[t][i][s][ii]));
		}

	      e->exrate_t[t][i][s][j] = e->ev_t[t][i][s][ii].n;
	      e->exitrate_t[t][i][s][j] = e->ev_t[t][i][s][ii].Fzm;

	      if(no_k_flag)
		{
		  e->l_t[t][i][s] += (p->lam_va[i][s]/e->MC_t[t][i][s]) *
		    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;
		  
		  for(r=0; r<NS; r++)
		    {
		      e->md_t[t][i][s][r] += (p->lam[i][s][r]/e->MC_t[t][i][s]) *
			(e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;
		    }
		}
	      else if(cobb_douglas_flag==0)
		{
		  e->kd_t[t][i][s] += (p->lam_va[i][s]/e->MC_t[t][i][s]) *
		    pow(e->rks_t[t][i][s]*(1.0-p->alpha[i][s])/e->ws_t[t][i][s]/p->alpha[i][s],p->alpha[i][s]-1.0) *
		    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;

		  e->l_t[t][i][s] += (p->lam_va[i][s]/e->MC_t[t][i][s]) *
		    pow(e->rks_t[t][i][s]*(1.0-p->alpha[i][s])/e->ws_t[t][i][s]/p->alpha[i][s],p->alpha[i][s]) *
		    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;

		  if(!noio_flag)
		    {
		      for(r=0; r<NS; r++)
			{
			  e->md_t[t][i][s][r] += (p->lam[i][s][r]/e->MC_t[t][i][s]) *
			    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;
			}
		    }
		}
	      else
		{
		  e->kd_t[t][i][s] += (p->lam_va[i][s]*p->alpha[i][s]/e->rks_t[t][i][s]) *
		    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;

		  e->l_t[t][i][s] += (p->lam_va[i][s]*(1.0-p->alpha[i][s])/e->ws_t[t][i][s]) *
		    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;

		  if(!noio_flag)
		    {
		      for(r=0; r<NS; r++)
			{
			  e->md_t[t][i][s][r] += (p->lam[i][s][r]/e->pm_t[t][i][r]) *
			    (e->Dhf_t[t][j][s][i] + e->Dhm_t[t][j][s][i]) * e->ev_t[t][i][s][ii].Z;
			}
		    }
		}

	      e->q2s_t[t][j][s][i] = p->H2[j][s][i] * pow(e->MC_t[t][i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		e->Df_t[t][j][s][i] * pow(1.0+p->ntb_f_ts[t][j][s][i],-p->eta) *
		pow(e->ev_t[t][i][s][ii].Z,p->eta/(p->eta-1.0));

	      if(!noio_flag)
		{
		  e->m2s_t[t][j][s][i] = p->M2[j][s][i] * pow(e->MC_t[t][i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
		    e->Dm_t[t][j][s][i] * pow(1.0+p->ntb_f_ts[t][j][s][i],-p->eta) *
		    pow(e->ev_t[t][i][s][ii].Z,p->eta/(p->eta-1.0));
		}
	      
	      e->y_t[t][i][s] += e->q2s_t[t][j][s][i]*(1.0+p->ntb_f_ts[t][j][s][i]) +
		e->m2s_t[t][j][s][i]*(1.0+p->ntb_m_ts[t][j][s][i]);
	      
	      e->va_t[t][i][s] += e->q2s_t[t][j][s][i]*(1.0+p->ntb_f_ts[t][j][s][i]) +
		e->m2s_t[t][j][s][i]*(1.0+p->ntb_m_ts[t][j][s][i]);

	      e->rgdp_t[t][i] += e->q2s_t[t][j][s][i]*(1.0+p->ntb_f_ts[t][j][s][i]) +
		e->m2s_t[t][j][s][i]*(1.0+p->ntb_m_ts[t][j][s][i]);

	      e->ngdp_t[t][i] += e->pq2_t[t][j][s][i]*e->q2s_t[t][j][s][i]*(1.0+p->ntb_f_ts[t][j][s][i]) +
		e->pm2_t[t][j][s][i]*e->m2s_t[t][j][s][i]*(1.0+p->ntb_m_ts[t][j][s][i]);
	      
	      if(!nokappa)
		{
		  e->lf_t[t][i][s] += p->kappa0[i][s][ii]*(1.0-n0)*(1.0-e->ev_t[t][i][s][ii].Fzp) +
		    p->kappa1[i][s][ii] * n0 * (1.0-e->ev_t[t][i][s][ii].Fzm);
		}
	    }	  
	}
    }

  // final aggregation
  for(i=0; i<NC; i++)
    {
      e->ll_t[t][i] = SUM(e->l_t[t][i],NS) + SUM( e->lf_t[t][i],NS);
      
      //if(k_adj_cost && t<(NT-1))
      if(t<(NT) && (!no_k_flag) && !bgp)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t>0)
		{
		  e->ll_t[t][i] += p->etaK * (e->k_t[t][i][s]/e->k_t[t-1][i][s]-1.0) * 
		    (e->k_t[t][i][s]/e->k_t[t-1][i][s]-1.0) * 
		    e->k_t[t-1][i][s];
		}
	      else
		{
		  e->ll_t[t][i] += p->etaK * (e->k_t[t][i][s]/p->k0[i][s]-1.0) * 
		    (e->k_t[t][i][s]/p->k0[i][s]-1.0) * 
		    p->k0[i][s];
		}
	    }
	}
      //if(l_adj_cost && t<(NT-1))
      if(t<(NT) && !bgp)
	{
	  for(s=0; s<NS; s++)
	    {
	      if(t>0)
		{
		  e->ll_t[t][i] += p->etaL * (e->l_t[t][i][s]/e->l_t[t-1][i][s]-1.0) * 
		    (e->l_t[t][i][s]/e->l_t[t-1][i][s]-1.0) * 
		    e->l_t[t-1][i][s];
		}
	      else
		{
		  e->ll_t[t][i] += p->etaL * (e->l_t[t][i][s]/p->l0[i][s]-1.0) * 
		    (e->l_t[t][i][s]/p->l0[i][s]-1.0) * 
		    p->l0[i][s];
		}
	    }
	}

      //if(!noio_flag && m_adj_cost && t<(NT-1))
      if(!noio_flag && t<(NT) && !bgp)
	{
	  for(s=0; s<NS; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  if(j!=i)
		    {
		      if(t>0)
			{
			  e->ll_t[t][i] +=
			    p->etaM * 
			    (e->m2_t[t][i][s][j]/e->m2_t[t-1][i][s][j]-1.0) * 
			    (e->m2_t[t][i][s][j]/e->m2_t[t-1][i][s][j]-1.0) * 
			    e->m2_t[t-1][i][s][j];
			}
		      else
			{
			  e->ll_t[t][i] +=
			    p->etaM * 
			    (e->m2_t[t][i][s][j]/p->m02[i][s][j]-1.0) * 
			    (e->m2_t[t][i][s][j]/p->m02[i][s][j]-1.0) * 
			    p->m02[i][s][j];
			}
		    }
		}
	    }
	}

      //if(f_adj_cost && t<(NT-1))
      if(t<(NT) && !bgp)
	{
	  for(s=0; s<NS; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  if(j!=i)
		    {
		      if(t>0)
			{
			  e->ll_t[t][i] +=
			    p->etaF * 
			    (e->q2_t[t][i][s][j]/e->q2_t[t-1][i][s][j]-1.0) * 
			    (e->q2_t[t][i][s][j]/e->q2_t[t-1][i][s][j]-1.0) * 
			    e->q2_t[t-1][i][s][j];
			}
		      else
			{
			  e->ll_t[t][i] +=
			    p->etaF * 
			    (e->q2_t[t][i][s][j]/p->q02[i][s][j]-1.0) * 
			    (e->q2_t[t][i][s][j]/p->q02[i][s][j]-1.0) * 
			    p->q02[i][s][j];
			}
		    }
		}
	    }
	}

	e->ngdp2_t[t][i] = 0.0;
	e->cc_t[t][i] = 0.0;
	for(s=0; s<NS; s++)
	  {
	    if(!no_k_flag)
	      {
		e->i_t[t][i][s] = e->pi_t[t][i] * p->eps[i][1][s] * e->ii_t[t][i]/e->p_t[t][i][s];
	      }
	    
	    //e->cc_t[t][i] += e->c_t[t][i][s];
	    e->cc_t[t][i] += p->eps[i][0][s] * pow(e->c_t[t][i][s],p->rho);
	    e->rgdp_t[t][i] -= e->m_t[t][i][s];
	    e->va_t[t][i][s] -= (e->md_t[t][i][s][0] + e->md_t[t][i][s][1] +
				 e->md_t[t][i][s][2] + e->md_t[t][i][s][3] + e->md_t[t][i][s][4]);
	    e->ngdp_t[t][i] -= e->pm_t[t][i][s]*e->m_t[t][i][s];

	    e->ngdp2_t[t][i] += e->p_t[t][i][s]*(e->c_t[t][i][s]+e->i_t[t][i][s]) + SUM(e->nxs_t[t][i][s],NC);
	  }
	e->cc_t[t][i] = pow(e->cc_t[t][i],1.0/p->rho);
	
	e->realloc_t[t][i] = 0.0;
	double wgt=0.0;
	for(s=0; s<NS; s++)
	{
	  e->realloc2_t[t][i][s] = ( (e->l_t[t][i][s]/e->ll_t[t][i])-(eee0[0].l_t[t][i][s]/eee0[0].ll_t[t][i]) ) /
	    ( (e->l_t[NT-1][i][s]/e->ll_t[NT-1][i])-(eee0[0].l_t[NT-1][i][s]/eee0[0].ll_t[NT-1][i]) );

	  double wgt2=p->l0[i][s];
	  e->realloc_t[t][i] = e->realloc_t[t][i] + wgt2*e->realloc2_t[t][i][s];
	  wgt=wgt+wgt2;
	}
      e->realloc_t[t][i] = e->realloc_t[t][i]/wgt;
    }

  return 0;
}

/*void test()
{
  double test=9;
  return;
  }*/

uint eval_bgp_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,s=0,j=0,t=NT,nx=0;
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);

  e->b_t[t][0] = bbgp[0];
  e->b_t[t][1] = bbgp[1];
  e->b_t[t][2] = bbgp[2];
  e->b_t[t][3] = bbgp[3];
  e->b_t[t][4] = bbgp[4];
  e->b_t[t][5] = bbgp[5];
  unstack_bgp_vars(e,myx);

  set_vars1(e,p,t,1);

  if(set_vars2(e,p,t,1))
    {
      fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
      return 1;
    }
  set_vars3(e,p,t,1);
  
  nx=0;

  myf[nx] = price_norm(e,t);
  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
    {
      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nPrice norm" RESET);
      return 1;
    }
  
  nx=nx+1;

  for(i=0; i<(NC-1); i++)
    {
      myf[nx] = bop(p,e,t,i);
      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	{
	  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nBop %d" RESET,i);
	  return 1;
	}

      nx = nx+1;
    }

  for(i=0; i<NC; i++)
    {
      myf[nx] = muc_mul(p,e,t,i);
      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	{
	  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMuc-Mul %d" RESET,i);
	  return 1;
	}

      nx=nx+1;

      for(s=0; s<NS; s++)
	{
	  if(s != (NS-1))
	    {
	      myf[nx] = mucs_mucr(p,e,t,i,s,NS-1);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMucs-Mucr %d %d\n" RESET,i,s);
		  return 1;
		}
	      nx=nx+1;
	    }

	  myf[nx] = mkt_clear_q(p,e,t,i,s);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear q %d %d\n" RESET,i,s);
	      printf("\n%0.4f %0.4f %0.4f\n",e->q_t[t][i][s],e->c_t[t][i][s],e->i_t[t][i][s]);
	      //printf("%0.4f %0.4f %0.4f\n",p->q0[0][s],p->c0[0][s],p->i0[0][s]);
	      printf("%0.4f %0.4f %0.4f\n",p->q0[i][s],p->c0[i][s],p->i0[i][s]);
	      printf("%0.4f %0.4f %0.4f %0.4f\n",e->pi_t[t][i],p->eps[i][1][s],e->ii_t[t][i],e->p_t[t][i][s]);
	      printf("%0.4f %0.4f %0.4f %0.4f %0.4f\n",e->k_t[t][i][0],e->k_t[t][i][1],e->k_t[t][i][2],e->k_t[t][i][3],e->k_t[t][i][4]);
	      printf("%0.4f %0.4f %0.4f %0.4f %0.4f\n",p->k0[i][0],p->k0[i][1],p->k0[i][2],p->k0[i][3],p->k0[i][4]);
	      printf("%0.4f %0.4f %0.4f %0.4f %0.4f\n",p->va0[i][0],p->va0[i][1],p->va0[i][2],p->va0[i][3],p->va0[i][4]);
	      printf("%0.4f\n",p->tauk[i]);

	      i=0;
	      printf("\n%0.4f %0.4f %0.4f\n",e->q_t[t][i][s],e->c_t[t][i][s],e->i_t[t][i][s]);
	      printf("%0.4f %0.4f %0.4f\n",p->q0[i][s],p->c0[i][s],p->i0[i][s]);
	      printf("%0.4f %0.4f %0.4f %0.4f\n",e->pi_t[t][i],p->eps[i][1][s],e->ii_t[t][i],e->p_t[t][i][s]);
	      printf("%0.4f %0.4f %0.4f %0.4f %0.4f\n",e->k_t[t][i][0],e->k_t[t][i][1],e->k_t[t][i][2],e->k_t[t][i][3],e->k_t[t][i][4]);
	      printf("%0.4f %0.4f %0.4f %0.4f %0.4f\n",p->k0[i][0],p->k0[i][1],p->k0[i][2],p->k0[i][3],p->k0[i][4]);
	      printf("%0.4f %0.4f %0.4f %0.4f %0.4f\n",p->va0[i][0],p->va0[i][1],p->va0[i][2],p->va0[i][3],p->va0[i][4]);
	      printf("%0.4f\n",p->tauk[i]);

	      return 1;
	    }
	  
	  nx=nx+1;

	  if(!noio_flag)
	    {
	      myf[nx] = mkt_clear_m(p,e,t,i,s);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear m %d %d\n" RESET,i,s);
		  return 1;
		}

	      nx=nx+1;
	    }

	  if(!no_k_flag)
	    {
	      myf[nx] = mkt_clear_k(p,e,t,i,s);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear k %d %d\n" RESET,i,s);
		  return 1;
		}
	      
	      nx=nx+1;
	    }

	  for(j=0; j<NC; j++)
	    {
	      myf[nx] = mkt_clear_q2(p,e,t,i,s,j);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear q2 %d %d %d\n" RESET,i,s,j);
		  return 1;
		}

	      nx=nx+1;

	      myf[nx] = foc_q2(p,e,t,i,s,j);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nFoc q2 %d %d %d\n" RESET,i,s,j);
		  return 1;
		}

	      nx=nx+1;
	      
	      if(!noio_flag)
		{
		  myf[nx] = mkt_clear_m2(p,e,t,i,s,j);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear m2 %d %d %d\n" RESET,i,s,j);
		      return 1;
		    }

		  nx=nx+1;

		  myf[nx]  = foc_m2(p,e,t,i,s,j);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nFoc m2 %d %d %d\n" RESET,i,s,j);
		      return 1;
		    }

		  nx=nx+1;
		}	      
	    }
	}
    }

  if(nx != solver_n)
    {
      fprintf(logfile,KRED "Wrong number of bgp eqns! nx = %d, nbgp = %d\n" RESET,nx,solver_n);
      return 1;
    }

  for(i=0; i<nx; i++)
    {
      if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	{
	  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected in position %d!\n" RESET,i);
	  return 1;
	}
    }
  
  return 0;
}

uint solve_bgp(double bb[NC])
{
  bbgp[0] = bb[0];
  bbgp[1] = bb[1];
  bbgp[2] = bb[2];
  bbgp[3] = bb[3];
  bbgp[4] = bb[4];
  bbgp[5] = bb[5];

  if(no_k_flag)
    {
      solver_n = nbgp_no_k;
    }
  else if(noio_flag)
    {
      solver_n = nbgp_noio;
    }
  else
    {
      solver_n = nbgp;
    }

  uint fixl_ = fixl;
  fixl=1;
  
  alloc_solver_mem();
  set_initial_bgp_guess();

  uint status=0;
  if(eval_bgp_once_flag)
    {
      status = bgp_func_f(solver_x,NULL,f0[0]);
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
      status=1;
    }
  else
    {
      gsl_multiroot_function_fdf f = {&bgp_func_f,&bgp_func_df,&bgp_func_fdf,solver_n,NULL};
      par=1;
      status = find_root_deriv_mkl(&f);
      if(status)
	{
	  fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
	  write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
	}
    }

  fixl=fixl_;

  free_solver_mem();
  return status;
}

uint eval_eqm_conds(const double * myx, double * myf, uint tn)
{
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);
  uint i=0,s=0,j=0,nx=0;
  int t = 0;
  int t0 = 0;

  if(scenario>=1)
    {
      e = &(eee1[tn]);
      p = &(ppp1[tn]);
      t0=TNAFTA;
    }

  unstack_eqm_vars(e,myx);

  e->b_t[0][0] = p->b0[0];
  e->b_t[0][1] = p->b0[1];
  e->b_t[0][2] = p->b0[2];
  e->b_t[0][3] = p->b0[3];
  e->b_t[0][4] = p->b0[4];
  e->b_t[0][5] = p->b0[5];

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  e->k_t[0][i][s] = p->k0[i][s];
	}
    }

  for(t=t0; t<(NT+1); t++)
    {
      if(set_vars1(e,p,t,0))
	{
	  fprintf(logfile,KRED "Error calling set_vars1!\n" RESET);
	  return 1;
	}
    }

  for(t=NT; t>=t0; t--)
    {
      if(set_vars2(e,p,t,0))
	{
	  fprintf(logfile,KRED "Error calling set_vars2!\n" RESET);
	  return 1;
	}
    }
  for(t=t0; t<(NT+1); t++)
    {
      set_vars3(e,p,t,0);
    }

  nx=0;
  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	{
	  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\n Price norm %d\n" RESET,t);
	  return 1;
	}

      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nBop %d %d\n" RESET,t,i);
	      return 1;
	    }
		      
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  myf[nx] = muc_mul(p,e,t,i);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMuc-Mul %d %d\n" RESET,t,i);
	      return 1;
	    }

	  
	  /*
	  printf("MUC = %0.6f\n,",
		 muc
		 (
		  e->c_t[t][i],
		  e->ll_t[t][i],
		  p->lbar[i],
		  p->eps[i][0],
		  p->rho,
		  p->phi[i],
		  p->psi,
		  2));
	  printf("MUL = %0.6f\n",
		 mul
		 (
		  e->c_t[t][i],
		  e->ll_t[t][i],
		  p->lbar[i],
		  p->eps[i][0],
		  p->rho,
		  p->phi[i],
		  p->psi));
	  printf("t = %d\n",t);
	  printf("l = %0.6f\n",e->ll_t[t][i]);
	  printf("c = %0.6f\n",e->c_t[t][i][2]);
	  printf("p = %0.6f\n",e->p_t[t][i][2]);
	  printf("w = %0.6f\n",e->w_t[t][i]);
		 
	  exit(1);


	  */
	  
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nEuler %d %d\n" RESET,t,i);
		  return 1;
		}

	      nx = nx+1;	       
	    }
	    

	  for(s=0; s<NS; s++)
	    {
	      if(s != (NS-1))
		{
		  myf[nx] = mucs_mucr(p,e,t,i,s,NS-1);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMucs-Mucr %d %d %d\n" RESET,t,i,s);
		  return 1;
		}

		  nx=nx+1;
		}

	      myf[nx] = mkt_clear_q(p,e,t,i,s);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear q %d %d %d\n" RESET,t,i,s);
		  return 1;
		}

	      nx=nx+1;

	      if(!noio_flag)
		{
		  myf[nx] = mkt_clear_m(p,e,t,i,s);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear m %d %d %d\n" RESET,t,i,s);
		      return 1;
		    }

		  nx=nx+1;
		}

	      if(t>t0 && (!no_k_flag))
		{
		  myf[nx] = noarb(p,e,t,i,s);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nNoarb %d %d %d\n" RESET,t,i,s);
		      return 1;
		    }

		  nx=nx+1;
		}
	      
	      myf[nx] = noarb_l(p,e,t,i,s);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nNoarb_l %d %d %d\n" RESET,t,i,s);
		  return 1;
		}

	      nx=nx+1;

	      if(!no_k_flag)
		{
		  myf[nx] = mkt_clear_k(p,e,t,i,s);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear k %d %d %d\n" RESET,t,i,s);
		      return 1;
		    }

		  nx=nx+1;
		}

	      for(j=0; j<NC; j++)
		{
		  myf[nx] = mkt_clear_q2(p,e,t,i,s,j);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear q2 %d %d %d %d\n" RESET,t,i,s,j);
		      printf("%0.4f %0.4f\n",e->q2_t[t][i][s][j],e->q2s_t[t][i][s][j]);

		      /*
		      e->q2s_t[t][j][s][i] = p->H2[j][s][i] * pow(e->MC_t[t][i][s] * (p->eta/(p->eta-1.0)),-p->eta) *
			e->Df_t[t][j][s][i] * pow(1.0+p->ntb_f_ts[t][j][s][i],-p->eta) *
			pow(e->ev_t[t][i][s][ii].Z,p->eta/(p->eta-1.0));
			j=p->Ji[i][ii];
			p->Ji[i][0]=1;
			p->Ji[i][1]=2;
			p->Ji[i][2]=3;
			p->Ji[i][3]=4;
			p->Ji[i][4]=5;
		      */

		      printf("%0.4f %0.4f %0.4f %0.4f\n",p->H2[i][s][j],e->MC_t[t][j][s],e->Df_t[t][i][s][j],e->ev_t[t][j][s][3].Zi);
		      return 1;
		    }

		  nx=nx+1;

		  myf[nx] = foc_q2(p,e,t,i,s,j);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nFoc q2 %d %d %d %d\n" RESET,t,i,s,j);
		      return 1;
		    }
		  
		  nx=nx+1;
	      
		  if(!noio_flag)
		    {
		      myf[nx] = mkt_clear_m2(p,e,t,i,s,j);
		      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
			{
			  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nMkt clear m2 %ds %d %d %d\n" RESET,t,i,s,j);
			  return 1;
			}

		      nx=nx+1;

		      myf[nx]  = foc_m2(p,e,t,i,s,j);
		      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
			{
			  fprintf(logfile,KRED "Error evaluating bgp eqns! NaN/Inf detected!\nFoc m2 %ds %d %d %d\n" RESET,t,i,s,j);
			  return 1;
			}

		      nx=nx+1;
		    }	      
		}
	    }
	}
    }

  if(nx != neqm)
    {
      fprintf(logfile,KRED "Error evaluating eqm eqns! nx = %d, neqm = %d\n" RESET,nx,neqm);
      return 1;
    }

  for(i=0; i<nx; i++)
    {
      if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	{
	  fprintf(logfile,KRED "Error evaluating equilibrium conditions! NaN/Inf detected!\n" RESET);
	  return 1;
	}
    }
  

  //double max = fabs_max(myf,neqm);
  //uint imax = fabs_imax(myf,neqm);

  return 0;
}

uint solve_eqm()
{
  char * sname;

  if(scenario==0)
    {
      sname = "output/seed0.bin";
      if(read_seed==1)
	{
	  free_solver_mem();
	  solver_n = neqm;
	  alloc_solver_mem();

	  if(read_vec_bin(solver_x->data, neqm, sname))
	    {
	      fprintf(logfile,KRED "Error loading equilibrium guess from seed file!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
      else
	{
	  if(set_initial_eqm_guess())
	    {
	      fprintf(logfile,KRED "Error constructing equilibrium guess!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
    }
  // otherwise we should use the solution from the previous exercise as the initial guess
  else
    {      
      sname = "output/seed1.bin";
      free_solver_mem();
      solver_n = neqm;
      alloc_solver_mem();

      eqm * e = &(eee0[0]);

      if(read_seed==1)
	{
	  if(read_vec_bin(solver_x->data, neqm, sname))
	    {
	      fprintf(logfile,KRED "Error loading equilibrium guess from seed file!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}
      else
	{
	  uint it;
	  for(it=0; it<NTH; it++)
	    {
	      copy_vars( &(eee1[it]), e  );
	    }
	  if(stack_eqm_vars(solver_x->data,e))
	    {
	      fprintf(logfile,KRED "Failed to stack variables from previous exercise!\n" RESET);
	      free_solver_mem();
	      return 1;
	    }
	}

    }

  uint status = 0;
  if(eval_eqm_once_flag)
    {
      status = eqm_func_f(solver_x,NULL,f0[0]);
      write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
      if(status)
	fprintf(logfile,KRED "Error evaluating equilibrium function!\n" RESET);
    }
  else
    {
      gsl_multiroot_function_fdf f = {&eqm_func_f,&eqm_func_df,&eqm_func_fdf,neqm,NULL};

      par=1;
      status = find_root_deriv_mkl(&f);
      if(status)
	{
	  fprintf(logfile,KRED "Error solving for equilibrium!\n" RESET);
	  write_vec_txt(f0[0]->data,solver_n,"output/F.txt");
	}
      
      if(write_seed==1 && !status)
	{
	  write_vec_bin(solver_x->data, neqm, sname);
	}

    }

  free_solver_mem();

  return status;
}

void calc_welfare(eqm * e, const params * p)
{
  int t, i;
  for(i=0; i<NC; i++)
    {
      t=NT;
      e->welfare_t[t][i] = (1.0/(1.0-p->beta[i]*pow(1.0,p->phi[i]*p->psi))) * 
	(pow(p->eps[i][0][0] * pow(e->c_t[t][i][0],p->rho) +
	     p->eps[i][0][1] * pow(e->c_t[t][i][1],p->rho) + 
	     p->eps[i][0][2] * pow(e->c_t[t][i][2],p->rho) +
	     p->eps[i][0][3] * pow(e->c_t[t][i][3],p->rho) +
	     p->eps[i][0][4] * pow(e->c_t[t][i][4],p->rho),
	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar[i]-e->ll_t[t][i]),(1.0-p->phi[i])*p->psi));	
      
      for(t=(NT-1); t>=0; t--)
	{
	  e->welfare_t[t][i] = p->beta[i] * e->welfare_t[t+1][i] + 
	    (pow(p->eps[i][0][0] * pow(e->c_t[t][i][0],p->rho) +
		 p->eps[i][0][1] * pow(e->c_t[t][i][1],p->rho) + 
		 p->eps[i][0][2] * pow(e->c_t[t][i][2],p->rho) +
		 p->eps[i][0][3] * pow(e->c_t[t][i][3],p->rho) +
		 p->eps[i][0][4] * pow(e->c_t[t][i][4],p->rho),
		 p->phi[i]*p->psi/p->rho) * 
	    pow((p->lbar[i]-e->ll_t[t][i]),(1.0-p->phi[i])*p->psi));
	  
	  e->welfare_t[t+1][i] = pow(e->welfare_t[t+1][i],1.0/p->psi);
	}
      e->welfare_t[0][i] = pow(e->welfare_t[0][i],1.0/p->psi);
    }

  if(scenario == 0)
    {
      for(i=0; i<NC; i++)
	{
	  t=NT;
	  e->welfare_cost_t[t][i] = (1.0/(1.0-e->pb_t[t])) * e->cpi_t[t][i] * e->cc_t[t][i];
	  
	  for(t=(NT-1); t>=0; t--)
	    {
	      e->welfare_cost_t[t][i] = e->cpi_t[t][i]*e->cc_t[t][i] + e->pb_t[t] * e->welfare_cost_t[t+1][i];
	    }
	  for(t=0; t<NT; t++)
	    {
	      e->welfare_cost_t[t][i] = e->pb_t[t] * e->welfare_cost_t[t+1][i];
	    }

 	}
    }
}

///////////////////////////////

int bgp_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
   //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_bgp_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int bgp_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&bgp_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int bgp_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(bgp_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&bgp_func_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
    }
}

int eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f)
{
  //fcnt = fcnt + 1;
  uint tn;
  if(data==NULL)
    {
      tn = 0;
    }
  else
    {
      tn = *((uint *)data);
    }
  if(eval_eqm_conds(x->data,f->data,tn))
    {
      return GSL_EBADFUNC;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J)
{
  if(jacobian(&eqm_func_f, x, J, 1))
    {
      return GSL_EFAILED;
    }
  else
    {
      return GSL_SUCCESS;
    }
}

int eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J)
{
  if(eqm_func_f(x,NULL,f))
    {
      return GSL_EFAILED;
    }
  else
    {
      gsl_vector_memcpy(f0[0],f);
      if(jacobian(&eqm_func_f, x, J, 0))
	{
	  return GSL_EFAILED;
	}
      else
	{
	  return GSL_SUCCESS;
	}
     }
}

#endif
