#ifndef __EQM_C__
#define __EQM_C__

#include "eqm.h"

// wage (NC), sectoral capital (NC*NS), sectoral labor (NC*NS), gross output prices (NC*NS),
// sectoral consumption (NC*(NS-1))
const uint nbgp = NC+3*NC*NS+NC*(NS-1);

void set_neqm()
{
  // key eqm vars used in solver: wage (NC), sectoral investment (NC*NS, but not in last period),
  // sectoral labor (NC*NS), gross output prices (NC*NS), sectoral consumption (NC*(NS-1)), 
  // bonds (NC-1, but not in first period), bond price (1, but not in last period)

  // deterministic no-Brexit counterfactual
  if(scenario==0)
    {
      neqm = (NT+1)*NC + 2*(NT+1)*NC*NS + (NT+1)*NC*(NS-1) + NT*NC*NS + NT*(NC-1) + NT;
      if(f_adj_cost)
	{
	  neqm = neqm + NT*NC*(NS-1)*NC + NT*NC*(NS-1);
	}
      if(m_adj_cost)
	{
	  neqm = neqm + NT*NC*(NS-1)*NC + NT*NC*(NS-1);
	}
    }
  else if(scenario==1)
    {
      uint nn = NT+1-TSHOCK;
      neqm = (nn)*NC + 2*(nn)*NC*NS + (nn)*NC*(NS-1) + (nn-1)*NC*NS + (nn-1)*(NC-1) + (nn-1);
      if(f_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*(NS-1)*NC + (nn-1)*NC*(NS-1);
	}
      if(m_adj_cost)
	{
	  neqm = neqm + (nn-1)*NC*(NS-1)*NC + (nn-1)*NC*(NS-1);
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
  SET_ALL_V(e->md_t,(NT+1)*NC*NS*(NS-1),0.0);
  SET_ALL_V(e->k_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->l_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->is_t,(NT+1)*NC*NS,0.0);
  SET_ALL_V(e->lp_t,(NT+1)*NC*NS,0.0);

  SET_ALL_V(e->exs_t,(NT+1)*NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->ims_t,(NT+1)*NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->nxs_t,(NT+1)*NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->rexs_t,(NT+1)*NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->rims_t,(NT+1)*NC*(NS-1)*NC,0.0);

  SET_ALL_V(e->pm_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->m_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->m2_t,(NT+1)*NC*(NS-1)*NC,0.0);

  SET_ALL_V(e->p_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->q_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->q2_t,(NT+1)*NC*(NS-1)*NC,0.0);

  SET_ALL_V(e->c_t,(NT+1)*NC*(NS-1),0.0);
  SET_ALL_V(e->i_t,(NT+1)*NC*NS,0.0);

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
  memcpy((double *)(e1->md_t),
	 (const double *)(e0->md_t),
	 sizeof(double)*(NT+1)*NC*NS*(NS-1));
  memcpy((double *)(e1->lp_t),
	 (const double *)(e0->lp_t),
	 sizeof(double)*(NT+1)*NC*NS);

  memcpy((double *)(e1->exs_t),
	 (const double *)(e0->exs_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);
  memcpy((double *)(e1->ims_t),
	 (const double *)(e0->ims_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);
  memcpy((double *)(e1->nxs_t),
	 (const double *)(e0->nxs_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);
  memcpy((double *)(e1->rexs_t),
	 (const double *)(e0->rexs_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);
  memcpy((double *)(e1->rims_t),
	 (const double *)(e0->rims_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);

  memcpy((double *)(e1->pm_t),
	 (const double *)(e0->pm_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->m_t),
	 (const double *)(e0->m_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->m2_t),
	 (const double *)(e0->m2_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);

  memcpy((double *)(e1->p_t),
	 (const double *)(e0->p_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->q_t),
	 (const double *)(e0->q_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
  memcpy((double *)(e1->q2_t),
	 (const double *)(e0->q2_t),
	 sizeof(double)*(NT+1)*NC*(NS-1)*NC);

  memcpy((double *)(e1->c_t),
	 (const double *)(e0->c_t),
	 sizeof(double)*(NT+1)*NC*(NS-1));
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
}

uint stack_bgp_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t = NT;
  
  COPY_SUBVECTOR_LOG(myx+nx,e->w_t[t],NC);
  nx=nx+NC;

  COPY_SUBVECTOR_LOG(myx+nx,e->k_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->l_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->py_t[t],NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,e->c_t[t],NC*(NS-1));
  nx=nx+NC*(NS-1);

  if(nx != nbgp)
    {
      fprintf(logfile,"Error stacking bgp vars! nx = %d, nbgp = %d\n",nx,nbgp);
      return 1;
    }

    return 0;
}

uint unstack_bgp_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t = NT;

  COPY_SUBVECTOR_EXP(e->w_t[t], myx+nx, NC);
  nx=nx+NC;

  COPY_SUBVECTOR_EXP(e->k_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->l_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->py_t[t],myx+nx,NC*NS);
  nx=nx+NC*NS;

  COPY_SUBVECTOR_EXP(e->c_t[t],myx+nx,NC*(NS-1));
  nx=nx+NC*(NS-1);

  if(nx != nbgp)
    {
      fprintf(logfile,"Error stacking bgp vars! nx = %d, nbgp = %d\n",nx,nbgp);
      return 1;
    }

    return 0;
}

uint stack_eqm_vars(double * myx, const eqm * e)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;
  
  if(scenario==1)
    {
      nn = NT+1-TSHOCK;
      t0=TSHOCK;
    }

  COPY_SUBVECTOR_LOG(myx+nx,&(e->w_t[t0]),(nn)*NC);
  nx = nx + (nn)*NC;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->is_t[t0]),(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->l_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->py_t[t0]),(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_LOG(myx+nx,&(e->c_t[t0]),(nn)*NC*(NS-1));
  nx = nx + (nn)*NC*(NS-1);

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

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->q2_t[t0]),(nn-1)*NC*(NS-1)*NC);
      nx = nx + (nn-1)*NC*(NS-1)*NC;

      COPY_SUBVECTOR_LOG(myx+nx,&(e->p_t[t0]),(nn-1)*NC*(NS-1));
      nx = nx + (nn-1)*NC*(NS-1);
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_LOG(myx+nx,&(e->m2_t[t0]),(nn-1)*NC*(NS-1)*NC);
      nx = nx + (nn-1)*NC*(NS-1)*NC;
      
      COPY_SUBVECTOR_LOG(myx+nx,&(e->pm_t[t0]),(nn-1)*NC*(NS-1));
      nx = nx + (nn-1)*NC*(NS-1);
    }

  if(nx != neqm)
    {
      fprintf(logfile,"Error stacking eqm vars! nx = %d, neqm = %d\n",nx,neqm);
      return 1;
    }

    return 0;
}

uint unstack_eqm_vars(eqm * e, const double * myx)
{
  uint nx = 0;
  uint t0 = 0;
  uint nn = NT+1;

  if(scenario==1)
    {
      nn = NT+1-TSHOCK;
      t0=TSHOCK;
    }

  COPY_SUBVECTOR_EXP(&(e->w_t[t0]),myx+nx,(nn)*NC);
  nx = nx + (nn)*NC;
  
  COPY_SUBVECTOR_EXP(&(e->is_t[t0]),myx+nx,(nn-1)*NC*NS);
  nx = nx + (nn-1)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->l_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->py_t[t0]),myx+nx,(nn)*NC*NS);
  nx = nx + (nn)*NC*NS;

  COPY_SUBVECTOR_EXP(&(e->c_t[t0]),myx+nx,(nn)*NC*(NS-1));
  nx = nx + (nn)*NC*(NS-1);

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

  if(f_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->q2_t[t0]),myx+nx,(nn-1)*NC*(NS-1)*NC);
      nx = nx + (nn-1)*NC*(NS-1)*NC;

      COPY_SUBVECTOR_EXP(&(e->p_t[t0]),myx+nx,(nn-1)*NC*(NS-1));
      nx = nx + (nn-1)*NC*(NS-1);
    }

  if(m_adj_cost)
    {
      COPY_SUBVECTOR_EXP(&(e->m2_t[t0]),myx+nx,(nn-1)*NC*(NS-1)*NC);
      nx = nx + (nn-1)*NC*(NS-1)*NC;
      
      COPY_SUBVECTOR_EXP(&(e->pm_t[t0]),myx+nx,(nn-1)*NC*(NS-1));
      nx = nx + (nn-1)*NC*(NS-1);
    }

  if(nx != neqm)
    {
      fprintf(logfile,"Error unstacking eqm vars! nx = %d, neqm = %d\n" ,nx,neqm);
      return 1;
    }

  return 0;
}

uint set_initial_bgp_guess()
{
  uint i, t, s;
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
	  e->py_t[t][i][s] = 1.0;
	}
    }

  if(stack_bgp_vars(solver_x->data,e))
    {
      fprintf(logfile,"Failed to create guess for balanced growth path!\n");
      return 1;
    }
  else
    {
      return 0;
    }
}

uint set_initial_eqm_guess()
{
  uint i,s,j,t;
  double bb[NC];

  eqm * e = &(eee0[0]);
  params * p = &(ppp0[0]);

  bb[0] = p->b0[0];
  bb[1] = p->b0[1];
  bb[2] = p->b0[2];

  free_solver_mem();
  solver_n = nbgp;  
  alloc_solver_mem();
  if(solve_bgp(bb))
    {
      fprintf(logfile, "Error solving for steady state!\n");
      return 1;
    }  
  free_solver_mem();
  solver_n = neqm;
  alloc_solver_mem();

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

  LOGSPACE_2D(p->k0,e->k_t[NT],NT+1,NC*NS,e->k_t);
  LOGSPACE_2D(p->l0,e->l_t[NT],NT+1,NC*NS,e->l_t);  
  LOGSPACE_2D(p->c0,e->c_t[NT],NT+1,NC*(NS-1),e->c_t);
  LOGSPACE_2D(p->q0,e->q_t[NT],NT+1,NC*(NS-1),e->q_t);  
  LOGSPACE_2D(p->m0,e->m_t[NT],NT+1,NC*(NS-1),e->m_t);
  LINSPACE_2D(tmpp,e->w_t[NT],NT+1,NC,e->w_t);
  LINSPACE_2D(tmpp2,e->py_t[NT],NT+1,NC*NS,e->py_t);
  LINSPACE_2D(tmpp2,e->p_t[NT],NT+1,NC*(NS-1),e->p_t);
  LINSPACE_2D(tmpp2,e->pm_t[NT],NT+1,NC*(NS-1),e->pm_t);

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

  for(t=0; t<(NT+1); t++)
    {
      for(i=0; i<NC; i++)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      for(j=0; j<NC; j++)
		{
		  e->q2_t[t][i][s][j] = (p->q02[i][s][j]/p->q0[i][s]) * e->q_t[t][i][s];
		  e->m2_t[t][i][s][j] = (p->m02[i][s][j]/p->m0[i][s]) * e->m_t[t][i][s];
		}
	    }
      
	}
    }

  if(stack_eqm_vars(solver_x->data,e))
    {
      fprintf(logfile,"Failed to create guess for transition equilibrium!\n");
      return 1;
    }
  else
    {      
      return 0;
    }
}

uint write_eqm_vars(const eqm * e, const params * p, char * fname, uint i)
{
  char fname2[128];

  if(scenario==1)
    {
      sprintf(fname2,"output/%s_t%d_s%d_c%d_r%d_d%d_a%d.csv",
	      fname,(int)(tariff*100),target_sector_flag,target_country_flag,retaliation_flag,duration_flag,adjustment_flag);
    }
  else
    {
      sprintf(fname2,"output/%s_a%d.csv",fname,adjustment_flag);
    }


  FILE * file = fopen(fname2,"w");

  if(file)
    {
      uint s,j,t;
      fprintf(file,"period,rgdp,ngdp,ii,ll,c,w,ae,te");
      for(s=0;s<NS;s++)
	{
	  fprintf(file,",y%d,va%d,lp%d,i%d,c%d,k%d,l%d,py%d,ae%d,te%d",s,s,s,s,s,s,s,s,s,s);
	}
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      fprintf(file,",rer%d,ex%d,im%d,nx%d,rex%d,rim%d",
		      j,j,j,j,j,j);
	      for(s=0;s<NS-1;s++)
		{
		  fprintf(file,
			  ",taus%d-%d,exs%d-%d,ims%d-%d,nxs%d-%d,rexs%d-%d,rims%d-%d,exsm%d-%d,exsf%d-%d,imsm%d-%d,imsf%d-%d,aes%d-%d,tes%d-%d",
			  s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j,s,j);
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
	  fprintf(file,"%0.16f,",e->w_t[t][i]);
	  fprintf(file,"%0.16f,",e->ae_t[t][i]);
	  fprintf(file,"%0.16f",e->te_t[t][i]);
	  for(s=0; s<NS; s++)
	    {
	      fprintf(file,",%0.16f,",e->y_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->rva_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->lp_t[t][i][s]);
	      fprintf(file,"%0.16f,",e->is_t[t][i][s]);

	      if(s!=CNS)
		fprintf(file,"%0.16f,",e->c_t[t][i][s]);
	      else
		fprintf(file,"%0.16f,",0.0);

	      fprintf(file,"%0.16f,",e->k_t[t][i][s]);
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
		  for(s=0; s<NS-1; s++)
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
		      fprintf(file,"%0.16f",e->tes2_t[t][i][s][j]);
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
      fprintf(logfile,"Error opening file to write equilibrium vars!\n");
      return 1;
    }
}

uint set_vars(eqm * e, const params * p, uint t, uint bgp)
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
  SET_ALL_V(e->exs_t[t],NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->ims_t[t],NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->nxs_t[t],NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->rex_t[t],NC*NC,0.0);
  SET_ALL_V(e->rim_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexf_t[t],NC*NC,0.0);
  SET_ALL_V(e->rimf_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexm_t[t],NC*NC,0.0);
  SET_ALL_V(e->rimm_t[t],NC*NC,0.0);
  SET_ALL_V(e->rexs_t[t],NC*(NS-1)*NC,0.0);
  SET_ALL_V(e->rims_t[t],NC*(NS-1)*NC,0.0);

  // bond market clearing
  e->b_t[t][2] = -(e->b_t[t][0]+e->b_t[t][1]);

  // compute sector-level aggregates
  for(i=0; i<NC; i++)
    {
      // value added, gross output, intermediate demand
      for(s=0; s<NS; s++)
	{
	  e->va_t[t][i][s] = (p->a_ts[t][i][s]) * prod_va(e->k_t[t][i][s],e->l_t[t][i][s],p->A[i][s],p->alpha[i][s]);
	  e->y_t[t][i][s] = e->va_t[t][i][s]/p->lam_va[i][s];
	  
	  if(t<(NT-1) && l_adj_cost==1)
	    {
	      if(t>0)
		{
		  e->y_t[t][i][s] = e->y_t[t][i][s] - 
		    p->etaL * (e->l_t[t][i][s]/e->l_t[t-1][i][s]-1.0) * 
		    (e->l_t[t][i][s]/e->l_t[t-1][i][s]-1.0) * 
		    e->l_t[t-1][i][s];
		}
	      else
		{
		  e->y_t[t][i][s] = e->y_t[t][i][s] - 
		    p->etaL * (e->l_t[t][i][s]/p->l0[i][s]-1.0) * 
		    (e->l_t[t][i][s]/p->l0[i][s]-1.0) * 
		    p->l0[i][s];
		}
	    }

	  e->ngdp_t[t][i] = e->ngdp_t[t][i] + e->py_t[t][i][s]*e->y_t[t][i][s];
	  e->rgdp_t[t][i] = e->rgdp_t[t][i] + e->y_t[t][i][s];
	  e->rva_t[t][i][s] = e->y_t[t][i][s];
	  

	  for(r=0; r<NS-1; r++)
	    {
	      e->md_t[t][i][s][r] = e->y_t[t][i][s]*p->lam[i][s][r];
	      e->rva_t[t][t][s] -= e->md_t[t][i][s][r];
	    }

	  e->lp_t[t][i][s] = e->rva_t[t][i][s]/e->l_t[t][i][s];
	}

      e->ll_t[t][i] = SUM(e->l_t[t][i], NS);

      if(t<NT)
	{
	  if(t==(NT-1) || k_adj_cost==0)
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
		  e->k_t[t+1][i][s] = (1.0-p->delta) * e->k_t[t][i][s] + 
		    phiK(e->is_t[t][i][s]/e->k_t[t][i][s],p->delta,p->etaK) * e->k_t[t][i][s];
		}
	    }
	}
      else
	{
	  for(s=0; s<NS; s++)
	    {
	      e->is_t[t][i][s] = p->delta * e->k_t[t][i][s];
	    }
	}

      e->kk_t[t][i] = SUM(e->k_t[t][i], NS);
      e->ii_t[t][i] = SUM(e->is_t[t][i], NS);

      // Armington composite prices
      if(!m_adj_cost || t==NT)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      e->pm_t[t][i][s] = 0.0;
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_m_ts[t][i][s][j];
		  e->pm_t[t][i][s] = e->pm_t[t][i][s] + 
		    pow(p->mu[i][s][j],1.0/(1.0-p->zeta[i][s])) * 
		    pow(tc*e->py_t[t][j][s],p->zeta[i][s]/(p->zeta[i][s]-1.0));
		}
	      e->pm_t[t][i][s] = (1.0/p->M[i][s]) * 
		pow(e->pm_t[t][i][s],(p->zeta[i][s]-1.0)/p->zeta[i][s]);
	    }
	}

      if(!f_adj_cost || t==NT)
	{
	  for(s=0; s<NS-1; s++)
	    {
	      e->p_t[t][i][s] = 0.0;
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_f_ts[t][i][s][j];
		  e->p_t[t][i][s] = e->p_t[t][i][s] + 
		    pow(p->theta[i][s][j],1.0/(1.0-p->sig[i][s])) * 
		    pow(tc*e->py_t[t][j][s],p->sig[i][s]/(p->sig[i][s]-1.0));
		}
	      e->p_t[t][i][s] = (1.0/p->H[i][s]) * 
		pow(e->p_t[t][i][s],(p->sig[i][s]-1.0)/p->sig[i][s]);
	    }
	}
      double pt_cns = e->py_t[t][i][CNS];

      // households' stochastic discount factor for dynamic firm's problem with adjustment costs
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

      // investment price
      e->pi_t[t][i] = 1.0/p->G[i];
      for(s=0;s<NS-1; s++)
	{
	  e->pi_t[t][i] = e->pi_t[t][i] * pow(e->p_t[t][i][s]/p->eps[i][1][s],p->eps[i][1][s]);
	}
      e->pi_t[t][i] = e->pi_t[t][i] * pow(pt_cns/p->eps[i][1][s],p->eps[i][1][s]);

      // demand for final goods and intermediates other than construction
      for(s=0; s<NS-1; s++)
	{
	  e->i_t[t][i][s] = e->pi_t[t][i] * p->eps[i][1][s] * e->ii_t[t][i]/e->p_t[t][i][s];
	  e->q_t[t][i][s] = e->c_t[t][i][s] + e->i_t[t][i][s];
	  e->m_t[t][i][s] = e->md_t[t][i][UPS][s] + e->md_t[t][i][DNS][s] + e->md_t[t][i][SVC][s] + e->md_t[t][i][CNS][s];

	  if(!m_adj_cost || t==NT)
	    {
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_m_ts[t][i][s][j];
		  e->m2_t[t][i][s][j] = e->m_t[t][i][s] * 
		    pow(tc*e->py_t[t][j][s],1.0/(p->zeta[i][s]-1.0)) * 
		    pow(e->pm_t[t][i][s]*p->mu[i][s][j]*pow(p->M[i][s],p->zeta[i][s]),1.0/(1.0-p->zeta[i][s]));
		}
	    }
	  
	  if(!f_adj_cost || t==NT)
	    {
	      for(j=0; j<NC; j++)
		{
		  double tc = 1.0+p->tau_f_ts[t][i][s][j];
		  e->q2_t[t][i][s][j] = e->q_t[t][i][s] * 
		    pow(tc*e->py_t[t][j][s],1.0/(p->sig[i][s]-1.0)) * 
		    pow(e->p_t[t][i][s]*p->theta[i][s][j]*pow(p->H[i][s],p->sig[i][s]),1.0/(1.0-p->sig[i][s]));		 
		}
	    }
	}

      // construction investment demand
      e->i_t[t][i][CNS] = e->pi_t[t][i] * p->eps[i][1][s] * e->ii_t[t][i]/pt_cns;

      e->cpi_t[t][i] = 0.0;
      e->cc_t[t][i] = 0.0;
      for(s=0; s<NS-1; s++)
	{
	  e->cpi_t[t][i] = e->cpi_t[t][i] + e->p_t[t][i][s]*p->c0[i][s];
	  e->cc_t[t][i] = e->cc_t[t][i]+e->c_t[t][i][s];
	}
      e->cpi_t[t][i] = e->cpi_t[t][i]/SUM(p->c0[i],NS-1);
      
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
	      e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->pb_t[t-1] - (1.0-p->delta)*e->pi_t[t][i];
	    }
	}
      else
	{
	  e->rk_t[t][i] = e->pi_t[t-1][i]*e->cpi_t[t][0]/e->pb_t[t-1] - (1.0-p->delta)*e->pi_t[t][i];
	}
    }

  for(i=0; i<NC; i++)
    {
      for(j=0; j<NC; j++)
	{
	  e->ngdp_t[t][i] = e->ngdp_t[t][i] - 
	    e->py_t[t][j][0]*e->m2_t[t][i][0][j] - 
	    e->py_t[t][j][1]*e->m2_t[t][i][1][j] - 
	    e->py_t[t][j][2]*e->m2_t[t][i][2][j];

	  e->rgdp_t[t][i] = e->rgdp_t[t][i] - 
	    e->m2_t[t][i][0][j] - 
	    e->m2_t[t][i][1][j] - 
	    e->m2_t[t][i][2][j];
	  
	  if(j!=i)
	    {
	      e->rer_t[t][i][j] = e->cpi_t[t][j]/e->cpi_t[t][i];
	      for(s=0; s<NS-1; s++)
		{
		  double m = e->py_t[t][i][s]*e->m2_t[t][j][s][i];
		  double f = e->py_t[t][i][s]*e->q2_t[t][j][s][i];
	
		  e->exs_t[t][i][s][j] = m+f;
		  e->exm_t[t][i][j] = e->exm_t[t][i][j] + m;
		  e->exf_t[t][i][j] = e->exf_t[t][i][j] + f;
		  e->ex_t[t][i][j] = e->ex_t[t][i][j] + m+f;

		  m = e->m2_t[t][j][s][i];
		  f = e->q2_t[t][j][s][i];
		  e->rexs_t[t][i][s][j] = m+f;
		  e->rexm_t[t][i][j] = e->rexm_t[t][i][j] + m;
		  e->rexf_t[t][i][j] = e->rexf_t[t][i][j] + f;
		  e->rex_t[t][i][j] = e->rex_t[t][i][j] + m+f;

		  m = e->py_t[t][j][s]*e->m2_t[t][i][s][j];
		  f = e->py_t[t][j][s]*e->q2_t[t][i][s][j];
		  e->ims_t[t][i][s][j] = m+f;
		  e->imm_t[t][i][j] = e->imm_t[t][i][j] + m;
		  e->imf_t[t][i][j] = e->imf_t[t][i][j] + f;
		  e->im_t[t][i][j] = e->im_t[t][i][j] + m+f;

		  m = e->m2_t[t][i][s][j];
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

      e->iy_t[t][i] = e->pi_t[t][i]*e->ii_t[t][i]/e->ngdp_t[t][i];
      e->lp_agg_t[t][i] = e->rgdp_t[t][i] / e->ll_t[t][i];
    }

  for(i=0; i<NC; i++)
    {
      e->ae_t[t][i] = 0.0;
      e->te_t[t][i] = 0.0;
      double w = 0.0;
      double wt = 0.0;
      
      for(s=0; s<NS-2; s++)
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
		  double xt = log((e->py_t[t][i][s]/((1.0+p->tau_m_ts[t][i][s][j])*e->py_t[t][j][s])));
		  double a0 = (1.0/(1.0-p->zeta[i][s]))*log(p->mu[i][s][j]/p->mu[i][s][i]);
		  e->aes2_m_t[t][i][s][j] = (yt-a0)/xt;

		  yt = log(e->q2_t[t][i][s][j] / e->q2_t[t][i][s][i]);
		  xt = log((e->py_t[t][i][s]/((1.0+p->tau_f_ts[t][i][s][j])*e->py_t[t][j][s])));
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

  for(i=0; i<NC; i++)
    {
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

uint eval_bgp_conds(const double * myx, double * myf, uint tn)
{
  uint i=0,s=0,t=NT,nx=0;
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);

  e->b_t[t][0] = bbgp[0];
  e->b_t[t][1] = bbgp[1];
  e->b_t[t][2] = bbgp[2];
  unstack_bgp_vars(e,myx);
  if(set_vars(e,p,t,1))
    {
      return 1;
    }
  
  nx=0;

  myf[nx] = price_norm(e,t);
  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
    {
      fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nPrice norm");
      return 1;
    }
  
  nx=nx+1;

  for(i=0; i<(NC-1); i++)
    {
      myf[nx] = bop(p,e,t,i);
      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	{
	  fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nBop %d",i);
	  return 1;
	}

      nx = nx+1;
    }

  for(i=0; i<NC; i++)
    {
      myf[nx] = muc_mul(p,e,t,i);
      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	{
	  fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nMuc-Mul %d",i);
	  return 1;
	}

      nx=nx+1;

      for(s=0; s<NS; s++)
	{
	  myf[nx] = mpk_rk(p,e,t,i,s);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nMPK = Rk %d %d",i,s);
	      return 1;
	    }

	  nx=nx+1;

	  myf[nx] = mpl_w(p,e,t,i,s);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nMPL = W %d %d",i,s);
	      return 1;
	    }

	  nx=nx+1;

	  myf[nx] = mkt_clear_y(p,e,t,i,s);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nMkt clearing for y %d %d",i,s);
	      return 1;
	    }

	  nx=nx+1;
	  
	  if(s!=SVC && s!=CNS)
	    {
	      myf[nx] = mucs_mucr(p,e,t,i,s,SVC);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected!\nMucs-Mucr %d %d\n",i,s);
		  return 1;
		}
	      nx=nx+1;
	    }

	}
    }

  if(nx != solver_n)
    {
      fprintf(logfile,"Wrong number of bgp eqns! nx = %d, nbgp = %d\n",nx,solver_n);
      return 1;
    }

  for(i=0; i<nx; i++)
    {
      if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	{
	  fprintf(logfile,"Error evaluating bgp eqns! NaN/Inf detected in position %d!\n",i);
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

  solver_n = nbgp;

  uint fixl_ = fixl;
  fixl=1;
  
  alloc_solver_mem();
  set_initial_bgp_guess();

  gsl_multiroot_function_fdf f = {&bgp_func_f,&bgp_func_df,&bgp_func_fdf,solver_n,NULL};
  par=1;
  uint status = find_root_deriv_mkl(&f);
  if(status)
    {
      fprintf(logfile,"Error evaluating equilibrium function!\n");
    }

  fixl=fixl_;

  free_solver_mem();
  return status;
}

uint eval_eqm_conds(const double * myx, double * myf, uint tn)
{
  eqm * e = &(eee0[tn]);
  params * p = &(ppp0[tn]);
  uint i=0,s=0,nx=0;
  int t = 0;
  int t0 = 0;

  if(scenario>=1)
    {
      e = &(eee1[tn]);
      p = &(ppp1[tn]);
      t0=TSHOCK;
    }

  unstack_eqm_vars(e,myx);

  e->b_t[0][0] = p->b0[0];
  e->b_t[0][1] = p->b0[1];
  e->b_t[0][2] = p->b0[2];

  for(i=0; i<NC; i++)
    {
      for(s=0; s<NS; s++)
	{
	  e->k_t[0][i][s] = p->k0[i][s];
	}
    }

  for(t=t0; t<(NT+1); t++)
    {
      if(set_vars(e,p,t,0))
	{
	  fprintf(logfile,"Error calling set_vars1!\n");
	  return 1;
	}
    }

  nx=0;
  for(t=t0; t<(NT+1); t++)
    {
      myf[nx] = price_norm(e,t);
      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	{
	  fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\n Price norm %d\n",t);
	  return 1;
	}
      nx=nx+1;

      for(i=0; i<(NC-1); i++)
	{
	  myf[nx] = bop(p,e,t,i);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nBop %d %d\n",t,i);
	      return 1;
	    }
	  nx = nx+1;
	}

      for(i=0; i<NC; i++)
	{
	  myf[nx] = muc_mul(p,e,t,i);
	  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
	    {
	      fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMuc-Mul %d %d\n",t,i);
	      return 1;
	    }
	  nx=nx+1;

	  if(t<NT)
	    {
	      myf[nx] = euler(p,e,t,i);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nEuler %d %d\n",t,i);
		  return 1;
		}
	      nx = nx+1;	       
	    }	    

	  for(s=0; s<NS; s++)
	    {
	      if(s!= SVC && s!=CNS)
		{
		  myf[nx] = mucs_mucr(p,e,t,i,s,SVC);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMucs-Mucr %d %d %d\n",t,i,s);
		      return 1;
		    }
		  nx=nx+1;
		}

	      if(t<NT)
		{
		  myf[nx] = mpk_rk(p,e,t+1,i,s);
		  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		    {
		      fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = Rk %d %d %d\n",t,i,s);
		      return 1;
		    }
		  nx=nx+1;
		}

	      myf[nx] = mpl_w(p,e,t,i,s);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = Rk %d %d %d\n",t,i,s);
		  return 1;
		}
	      nx=nx+1;
		  
	      myf[nx] = mkt_clear_y(p,e,t,i,s);
	      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
		{
		  fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = Mkt clearing for y %d %d %d\n",t,i,s);
		  return 1;
		}
	      nx=nx+1;

	      if(m_adj_cost && t<NT)
		{
		  if(s!=CNS)
		    {
		      myf[nx] = prod_m_chk(p,e,t,i,s);
		      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
			{
			  fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = Prod m chk %d %d %d\n",t,i,s);
			  return 1;
			}			  
		      nx = nx+1;
			  
		      uint j;
		      for(j=0; j<NC; j++)
			{
			  myf[nx] = foc_m2(p,e,t,i,s,j);
			  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
			    {
			      fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = foc q2 %d %d %d\n",t,i,s);
			      return 1;
			    }		  
			  nx = nx+1;
			}
		    }
		}
		  
	      if(f_adj_cost && t<NT)
		{
		  if(s!=CNS)
		    {
		      myf[nx] = prod_q_chk(p,e,t,i,s);
		      if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
			{
			  fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = Prod q chk %d %d %d\n",t,i,s);
			  return 1;
			}			  
		      nx = nx+1;
			  
		      uint j;
		      for(j=0; j<NC; j++)
			{
			  myf[nx] = foc_q2(p,e,t,i,s,j);
			  if(gsl_isnan(myf[nx]) || gsl_isinf(myf[nx]))
			    {
			      fprintf(logfile,"Error evaluating eqm eqns! NaN/Inf detected!\nMPK = foc q2 %d %d %d\n",t,i,s);
			      return 1;
			    }			  
			  nx = nx+1;
			}
		    }
		}
	    }
	}
    }


  if(nx != neqm)
    {
      fprintf(logfile,"Error evaluating eqm eqns! nx = %d, neqm = %d\n",nx,neqm);
      return 1;
    }

  for(i=0; i<nx; i++)
    {
      if(gsl_isnan(myf[i]) || gsl_isinf(myf[i]))
	{
	  fprintf(logfile,"Error evaluating equilibrium conditions! NaN/Inf detected!\n");
	  return 1;
	}
    }
  
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
	      fprintf(logfile,"Error loading equilibrium guess from seed file!\n");
	      free_solver_mem();
	      return 1;
	    }
	}
      else
	{
	  if(set_initial_eqm_guess())
	    {
	      fprintf(logfile,"Error constructing equilibrium guess!\n");
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
	      fprintf(logfile,"Error loading equilibrium guess from seed file!\n");
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
	      fprintf(logfile,"Failed to stack variables from previous exercise!\n");
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
	fprintf(logfile,"Error evaluating equilibrium function!\n");
    }
  else
    {
      gsl_multiroot_function_fdf f = {&eqm_func_f,&eqm_func_df,&eqm_func_fdf,neqm,NULL};

      par=1;
      status = find_root_deriv_mkl(&f);
      if(status)
	{
	  fprintf(logfile,"Error solving for equilibrium!\n");
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
	     p->eps[i][0][2] * pow(e->c_t[t][i][2],p->rho),

	     p->phi[i]*p->psi/p->rho) * 
	 pow((p->lbar[i]-e->ll_t[t][i]),(1.0-p->phi[i])*p->psi));	
      
      for(t=(NT-1); t>=0; t--)
	{
	  e->welfare_t[t][i] = p->beta[i] * e->welfare_t[t+1][i] + 
	    (pow(p->eps[i][0][0] * pow(e->c_t[t][i][0],p->rho) +
		 p->eps[i][0][1] * pow(e->c_t[t][i][1],p->rho) + 
		 p->eps[i][0][2] * pow(e->c_t[t][i][2],p->rho),
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
