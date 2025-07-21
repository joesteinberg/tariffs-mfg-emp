#ifndef __SOLVER_C__
#define __SOLVER_C_

#include "solver.h"

//const double epsjac = GSL_SQRT_DBL_EPSILON;
const double epsjac = 1.0e-9;
const double root_tol = 1.0e-9;
const uint max_root_iter = 250;

int jacobian(
	     int (* F) (const gsl_vector * x, void * data, gsl_vector * f),
	     const gsl_vector * x,
	     gsl_matrix * J,
	     uint set_f0
	     )
{
  jcnt = jcnt + 1;

  if(set_f0)
    {
      (*F)(x,NULL,f0[0]);
    }

  if(par)
    {
      uint k;
      for(k=0;k<NTH;k++)
	{
	  if(k>0)
	    {
	      gsl_vector_memcpy(f0[k],f0[0]);
	    }
	}
    }

  uint abort=0;
  uint i;

#pragma omp parallel for private(i) schedule(static) if(par)
  for(i=0; i<solver_n; i++)
    {
#pragma omp flush(abort)
      if(!abort)
	{
	  uint it=0;

#ifdef _OPENMP
	  if(par)
	  {
	      it = omp_get_thread_num();
	  }
#endif
	  
	  gsl_vector_memcpy(xh[it],x);
	  //double h = fabs(gsl_vector_get(xh[it],i))*epsjac + epsjac;
	  double h = epsjac;
	  gsl_vector_set(xh[it],i,gsl_vector_get(xh[it],i)+h);

	  (*F)(xh[it],&it,fh[it]);
	  uint j;
	  for(j=0; j<solver_n; j++)
	    {	      
	      double tmp = (gsl_vector_get(fh[it],j)-gsl_vector_get(f0[it],j))/h;
	      gsl_matrix_set(J,j,i,tmp);

	      if(gsl_isnan(tmp))
		{
		  fprintf(logfile,"\t\tNaN detected in Jacobian processing...\n");
		  fprintf(logfile,"\t\tThread: %d\n",it);
		  fprintf(logfile,"\t\ti, x[i], xh[i]: %d, %0.4f, %0.4f\n",i,x->data[i],xh[it]->data[i]);
		  fprintf(logfile,"\t\tj, f[j], fh[j]: %d, %0.4f, %0.4f\n",j,f0[it]->data[j],fh[it]->data[j]);
		  abort = 1;
#pragma omp flush(abort)
		  break;
		}
	    }
	}
    }

  double col_chk[solver_n];
  double row_chk[solver_n];
  set_all_v(col_chk,solver_n,0.0);
  set_all_v(row_chk,solver_n,0.0);
  for(i=0; i<solver_n; i++)
    {
      uint j;
      for(j=0; j<solver_n; j++)
	{
	  if(gsl_isnan(gsl_matrix_get(J,j,i)))
	    {
	      fprintf(logfile,"\t\tElement (%d,%d) of Jacobian matrix is NaN!\n",i,j);
	      return GSL_EFAILED;
	    }
	  col_chk[i] = col_chk[i] + fabs(gsl_matrix_get(J,j,i));
	  row_chk[j] = row_chk[j] + fabs(gsl_matrix_get(J,j,i));
	}
    }
  
  for(i=0; i<solver_n; i++)
    {
      if(fabs(col_chk[i]<1.0e-14))
	{
	  fprintf(logfile,"Column (variable) %d of Jacobian matrix is all zeros!\n",i);
	  return GSL_EFAILED;
	}
      else if(fabs(row_chk[i]<1.0e-14))
	{
	  fprintf(logfile,"Row (function) %d of Jacobian matrix is all zeros!\n",i);
	  return GSL_EFAILED;
	}

    }

  return GSL_SUCCESS;
}

uint find_root_deriv_mkl(gsl_multiroot_function_fdf * f)
{
  fprintf(logfile,"\n\tInitializing solver on %zu-variable system...\n",f->n);

  uint status = 0, iter = 0;
  time_t start1 = 0, stop1 = 0, start2 = 0, stop2 = 0;
  time(&start1);

  gnewton_state_t * s = gnewton_alloc(f);

  jcnt=0;
  time(&start2);
  status = gnewton_set(s,solver_x);
  time(&stop2);

  if(status)
  {
    fprintf(logfile,"\t\tInitialization failed!\n");
  }
  else
    {
      status = GSL_CONTINUE;

      fprintf(logfile,"\t\tInitialization complete. Time = %0.2f\n",difftime(stop2,start2));
      double fmax = gsl_vector_max(s->f);
      double fmin = gsl_vector_min(s->f);
      if(fabs(fmax)<root_tol && fabs(fmin)<root_tol)
	{
	  fprintf(logfile,"\t\tInitial guess satisfies all equilibrium conditions!\n");	  
	  status = GSL_SUCCESS;
	}

      if(status == GSL_CONTINUE)
	{
	  do
	    {
	      iter++;
	      jcnt=0;
	      time(&start2);
	      status = gnewton_iterate(s);
	      time(&stop2);

	      if (status)   /* check if solver is stuck */
		break;


	      fmax = gsl_vector_max(s->f);
	      fmin = gsl_vector_min(s->f);
	      if(gsl_isnan(fmax) || gsl_isnan(fmin))
		{
		  fprintf(logfile,"\t\tNaN detected in function value!\n");
		  status = GSL_EFAILED;
		  break;
		}

	      if(gsl_isinf(fmax) || gsl_isinf(fmin))
		{
		  fprintf(logfile,"\t\tInf detected in function value!\n");
		  status = GSL_EFAILED;
		  break;
		}

	      //status = gsl_multiroot_test_residual (s->f, sqrt(sqrt(f->n)) * root_tol);
	      if(fabs(fmax)<root_tol && fabs(fmin)<root_tol)
		{
		  status = GSL_SUCCESS;
		}
	      else
		{
		  status = GSL_CONTINUE;
		}


	      if(status != GSL_SUCCESS && gsl_vector_max(s->dx) < 1.0e-14 && gsl_vector_min(s->dx) > -1.0e-14)
		{
		  status = GSL_EFAILED;
		  fprintf(logfile,"\t\tdx is zero vector but equations not solved!\n");
		  break;
		}

	      if(solver_verbose)
		{
		  if(f->n <= 4)
		    {
		      fprintf(logfile,"\t\t%d, %0.2f secs \tf = (",iter,difftime(stop2,start2));
		      uint i;
		      for(i=0; i<f->n; i++)
			{
			  fprintf(logfile,"\t%0.6g",gsl_vector_get(s->f,i));
			}
		      fprintf(logfile,"\t)\n");
		    }
		  else
		    {
		      fprintf(logfile,"\t\t%d, %0.0f secs\tmax = %0.2e, %d\tmin = %0.2e, %d\n",
			      iter,
			      difftime(stop2,start2),
			      gsl_vector_max(s->f),
			      (int)gsl_vector_max_index(s->f),
			      gsl_vector_min(s->f),
			      (int)gsl_vector_min_index(s->f));
		    }
		}
	    }
	  while (status == GSL_CONTINUE && iter < max_root_iter);
	}
      time(&stop1);

      double sum = 0.0;
      sum = vsum(s->f);
      fprintf(logfile,"\tFinished: %s.",gsl_strerror (status));
      if(status==0)
	{
	  fprintf(logfile," Time = %0.2f, iter = %u, sum |f| = %0.6g",
		  difftime(stop1,start1),iter,sum);
	}
      fprintf(logfile,"\n");
    }

  gsl_vector_memcpy(solver_x,s->x);
  gnewton_free(s);

  if(status==0)
    {
      return 0;
    }
  else
    {
      fprintf(logfile,"\tError solving system of equations!\n");
      return 1;
    }
}

#endif
