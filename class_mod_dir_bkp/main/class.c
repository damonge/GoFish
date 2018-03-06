/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }
  
  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }

  /****** all calculations done, now free the structures ******/

  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}


// int tune_scalar_field_parameters(
// 		                 struct precision * ppr,
// 		                 struct background * pba
// 		                 )  {
//   if (!pba->has_smg) {
//     printf(" no scalar field, skipping tune algorithm \n");
//     return _SUCCESS_; 
//   }
//   
//   /** The parameter to be tuned can be defined as a preprocessor variable */
// 
// //For galileon
// // #define _TUNE_PARAM_ pba->smg_Lambda_2
// //   int iter=0, maximum_iter = 1000; //maximum number of iterations
// //   double smg_TUNE_PARAM_min = -1000; //-2/pow(pba->H0,2); //_TUNE_PARAM_
// //   double smg_TUNE_PARAM_max = 100;//2.5/pow(pba->H0,2);
// //   double Omega0_smg_try = 0.;//absurd value
// //   double tolerance = 1e-4; //Need to pass a preccision argument
// 
// //For quintessence
// #define _TUNE_PARAM_ pba->smg_lambda
//   int iter=0, maximum_iter = 1000; //maximum number of iterations
//   double smg_TUNE_PARAM_min = 1; //-2/pow(pba->H0,2); //_TUNE_PARAM_
//   double smg_TUNE_PARAM_max = 100;//2.5/pow(pba->H0,2);
//   double Omega0_smg_try = 0.;//absurd value
//   double tolerance = 1e-4; //Need to pass a preccision argument
// 
//  int wish_background_verbose = pba->background_verbose;
//   pba->background_verbose = 0; /*remove unnecessary output while tuning*/  
//   
//   
//   if(pba->tuning_smg == _TRUE_) { 
// 
//     while (fabs(pba->Omega0_smg-Omega0_smg_try) > tolerance){
//       
//       _TUNE_PARAM_=0.5*(smg_TUNE_PARAM_max + smg_TUNE_PARAM_min);
//       
//       //NOTE: it would be nice if the string could be updated with the preprocessor _TUNE_PARAM_
//       if(wish_background_verbose>0)
// 	printf(" _TUNE_PARAM_ = %.3e, Omega_0_try = %e \n",_TUNE_PARAM_, Omega0_smg_try);
//     
//       if (background_init(ppr,pba) == _FAILURE_) {
// 	printf("\n\nError running background_init with _TUNE_PARAM_ %g \n=>%s\n",_TUNE_PARAM_,pba->error_message);
// 	return _FAILURE_;
//       }
//     
//       Omega0_smg_try = 
//       pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_rho_smg]/
//       pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_rho_crit];
//       /*this line returns the real value of Omega0_smg infered from the real evolution of the fields */
//     
//       if (Omega0_smg_try > pba->Omega0_smg) // might need to change < for >
// 	smg_TUNE_PARAM_min = _TUNE_PARAM_;
//       else
// 	smg_TUNE_PARAM_max = _TUNE_PARAM_;
//       /* note: the above four lines will find the correct value by bisection only if Omega0_smg_try is a monotonous function of sc_lambda. I hope that this is the case for the model we are after. It might not be the case in general. Then, one needs to think more. E.g. maybe one should pass values of smg_TUNE_PARAM_min and smg_TUNE_PARAM_max chosen not by chance, but knowning that they define the interval inside which the function Omega0_smg(smg_lambda) is monotonous. Moreover, depending on the fact that the function Omega0_smg(smg_lambda) is growing or decreasing, you may need to invert "min" and "max" inside the four lines above. */
//     
//       if (background_free(pba) == _FAILURE_) {
// 	printf("\n\nError in background_free \n=>%s\n",pba->error_message);
// 	return _FAILURE_;
//       }
//       /* note: it is important to call background_free now, otherwise we could not call background_init again */
//     
//       if (iter > maximum_iter) {
// 	printf("\n\nError in tune_parameters_smg, no tuning after %u iterations \n=>%s\n",iter, pba->error_message);
// 	return _FAILURE_;
//       }
//       iter ++;
//     
//     } /* end of bisection loop */
//   
//   }
//   
//   // in case that no tuning is wanted, run a background and overwrite the values of Omega0 for each component
//   else {
//     
//     if (background_init(ppr,pba) == _FAILURE_) {
//       //NOTE: it would be nice if the string could be updated with the preprocessor _TUNE_PARAM_
//       printf("\n\nError running background_init with _TUNE_PARAM_ %g \n=>%s\n",_TUNE_PARAM_,pba->error_message);
//       return _FAILURE_;
//     }
//    double Omega0_notune = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_rho_smg]/
//    pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_rho_crit];
// 
//    if(wish_background_verbose > 1)
//      printf("Scalar not tuned: Omega0_smg = %f, initial input %f \n",Omega0_notune, pba->Omega0_smg);
//    
//    /* Then the values introduced by the user are not the real ones, e.g. rho_m(a_0) = Omega_m(input)*H_0^2 (in CLASS units)
//     * This might be useful to run MCMC's without tuning parameters: 
//     * everything is computed and then the right values of Omega0_cdm, etc... 
//     * are stored as derived parameters (might introduce a bias though)
//     * TODO: consider updating the values of all the omegas after an non-tuned model is run.
//     */
// 
//   }
//   
//   // put the background verbose back
//   pba->background_verbose = wish_background_verbose;
//   
//   return _SUCCESS_;
// 
// }
