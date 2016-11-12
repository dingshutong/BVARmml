module BMA


INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use structures
use matutil
use VARutils
use Bayesutils

!================================================================================================================
contains 


!==================================== Model Bayesian Averaging =====================================!

subroutine bma_weights(N,H, m, q, p, I, ms, mc_burn,mc_reps,DATAMAT, dim_max, Models, priorspec, NUM, EVALvec,posteval,methods, &
!--------------- output 
lnmarglik,lnpredlik,bayes_fcasts,non_deterfcasts)

implicit none

!HYPER,LEVEL,ENDOEXO)  
			
! ------------------------- for preditive likelihood 
! S, S_eval, Burn, NUM_EVAL,RESMAT_TEMP, X_pred,&
!--------------------------- for forecasting 
! H,CRED_PRED, EVALvec, &					    &
					 

!This subroutine covers all combinations of VAR models from a data set
!The routine calls the Least squares routine 'freqanal', the prior routine
!'make_prior' and performs a posterior distribution Gibbs sampling 'gibbs_samp'
!The routine returns the predictive density of each model.  
!
!Written by M.K. Andersson, Sveriges Riksbank, 2004-09-22

! N				= number of obs
! m				= number of potentially endogenous vars
! q				= number of deterministic vars
! p				= number of lags in estimated model
! I				= number of endogenous variables in model
! H				= Forecast horizon
! S				= MCMC reps
! S_eval		= MCMC reps for predictive likelihood (if density estimation)
! Burn			= Burn-in for MCMC
! dim_max		= Max dimension of VAR
! NUM_EVAL		= Number of origins to evaluate the predictive likelihood at
! EVALvec		= Horizon(s) to evaluate predictive likelihood at
! CRED_PRED		= Level for credible intervals for forecasts
! DATAMAT		= Data matrix m modelled vars first, deterministic vars
! LEVEL			= Indicator if variable is used if levels or difference LEV/DIF
! ENDOEXO		= Indicator if variable is endogenous or exogenous END/EXO
! HYPER			= vector of hyper parameters for prior
! Models		= vector with indexes of variables to include in model
! X_pred		= deterministic variables for forecast period
! RESMAT_TEMP	= log predictive/marginal likelihood and log model prior (output)
! pred_temp		= (output)
! perc_bymodel	= (output)
! perc4_bymodel = (output)
! prior_val		= prior inclusion probabilities for variables
! posteval		= How to evaluate predictive likelihood
!				  UPDATE = update parameter prior as forecast origin is updated 
!				  BEFORE = use parameter prior for first forecast origin
!				  ALLDATA = use parameter prior based on all N obs
! priortype     = NORMDIFF or NORMWISH
! priordf		= prior degrees of freedom for Normal-Wishart

INTEGER, INTENT(IN)                              :: N,h,m,ms,q,p,I,dim_max, models(dim_max+1), evalvec(:), mc_burn,mc_reps, NUM

DOUBLE PRECISION, DIMENSION(N,m+q), INTENT(IN)   :: DATAMAT
character*20, intent(in)                         :: methods
!double precision, dimension(:,:), intent(in)     :: X_pred
type(PriorInfo),  intent(in)                      :: priorspec
type(VARPrior)                                    :: tPrior
character*12, intent(in)                         :: posteval
!type(modelinfo), intent(in)                      :: tmodel
real*8, intent(out)                               ::  lnpredlik,lnmarglik, bayes_fcasts(:,:)
integer,optional, intent(in)                     :: non_deterfcasts

!Double precision, intent(in) ,  optional         :: HYPER(:)
!Character*12, DIMENSION(m+q),INTENT(IN),optional :: LEVEL,ENDOEXO


integer k,T
integer varoutmod(m-I), vars(I+q)							
type(VARData) tData
!type(NDsampler):: tSampler

k = p*I + q
T = N-p
vars(1:I) = models(1:I)
vars(I+1:I+q) = (/ m+1:m+q /)


!**********************************************************************
!************** Transform to multivariate regression form *************
!**************                  Y = ZB + U               ************* 
!**********************************************************************

call make_data( Datamat(1:N,vars), I, p, q, 'NW', tData )

! Check on alignment of X_pred
!do jj = 1, q
!	do ii = T-evalh-num_eval+1, T
!		if ( tData%Z(ii,jj) /= X_pred(ii-T+evalh+num_eval,jj) ) then
!			print *, "Misalingment in X_pred, element,", ii, ",", jj, " in Z does not match"
!!			stop
!		endif
!	end do
!end do

call OLS( tData )

call make_Littermanprior( tData, priorspec, tPrior)


if ( priorspec%priortype == 'NW' ) then

    call bayes_marglik_NW( tdata%Y, tdata%ZPZ, tdata%Z, tprior%NW, lnmarglik )
    call bayes_predlik_NW(tdata,tprior%NW, NUM, EVALvec,mc_reps,posteval,lnpredlik )
    call bayes_forecasts_NW(tdata,tprior%NW,h,mc_reps,bayes_fcasts,non_deterfcasts)

else if (priorspec%priortype == 'ND') then

   call bayes_lnpredlik_NIW(tData,tPrior%ND,NUM,EVALvec,mc_reps,mc_burn,posteval,lnpredlik)
   call bayes_margins_NIW(tdata,tprior%ND,ms,mc_burn,mc_reps,mc_burn,mc_reps,methods,lnmarglik)
   call bayes_forecasts_NIW(tdata,tprior%ND,h,mc_burn,mc_reps,bayes_fcasts,non_deterfcasts)
     
end if 

end subroutine


!==================================== Bayesian Forecasting Combinations =====================================!


!==================================================================================================
subroutine bma_forecasts( mod_post, res_preds, no_mods_adj, h, preds_bma)

!***************************************************
!******* Bayes Factor Weighted Forecasts ***********
!***************************************************

implicit none 

double precision, intent(in) :: mod_post(:), res_preds(:,:)
integer, intent(in) :: no_mods_adj, h
double precision, intent(out) :: preds_bma(:)

integer kk, i

PREDS_BMA = 0.0d0
do kk = 1,no_mods_adj
	PREDS_BMA(1:h) = PREDS_BMA(1:h) + mod_post(kk)*res_preds(kk,1:h)
enddo

end subroutine


end module 