module simutils

INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use structures 
use times
use datautils
use priorspost
use bma
use modelutils
use postproc
use VARutils


contains


!===================================================================simulate one ==================================================================================

subroutine sim_result1(ms,method,tDgp,datamat,fcast_vec,prior,mod,varprob,dimprob,modprob,modselect,rmse)

integer,         intent(in) :: ms
type(PriorInfo), intent(in) :: prior
type(dgp), intent(in)       :: tdgp 
type(mods), intent(in)      :: mod
character*20, intent(in)    :: method
real*8,      intent(in)     :: datamat(:,:,:),fcast_vec(:,:)
real*8,      intent(out)    :: varprob(:,:,:),dimprob(:,:,:),modprob(:,:),modselect(:,:)
type(rmses), intent(out)    :: rmse

type(mats) :: mat

integer H,sim_i,s,i,mod_endo,mod_k,num_pr,median_mod

real*8, allocatable :: pred_temp(:,:),lnmarglik(:),lnpredlik(:),ftemp(:,:), &
                       mod_post(:,:,:),sim_vprob(:,:,:,:),sim_dprob(:,:,:,:),pred_bma(:),pred_sort(:,:), &
                       sum_select(:,:),vprob_mean(:),dprob_mean(:),sum_modpost(:,:)
! initialize 


H=tdgp%H
s=tdgp%simreps
m=tdgp%m
num_pr=size(mod%lnprior,2)

! allocation

allocate(lnmarglik(mod%n),lnpredlik(mod%n),ftemp(mod%n,tdgp%H),mod_post(mod%n,num_pr,2),pred_sort(mod%n,tdgp%H))
allocate(sim_vprob(tdgp%m,tdgp%simreps,num_pr,2),sim_dprob(tdgp%dims(2),tdgp%simreps,num_pr,2),pred_bma(tdgp%H))
allocate(sum_select(num_pr,2),sum_modpost(num_pr,2),vprob_mean(tdgp%m),dprob_mean(tdgp%dims(2)))


call clear_rmses(rmse)

call init_rmses(H,2,num_pr,2,rmse)


call clear_mats(mat)

allocate(mat%bma(tdgp%H,tdgp%simreps,num_pr,2),mat%top(tdgp%H,tdgp%simreps,num_pr,2),mat%med(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%rw(tdgp%H,tdgp%simreps,num_pr,2),mat%ar(tdgp%H,tdgp%simreps,num_pr,2),mat%rmean(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%main(tdgp%H,tdgp%simreps,num_pr,2))

! run

sum_modpost = 0d0
sum_select =0d0


do sim_i = 1,s 

Write(*,*)  '  '
write(*,*)  '**********************************  ' 
write(*,*)  'Monte Carlo Peplicate  ', sim_i
write(*,*)  '**********************************  '


Do mod_i=1,mod%n 
     mod_endo=mod%models(mod_i,tdgp%dims(2)+1)
     mod_k   =mod_endo*mod%p+mod%q
     allocate(pred_temp(tdgp%H,mod_endo))
     call bma_weights(tdgp%N,tdgp%H, tdgp%m, mod%q, mod%p, mod_endo,ms,tdgp%mcburn,tdgp%mcreps,DATAMAT(:,:,sim_i), tdgp%dims(2), mod%Models(mod_i,:), prior, mod%NUM, mod%EVALvec,mod%posteval,method, &
lnmarglik(mod_i),lnpredlik(mod_i),pred_temp)  
     ! first endogeneous variable y is our variable of interest
    ftemp(mod_i,:)=pred_temp(:,1)
    deallocate(pred_temp)
end do




! For the model Bayes weights

do i=1,num_pr

   mod_post(:,i,1) = lnpredlik
   mod_post(:,i,2) = lnmarglik
        
  do j=1,2
    
    median_mod=1
        
  ! for marginal and predictive likelihood weights
    mod_post(:,i,j)=mod_post(:,i,j) + mod%lnprior(:,i) 
    mod_post(:,i,j)=exp(mod_post(:,i,j)-maxval(mod_post(:,i,j)))
    mod_post(:,i,j)=mod_post(:,i,j)/sum(mod_post(:,i,j))
    call calc_varprobs(mod%n,tdgp%dims(2),tdgp%m,mod_post(:,i,j),mod%models,sim_vprob(:,sim_i,i,j),sim_dprob(:,sim_i,i,j),median_mod)
    call bma_forecasts( mod_post(:,i,j),ftemp,mod%n,tdgp%H,pred_bma)
    call mod_summary_new( mod_post(:,i,j), ftemp, tdgp%dims, mod%models, mod%true, mod%n, mod%no_mods_del, &
					  sum_modpost(i,j), pred_sort, sum_select(i,j))
    call fcast4mc( pred_bma, fcast_vec(:,sim_i), pred_sort, ftemp, median_mod, tdgp%H, mat%bma(:,sim_i,i,j), mat%top(:,sim_i,i,j), mat%ar(:,sim_i,i,j), mat%med(:,sim_i,i,j), mat%rw(:,sim_i,i,j), mat%rmean(:,sim_i,i,j), mat%main(:,sim_i,i,j))
  end do
end do


end do


!*********************************************
! post simulation calculations
!*********************************************
! Variable Probabilities

do i=1,num_pr
  
  do j=1,2
   
    vprob_mean = 0.0d0
    dprob_mean = 0.0d0
    do sim_i = 1,tdgp%simreps
	    vprob_mean = vprob_mean + sim_vprob(:,sim_i,i,j)
	    dprob_mean = dprob_mean + sim_dprob(:,sim_i,i,j)
    enddo
    varprob(:,i,j)= vprob_mean/tdgp%simreps
    dimprob(:,i,j)= dprob_mean/tdgp%simreps
    modprob(i,j)  = sum_modpost(i,j)/tdgp%simreps
    modselect(i,j)= sum_select(i,j)/tdgp%simreps
    
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%bma(:,:,i,j),mat%Main(:,:,i,j),rmse%bma(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%top(:,:,i,j),mat%Main(:,:,i,j),rmse%top(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%ar(:,:,i,j),mat%Main(:,:,i,j),rmse%ar(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%med(:,:,i,j),mat%Main(:,:,i,j),rmse%med(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rw(:,:,i,j),mat%Main(:,:,i,j),rmse%rw(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rmean(:,:,i,j),mat%Main(:,:,i,j),rmse%rmean(:,:,i,j))
    
  end do
end do

deallocate(lnmarglik,lnpredlik,ftemp,mod_post,pred_sort)
deallocate(sim_vprob,sim_dprob,pred_bma)
deallocate(sum_select,sum_modpost,vprob_mean,dprob_mean)


end subroutine sim_result1











!===================================================================simulate two ==================================================================================

subroutine sim_result2(ms,method,tDgp,datamat,fcast_vec,prior,mod,vartable,dimtable,modtable,rmse)

implicit none

integer,         intent(in) :: ms
type(PriorInfo), intent(in) :: prior
type(dgp), intent(in)       :: tdgp 
type(mods), intent(in)      :: mod
character*20, intent(in)    :: method
real*8,      intent(in)     :: datamat(:,:,:),fcast_vec(:,:)
real*8,      intent(out)    :: vartable(:,:),dimtable(:,:),modtable(:,:)
type(rmses), intent(out)    :: rmse

type(mats) :: mat

integer H,sim_i,s,i,j,mod_i,mod_endo,mod_k,num_pr,median_mod,kk,ll,jj,j1,l1,dd,j2,d2,m

real*8, allocatable :: pred_temp(:,:),lnmarglik(:),lnpredlik(:),ftemp(:,:), &
                       mod_post(:,:,:),sim_vprob(:,:,:,:),sim_dprob(:,:,:,:),pred_bma(:),pred_sort(:,:), &
                       sum_select(:,:),sum_modpost(:,:),vtemp(:),vsort(:),  & 
                       varprob(:,:,:),dimprob(:,:,:),modprob(:,:),modselect(:,:)
! initialize 


H=tdgp%H
s=tdgp%simreps
m=tdgp%m
d2=tdgp%dims(2)
num_pr=size(mod%lnprior,2)

! allocation

allocate(lnmarglik(mod%n),lnpredlik(mod%n),ftemp(mod%n,tdgp%H),mod_post(mod%n,num_pr,2),pred_sort(mod%n,tdgp%H))
allocate(sim_vprob(tdgp%m,tdgp%simreps,num_pr,2),sim_dprob(tdgp%dims(2),tdgp%simreps,num_pr,2),pred_bma(tdgp%H))
allocate(sum_select(num_pr,2),sum_modpost(num_pr,2))
allocate(varprob(m,num_pr,2),dimprob(D2,num_pr,2),modprob(num_pr,2),modselect(num_pr,2))
allocate(vtemp(tdgp%m),vsort((tdgp%m_syst-1)*2+1))

!
call clear_rmses(rmse)

call init_rmses(H,2,num_pr,2,rmse)


call clear_mats(mat)

allocate(mat%bma(tdgp%H,tdgp%simreps,num_pr,2),mat%top(tdgp%H,tdgp%simreps,num_pr,2),mat%med(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%rw(tdgp%H,tdgp%simreps,num_pr,2),mat%ar(tdgp%H,tdgp%simreps,num_pr,2),mat%rmean(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%main(tdgp%H,tdgp%simreps,num_pr,2))

! run

sum_modpost = 0d0
sum_select =0d0
kk = (tdgp%m_syst-1)*2+1
dd = tdgp%dims(2)

do sim_i = 1,s 

Write(*,*)  '  '
write(*,*)  '**********************************  ' 
write(*,*)  'Monte Carlo Peplicate  ', sim_i
write(*,*)  '**********************************  '

Do mod_i=1,mod%n 
     mod_endo=mod%models(mod_i,tdgp%dims(2)+1)
     mod_k   =mod_endo*mod%p+mod%q
     allocate(pred_temp(tdgp%H,mod_endo))
     call bma_weights(tdgp%N,tdgp%H, tdgp%m, mod%q, mod%p, mod_endo,ms,tdgp%mcburn,tdgp%mcreps,DATAMAT(:,:,sim_i), tdgp%dims(2), mod%Models(mod_i,:), prior, mod%NUM, mod%EVALvec,mod%posteval,method, &
lnmarglik(mod_i),lnpredlik(mod_i),pred_temp)  
     ! first endogeneous variable y is our variable of interest
    ftemp(mod_i,:)=pred_temp(:,1)
    deallocate(pred_temp)
end do




! For the model Bayes weights

do i=1,num_pr

   mod_post(:,i,1) = lnpredlik
   mod_post(:,i,2) = lnmarglik
        
  do j=1,2
    
    median_mod=1
        
  ! for marginal and predictive likelihood weights
    mod_post(:,i,j)=mod_post(:,i,j) + mod%lnprior(:,i) 
    mod_post(:,i,j)=exp(mod_post(:,i,j)-maxval(mod_post(:,i,j)))
    mod_post(:,i,j)=mod_post(:,i,j)/sum(mod_post(:,i,j))
    call calc_varprobs(mod%n,tdgp%dims(2),tdgp%m,mod_post(:,i,j),mod%models,sim_vprob(:,sim_i,i,j),sim_dprob(:,sim_i,i,j),median_mod)
    call bma_forecasts( mod_post(:,i,j),ftemp,mod%n,tdgp%H,pred_bma)
    call mod_summary_new( mod_post(:,i,j), ftemp, tdgp%dims, mod%models, mod%true, mod%n, mod%no_mods_del, &
					  sum_modpost(i,j), pred_sort, sum_select(i,j))
    call fcast4mc( pred_bma, fcast_vec(:,sim_i), pred_sort, ftemp, median_mod, tdgp%H, mat%bma(:,sim_i,i,j), mat%top(:,sim_i,i,j), mat%ar(:,sim_i,i,j), mat%med(:,sim_i,i,j), mat%rw(:,sim_i,i,j), mat%rmean(:,sim_i,i,j), mat%main(:,sim_i,i,j))
  end do
end do

end do


!*********************************************
! post simulation calculations
!*********************************************
! Variable Probabilities

do i=1,num_pr
  
  do j=1,2
   
    varprob(1:m,i,j) = 0.0d0
    dimprob(1:d2,i,j) = 0.0d0
    do sim_i = 1,tdgp%simreps
	    varprob(1:m,i,j) = varprob(1:m,i,j) + sim_vprob(1:m,sim_i,i,j)
	    dimprob(1:d2,i,j) = dimprob(1:d2,i,j) + sim_dprob(1:d2,sim_i,i,j)
    enddo

    modprob(i,j)  = sum_modpost(i,j)/tdgp%simreps
    modselect(i,j)= sum_select(i,j)/tdgp%simreps
    
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%bma(:,:,i,j),mat%Main(:,:,i,j),rmse%bma(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%top(:,:,i,j),mat%Main(:,:,i,j),rmse%top(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%ar(:,:,i,j),mat%Main(:,:,i,j),rmse%ar(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%med(:,:,i,j),mat%Main(:,:,i,j),rmse%med(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rw(:,:,i,j),mat%Main(:,:,i,j),rmse%rw(:,:,i,j))
    call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rmean(:,:,i,j),mat%Main(:,:,i,j),rmse%rmean(:,:,i,j))
    
  end do
end do

 varprob= varprob/tdgp%simreps
 dimprob= dimprob/tdgp%simreps


do i = 1, num_pr
    
    jj=1;ll=kk
    j1=1;l1=dd
    j2=0
    
    do j =1,2    
      vtemp =  varprob(1:m,i,j)
      vsort(1:(tdgp%m_syst-1))= vtemp(2:tdgp%m_syst) 
      vsort(tdgp%m_syst) = maxval(vtemp((tdgp%m_syst+1):))
      vsort((tdgp%m_syst+1):) = vsort(1:(tdgp%m_syst-1))/vsort(tdgp%m_syst)
      vartable(i,jj:ll)=vsort
      dimtable(i,j1:l1)= dimprob(1:d2,i,j)
      modtable(i,1+j2) = modprob(i,j)
      modtable(i,2+j2) = modselect(i,j)
      jj=jj+kk ; ll = ll+kk
      j1=j1+dd ; l1 = l1+dd
      j2=j2+2
    end do
    
end do     


deallocate(lnmarglik,lnpredlik,ftemp,mod_post,pred_sort)
deallocate(sim_vprob,sim_dprob,pred_bma)
deallocate(sum_select,sum_modpost)
deallocate(varprob,dimprob,modprob,modselect)
deallocate(vtemp,vsort)



end subroutine sim_result2



!===================================================================simulate three ==================================================================================

subroutine sim_result3(ms,ar_p,method,tDgp,datamat,fcast_vec,prior,mod,vartable,dimtable,modtable,rmse,non_deterfcasts)

implicit none

integer,         intent(in) :: ms,ar_p
type(PriorInfo), intent(in) :: prior
type(dgp), intent(in)       :: tdgp 
type(mods), intent(in)      :: mod
character*20, intent(in)    :: method
real*8,      intent(in)     :: datamat(:,:,:),fcast_vec(:,:)
real*8,      intent(out)    :: vartable(:,:),dimtable(:,:),modtable(:,:)
type(rmses), intent(out)    :: rmse

type(mats) :: mat

integer H,sim_i,s,i,j,mod_i,mod_endo,mod_k,num_pr,median_mod,kk,ll,jj,j1,l1,dd,j2,d2,m

real*8, allocatable :: pred_temp(:,:),lnmarglik(:),lnpredlik(:),ftemp(:,:), &
                       mod_post(:,:,:),sim_vprob(:,:,:,:),sim_dprob(:,:,:,:),pred_bma(:),pred_sort(:,:), &
                       sum_select(:,:),sum_modpost(:,:),vtemp(:),vsort(:),  & 
                       varprob(:,:,:),dimprob(:,:,:),modprob(:,:),modselect(:,:),pred_ar(:)

integer,optional, intent(in)  :: non_deterfcasts                       
                       
! initialize 


H=tdgp%H
s=tdgp%simreps
m=tdgp%m
d2=tdgp%dims(2)
num_pr=size(mod%lnprior,2)

! allocation

allocate(lnmarglik(mod%n),lnpredlik(mod%n),ftemp(mod%n,tdgp%H),mod_post(mod%n,num_pr,2),pred_sort(mod%n,tdgp%H))
allocate(sim_vprob(tdgp%m,tdgp%simreps,num_pr,2),sim_dprob(tdgp%dims(2),tdgp%simreps,num_pr,2),pred_bma(tdgp%H))
allocate(sum_select(num_pr,2),sum_modpost(num_pr,2))
allocate(varprob(m,num_pr,2),dimprob(D2,num_pr,2),modprob(num_pr,2),modselect(num_pr,2))
allocate(vtemp(tdgp%m),vsort((tdgp%m_syst-1)*2+1))
allocate(pred_ar(H))
!
call clear_rmses(rmse)

call init_rmses(H,2,num_pr,2,rmse)


call clear_mats(mat)

allocate(mat%bma(tdgp%H,tdgp%simreps,num_pr,2),mat%top(tdgp%H,tdgp%simreps,num_pr,2),mat%med(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%rw(tdgp%H,tdgp%simreps,num_pr,2),mat%ar(tdgp%H,tdgp%simreps,num_pr,2),mat%rmean(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%main(tdgp%H,tdgp%simreps,num_pr,2))

! run

sum_modpost = 0d0
sum_select =0d0
kk = (tdgp%m_syst-1)*2+1
dd = tdgp%dims(2)

do sim_i = 1,s 

Write(*,*)  '  '
write(*,*)  '**********************************  ' 
write(*,*)  'Monte Carlo Peplicate  ', sim_i
write(*,*)  '**********************************  '

Do mod_i=1,mod%n 
     mod_endo=mod%models(mod_i,tdgp%dims(2)+1)
     mod_k   =mod_endo*mod%p+mod%q
     allocate(pred_temp(tdgp%H,mod_endo))
     call bma_weights(tdgp%N,tdgp%H, tdgp%m, mod%q, mod%p, mod_endo,ms,tdgp%mcburn,tdgp%mcreps,DATAMAT(:,:,sim_i), tdgp%dims(2), mod%Models(mod_i,:), prior, mod%NUM, mod%EVALvec,mod%posteval,method, &
lnmarglik(mod_i),lnpredlik(mod_i),pred_temp,non_deterfcasts)  
     ! first endogeneous variable y is our variable of interest
    ftemp(mod_i,:)=pred_temp(:,1)
    deallocate(pred_temp)
end do


   if (ar_p .ne.  mod%p)  then 
        call   fcast_uar(ar_p,datamat(:,:,sim_i),tdgp%N,tdgp%m,mod%q,tdgp%H,prior,tdgp%mcreps,tdgp%mcburn,pred_ar,non_deterfcasts)
    end if 

! For the model Bayes weights

do i=1,num_pr

   mod_post(:,i,1) = lnpredlik
   mod_post(:,i,2) = lnmarglik
        
  do j=1,2
    
    median_mod=1
        
  ! for marginal and predictive likelihood weights
    mod_post(:,i,j)=mod_post(:,i,j) + mod%lnprior(:,i) 
    mod_post(:,i,j)=exp(mod_post(:,i,j)-maxval(mod_post(:,i,j)))
    mod_post(:,i,j)=mod_post(:,i,j)/sum(mod_post(:,i,j))
    call calc_varprobs(mod%n,tdgp%dims(2),tdgp%m,mod_post(:,i,j),mod%models,sim_vprob(:,sim_i,i,j),sim_dprob(:,sim_i,i,j),median_mod)
    call bma_forecasts( mod_post(:,i,j),ftemp,mod%n,tdgp%H,pred_bma)
    call mod_summary_new( mod_post(:,i,j), ftemp, tdgp%dims, mod%models, mod%true, mod%n, mod%no_mods_del, &
					  sum_modpost(i,j), pred_sort, sum_select(i,j))
    call fcast4mc( pred_bma, fcast_vec(:,sim_i), pred_sort, ftemp, median_mod, tdgp%H, mat%bma(:,sim_i,i,j), mat%top(:,sim_i,i,j), mat%ar(:,sim_i,i,j), mat%med(:,sim_i,i,j), mat%rw(:,sim_i,i,j), mat%rmean(:,sim_i,i,j), mat%main(:,sim_i,i,j))
    if (ar_p .ne. mod%p) then
      mat%ar(:,sim_i,i,j) = pred_ar
    end if 
  end do
end do

end do


!*********************************************
! post simulation calculations
!*********************************************
! Variable Probabilities

do i=1,num_pr
  
  do j=1,2
   
    varprob(1:m,i,j) = 0.0d0
    dimprob(1:d2,i,j) = 0.0d0
    do sim_i = 1,tdgp%simreps
	    varprob(1:m,i,j) = varprob(1:m,i,j) + sim_vprob(1:m,sim_i,i,j)
	    dimprob(1:d2,i,j) = dimprob(1:d2,i,j) + sim_dprob(1:d2,sim_i,i,j)
    enddo

    modprob(i,j)  = sum_modpost(i,j)/tdgp%simreps
    modselect(i,j)= sum_select(i,j)/tdgp%simreps
    
    call calc_mse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%bma(:,:,i,j),mat%Main(:,:,i,j),rmse%bma(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%top(:,:,i,j),mat%Main(:,:,i,j),rmse%top(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%ar(:,:,i,j),mat%Main(:,:,i,j),rmse%ar(:,:,i,j))
    !call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%med(:,:,i,j),mat%Main(:,:,i,j),rmse%med(:,:,i,j))
    !call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rw(:,:,i,j),mat%Main(:,:,i,j),rmse%rw(:,:,i,j))
    !call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rmean(:,:,i,j),mat%Main(:,:,i,j),rmse%rmean(:,:,i,j))
    
  end do
end do

 varprob= varprob/tdgp%simreps
 dimprob= dimprob/tdgp%simreps


do i = 1, num_pr
    
    jj=1;ll=kk
    j1=1;l1=dd
    j2=0
    
    do j =1,2    
      vtemp =  varprob(1:m,i,j)
      vsort(1:(tdgp%m_syst-1))= vtemp(2:tdgp%m_syst) 
      vsort(tdgp%m_syst) = maxval(vtemp((tdgp%m_syst+1):))
      vsort((tdgp%m_syst+1):) = vsort(1:(tdgp%m_syst-1))/vsort(tdgp%m_syst)
      vartable(i,jj:ll)=vsort
      dimtable(i,j1:l1)= dimprob(1:d2,i,j)
      modtable(i,1+j2) = modprob(i,j)
      modtable(i,2+j2) = modselect(i,j)
      jj=jj+kk ; ll = ll+kk
      j1=j1+dd ; l1 = l1+dd
      j2=j2+2
    end do
    
end do     


deallocate(lnmarglik,lnpredlik,ftemp,mod_post,pred_sort)
deallocate(sim_vprob,sim_dprob,pred_bma)
deallocate(sum_select,sum_modpost)
deallocate(varprob,dimprob,modprob,modselect)
deallocate(vtemp,vsort)



end subroutine sim_result3


!===================================================================simulate four ==================================================================================

subroutine sim_result4(ms,ar_p,method,tDgp,datamat,fcast_vec,prior,mod,vartable,dimtable,modtable,rmse,non_deterfcasts,simrange)

implicit none

integer,         intent(in) :: ms,ar_p
type(PriorInfo), intent(in) :: prior
type(dgp), intent(in)       :: tdgp 
type(mods), intent(in)      :: mod
character*20, intent(in)    :: method
real*8,      intent(in)     :: datamat(:,:,:),fcast_vec(:,:)
real*8,      intent(out)    :: vartable(:,:),dimtable(:,:),modtable(:,:)
type(rmses), intent(out)    :: rmse

type(mats) :: mat

integer H,sim_i,s,i,j,mod_i,mod_endo,mod_k,num_pr,median_mod,kk,ll,jj,j1,l1,dd,j2,d2,m,sim_ii

real*8, allocatable :: pred_temp(:,:),lnmarglik(:),lnpredlik(:),ftemp(:,:), &
                       mod_post(:,:,:),sim_vprob(:,:,:,:),sim_dprob(:,:,:,:),pred_bma(:),pred_sort(:,:), &
                       sum_select(:,:),sum_modpost(:,:),vtemp(:),vsort(:),  & 
                       varprob(:,:,:),dimprob(:,:,:),modprob(:,:),modselect(:,:),pred_ar(:)

integer,optional, intent(in)  :: non_deterfcasts, simrange(:)                       
                       
! initialize 


H=tdgp%H
s=tdgp%simreps
m=tdgp%m
d2=tdgp%dims(2)
num_pr=size(mod%lnprior,2)

! allocation

allocate(lnmarglik(mod%n),lnpredlik(mod%n),ftemp(mod%n,tdgp%H),mod_post(mod%n,num_pr,2),pred_sort(mod%n,tdgp%H))
allocate(sim_vprob(tdgp%m,tdgp%simreps,num_pr,2),sim_dprob(tdgp%dims(2),tdgp%simreps,num_pr,2),pred_bma(tdgp%H))
allocate(sum_select(num_pr,2),sum_modpost(num_pr,2))
allocate(varprob(m,num_pr,2),dimprob(D2,num_pr,2),modprob(num_pr,2),modselect(num_pr,2))
allocate(vtemp(tdgp%m),vsort((tdgp%m_syst-1)*2+1))
allocate(pred_ar(H))
!
call clear_rmses(rmse)

call init_rmses(H,2,num_pr,2,rmse)


call clear_mats(mat)

allocate(mat%bma(tdgp%H,tdgp%simreps,num_pr,2),mat%top(tdgp%H,tdgp%simreps,num_pr,2),mat%med(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%rw(tdgp%H,tdgp%simreps,num_pr,2),mat%ar(tdgp%H,tdgp%simreps,num_pr,2),mat%rmean(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%main(tdgp%H,tdgp%simreps,num_pr,2))

! run

sum_modpost = 0d0
sum_select =0d0
kk = (tdgp%m_syst-1)*2+1
dd = tdgp%dims(2)

do sim_i = simrange(1),simrange(2) 

sim_ii = sim_i - (simrange(1) - 1)

Write(*,*)  '  '
write(*,*)  '**********************************  ' 
write(*,*)  'Monte Carlo Peplicate  ', sim_i
write(*,*)  '**********************************  '

Do mod_i=1,mod%n 
     mod_endo=mod%models(mod_i,tdgp%dims(2)+1)
     mod_k   =mod_endo*mod%p+mod%q
     allocate(pred_temp(tdgp%H,mod_endo))
     call bma_weights(tdgp%N,tdgp%H, tdgp%m, mod%q, mod%p, mod_endo,ms,tdgp%mcburn,tdgp%mcreps,DATAMAT(:,:,sim_i), tdgp%dims(2), mod%Models(mod_i,:), prior, mod%NUM, mod%EVALvec,mod%posteval,method, &
lnmarglik(mod_i),lnpredlik(mod_i),pred_temp,non_deterfcasts)  
     ! first endogeneous variable y is our variable of interest
    ftemp(mod_i,:)=pred_temp(:,1)
    deallocate(pred_temp)
end do


   if (ar_p .ne.  mod%p)  then 
        call   fcast_uar(ar_p,datamat(:,:,sim_i),tdgp%N,tdgp%m,mod%q,tdgp%H,prior,tdgp%mcreps,tdgp%mcburn,pred_ar,non_deterfcasts)
    end if 

! For the model Bayes weights

do i=1,num_pr

   mod_post(:,i,1) = lnpredlik
   mod_post(:,i,2) = lnmarglik
        
  do j=1,2
    
    median_mod=1
        
  ! for marginal and predictive likelihood weights
    mod_post(:,i,j)=mod_post(:,i,j) + mod%lnprior(:,i) 
    mod_post(:,i,j)=exp(mod_post(:,i,j)-maxval(mod_post(:,i,j)))
    mod_post(:,i,j)=mod_post(:,i,j)/sum(mod_post(:,i,j))
    call calc_varprobs(mod%n,tdgp%dims(2),tdgp%m,mod_post(:,i,j),mod%models,sim_vprob(:,sim_ii,i,j),sim_dprob(:,sim_ii,i,j),median_mod)
    call bma_forecasts( mod_post(:,i,j),ftemp,mod%n,tdgp%H,pred_bma)
    call mod_summary_new( mod_post(:,i,j), ftemp, tdgp%dims, mod%models, mod%true, mod%n, mod%no_mods_del, &
					  sum_modpost(i,j), pred_sort, sum_select(i,j))
    call fcast4mc( pred_bma, fcast_vec(:,sim_i), pred_sort, ftemp, median_mod, tdgp%H, mat%bma(:,sim_ii,i,j), mat%top(:,sim_ii,i,j), mat%ar(:,sim_ii,i,j), mat%med(:,sim_ii,i,j), mat%rw(:,sim_ii,i,j), mat%rmean(:,sim_ii,i,j), mat%main(:,sim_ii,i,j))
    if (ar_p .ne. mod%p) then
      mat%ar(:,sim_ii,i,j) = pred_ar
    end if 
  end do
end do

end do


!*********************************************
! post simulation calculations
!*********************************************
! Variable Probabilities

do i=1,num_pr
  
  do j=1,2
   
    varprob(1:m,i,j) = 0.0d0
    dimprob(1:d2,i,j) = 0.0d0
    do sim_i = 1,tdgp%simreps
	    varprob(1:m,i,j) = varprob(1:m,i,j) + sim_vprob(1:m,sim_i,i,j)
	    dimprob(1:d2,i,j) = dimprob(1:d2,i,j) + sim_dprob(1:d2,sim_i,i,j)
    enddo

    modprob(i,j)  = sum_modpost(i,j)/tdgp%simreps
    modselect(i,j)= sum_select(i,j)/tdgp%simreps
    
    call calc_mse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%bma(:,:,i,j),mat%Main(:,:,i,j),rmse%bma(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%top(:,:,i,j),mat%Main(:,:,i,j),rmse%top(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%ar(:,:,i,j),mat%Main(:,:,i,j),rmse%ar(:,:,i,j))
    !call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%med(:,:,i,j),mat%Main(:,:,i,j),rmse%med(:,:,i,j))
    !call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rw(:,:,i,j),mat%Main(:,:,i,j),rmse%rw(:,:,i,j))
    !call calc_rmse(tdgp%H,tdgp%simreps,mod%rmse_obs,mat%rmean(:,:,i,j),mat%Main(:,:,i,j),rmse%rmean(:,:,i,j))
    
  end do
end do

 varprob= varprob/tdgp%simreps
 dimprob= dimprob/tdgp%simreps


do i = 1, num_pr
    
    jj=1;ll=kk
    j1=1;l1=dd
    j2=0
    
    do j =1,2    
      vtemp =  varprob(1:m,i,j)
      vsort(1:(tdgp%m_syst-1))= vtemp(2:tdgp%m_syst) 
      vsort(tdgp%m_syst) = maxval(vtemp((tdgp%m_syst+1):))
      vsort((tdgp%m_syst+1):) = vsort(1:(tdgp%m_syst-1))/vsort(tdgp%m_syst)
      vartable(i,jj:ll)=vsort
      dimtable(i,j1:l1)= dimprob(1:d2,i,j)
      modtable(i,1+j2) = modprob(i,j)
      modtable(i,2+j2) = modselect(i,j)
      jj=jj+kk ; ll = ll+kk
      j1=j1+dd ; l1 = l1+dd
      j2=j2+2
    end do
    
end do     


deallocate(lnmarglik,lnpredlik,ftemp,mod_post,pred_sort)
deallocate(sim_vprob,sim_dprob,pred_bma)
deallocate(sum_select,sum_modpost)
deallocate(varprob,dimprob,modprob,modselect)
deallocate(vtemp,vsort)



end subroutine sim_result4





!===================================================================Real data sets==================================================================================


subroutine reals_result1(data_N,ms,ar_p,method,tDgp,datamat,prior,mod,vartable,dimtable,rmse,mat,non_deterfcasts,simvprob,simdprob,rmse_obs)

implicit none

integer,         intent(in) :: ms,ar_p,data_n
type(PriorInfo), intent(in) :: prior
type(dgp), intent(in)       :: tdgp 
type(mods), intent(in)      :: mod
character*20, intent(in)    :: method
real*8,      intent(in)     :: datamat(:,:)
real*8,      intent(out)     :: vartable(:,:),dimtable(:,:)
integer,optional, intent(out) :: rmse_obs(:)  
real*8, optional, intent(out):: simvprob(:,:,:,:),simdprob(:,:,:,:)
type(rmses), intent(out)    :: rmse
type(mats),  intent(out)    :: mat
integer,optional, intent(in)  :: non_deterfcasts    
integer H,sim_i,s,i,j,mod_i,mod_endo,mod_k,num_pr,median_mod,kk,ll,jj,j1,l1,dd,j2,d2,m,n_est,h_i

real*8, allocatable :: pred_temp(:,:),lnmarglik(:),lnpredlik(:),ftemp(:,:), &
                       mod_post(:,:,:),sim_vprob(:,:,:,:),sim_dprob(:,:,:,:),pred_bma(:),pred_sort(:,:), &
                       sum_select(:,:),sum_modpost(:,:), varprob(:,:,:),dimprob(:,:,:), & 
                       modprob(:,:),modselect(:,:),pred_ar(:),fcast_vec(:)

integer, allocatable :: rmseobs(:)
                                        
! initialize 


H=tdgp%H
s=tdgp%simreps
m=tdgp%m
d2=tdgp%dims(2)
num_pr=size(mod%lnprior,2)

! allocation

allocate(lnmarglik(mod%n),lnpredlik(mod%n),ftemp(mod%n,tdgp%H),mod_post(mod%n,num_pr,2),pred_sort(mod%n,tdgp%H))
allocate(sim_vprob(tdgp%m,tdgp%simreps,num_pr,2),sim_dprob(tdgp%dims(2),tdgp%simreps,num_pr,2),pred_bma(tdgp%H))
allocate(sum_select(num_pr,2),sum_modpost(num_pr,2))
allocate(varprob(m,num_pr,2),dimprob(D2,num_pr,2),modprob(num_pr,2),modselect(num_pr,2))
allocate(pred_ar(H),rmseobs(H))
!
call clear_rmses(rmse)

call init_rmses(H,2,num_pr,2,rmse)


call clear_mats(mat)

allocate(mat%bma(tdgp%H,tdgp%simreps,num_pr,2),mat%top(tdgp%H,tdgp%simreps,num_pr,2),mat%med(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%rw(tdgp%H,tdgp%simreps,num_pr,2),mat%ar(tdgp%H,tdgp%simreps,num_pr,2),mat%rmean(tdgp%H,tdgp%simreps,num_pr,2))
allocate(mat%main(tdgp%H,tdgp%simreps,num_pr,2))

! run

sum_modpost = 0d0
sum_select =0d0
kk = (tdgp%m_syst-1)*2+1
dd = tdgp%dims(2)

allocate(fcast_vec(-10:tdgp%H))

rmseobs = mod%rmse_obs 

do sim_i = 1,s 

Write(*,*)  '  '
write(*,*)  '**********************************  ' 
write(*,*)  'Forecasting origin at  ', sim_i
write(*,*)  '**********************************  '


fcast_vec = 0d0
n_est = tdgp%N+sim_i-1
fcast_vec = datamat((n_est-10):(n_est+tdgp%H),1)

!---- Accunting for losing observations due to no real data for comparsion ----

if (n_est .gt. (data_N -tdgp%H)) then
      Do h_i = (data_n - n_est + 1), H 
             rmseobs(h_i) = rmseobs(h_i) - 1
      end do
end if 


Do mod_i=1,mod%n 
     mod_endo=mod%models(mod_i,tdgp%dims(2)+1)
     mod_k   =mod_endo*mod%p+mod%q
     allocate(pred_temp(tdgp%H,mod_endo))
     call bma_weights(n_est,tdgp%H, tdgp%m, mod%q, mod%p, mod_endo,ms,tdgp%mcburn,tdgp%mcreps,DATAMAT, tdgp%dims(2), mod%Models(mod_i,:), prior, mod%NUM, mod%EVALvec,mod%posteval,method, &
lnmarglik(mod_i),lnpredlik(mod_i),pred_temp,non_deterfcasts)  
     ! first endogeneous variable y is our variable of interest
    ftemp(mod_i,:)=pred_temp(:,1)
    deallocate(pred_temp)
end do


   if (ar_p .ne.  mod%p)  then 
        call   fcast_uar(ar_p,datamat(1:n_est,:),n_est,tdgp%m,mod%q,tdgp%H,prior,tdgp%mcreps,tdgp%mcburn,pred_ar,non_deterfcasts)
    end if 

! For the model Bayes weights

do i=1,num_pr

   mod_post(:,i,1) = lnpredlik
   mod_post(:,i,2) = lnmarglik
        
  do j=1,2
    
    median_mod=1
        
  ! for marginal and predictive likelihood weights
    mod_post(:,i,j)=mod_post(:,i,j) + mod%lnprior(:,i) 
    mod_post(:,i,j)=exp(mod_post(:,i,j)-maxval(mod_post(:,i,j)))
    mod_post(:,i,j)=mod_post(:,i,j)/sum(mod_post(:,i,j))
    call calc_varprobs(mod%n,tdgp%dims(2),tdgp%m,mod_post(:,i,j),mod%models,sim_vprob(:,sim_i,i,j),sim_dprob(:,sim_i,i,j),median_mod)
    call bma_forecasts( mod_post(:,i,j),ftemp,mod%n,tdgp%H,pred_bma)
    call mod_summary_new( mod_post(:,i,j), ftemp, tdgp%dims, mod%models, mod%true, mod%n, mod%no_mods_del, &
					  sum_modpost(i,j), pred_sort, sum_select(i,j))
    call fcast4mc( pred_bma, fcast_vec, pred_sort, ftemp, median_mod, tdgp%H, mat%bma(:,sim_i,i,j), mat%top(:,sim_i,i,j), mat%ar(:,sim_i,i,j), mat%med(:,sim_i,i,j), mat%rw(:,sim_i,i,j), mat%rmean(:,sim_i,i,j), mat%main(:,sim_i,i,j))
    if (ar_p .ne. mod%p) then
      mat%ar(:,sim_i,i,j) = pred_ar
    end if 
  end do
end do

end do

if (present(simvprob)) then
   simvprob = sim_vprob
end if 

if (present(simdprob)) then
   simdprob = sim_dprob
end if 

!*********************************************
! post simulation calculations
!*********************************************
! Variable Probabilities

do i=1,num_pr
  
  do j=1,2
   
    varprob(1:m,i,j) = 0.0d0
    dimprob(1:d2,i,j) = 0.0d0
    do sim_i = 1,tdgp%simreps
	    varprob(1:m,i,j) = varprob(1:m,i,j) + sim_vprob(1:m,sim_i,i,j)
	    dimprob(1:d2,i,j) = dimprob(1:d2,i,j) + sim_dprob(1:d2,sim_i,i,j)
    enddo
    
    
    call calc_mse(tdgp%H,tdgp%simreps,rmseobs,mat%bma(:,:,i,j),mat%Main(:,:,i,j),rmse%bma(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,rmseobs,mat%top(:,:,i,j),mat%Main(:,:,i,j),rmse%top(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,rmseobs,mat%ar(:,:,i,j),mat%Main(:,:,i,j),rmse%ar(:,:,i,j))  
    call calc_mse(tdgp%H,tdgp%simreps,rmseobs,mat%med(:,:,i,j),mat%Main(:,:,i,j),rmse%med(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,rmseobs,mat%rw(:,:,i,j),mat%Main(:,:,i,j),rmse%rw(:,:,i,j))
    call calc_mse(tdgp%H,tdgp%simreps,rmseobs,mat%rmean(:,:,i,j),mat%Main(:,:,i,j),rmse%rmean(:,:,i,j))
    
  end do
end do

 varprob= varprob/tdgp%simreps
 dimprob= dimprob/tdgp%simreps

do i =1,num_pr
  
 vartable(:,(i*2-1)) = varprob(:,i,1)
 vartable(:,(i*2)) = varprob(:,i,2)
 dimtable(:,(i*2-1)) =dimprob(:,i,1)
 dimtable(:,(i*2)) = dimprob(:,i,2)
 
end do

if (present(rmse_obs)) then 
    rmse_obs = rmseobs
end if    
 

deallocate(lnmarglik,lnpredlik,ftemp,mod_post,pred_sort)
deallocate(sim_vprob,sim_dprob,pred_bma,rmseobs)
deallocate(sum_select,sum_modpost)
deallocate(varprob,dimprob,modprob,modselect,fcast_vec)

end subroutine reals_result1








!===================================================================Real data sets==================================================================================

subroutine reals_result2(ar_p,N,H,m,ms,mcburn,mcreps,present_top,dims,datamat,method,prior,mod,varprob,dimprob,mat,modpost,top_mods,non_deterfcasts)


implicit none

integer, intent(in)         :: mcburn,mcreps,dims(:),N,H,m,ms,ar_p ,present_top
real*8, intent(in)          ::datamat(:,:)
type(PriorInfo), intent(in) :: prior
character*20, intent(in)    :: method
type(mods), intent(in)      :: mod
real*8,      intent(out)    :: varprob(:,:),dimprob(:,:),modpost(:,:,:)
type(mats),  intent(out)    :: mat
integer,optional, intent(in)  :: non_deterfcasts 
integer, intent(out)        :: top_mods(:,:,:)

integer :: mod_i, mod_endo,mod_k,num_pr,i,j,median_mod
real*8, allocatable :: ftemp(:,:), pred_temp(:,:),lnmarglik(:), lnpredlik(:),pred_ar(:),mod_post(:,:,:),pred_bma(:),pred_sort(:,:),fcast_vec(:)


num_pr=size(mod%lnprior,2)

allocate(ftemp(mod%n,H),lnmarglik(mod%n),lnpredlik(mod%n))
allocate(pred_ar(H),mod_post(mod%n,num_pr,2))

call clear_mats(mat)

allocate(mat%bma(H,1,num_pr,2),mat%top(H,1,num_pr,2),mat%med(H,1,num_pr,2))
allocate(mat%rw(H,1,num_pr,2),mat%ar(H,1,num_pr,2),mat%rmean(H,1,num_pr,2))
allocate(mat%main(H,1,num_pr,2))

allocate(pred_bma(H),pred_sort(mod%n,H))

allocate(fcast_vec(-10:H))

fcast_vec = 0d0

fcast_vec = datamat(N-10:N+H,1)

Do mod_i=1,mod%n 
     mod_endo=mod%models(mod_i,dims(2)+1)
     mod_k   =mod_endo*mod%p+mod%q
     allocate(pred_temp(H,mod_endo))
     call bma_weights(N,H, m, mod%q, mod%p, mod_endo,ms,mcburn,mcreps,datamat, dims(2), mod%Models(mod_i,:), prior, mod%NUM, mod%EVALvec,mod%posteval,method, &
lnmarglik(mod_i),lnpredlik(mod_i),pred_temp,non_deterfcasts)  
     ! first endogeneous variable y is our variable of interest
    ftemp(mod_i,:)=pred_temp(:,1)
    deallocate(pred_temp)
end do
   
    call write_data(ftemp,"ftemp")
    call wrrrn("models",dble(mod%models))

   if (ar_p .ne.  mod%p)  then 
        call   fcast_uar(ar_p,datamat,N,m,mod%q,H,prior,mcreps,mcburn,pred_ar,non_deterfcasts)
    end if 


do i=1,num_pr

   mod_post(:,i,1) = lnpredlik
   mod_post(:,i,2) = lnmarglik
        
  do j=1,2
  
    median_mod=1
    mod_post(:,i,j)=mod_post(:,i,j) + mod%lnprior(:,i) 
    mod_post(:,i,j)=exp(mod_post(:,i,j)-maxval(mod_post(:,i,j)))
    mod_post(:,i,j)=mod_post(:,i,j)/sum(mod_post(:,i,j))
    call calc_varprobs(mod%n,dims(2),m,mod_post(:,i,j),mod%models,varprob(:,((i-1)*2+j)),dimprob(:,((i-1)*2+j)),median_mod)
    call bma_forecasts( mod_post(:,i,j),ftemp,mod%n,H,pred_bma)
	call mod_summary_real( mod_post(:,i,j), ftemp, dims, mod%models, mod%n, mod%no_mods_del, present_top, &
					  modpost(:,i,j), top_mods(:,i,j), pred_sort)
    call fcast4mc( pred_bma, fcast_vec, pred_sort, ftemp, median_mod, H, mat%bma(:,1,i,j), mat%top(:,1,i,j), mat%ar(:,1,i,j), mat%med(:,1,i,j), mat%rw(:,1,i,j), mat%rmean(:,1,i,j), mat%main(:,1,i,j))
    if (ar_p .ne. mod%p) then
      mat%ar(:,1,i,j) = pred_ar
    end if 

 end do

end do



deallocate(ftemp,lnmarglik,lnpredlik,pred_ar,mod_post,fcast_vec)
deallocate(pred_bma,pred_sort)




end subroutine 



















































! ================================================================
subroutine clear_mats(tmat)

type(mats), intent(inout) :: tmat

if (allocated(tmat%bma)) deallocate(tmat%bma)
if (allocated(tmat%top)) deallocate(tmat%top) 
if (allocated(tmat%med)) deallocate(tmat%med)
if (allocated(tmat%ar)) deallocate(tmat%ar)
if (allocated(tmat%rw)) deallocate(tmat%rw)
if (allocated(tmat%rmean)) deallocate(tmat%rmean)
if (allocated(tmat%main)) deallocate(tmat%main)

end subroutine

! ================================================================
subroutine clear_weights(ot)

type(weights), intent(inout) :: ot

if (allocated(ot%priorval)) deallocate(ot%priorval)
if (allocated(ot%dims)) deallocate(ot%dims) 
if (allocated(ot%var)) deallocate(ot%var)
if (allocated(ot%sumposts)) deallocate(ot%sumposts)
if (allocated(ot%topmod)) deallocate(ot%topmod)
if (allocated(ot%ranksum)) deallocate(ot%ranksum)

end subroutine

! ================================================================
subroutine clear_rmses(otr)

type(rmses), intent(inout) :: otr
if (allocated(otr%bma)) deallocate(otr%bma)
if (allocated(otr%top)) deallocate(otr%top) 
if (allocated(otr%med)) deallocate(otr%med)
if (allocated(otr%ar)) deallocate(otr%ar) 
if (allocated(otr%rw)) deallocate(otr%rw)
if (allocated(otr%rmean)) deallocate(otr%rmean) 

end subroutine



subroutine init_weights(m,t,dims2,mn,top,s,ot)
integer, intent(in) :: m,t,dims2,mn,top,s
  type(weights), intent(inout) :: ot
     allocate(ot%priorval(m,t),ot%dims(dims2,t),ot%var(m,t),ot%sumposts(mn,t),ot%topmod(s,top,t),ot%ranksum(mn,t))
end subroutine


subroutine init_rmses(h,t,t1,t2,otr)
integer, intent(in) :: h,t,t1,t2
  type(rmses), intent(inout) :: otr
     allocate(otr%bma(h,t,t1,t2),otr%top(h,t,t1,t2),otr%med(h,t,t1,t2),otr%ar(h,t,t1,t2),otr%rw(h,t,t1,t2),otr%rmean(h,t,t1,t2))
end subroutine





end module simutils
