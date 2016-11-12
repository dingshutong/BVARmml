module postproc

INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use structures

implicit none 

contains

!=============================================================================
subroutine calc_rmse(H,num_preds,rmse_obs,F,X,rmse)

implicit none
integer, intent(in) :: H, num_preds
integer, dimension(H), intent(in) :: rmse_obs
double precision, dimension(H,num_preds) :: F,X				!F is forecasts, X is data
double precision, dimension(H,2) :: rmse

integer :: i,j,HH

rmse = 0.0d0

HH = H
do i = 1,H
	if(rmse_obs(i) == 0) then
		HH = i-1
		exit
	endif
enddo

do i = 1,HH
	do j = 1,rmse_obs(i)
		rmse(i,1) = rmse(i,1) + (X(i,j)-F(i,j))*(X(i,j)-F(i,j))
	enddo
	rmse(i,1) = rmse(i,1)/dble(rmse_obs(i))
	rmse(i,1) = dsqrt(rmse(i,1))
enddo

do i = 1,HH
	do j = 1,rmse_obs(i)
		rmse(i,2) = rmse(i,2) + ( SQRT((X(i,j)-F(i,j))*(X(i,j)-F(i,j))) - rmse(i,1))**2
	enddo
	rmse(i,2) = sqrt(rmse(i,2)/dble(rmse_obs(i)-1))
enddo

end subroutine

!=============================================================================
subroutine calc_mse(H,num_preds,rmse_obs,F,X,mse)

implicit none
integer, intent(in) :: H, num_preds
integer, dimension(H), intent(in) :: rmse_obs
double precision, dimension(H,num_preds) :: F,X				!F is forecasts, X is data
double precision, dimension(H,2) :: mse

integer :: i,j,HH

mse = 0.0d0

HH = H
do i = 1,H
	if(rmse_obs(i) == 0) then
		HH = i-1
		exit
	endif
enddo

do i = 1,HH
	do j = 1,rmse_obs(i)
		mse(i,1) = mse(i,1) + (X(i,j)-F(i,j))*(X(i,j)-F(i,j))
	enddo
	mse(i,1) = mse(i,1)/dble(rmse_obs(i))
enddo

do i = 1,HH
	do j = 1,rmse_obs(i)
		mse(i,2) = mse(i,2) + ( (X(i,j)-F(i,j))*(X(i,j)-F(i,j)) - mse(i,1))**2
	enddo
	mse(i,2) = sqrt(mse(i,2)/dble(rmse_obs(i)-1))
enddo

end subroutine


!------------------------------------------------------

subroutine calc_varprobs(no_mods_adj,dim,nvar,mod_post,models,var_probs,dim_probs,median_mod)

implicit none 
integer, intent(in) :: no_mods_adj, dim, nvar
double precision, intent(in) :: mod_post(no_mods_adj)
integer, intent(in) ::  models(no_mods_adj,dim+1)
double precision, intent(out) :: var_probs(nvar)
double precision, intent(out) :: dim_probs(dim)

integer, intent(inout) :: median_mod ! set positive to return index of median model in RESULTS

integer i, j, median_vars(nvar)

var_probs = 0.0d0
dim_probs = 0.0d0


! variable and dimension probabilities
do i = 1, no_mods_adj
	! variable probabilities
	do j = 1, nvar
		if ( any( models(i,1:models(i,dim+1)) == j ) ) then
			var_probs(j) = var_probs(j) + mod_post(i)
		endif
	enddo
	! dimension probabilities
	do j = 1,dim
		if(models(i,dim+1) == j) then
			dim_probs(j) = dim_probs(j) + mod_post(i)
		endif
	enddo
enddo



! find median model
if ( median_mod > 0 ) then
	median_mod = 0
	j = 0
	do i = 1, nvar
		if ( var_probs(i) >= 0.5d0 ) then
			j = j+1
			median_vars(j) = i
		endif
	enddo
		
	do i = 1, no_mods_adj
		if ( all( models(i,1:j) == median_vars(1:j) ) ) then
			median_mod = i
			goto 1000
		endif
	enddo

end if

1000 continue

end subroutine

!------------------

subroutine mod_summary( mod_post, res_preds, bma_comb, dims, models, no_mods_adj, no_mods_del, present, &
					  mod_ranksum, sum_modpost, top_mods, pred_sort, class_bma )

implicit none

double precision, intent(in) :: mod_post(:), res_preds(:,:)
integer, intent(in) :: bma_comb(:), dims(:), models(:,:), no_mods_adj, no_mods_del, present
integer, intent(inout) :: mod_ranksum(:)
double precision, intent(inout) :: sum_modpost(:)
integer, intent(out) :: top_mods(:), class_bma
double precision, intent(out) :: pred_sort(:,:)

integer, allocatable :: indx(:), mod_sort(:,:)
integer kk, i

! sort by model prob
ALLOCATE(INDX(1:NO_MODS_adj))
allocate(mod_sort(no_mods_adj,1+dims(2)))

CALL INDEXX(NO_MODS_adj,mod_post,INDX)
			
DO KK = 1,NO_MODS_adj
	PRED_SORT(KK,:) = RES_PREDS(INDX(no_mods_adj-KK+1)+no_mods_del,:)
	mod_sort(kk,:) = models(INDX(no_mods_adj-KK+1)+no_mods_del,:)
ENDDO

mod_ranksum(indx) = mod_ranksum(indx) + (/ no_mods_adj:1:-1 /)
top_mods = indx(no_mods_adj:no_mods_adj-present+1:-1)
sum_modpost = sum_modpost + mod_post

do i = 1,no_mods_adj
	if (all(bma_comb(1:dims(2)) == mod_sort(i,1:dims(2)))) then
		class_bma = i
		exit
	endif
enddo  

deallocate( indx, mod_sort )

end subroutine


!------------------

subroutine mod_summary_real( mod_post, res_preds, dims, models, no_mods_adj, no_mods_del, present, &
					  top_modpost, top_mods, pred_sort)

implicit none

double precision, intent(in) :: mod_post(:), res_preds(:,:)
integer, intent(in) :: dims(:), models(:,:), no_mods_adj, no_mods_del, present
double precision, intent(out) :: top_modpost(:)
integer, intent(out) :: top_mods(:)
double precision, intent(out) :: pred_sort(:,:)

integer, allocatable :: indx(:), mod_sort(:,:),topmod(:)
integer kk, i

! sort by model prob
ALLOCATE(INDX(1:NO_MODS_adj))
allocate(mod_sort(no_mods_adj,1+dims(2)),topmod(present))

CALL INDEXX(NO_MODS_adj,mod_post,INDX)
			
DO KK = 1,NO_MODS_adj
	PRED_SORT(KK,:) = RES_PREDS(INDX(no_mods_adj-KK+1)+no_mods_del,:)
	mod_sort(kk,:) = models(INDX(no_mods_adj-KK+1)+no_mods_del,:)
ENDDO

topmod = indx(no_mods_adj:no_mods_adj-present+1:-1)
top_mods    = topmod
top_modpost = mod_post(topmod)

deallocate( indx, mod_sort,topmod )

end subroutine
!-----------------

subroutine fcast4mc( preds_bma, fcast_vec, pred_sort, res_preds, median_mod, h, &
               bma_mat, top_mat, ar_mat, med_mat, rw_mat, rmean_mat, main_mat )

!	*************** storing the forecast for RMSE ************************

implicit none 

double precision, intent(in) :: preds_bma(:), fcast_vec(-10:), &
								pred_sort(:,:), res_preds(:,:)
integer, intent(in) :: h, median_mod
double precision, intent(out) :: bma_mat(:), top_mat(:), ar_mat(:), med_mat(:), rw_mat(:), rmean_mat(:), main_mat(:)

integer i
double precision rmean

	! BMA forecasts
	BMA_mat(1:H) = preds_bma(1:H)
	! Top model forecasts
	TOP_mat(1:H) = Pred_sort(1,1:H)
	! Autoregressive forecasts
	AR_mat(1:H) = res_preds(1,1:H)
	! Median model forecasts
	MED_mat(1:H) = res_preds(median_mod,1:H)
	! Random Walk Forecasts
	do i = 1,H
		RW_mat(i) = fcast_vec(0)
	enddo
	! Recent Mean Forecasts
	Rmean = sum(fcast_vec(-7:0))/8.0d0
	do i = 1,H
		RMean_mat(i) = Rmean
	enddo
	! Data for forecast period
	Main_mat(1:H) = fcast_vec(1:H)

end subroutine

!========================

Subroutine fcast_uar(ar_p,datamat,n,m,q,H,priorspec,mc_reps,mc_burn,pred_ar,non_deterfcasts)

use VARutils
use priorspost
use bayesutils

IMPlicit none



integer,          intent(in)      :: ar_p,n,m,q,H,mc_reps,mc_burn
type(PriorInfo),  intent(in)      :: priorspec
real*8,           intent(in)      :: datamat(:,:)
real*8,           intent(out)     :: pred_ar(:)
integer,optional, intent(in)      :: non_deterfcasts

type(VARData) tData
type(VARprior) tprior
real*8, allocatable :: bayes_fcasts(:,:)


integer vars(1+q)

vars(1) = 1
vars(2:1+q) = (/ m+1:m+q /)
allocate(bayes_fcasts(h,1))


call make_data( Datamat(1:N,vars), 1, ar_p, q, 'NW', tData )

call OLS( tData )

call make_Littermanprior( tData, priorspec, tPrior)

if ( priorspec%priortype == 'NW' ) then

    call bayes_forecasts_NW(tdata,tprior%NW,h,mc_reps,bayes_fcasts,non_deterfcasts)

else if (priorspec%priortype == 'ND') then

   call bayes_forecasts_NIW(tdata,tprior%ND,h,mc_burn,mc_reps,bayes_fcasts,non_deterfcasts)
     
end if 

pred_ar=bayes_fcasts(:,1)

deallocate(bayes_fcasts)

end subroutine 
!-----------------

subroutine outputx( thread,title, ot,otr,models, ser_name, bma_comb, dims, no_mods_adj, m_syst, mc_rep, present, h, m )

implicit none

type(weights), intent(in)  ::  ot
type(rmsess),  intent(in)   ::  otr

character*(*), intent(in) :: title
integer, intent(in) :: bma_comb(:), dims(:),thread
integer, intent(in) :: models(:,:)

character*(12), intent(in) :: ser_name(:)
integer, intent(in) :: no_mods_adj, m_syst, mc_rep, present, h, m

integer, allocatable :: INDX(:)
integer i, k, j, ii


write(20,*) '  '
write(20,*) '==================================================='
write(20,*) '  '//title
write(20,*) '==================================================='
write(20,*) 'Variable  Prior Probability  Posterior Probability'
do i = 1,m	
	write(20,201) ser_name(i), ot%priorval(i,thread),ot%var(i,thread)
enddo
201 format(tr1, A12,tr2,2(F8.5,tr1))

write(20,*) '==================================================='
WRITE(20,*) 'Dim		Post Prob for Dim'
do i = dims(1),dims(2)
	write(20,202) i,ot%dims(i,thread)
enddo	
202 format(tr1, I3, tr2, F8.5)

ALLOCATE(INDX(1:NO_MODS_adj))
CALL INDEXX(NO_MODS_adj,ot%sumposts(:,thread),INDX)
write(20,*) '==================================================='
write(20,*) 'Classificatíon of Models'
write(20,*) '==================================================='
write(20,*) 'DGP is ', ( ser_name(bma_comb(i)), i=1,m_syst)
write(20,*) ' '


k = dims(2) 
write(20,203,advance='no') ( ' ', j = 1, k )
write(20,*) '#Top	 #Top5	#Top10		Mean rank    av.prop'
do i = 1, no_mods_adj
	k = models(i,dims(2)+1)
	write(20,203,advance='no') ( ser_name(models(i,j)), j=1,models(i,dims(2)+1) ) 
	if ( k < dims(2) ) then
		k = dims(2) - k
		write(20,203,advance='no') ( ' ', j = 1, k )
	endif
	write(20,204) dble(count( ot%topmod(:,1,thread) == i ))/dble(mc_rep), dble(count( ot%topmod(:,1:5,thread) == i ))/dble(mc_rep),&
				  dble(count( ot%topmod(:,1:present,thread) == i ))/dble(mc_rep), &
	  			  dble(ot%ranksum(i,thread))/dble(mc_rep), ot%sumposts(i,thread)/dble(mc_rep)
enddo
write(20,*) '==================================================='

write(20,*) '==================================================='
write(20,*) 'Classificatíon of Models, sorted'
write(20,*) '==================================================='
write(20,*) 'DGP is ', ( ser_name(bma_comb(i)), i=1,m_syst)
write(20,*) ' '
k = dims(2) 
write(20,203,advance='no') ( ' ', j = 1, k )
write(20,*) '#Top	 #Top5	#Top10		Mean rank    av.prop'
do ii = 1, no_mods_adj
	i = indx(no_mods_adj-ii+1)
	k = models(i,dims(2)+1)
	write(20,203,advance='no') ( ser_name(models(i,j)), j=1,models(i,dims(2)+1) ) 
	if ( k < dims(2) ) then
		k = dims(2) - k
		write(20,203,advance='no') ( ' ', j = 1, k )
	endif
	write(20,204) dble(count( ot%topmod(:,1,thread) == i ))/dble(mc_rep), dble(count( ot%topmod(:,1:5,thread) == i ))/dble(mc_rep),&
				  dble(count( ot%topmod(:,1:present,thread) == i ))/dble(mc_rep), &
	  			  dble(ot%ranksum(i,thread))/dble(mc_rep), ot%sumposts(i,thread)/dble(mc_rep)
enddo
write(20,*) '==================================================='


203 format( <k>(tr1,A12) )
204 format( tr1, 3(f5.3,tr1), f10.3, tr2, f5.3 )

write(20,*) ' ' 
write(20,*) '           RMSE           '
write(20,*) ' Hor      BMA      TOP    MEDIAN     AR       RW     Rec.Mean'
do i = 1,H
	write(20,205) i,otr%bma(i,thread), otr%top(i,thread), otr%med(i,thread),otr%ar(i,thread), otr%rw(i,thread), otr%rmean(i,thread)
enddo
write(20,*) '==================================================='
write(20,*) ' ' 
205 format(tr1,I3,tr2,6(F8.3,tr1)) 
	   


!WRITE(20,*) '========================================='
!WRITE(20,*) ' FORECASTS START AT T=',PRED_DATE
!WRITE(20,*) ' RANK BASED ON:'
!WRITE(20,114) ' Evaluation Horizon',EVALH
!WRITE(20,*) '========================================='
!WRITE(20,*) ''
!write(20,*) '=============================================================='
!WRITE(20,*) 'BMA Forecasts 1:st differences or levels:    ',ser_name(1),endoexo(1)
!write(20,*) '=============================================================='
!WRITE(20,*) ''
!WRITE(20,*) '                     Credible Regions'
!WRITE(20,110) 'Date medel   Median         ',cred_pred(1),'          ', cred_pred(2)

!	do i = 1, h
!		write(20,133) pred_date+i,preds_bma(i,1),median(i), &
!				((lower(j,i), upper(j,i)), j=1,size(cred_pred))
!	enddo

!IF (LEV(1:3) == 'DIF') THEN
!	write(20,*) '=============================================================='
!	WRITE(20,*) 'BMA Forecasts 4:th differences:              ',ser_name(1),endoexo(1)
!	write(20,*) '=============================================================='
!	WRITE(20,*) ''
!	WRITE(20,*) '                     Credible Regions'
!	WRITE(20,110) 'Date  medel Median         ',cred_pred(1),'          ', cred_pred(2)

!	do i = 1, h
!		write(20,133) pred_date+i, preds4_bma(i,1),median4(i), &
!				((lower4(j,i), upper4(j,i)), j=1,size(cred_pred))
!	enddo

!ENDIF

!WRITE(20,*) '================================================'

!if (present .gt. no_mods_adj) then
!	present = no_mods_adj
!endif

!DO I = 1,present
!	write(20,108) 'Rank: ',I,' Model probability: ',mod_post(indx(no_mods_adj-i+1))
!	write(20,*)	'            ',ser_name(1), &
!							ser_name(models(indx(no_mods_adj-i+1),1:models(indx(no_mods_adj-i+1),dims(2)+1)-2))
!	write(20,*) ' '
!ENDDO
!WRITE(20,*) '================================================'
!write(20,*) ' '
!write(20,*) ' '


!***************

deallocate( indx )

108 format(tr1,A6,I3,TR2,A20,F4.3)
110 format(tr1,A20,F4.2,A10, F4.2)
111 format(tr1,i3,tr2,F4.3,tr2,4(I3,tr1))
!112 FORMAT(TR1,I2,TR3,5(F15.3,TR3))
113 FORMAT(TR1,I2,TR3,4(F15.3,TR1))
114 Format(TR1,A19,TR2,I3)
133 FORMAT(TR1,I3,TR2,6(F7.3,TR1))

end subroutine

!=======================================================================================================

function comb(n,k)

!This function returns the value of the combination operator
!
!Written by M.K. Andersson, Sveriges Riksbank, 2004-09-22



integer, intent(in) :: n,k
integer :: comb

integer :: j

double precision :: c,q1,q2,q3


q1 = 1
q2 = 1
q3 = 1

do j=2,n
	q1=q1*j
enddo

do j = 2,k
	q2=q2*j
enddo

do j = 2,n-k
	q3=q3*j
enddo

c = q1/(q2*q3)

comb = int(c)

end function

!*****************************************************************

SUBROUTINE indexx(n,arr,indx)

! This subroutine is taken from 'Numerical Recipes in Fortran'
INTEGER n,indx(n),M,NSTACK
DOUBLE PRECISION arr(n)
PARAMETER (M=7,NSTACK=50)
!Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
!is in ascending order for j = 1;2; : : :;N. The input quantities n and arr are not changed.

INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
DOUBLE PRECISION a
	do j=1,n
		indx(j)=j
	enddo
	jstack=0
	l=1
	ir=n
1	if(ir-l.lt.M)then
		do j=l+1,ir
			indxt=indx(j)
			a=arr(indxt)
			do i=j-1,l,-1
				if(arr(indx(i)).le.a)goto 2
				indx(i+1)=indx(i)
			enddo 
			i=l-1
2			indx(i+1)=indxt
		enddo 
		if(jstack.eq.0)return
		ir=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+ir)/2
		itemp=indx(k)
		indx(k)=indx(l+1)
		indx(l+1)=itemp
		if(arr(indx(l)).gt.arr(indx(ir)))then
			itemp=indx(l)
			indx(l)=indx(ir)
			indx(ir)=itemp
		endif
		if(arr(indx(l+1)).gt.arr(indx(ir)))then
			itemp=indx(l+1)
			indx(l+1)=indx(ir)
			indx(ir)=itemp
		endif
		if(arr(indx(l)).gt.arr(indx(l+1)))then
			itemp=indx(l)
			indx(l)=indx(l+1)
			indx(l+1)=itemp
		endif
		i=l+1
		j=ir
		indxt=indx(l+1)
		a=arr(indxt)
3		continue
			i=i+1
		if(arr(indx(i)).lt.a)goto 3
4		continue
			j=j-1
		if(arr(indx(j)).gt.a)goto 4
		if(j.lt.i)goto 5
		itemp=indx(i)
		indx(i)=indx(j)
		indx(j)=itemp
		goto 3
5		indx(l+1)=indx(j)
		indx(j)=indxt
		jstack=jstack+2
		if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
		if(ir-i+1.ge.j-l)then
			istack(jstack)=ir
			istack(jstack-1)=i
			ir=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i	
		endif
	endif
	goto 1
END subroutine


!------------------

subroutine mod_summary_new( mod_post, res_preds, dims, models, mod_true, no_mods_adj, no_mods_del, &
					  sum_modpost, pred_sort, sum_modselect)

implicit none

double precision, intent(in) :: mod_post(:), res_preds(:,:)
integer, intent(in) :: dims(:), models(:,:), no_mods_adj, no_mods_del, mod_true

double precision, intent(inout) :: sum_modpost ,sum_modselect
double precision, intent(out)   :: pred_sort(:,:)


real*8 :: mod_select

integer, allocatable :: indx(:)
integer kk, i

! sort by model prob
ALLOCATE(INDX(1:NO_MODS_adj))


CALL INDEXX(NO_MODS_adj,mod_post,INDX)
			
DO KK = 1,NO_MODS_adj
	PRED_SORT(KK,:) = RES_PREDS(INDX(no_mods_adj-KK+1)+no_mods_del,:)
ENDDO

if (INDX(no_mods_adj) .eq. mod_true) then
         sum_modselect = sum_modselect + 1d0
end if 

sum_modpost = sum_modpost + mod_post(mod_true)

!my = models(INDX(no_mods_adj)+no_mods_del,:)

!if (my(dims(2)+1) .ne. m_syst) then
!    mod_select = 0d0
!else
!    mod_select = 1d0 
!    do i =1,m_syst
!       if (my(i) .ne. i) then
!           my(i)=0
!       else 
!           my(i)=1
!       end if
!       mod_select=mod_select*dble(my(i))  
!    end do
!end if





deallocate( indx)

end subroutine



!-----------------------------------------------------------------

FUNCTION select(k,n,arr)
      INTEGER :: k,n
      DOUBLE PRECISION :: select,arr(n)
      INTEGER i,ir,j,l,mid
      DOUBLE PRECISION :: a,temp
      l=1
      ir=n
1     if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
        endif
        select=arr(k)
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)

          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i

      endif
      goto 1
END function



end module