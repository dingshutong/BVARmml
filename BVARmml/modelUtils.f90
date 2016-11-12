module modelutils


INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use matrix

implicit none 

contains





!====================================Returns the total number of models and the model combinations=====================================!
subroutine mod_gen(m,s,dims,mod_nums,mod_combs)

! INTENT(IN)
! m      ---------- total number of variables (include study variables)
! s      ---------- total number of study variables (default is 1)
! dims(1)---------- minimum number of variables in the model
! dims(2)---------- maximum number of variables in the model
! INTENT(OUT)
!mod_nums     ---------- The total number of possiable combinations of models 
!mod_combs    ---------- The all possiable combinations of models 

integer, intent(in)  :: dims(2),m,s
integer, intent(out) :: mod_nums,mod_combs(:,:)
integer mn

call mod_num(m,s,dims,mn)
call mod_search(m,s,mn,dims(2),mod_combs)
mod_nums=mn


end subroutine


!====================================TOtal number of possiable combinations of models=====================================!

subroutine mod_num(m,s,dims,nums,no_mods_dels)

! INTENT(IN)
! m      ---------- total number of variables (include study variables)
! s      ---------- total number of study variables (default is 1)
! dims(1)---------- minimum number of variables in the model
! dims(2)---------- maximum number of variables in the model
! INTENT(OUT)
!nums    ---------- The total number of possiable combinations of models 

implicit none

integer, intent(in)  :: dims(2),m,s
integer, intent(out) :: nums
integer, optional, intent(out) :: no_mods_dels
integer i,num_vec(dims(2)),no_mods,no_mods_del,no_mods_adj

if (dims(1)<s) then
  write (*,*) 'minimum dim should not smaller than the num of study variables'
  pause
end if




NO_MODS = 0
NUM_VEC = 0

if ( dims(1) == s ) then
	num_vec(1) = 1
	no_mods = 1
endif



DO I = 2, (DIMS(2)-s+1)
	num_vec(i) = combs(m-s,I-1)
	NO_MODS = NO_MODS + NUM_VEC(I)
ENDDO

if (dims(1) .gt. dims(2)) then
	write(*,*) 'dim_min cannot be greater than dim_max'
else 
    if (DIMS(1) .gt. s) then
	   no_mods_del = sum(num_vec(1:dims(1)-s))
	   no_mods_adj = no_mods - no_mods_del
     else
	   no_mods_del = 0
	   no_mods_adj = no_mods
    endif
end if
nums = no_mods_adj
if (present(no_mods_dels)) then
 no_mods_dels = no_mods_del
end if  
end subroutine 


!====================================Finding Combinations=====================================!


subroutine mod_search(m,s,num,dim,mod_combs)

implicit none
integer, intent(in) :: m,num,dim,s
integer, dimension(num,dim+1), intent(out) :: mod_combs

integer p_start, p_end, last, mods2add, csize, i, allvars(m), mm

mod_combs = 0
allvars = (/ (i,i=1,m) /)

!univariat

do i=1,num 
   mod_combs(i,1:s) = allvars(1:s) 
end do

mod_combs(1,dim+1) = s

p_start = 1
p_end = 1
last = 1
mods2add = m

do csize = s+1, dim
	mm=m-csize+2
	do i = p_start, p_end
		mods2add = m - mod_combs(i,csize-1)
		call mod_search2( mods2add, dim, csize-1, mod_combs(i,1:csize-1), mod_combs(last+1:last+mods2add,:) )
		mod_combs(last+1:last+mods2add,dim+1) = csize
		last = last+mods2add
	enddo
	p_start = p_end + 1
	p_end = last

enddo

!open(75,file='test.ut')
!do i = 1,num
!	write(75,*) mod_combs(i,1:dim)
!enddo

end subroutine

!--------------------------------------------------

subroutine mod_search2( mods2add, dim, psize, base, models  )

integer, intent(in) :: mods2add, dim, psize, base(psize)
integer, intent(out) :: models(mods2add,dim)

integer i

do i = 1, mods2add
	models(i,1:psize) = base
!	models(i,psize+1) = vars2add(i)
	models(i,psize+1) = base(psize) + i
enddo

end subroutine


!===============================================================

subroutine mod_prior(m,prior_val,mod_endo,model,mod_lnprior) 

real*8, intent(in)  :: prior_val(:)
integer, intent(in) :: model(:),mod_endo,m
real*8, intent(out) :: mod_lnprior

integer :: ii, jj, ll,varoutmod(m-mod_endo),indic
real*8  :: priorin,priorout

priorin = 1.0d0

do ii = 1, mod_endo
	priorin = priorin*prior_val(model(ii))
enddo

priorout = 1.0d0
ll =1
do ii = 1,m
	indic = 0
	do jj = 1,mod_endo
		if(ii .eq. model(jj)) then
			indic = 1
		endif
	enddo
	if (indic==0) then
		varoutmod(ll) = ii
		ll=ll+1
	endif
enddo
do ii = 1,(m-mod_endo)
	priorout = priorout*(1-prior_val(varoutmod(ii)))
enddo


mod_lnprior = log( priorin*priorout )

end subroutine




! =============================================================
SUBROUTINE mod_priorsize(mod_lnprior,MOD_N,dim,models,priorsize)

integer, intent(in):: mod_n,dim,models(:,:)
real*8, intent(in) :: mod_lnprior(:)
real*8, intent(out):: priorsize

integer i,j
real*8, allocatable :: dim_probs(:)


allocate(dim_probs(dim))

dim_probs =0d0

do i = 1, mod_N
	! dimension probabilities
	do j = 1,dim
		if(models(i,dim+1) == j) then
			dim_probs(j) = dim_probs(j) + exp(mod_lnprior(i))
		endif
	enddo
enddo

do i =1,dim
   dim_probs(i) = dim_probs(i)*dble(i)
end do

priorsize = sum(dim_probs)

end subroutine 

end module modelutils