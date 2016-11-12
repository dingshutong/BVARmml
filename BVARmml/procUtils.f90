module procutils


INCLUDE 'link_fnl_static.h'!DEC$ OBJCOMMENT LIB:'libiomp5md.lib'
use imsl_libraries
use structures

implicit none 

contains

!-----------------

subroutine table_var( vartable, ser_name, m_syst)

implicit none

integer, intent(in) :: m_syst
real*8,  intent(in) :: vartable(:,:)
character*(12), intent(in) :: ser_name(:)

integer i, tabrow

tabrow=size(vartable,1)

write(20,*) '  '
write(20,*) '======================================================'
write(20,*) 'Predictive weights and marginal weights for variables'
write(20,*) '======================================================'
write(20,*) "Predictive                      ",  "Marginal"

if (m_syst .eq. 2) then 

write(20,*) ser_name(2), "max(x_i)  ", trim(adjustl(ser_name(2)))//"/max(x_i)  ", ser_name(2), "max(x_i)   ", trim(adjustl(ser_name(2)))//"/max(x_i)"

do i = 1,tabrow	
	   write(20,501) vartable(i,:)        
enddo
else if (m_syst .eq. 3)	then

write(20,*) ser_name(2:m_syst), "max(x_i)  ", trim(adjustl(ser_name(2)))//"/max(x_i)  ", trim(adjustl(ser_name(3)))//"/max(x_i)  ", &
ser_name(2:m_syst), "max(x_i)  ", trim(adjustl(ser_name(2)))//"/max(x_i)  ", trim(adjustl(ser_name(3)))//"/max(x_i)  "

do i = 1,tabrow	
	   write(20,502) vartable(i,:)        
enddo
end if 
501 format(tr1,6(F8.5,tr1))
502 format(tr1,10(F8.5,tr1))

write(20,*) '=================================================================='


end subroutine



!-----------------

subroutine table_dim( dimtable)

implicit none

real*8,  intent(in) :: dimtable(:,:)

integer i, tabrow

tabrow=size(dimtable,1)

write(20,*) '  '
write(20,*) '======================================================'
write(20,*) 'Predictive weights and marginal weights for dimensions'
write(20,*) '======================================================'
write(20,*) "Predictive              ",  "Marginal"
write(20,*) "1   ", "2   ", "3   ", "4   ", "1   ", "2   ", "3   ", "4   "

do i = 1,tabrow	
	   write(20,503) dimtable(i,:)        
enddo

503 format(tr1,8(F8.5,tr1))
write(20,*) '==================================================================='


end subroutine


!-----------------

subroutine table_mod( modtable)

implicit none

real*8,  intent(in) :: modtable(:,:)

integer i, tabrow

tabrow=size(modtable,1)

write(20,*) '  '
write(20,*) '==================================================================== '
write(20,*) 'Averaged posterior probability and propotion selected for true model'
write(20,*) '====================================================================='
write(20,*) "Predictive          ",  "Marginal"
write(20,*) "Weighted  ", "Selected  ", "Weighted  ", "Selected  "

do i = 1,tabrow	
	   write(20,504) modtable(i,:)        
enddo

504 format(tr1,8(F8.3,tr1))

write(20,*) '==================================================='


end subroutine


! =====

SUBROUTINE table_rmse(title,p,otr,modprior)  

integer, intent(in)      :: p
type(rmses), intent(in)  :: otr
real*8,      intent(in)  :: modprior(:)
character*(*), intent(in):: title     

integer ::i,j,k,l,H,nprior

H=size(otr%bma,1)
nprior=size(modprior)

write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) '  '//title, p
write(20,*) '====================================================================='


write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) 'MSE and SE(MSE) using predictive weights and marginal weights'
write(20,*) '====================================================================='

write(20,*) ' ' 

do j = 1,nprior 
write(20,*) '===================================================================================================================================='
write(20,*) ' ' 
write(20,*) '    FOR the variable inclusion probability = ', modprior(j)
    do l =1,2    

      write(20,*) ' ' 

       IF (l .eq. 1) then
          write(20,*) '           MSE           '
       else if (l .eq. 2) then
          write(20,*) '         se(MSE)           '
       end if

       write(20,*) ' Hor      BMA(p)     BMA(m)     TOP(p)      TOP(m)      MED(p)     MED(m)     RW(p)      RW(m)      RMEAN(p)     RMEAN(m)     AR(p)      AR(m)'
       do i = 1,H
	      write(20,505) i,otr%bma(i,l,j,1), otr%bma(i,l,j,2), otr%top(i,l,j,1), otr%top(i,l,j,2), otr%med(i,l,j,1), otr%med(i,l,j,2), &
	        otr%rw(i,l,j,1), otr%rw(i,l,j,2), otr%rmean(i,l,j,1), otr%rmean(i,l,j,2),otr%ar(i,l,j,1), otr%ar(i,l,j,2)
       enddo
        write(20,*) '===================================================================================================================================='
        write(20,*) ' ' 
        505 format(tr1,I3,tr2,12(F8.3,tr1)) 
     end do
     
end do

end subroutine 



! =====

SUBROUTINE table_rmse1(title,p,otr,modprior)  

integer, intent(in)      :: p
type(rmses), intent(in)  :: otr
real*8,      intent(in)  :: modprior(:)
character*(*), intent(in):: title     

integer ::i,j,k,l,H,nprior

H=size(otr%bma,1)
nprior=size(modprior)

write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) '  '//title, p
write(20,*) '====================================================================='


write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) 'RMSE and SE(RMSE) using predictive weights and marginal weights'
write(20,*) '====================================================================='

write(20,*) ' ' 

do j = 1,nprior 
write(20,*) '===================================================================================================================================='
write(20,*) ' ' 
write(20,*) '    FOR the variable inclusion probability = ', modprior(j)
    do l =1,2    

      write(20,*) ' ' 

       IF (l .eq. 1) then
          write(20,*) '           MSE           '
       else if (l .eq. 2) then
          write(20,*) '         se(MSE)           '
       end if

       write(20,*) ' Hor      BMA(p)     BMA(m)     TOP(p)      TOP(m)    AR(p)      AR(m)'
       do i = 1,H
	      write(20,507) i,otr%bma(i,l,j,1), otr%bma(i,l,j,2), otr%top(i,l,j,1), otr%top(i,l,j,2), &
	       otr%ar(i,l,j,1), otr%ar(i,l,j,2)
       enddo
        write(20,*) '===================================================================================================================================='
        write(20,*) ' ' 
        507 format(tr1,I3,tr2,6(F8.3,tr1)) 
     end do
     
end do

end subroutine 



!-----------------

subroutine table_var_real( vartable, ser_name)

implicit none

real*8,  intent(in) :: vartable(:,:)
character*(12), intent(in) :: ser_name(:)

integer i, tabrow

tabrow=size(vartable,1)

write(20,*) '  '
write(20,*) '======================================================'
write(20,*) 'Predictive weights and marginal weights for variables'
write(20,*) '======================================================'
write(20,*) "Varprior=0.2                      ",  "Varprior=0.5"



write(20,*) "ser_name", "Predictive    ", "Marginal  ", "Predictive   ", "Marginal"

do i = 1,tabrow	
	   write(20,801) ser_name(i), vartable(i,:)        
enddo

801 format(tr1, A12,tr2,4(F8.5,tr1))

write(20,*) '=================================================================='


end subroutine




!-----------------

subroutine table_var_real_seq( vartable, ser_name,modprior)

implicit none

real*8,  intent(in) :: vartable(:,:,:,:)
character*(12), intent(in) :: ser_name(:)
real*8,      intent(in)    :: modprior(:)

integer i,j,k,l, tabrow,tabcol,nprior

tabrow=size(vartable,2)
tabcol=size(vartable,1)
nprior=size(modprior)

do j = 1,nprior
 do   k =1,2

write(20,*) '  '
if (k .eq. 1) then
write(20,*) '======================================================'
write(20,*) 'Predictive weights for variables in each origin'
write(20,*) '======================================================'
 else if (k .eq. 2) then
write(20,*) '======================================================'
write(20,*) 'Marginal weights for variables in each origin'
write(20,*) '======================================================'
end if

write(20,*) "With Variable prior inclusion probability =  ", modprior(j)  

write(20,802) "Origin" , ( ser_name(l), l=1,tabcol ) 

do i = 1,tabrow	
	   write(20,803) i , vartable(:,i,j,k)        
enddo

write(20,*) '=================================================================='

end do
end do

803 format(tr1, i3,tr2,<tabcol>(F8.5,tr1))
802 format( tr1,A12, <tabcol>(tr1,A12))

end subroutine





!-----------------

subroutine table_dim_real( dimtable)

implicit none

real*8,  intent(in) :: dimtable(:,:)

integer i, tabrow

tabrow=size(dimtable,1)

write(20,*) '  '
write(20,*) '======================================================'
write(20,*) 'Predictive weights and marginal weights for dimensions'
write(20,*) '======================================================'
write(20,*) "Varprior=0.2                      ",  "Varprior=0.5"



write(20,*) "Dim  ", "Predictive    ", "Marginal  ", "Predictive   ", "Marginal"

do i = 1,tabrow	
	   write(20,902) i, dimtable(i,:)        
enddo

902 format(tr1, i3,tr2,4(F8.5,tr1))

write(20,*) '=================================================================='


end subroutine



!-----------------

subroutine table_dim_real_seq( dimtable, ser_name,modprior)

implicit none

real*8,  intent(in) :: dimtable(:,:,:,:)
character*(12), intent(in) :: ser_name(:)
real*8,      intent(in)    :: modprior(:)

integer i,j,k,l, tabrow,tabcol,nprior

tabrow=size(dimtable,2)
tabcol=size(dimtable,1)
nprior=size(modprior)

do j = 1,nprior
 do   k =1,2

write(20,*) '  '
if (k .eq. 1) then
write(20,*) '======================================================'
write(20,*) 'Predictive weights for variables in each origin'
write(20,*) '======================================================'
 else if (k .eq. 2) then
write(20,*) '======================================================'
write(20,*) 'Marginal weights for variables in each origin'
write(20,*) '======================================================'
end if
write(20,*) "With Variable prior inclusion probability =  ", modprior(j)  

write(20,*) "Origin   " , "1     ", "2     ", "3     ", "4    " 

do i = 1,tabrow	
	   write(20,903) i , dimtable(:,i,j,k)        
enddo

write(20,*) '=================================================================='

end do
end do

903 format(tr1, i3,tr2,4(F8.5,tr1))

end subroutine



! =====

SUBROUTINE table_real_rmse(title,otr,modprior,rmse_obs)  

type(rmses), intent(in)  :: otr
real*8,      intent(in)  :: modprior(:)
character*(*), intent(in):: title  
integer, intent(in)      :: rmse_obs(:)  

integer ::i,j,k,l,H,nprior

H=size(otr%bma,1)
nprior=size(modprior)

write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) '  '//title
write(20,*) '====================================================================='


write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) 'MSE and SE(MSE) using predictive weights and marginal weights'
write(20,*) '====================================================================='

write(20,*) ' ' 

do j = 1,nprior 
write(20,*) '===================================================================================================================================='
write(20,*) ' ' 
write(20,*) '    FOR the variable inclusion probability = ', modprior(j)
    do l =1,2    

      write(20,*) ' ' 

       IF (l .eq. 1) then
          write(20,*) '           RMSE           '
       else if (l .eq. 2) then
          write(20,*) '         se(RMSE)           '
       end if

       write(20,*) ' Hor      NoFcast      BMA(p)     BMA(m)     TOP(p)      TOP(m)    MED(p)      MED(m)        AR(p)      AR(m)'
       do i = 1,H
	      write(20,605) i, rmse_obs(i), otr%bma(i,l,j,1), otr%bma(i,l,j,2), otr%top(i,l,j,1), otr%top(i,l,j,2), otr%MED(i,l,j,1), otr%med(i,l,j,2), &
	       otr%ar(i,l,j,1), otr%ar(i,l,j,2)
       enddo
        write(20,*) '===================================================================================================================================='
        write(20,*) ' ' 
        605 format(tr1,I3,tr2,I3,tr3,8(F8.3,tr1)) 
     end do
     
end do

end subroutine 



! =====

SUBROUTINE table_real_mat1(title,otr,modprior)  

type(mats), intent(in)   :: otr
real*8,      intent(in)  :: modprior(:)
character*(*), intent(in):: title     

integer ::i,j,k,l,H,nprior,seq


H=size(otr%bma,1)
seq=size(otr%bma,2)
nprior=size(modprior)


write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) '  '//title
write(20,*) '====================================================================='


write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) 'Sequential forecasts using predictive weights and marginal weights'
write(20,*) '====================================================================='

write(20,*) ' ' 

do j = 1,nprior 
write(20,*) '===================================================================================================================================='
write(20,*) ' ' 
write(20,*) '    FOR the variable inclusion probability = ', modprior(j)
    do l =1,3   

      write(20,*) ' ' 
             
       select case (l)
    
       case (1)
          write(20,*) '           BMA(p)                                                              BMA(m)                                           '
          write(20,*) ' SEQ      1     2     3      4       5      6     7      8       1     2     3      4       5      6     7      8'
          do i = 1,seq
	          write(20,607) i,otr%bma(:,i,j,1),otr%bma(:,i,j,2)
          enddo
       case (2)
          write(20,*) '           TOP(p)                                                              TOP(m)                                           '
          write(20,*) ' SEQ      1     2     3      4       5      6     7      8       1     2     3      4       5      6     7      8'
          do i = 1,seq
	          write(20,607) i,otr%top(:,i,j,1),otr%top(:,i,j,2)
          enddo
       case (3)
          write(20,*) '           MED(p)                                                              MED(m)                                           '
          write(20,*) ' SEQ      1     2     3      4       5      6     7      8       1     2     3      4       5      6     7      8'
          do i = 1,seq
	          write(20,607) i,otr%MED(:,i,j,1),otr%MED(:,i,j,2)
          enddo
       
       end select
       
        write(20,*) '===================================================================================================================================='
        write(20,*) ' ' 
        607 format(tr1,I3,tr2,16(F8.3,tr1)) 
     end do
     
end do

end subroutine 


! =====

SUBROUTINE table_real_mat2(title,otr,modprior)  

type(mats), intent(in)  :: otr
real*8,      intent(in)  :: modprior(:)
character*(*), intent(in):: title     

integer ::i,j,k,l,H,nprior

H=size(otr%bma,1)
nprior=size(modprior)

write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) '  '//title
write(20,*) '====================================================================='


write(20,*) ' ' 
write(20,*) '==================================================================== '
write(20,*) 'Terminal forecasts  using predictive weights and marginal weights'
write(20,*) '====================================================================='

write(20,*) ' ' 

do j = 1,nprior 
write(20,*) '===================================================================================================================================='
write(20,*) ' ' 
write(20,*) '    FOR the variable inclusion probability = ', modprior(j)
   
      write(20,*) ' ' 

       write(20,*) ' Hor      BMA(p)     BMA(m)     TOP(p)      TOP(m)    MED(p)      MED(m)        AR(p)      AR(m)'
       do i = 1,H
	      write(20,608) i,otr%bma(i,1,j,1), otr%bma(i,1,j,2), otr%top(i,1,j,1), otr%top(i,1,j,2), otr%MED(i,1,j,1), otr%med(i,1,j,2), &
	       otr%ar(i,1,j,1), otr%ar(i,1,j,2)
       enddo
        write(20,*) '===================================================================================================================================='
        write(20,*) ' ' 
        608 format(tr1,I3,tr2,8(F8.3,tr1)) 
     
end do

end subroutine 


! ==============

subroutine present_topmod(models,modprior,ser_name,present_top,dimmax,top_modpost,top_mods)

integer, intent(in)   :: models(:,:),present_top, dimmax,top_mods(:,:,:)
real*8,  intent(in)   :: top_modpost(:,:,:)
real*8,      intent(in)  :: modprior(:)
character*(12), intent(in) :: ser_name(:)
integer :: mt,nprior,ii,i,k,j,r1,r2

mt=present_top
nprior = size(modprior)

write(20,*) ' ' 
write(20,*) '==================================================='
write(20,*) 'Classificatíon of Models, sorted top ', mt
write(20,*) '==================================================='
write(20,*) ' '
k = dimmax 

do r1=1,nprior
write(20,*) '============================================================='
write(20,*) ' ' 
write(20,*) '    FOR the variable inclusion probability = ', modprior(r1)

   do r2 = 1,2
write(20,*) '============================================================='
write(20,*) ' ' 
if (r2 .eq. 1) then 
    write(20,*) '    Using predictive likelihood weights        '
else if  (r2 .eq. 2) then 
    write(20,*) '    Using marginal likelihood weights        '
end if
write(20,*) ' ' 

write(20,703,advance='no') ( ' ', j = 1, k )
write(20,*) '   mod.weights'
do ii = 1,mt 
	i = top_mods(ii,r1,r2)
	k = models(i,dimmax+1)
	write(20,703,advance='no') ( ser_name(models(i,j)), j=1,models(i,dimmax+1) ) 
	if ( k < dimmax ) then
		k = dimmax - k
		write(20,703,advance='no') ( ' ', j = 1, k )
	endif
	write(20,704) top_modpost(ii,r1,r2)
enddo
write(20,*) '==================================================='
 
end do
end do

703 format( <k>(tr1,A12) )
704 format( tr1, 1(f5.3,tr1))


end subroutine 


end module procutils