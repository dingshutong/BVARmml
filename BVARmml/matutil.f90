module matutil

implicit none 

! Usage: compact_to_square( vector, squaremat, upper, triang )
! vector is vech( squaremat ) returns squaremat
! upper = .true. => vech is upper triangle of squaremat, otherwise lower triangle
! triang = .true. => squaremat is triangular, otherwise symmetric
interface compact_to_square
    module procedure compact_to_square_2dim, compact_to_square_3dim
end interface

! Usage: save_vech( tosave, savein, row, startcol, upper )
interface save_vech
    module procedure save_vech_2dim, save_vech_3dim
end interface

!usage: KronAB = kron(A,B)
!Kronecker routines written by M. Salabasis
	INTERFACE KRON
		MODULE PROCEDURE KRON_MM,KRON_VM,KRON_MV,KRON_SM,KRON_MS,KRON_VS,KRON_SV,KRON_SS,KRON_VV
	END INTERFACE

contains

subroutine trimult( A, B, C, Aup, Bup, Atran, Btran )

double precision, intent(in) :: A(:,:), B(:,:)
double precision, intent(out) :: C(:,:)
integer, intent(in), optional :: Aup, Bup, Atran, Btran

! Multiply triangular matrices
! Default is that A and B matrices are lower triangular and not transposed
! C can NOT share memory locations with A and B
! set optional arguments Aup, Bup, Atran and Btran to a non-zero value to indicate that
! matrices are upper triangular or should be transposed

integer iAup, iBup, iAtran, iBtran, i, j, k, r
double precision X

! defaults
iAup = 0; iBup = 0; iAtran = 0; iBtran = 0

! optional arguments

if ( present(Aup) ) then
  iAup = Aup
endif
if ( present(Bup) ) then
  iBup = Bup
endif
if ( present(Atran) ) then
  iAtran = Atran
endif
if ( present(Btran) ) then
  iBtran = Btran
endif

! check dimensions
r = size(A,1)
if ( r /= size(A,2) ) then
  stop 'A not square in trimult'
endif
if ( r /= size(B,1) ) then
  stop 'A and B inconsisent in trimult'
endif
if ( r /= size(B,2) ) then
  stop 'A and B inconsisent in trimult'
endif
if ( r /= size(C,1) ) then
  stop 'A and C inconsisent in trimult'
endif
if ( r /= size(C,2) ) then
  stop 'A and C inconsisent in trimult'
endif

if ( iAup == 0 .and. iBup == 0 ) then
  ! both lower triangular

  if ( iAtran == 0 .and. iBtran == 0 ) then
    ! No transposes, result is LT (lower*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = j, i
          X = X + A(i,k)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran == 0 .and. iBtran /= 0 ) then
    ! B transposed, result is full matrix (lower*upper)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = 1, min(i,j)
          X = X + A(i,k)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran /= 0 .and. iBtran == 0 ) then
    ! A transposed, result is full matrix (upper*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = max(i,j), r
          X = X + A(k,i)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  else
    ! Both transposed, result is UT (upper*upper)

    do j = 1, size(C,2)
      do i = 1, size(C,1)
        X = 0d0
        do k = i, j
          X = X + A(k,i)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  endif

elseif ( iAup == 0 .and. iBup /= 0 ) then
  ! B upper triangular

  if ( iAtran == 0 .and. iBtran == 0 ) then
    ! No transposes, result is full matrix (lower*upper)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = 1, min(i,j)
          X = X + A(i,k)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran == 0 .and. iBtran /= 0 ) then
    ! B transposed, result is LT (lower*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = j, i
          X = X + A(i,k)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran /= 0 .and. iBtran == 0 ) then
    ! A transposed, result is UT (upper*upper)

    do j = 1, size(C,2)
      do i = 1, size(C,1)
        X = 0d0
        do k = i, j
          X = X + A(k,i)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  else
    ! Both transposed, result is full matrix (upper*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = max(i,j), r
          X = X + A(k,i)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  endif

elseif ( iAup /= 0 .and. iBup == 0 ) then
  ! A upper triangular

  if ( iAtran == 0 .and. iBtran == 0 ) then
    ! No transposes, result is full matrix (upper*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = max(i,j), r
          X = X + A(i,k)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran == 0 .and. iBtran /= 0 ) then
    ! B transposed, result is UT (upper*upper)

    do j = 1, size(C,2)
      do i = 1, size(C,1)
        X = 0d0
        do k = i, j
          X = X + A(i,k)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran /= 0 .and. iBtran == 0 ) then
    ! A transposed, result is LT (lower*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = j, i
          X = X + A(k,i)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  else
    ! Both transposed, result is full matrix (lower*upper)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = 1, min(i,j)
          X = X + A(k,i)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  endif

else
  ! Both A and B upper triangular

  if ( iAtran == 0 .and. iBtran == 0 ) then
    ! No transposes, result is UT (upper*upper)

    do j = 1, size(C,2)
      do i = 1, size(C,1)
        X = 0d0
        do k = i, j
          X = X + A(i,k)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran == 0 .and. iBtran /= 0 ) then
    ! B transposed, result is full matrix (upper*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = max(i,j), r
          X = X + A(i,k)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo

  elseif ( iAtran /= 0 .and. iBtran == 0 ) then
    ! A transposed, result is full matrix (lower*upper)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = 1, min(i,j)
          X = X + A(k,i)*B(k,j)
        enddo
        C(i,j) = X
      enddo
    enddo

  else
    ! Both transposed, result is LT (lower*lower)

    do j = 1, r
      do i = 1, r
        X = 0d0
        do k = j, i
          X = X + A(k,i)*B(j,k)
        enddo
        C(i,j) = X
      enddo
    enddo
  
  endif

endif

end subroutine


!*****************************************************************

function vector( x_mat )

double precision, intent(in) :: X_mat(:,:)
double precision vector(size(x_mat))

call vectorize( x_mat, vector )

end function

subroutine vectorize(X_mat,X_vec)

! This routine vectorizes a (nrow x ncol) matrix:
!			X1	X2
! X_mat =	X3	X4		X_vec = [X1 X3 X5 X2 X4 X6]'	
!			X5	X6
!
!Written by M.K. Andersson, Sveriges Riksbank, 2004-09-22

implicit none
double precision, dimension(:,:), intent(in) :: X_mat
double precision, dimension(:), intent(out) :: X_vec

integer nrow, ncol, i

nrow = size( X_mat, 1 )
ncol = size( X_mat, 2 )
i = size( X_vec )

if ( i /= nrow*ncol ) then
    print *, "Error in dimensions for vectorize"
    stop
endif

do i = 1,ncol
	X_vec(nrow*(i-1)+1:nrow*i) = X_mat(:,i)
enddo

end subroutine

!*****************************************************************

subroutine save_vech_2dim( tosave, savein, row, startcol, upper )

! take vech( tosave ) and put in row 'row' of savein, 
! starting at startcol.
! tosave must by square matrix

double precision, intent(in) :: tosave(:,:)
double precision, intent(inout) :: savein(:,:)
integer, intent(in) :: row, startcol
logical, optional, intent(in) :: upper

integer m, i, j, k
logical up

up = .false.
if ( present(upper) ) up = upper

m = size( tosave, 1 )
if ( m /= size( tosave, 2 ) ) then
    stop 'Input matrix for save_vech must be square'
endif

if ( startcol + m*(m+1)/2 - 1 > size(savein,2) ) then
    stop 'To few columns in output matrix for save_vech'
endif

if ( up ) then

    k = startcol
    do j = 1, m
        do i = 1, j
            savein(row,k) = tosave(i,j)
            k = k+1
        end do
    end do

else

    k = startcol
    do j = 1, m
        do i = j, m
            savein(row,k) = tosave(i,j)
            k = k+1
        end do
    end do

endif

end subroutine

subroutine save_vech_3dim( tosave, savein, row, startcol, upper )

! take vech( tosave ) and put in row 'row' of savein, 
! starting at startcol.
! tosave must by square matrix

double precision, intent(in) :: tosave(:,:,:)
double precision, intent(inout) :: savein(:,:)
integer, intent(in) :: row, startcol
logical, optional, intent(in) :: upper

integer m, n, i, j, p, k
logical up

up = .false.
if ( present(upper) ) up = upper

m = size( tosave, 1 )
if ( m /= size( tosave, 2 ) ) then
    stop 'Input matrix for save_vech must be square'
endif
n = size( tosave, 3 )

j = size(savein,2)
if ( startcol + n*m*(m+1)/2 - 1 > size(savein,2) ) then
    stop 'To few columns in output matrix for save_vech'
endif

k = startcol
if ( up ) then

    do p = 1, n
        do j = 1, m
            do i = 1, j
                savein(row,k) = tosave(i,j,p)
                k = k+1
            end do
        end do
    end do

else

    do p = 1, n
        do j = 1, m
            do i = j, m
                savein(row,k) = tosave(i,j,p)
                k = k+1
            end do
        end do
    end do

endif

end subroutine

!*****************************************************************

subroutine compact_to_square_2dim( vector, trimat, upper, triang )

double precision, intent(in) :: vector(:)
double precision, intent(out) :: trimat(:,:)
logical, optional, intent(in) :: upper
logical, optional, intent(in) :: triang

integer m, i, j, k
logical up, tri

m = size(trimat,1)
if ( m /= size(trimat,2) ) stop 'Output matrix for compact_to_square must be square'
if ( m*(m+1)/2 /= size(vector) ) stop 'Size of input and output to compact_to_square must be same'

tri = .false.
if ( present(triang) ) tri = triang

up = .false.
if ( present(upper) ) up = upper
trimat = 0.0d0

k = 1
if ( up ) then

    do j = 1, m
        do i = 1, j
            trimat(i,j) = vector(k)
            k = k+1
            if ( .not. tri ) trimat(j,i) = trimat(i,j)
        end do
    end do

else

    do j = 1, m
        do i = j, m
            trimat(i,j) = vector(k)
            k = k+1
            if ( .not. tri ) trimat(j,i) = trimat(i,j)
        end do
    end do

endif

end subroutine

subroutine compact_to_square_3dim( vector, trimat, upper, triang )

double precision, intent(in) :: vector(:)
double precision, intent(out) :: trimat(:,:,:)
logical, optional, intent(in) :: upper
logical, optional, intent(in) :: triang

integer m, n, i, j, p, k
logical up, tri

m = size(trimat,1)
if ( m /= size(trimat,2) ) stop 'Output matrix for compact_to_square must be square'
n = size(trimat,3)
if ( n*m*(m+1)/2 /= size(vector) ) stop 'Size of input and output to compact_to_square must be same'

tri = .false.
if ( present(triang) ) tri = triang

up = .false.
if ( present(upper) ) up = upper
trimat = 0.0d0

k = 1
if ( up ) then

    do p = 1, n
        do j = 1, m
            do i = 1, j
                trimat(i,j,p) = vector(k)
                k = k+1
                if ( .not. tri ) trimat(j,i,p) = trimat(i,j,p)
            end do
        end do
    end do

else

    do p = 1, n
        do j = 1, m
            do i = j, m
                trimat(i,j,p) = vector(k)
                k = k+1
                if ( .not. tri ) trimat(j,i,p) = trimat(i,j,p)
            end do
        end do
    end do

endif

end subroutine

!*****************************************************************

!Kronecker routines written by M. Salabasis

FUNCTION KRON_MM(A,B)
REAL(KIND(1D0)),DIMENSION(:,:) ::A,B
REAL(KIND(1D0)), DIMENSION(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2)) ::KRON_MM
INTEGER :: I,J
	DO I=1,SIZE(A,1)
		DO J=1,SIZE(A,2)
			KRON_MM((I-1)*SIZE(B,1)+1:I*SIZE(B,1),(J-1)*SIZE(B,2)+1:J*SIZE(B,2))=A(I,J)*B(:,:)
		END DO
	END DO
END FUNCTION KRON_MM

FUNCTION KRON_VM(A,B)
REAL(KIND(1D0)),DIMENSION(:) ::A
REAL(KIND(1D0)), DIMENSION(:,:) ::B
REAL(KIND(1D0)), DIMENSION(SIZE(B,1)*SIZE(A,1),SIZE(B,2)) ::KRON_VM
INTEGER :: I
	DO I=1,SIZE(A,1)
		KRON_VM((I-1)*SIZE(B,1)+1:I*SIZE(B,1),:)=A(I)*B(:,:)
	END DO
END FUNCTION KRON_VM

FUNCTION KRON_MV(A,B)
REAL(KIND(1D0)), DIMENSION(:,:) ::A
REAL(KIND(1D0)),DIMENSION(:) ::B
REAL(KIND(1D0)), DIMENSION(SIZE(A,1)*SIZE(B,1),SIZE(A,2)) ::KRON_MV
INTEGER :: I,J
	DO I=1,SIZE(A,1)
		DO J=1,SIZE(A,2)
			KRON_MV((I-1)*SIZE(B,1)+1:I*SIZE(B,1),J)=A(I,J)*B(:)
		END DO
	END DO
END FUNCTION KRON_MV

FUNCTION KRON_SM(A,B)
REAL(KIND(1D0)) ::A
REAL(KIND(1D0)), DIMENSION(:,:) ::B
REAL(KIND(1D0)), DIMENSION(SIZE(B,1),SIZE(B,2)) ::KRON_SM
	KRON_SM(:,:)=A*B(:,:)
END FUNCTION KRON_SM

FUNCTION KRON_MS(A,B)
REAL(KIND(1D0)), DIMENSION(:,:) ::A
REAL(KIND(1D0)) ::B
REAL(KIND(1D0)), DIMENSION(SIZE(A,1),SIZE(A,2)) ::KRON_MS
	KRON_MS(:,:)=A(:,:)*B
END FUNCTION KRON_MS

FUNCTION KRON_SV(A,B)
REAL(KIND(1D0)) ::A
REAL(KIND(1D0)), DIMENSION(:) ::B
REAL(KIND(1D0)), DIMENSION(SIZE(B,1),1) ::KRON_SV
	KRON_SV(:,1) = A*B(:)
END FUNCTION KRON_SV

FUNCTION KRON_VS(A,B)
REAL(KIND(1D0)), DIMENSION(:) ::A
REAL(KIND(1D0)) ::B
REAL(KIND(1D0)), DIMENSION(SIZE(A,1),1) ::KRON_VS
	KRON_VS(:,1) = A(:)*B
END FUNCTION KRON_VS

FUNCTION KRON_SS(A,B)
REAL(KIND(1D0)) ::A,B
REAL(KIND(1D0)), DIMENSION(1,1) ::KRON_SS
	KRON_SS(1,1) = A*B
END FUNCTION KRON_SS

FUNCTION KRON_VV(A,B)
REAL(KIND(1D0)), DIMENSION(:) ::A,B
REAL(KIND(1D0)), DIMENSION(SIZE(A,1)*SIZE(B,1),1) ::KRON_VV
INTEGER ::I
	DO I=1,SIZE(A,1)
		KRON_VV((I-1)*SIZE(B,1)+1:I*SIZE(B,1),1) = A(I)*B(:)
	END DO
END FUNCTION KRON_VV

!-----------------------------------------------------------------

subroutine cholesky(U,n,p) !copied from Numerical Recipes, ch. 2.9
integer, intent(in):: n
double precision, intent(inout):: U(:,:),p(:)
! Given U positive-definite symmetric matrix U(1:n,1:n),
! this routine constructs its Cholesky decomposition, U = L·L'. 
! On input, only the upper triangle of U need be given; it is not modified. 
! The Cholesky factor L is returned in the lower triangle
! of U, except for its diagonal elements which are returned in p(1:n).
integer i,j,k
double precision sum

do i=1,n !is for rows

	do j= i,n ! for cols
	
		sum=U(i,j)
		do  k=i-1,1,-1
			sum=sum-U(i,k)*U(j,k)
		enddo 

		if(i == j)then
			if(sum .le. 0d0) then
				write(*,*) 'choldc failed'		!U, with rounding errors, is not
			endif										!positive definite.
			p(i)=dsqrt(sum)						
	
		else
		
			U(j,i)=sum/p(i)
		
		endif

	enddo

enddo

END subroutine

!------------------------------------------------------------

SUBROUTINE solve_chol(U,n,p,b, x, ssr) !copied from Numerical Recipes, ch. 2.9
integer, intent(in) :: n
double precision, intent(in):: U(:,:), p(:), b(:)
double precision, intent(out) ::x(:), ssr

! Solves the set of n linear equations U.x = b, 
! where U is U positive-definite symmetric matrix.
! U and p are input as the output of the routine choldc.
! Only the lower triangle of U is accessed. 
! b(1:n) is input as the right-hand side vector. 
! The solution vector is returned in x(1:n). 
! U, n, and p are not modified and can be left
! in place for successive calls with different right-hand sides b.
! b is not modified unless you identify b and x in the calling sequence, 
! which is allowed.

integer i,k
double precision:: sum
do i=1,n			!Solve L ?y = b, storing y in x.

	sum=b(i)

	do k=i-1,1,-1
		sum=sum-U(i,k)*x(k)
	
	enddo

	x(i)=sum/p(i)

enddo
ssr = 0d0
do i=n,1,-1
	sum=x(i)
	do k=i+1,n
		sum=sum-U(k,i)*x(k)
	
	enddo
	x(i)=sum/p(i)
	ssr = ssr +x(i)*b(i)
enddo

end subroutine



end module