module prettyprint

implicit none

integer, parameter :: outputwidth = 72

contains

subroutine write_table( matrix, title, rowlabels, collabels, sigdigits, unit )

double precision, intent(in) :: matrix(:,:)
character(len=*), optional, intent(in) :: title, rowlabels(:), collabels(:)
integer, optional, intent(in) :: sigdigits(:)
integer, optional, intent(in) :: unit

integer titlewidth, rowlabwidth, totalwidth, collabwidth(size(matrix,2)), &
        colwidth(size(matrix,2)), coldigits(size(matrix,2)), &
        colnumwidth(size(matrix,2))
integer rows, cols, i, j, k, skip, after, dd
character(len=60) frmt(size(matrix,2)), afrmt
double precision compvals(size(matrix,1))

rows = size( matrix, 1 )
cols = size( matrix, 2 )

if ( present(title) ) then
    titlewidth = len_trim(adjustl(title))
else
    titlewidth = 0
endif
rowlabwidth = 0
if ( present(rowlabels) ) then
!    if ( size(rowlabels) /= rows ) stop "Dim of row labels does not match matrix dim"
    do i = 1, max( rows, size(rowlabels) )
        rowlabwidth = max( len_trim(adjustl(rowlabels(i))), rowlabwidth )
    end do
endif
collabwidth = 0
if ( present(collabels) ) then
!    if ( size(collabels) /= cols ) stop "Dim of column labels does not match matrix dim"
    do i = 1, cols
     k = mod( i+size(collabels)-1, size(collabels) ) + 1
     collabwidth(i) = len_trim(adjustl(collabels(k)))
    end do
endif
coldigits = 4 ! four significant digits 
if ( present(sigdigits) ) then
    if ( size(sigdigits) /= cols ) stop "Dim of significant digits does not match matrix dim"
    coldigits = sigdigits
endif

do j = 1, cols
    if ( coldigits(j) < 0 ) then
        ! fixed format with this many digits
        frmt(j) = '(tr1,f  .  )'
        write( frmt(j)(7:8), '(i2.2)' ) abs(coldigits(j))+2
        write( frmt(j)(10:11), '(i2.2)' ) abs(coldigits(j))-1
        colnumwidth(j) = abs(coldigits(j))+2
    else
       k = 0
       do i = 1, rows
            if ( matrix(i,j) /= 0.0d0 ) then
                k = k+1
                compvals(k) = abs(matrix(i,j))
            endif
        enddo
        if ( log10( maxval(compvals(1:k)) ) < abs(coldigits(j)) .and. &
             log10( minval(compvals(1:k)) ) > -abs(coldigits(j)) ) then
            frmt(j) = '(tr1,f  .  )'
            write( frmt(j)(7:8), '(i2.2)' ) 2*abs(coldigits(j))+1
            write( frmt(j)(10:11), '(i2.2)' ) abs(coldigits(j))
            colnumwidth(j) = 2*abs(coldigits(j))+1
        else
            frmt(j) = '(tr1,g  .  )'
            write( frmt(j)(7:8), '(i2.2)' ) abs(coldigits(j))+7
            write( frmt(j)(10:11), '(i2.2)' ) abs(coldigits(j))
            colnumwidth(j) = abs(coldigits(j))+7
        endif
    endif
    colwidth(j) = max( collabwidth(j), colnumwidth(j) )
end do

totalwidth = rowlabwidth + 1 + sum(colwidth) + 1*(cols - 1)
totalwidth = max(totalwidth,titlewidth)

! write '='
call write_dline( totalwidth, unit )
if ( present(title) ) then
    call write_centered( title, totalwidth, unit )
    call write_sline( totalwidth, unit )
endif

if ( present(collabels) ) then
    after = colwidth(1) - collabwidth(1) + 1
    skip = rowlabwidth + 1 + after/2
    after = after - after/2
    if ( present(unit) ) then
        write(unit,104,advance='no') trim(adjustl(collabels(1)))
    else
        write(*,104,advance='no') trim(adjustl(collabels(1)))
    endif
    skip = 0
    do j = 2, cols
        after = colwidth(j) - collabwidth(j) + 1
        skip = after/2
        after = after - after/2
        k = mod( j+size(collabels)-1, size(collabels) ) + 1
        if ( present(unit) ) then
            write(unit,104,advance='no') trim(adjustl(collabels(k)))
        else
            write(*,104,advance='no') trim(adjustl(collabels(k)))
        endif
    end do
    if ( present(unit) ) then
        write(unit,*) 
    else
        write(*,*)
    endif
    ! write '-'
    call write_sline( totalwidth, unit )
endif

! write data
do i = 1, rows
    if ( present(rowlabels) ) then
        k = mod( i, size(rowlabels) ) + 1
        if ( present(unit) ) then
            write(unit,106,advance='no') rowlabels(k)
        else
            write(*,106,advance='no') rowlabels(k)
        endif
    endif
    do j = 1, cols
        afrmt = frmt(j)
        dd = abs(coldigits(j))
        if ( present(unit) ) then
            write(unit,afrmt,advance='no') matrix(i,j)
        else
            write(*,afrmt,advance='no') matrix(i,j)
        endif
    end do
    if ( present(unit) ) then
        write(unit,*) 
    else
        write(*,*)
    endif
end do

call write_dline( totalwidth, unit )

!103 format (tr<skip>,A)
104 format (tr<skip>,A,tr<after>)
105 format (A,tr<after>)
106 format (tr1,A<rowlabwidth>)

end subroutine

!---------------------------------------------------------------

subroutine write_dline( length, unit )

integer, optional, intent(in) :: length
integer, optional, intent(in) :: unit

integer linelength

if ( present( length ) ) then
    linelength = length
else
    linelength = outputwidth
endif

if ( present(unit) ) then
    write(unit,101)
else
    write(*,101)
endif

101 format (tr1,<linelength>('='))

end subroutine

!---------------------------------------------------------------

subroutine write_sline( length, unit )

integer, optional, intent(in) :: length
integer, optional, intent(in) :: unit

integer linelength

if ( present( length ) ) then
    linelength = length
else
    linelength = outputwidth
endif

if ( present(unit) ) then
    write(unit,101)
else
    write(*,101)
endif

101 format (tr1,<linelength>('-'))

end subroutine

!---------------------------------------------------------------

subroutine write_centered( string, length, unit )

character(len=*), intent(in) :: string
integer, optional, intent(in) :: length
integer, optional, intent(in) :: unit

integer i, skip, linelength

if ( present( length ) ) then
    linelength = length
else
    linelength = outputwidth
endif

i = len_trim(adjustl(string))

skip = (linelength-i)/2
if ( present(unit) ) then
    write(unit,103) trim(adjustl(string))
else
    write(*,103) trim(adjustl(string))
endif

103 format (tr<skip+1>,A)

end subroutine

!---------------------------------------------------------------

subroutine write_string_dble( string, dble, unit, width, decimal )

character(len=*), intent(in) :: string
double precision, intent(in) :: dble
integer, optional, intent(in) :: unit, width, decimal

integer w, d, chk
character(len=30) frmt

frmt = '(tr1,a,f  .  )'
w = 8
d = 3
if ( present(width) ) w = width
if ( present(decimal) ) d = decimal
if ( abs(dble) > 1d1*tiny(0.0d0) ) then
    chk = int(log10(abs(dble))) + 1
    if ( chk > 0 .and. dble < 0 ) chk = chk+1
    if ( chk > w-d-1 .or. chk < 1-d ) then
        frmt = '(tr1,a,g  .  )'
        w = 11
        d = 4
    endif
endif
write( frmt(9:10), '(i2.2)' ) w
write( frmt(12:13), '(i2.2)' ) d
    
if ( present(unit) ) then
    write(unit,frmt) string, dble
else
    write(*,frmt) string, dble
endif

103 format (tr1,a,f<w>.<d>)

end subroutine

!---------------------------------------------------------------

subroutine write_string_int( string, integ, unit, width, string2, integ2, string3, integ3 )

character(len=*), intent(in) :: string
integer, intent(in) :: integ
integer, optional, intent(in) :: unit, width
character(len=*), optional, intent(in) :: string2, string3
integer, optional, intent(in) :: integ2, integ3

integer w
character*3 adv

w = get_width_for_int( integ, width )

if ( present(string2) .or. present(string3) .or. present(integ2) .or. present(integ3) ) then
    adv = 'no'
else
    adv = 'yes'
endif
if ( present(unit) ) then
    write(unit,103,advance=adv) string, integ
else
    write(*,103,advance=adv) string, integ
endif

if ( present(string2) ) then

    if ( present(string3) .or. present(integ2) .or. present(integ3) ) then
        adv = 'no'
    else
        adv = 'yes'
    endif
    if ( present(unit) ) then
        write(unit,104,advance=adv) string2
    else
        write(*,104,advance=adv) string2
    endif

endif

if ( present(integ2) ) then

    w = get_width_for_int( integ2, width )
    if ( present(string3) .or. present(integ3) ) then
        adv = 'no'
    else
        adv = 'yes'
    endif
    if ( present(unit) ) then
        write(unit,105,advance=adv) integ2
    else
        write(*,105,advance=adv) integ2
    endif

endif

if ( present(string3) ) then

    if ( present(integ3) ) then
        adv = 'no'
    else
        adv = 'yes'
    endif
    if ( present(unit) ) then
        write(unit,104,advance=adv) string3
    else
        write(*,104,advance=adv) string3
    endif

endif

if ( present(integ3) ) then

    w = get_width_for_int( integ3, width )
    if ( present(unit) ) then
        write(unit,105) integ3
    else
        write(*,105) integ3
    endif

endif

103 format (tr1,a,tr1,i<w>)
104 format (tr1,a)
105 format (tr1,i<w>)

end subroutine

!---------------------------------------------------------------

function get_width_for_int( integ, width )

integer get_width_for_int

integer, intent(in) :: integ
integer, optional, intent(in) :: width

integer w

if ( present(width) ) then
    w = width
elseif ( integ == 0 ) then
    w = 1
else
    w = int(log10(dble(abs(integ)))) + 1
    if ( integ < 0 ) w = w + 1
endif

get_width_for_int = w

end function

!---------------------------------------------------------------

subroutine write_string( string, unit )

character(len=*), intent(in) :: string
integer, optional, intent(in) :: unit

if ( present(unit) ) then
    write(unit,103) string
else
    write(*,103) string
endif

103 format (tr1,a)

end subroutine

!---------------------------------------------------------------

subroutine write_stringarray( string, strarr, unit, width )

character(len=*), intent(in) :: string, strarr(:)
integer, optional, intent(in) :: unit, width

integer i, k, w

w = outputwidth
if ( present(width) ) w = width

k = len(string)
if ( present(unit) ) then
    write(unit,103,advance='no') string
else
    write(*,103,advance='no') string
endif

do i = 1, size(strarr)
    k = k + len_trim(adjustl(strarr(i)))
    if ( k > w ) then
        if ( present(unit) ) then
            write(unit,*) 
            write(unit,103,advance='no') '         '
        else
            write(*,*)
            write(*,103,advance='no') '         '
        endif
        k = 10 + len_trim(adjustl(strarr(i)))
    endif
    if ( present(unit) ) then
        write(unit,103,advance='no') trim(adjustl(strarr(i)))
    else
        write(*,103,advance='no') trim(adjustl(strarr(i)))
    endif
end do

if ( present(unit) ) then
    write(unit,*)
else
    write(*,*)
endif

103 format (tr1,a)

end subroutine

!---------------------------------------------------------------

subroutine make_array_label( labelarr, prefix, m, n, dim3, vech )

character(len=*), intent(out) :: labelarr(:)
character(len=*), intent(in) :: prefix
integer, intent(in) :: m, n
integer, optional, intent(in) :: dim3
character(len=*), optional, intent(in) :: vech

integer d3, i, j, k, jj, w, ww, www, ii

if ( present(dim3) ) then
    d3 = dim3
    if ( d3 < 1 ) stop "Third dimension in make_array_label must be positive"
else
    d3 = 1
endif

! check length of labelarr
if ( present(vech) ) then
    if ( n /= m ) then
        stop "Array must be square with vech in make_array_label"
    endif
    if ( size( labelarr ) /= d3*m*(m+1)/2 ) then
        stop "Size mismatch in make_array_label"
    endif
elseif ( size( labelarr ) /= d3*m*n ) then
    stop "Size mismatch in make_array_label"
endif

ii = len_trim( prefix )
jj = int(log10(dble(m))) + int(log10(dble(n))) + 2 + 3
if ( d3 > 1 ) then
    jj = jj + int(log10(dble(m))) + 2
endif
if ( ii+jj > len(labelarr(1)) ) then
    print *, "Array labels to short in make_array_labels, ", ii+jj, " needed"
    stop
endif

if ( present( vech ) ) then
    ! only do labels for upper or lower triangle

    if ( vech(1:1) == 'u' .or. vech(1:1) == 'U' ) then
        ! upper triangle
        jj = 1
        do k = 1, d3
            do j = 1, n
                do i = 1, j
                    w = int(log10(dble(i))) + 1
                    ww = int(log10(dble(j))) + 1
                    if ( d3 == 1 ) then
                        write(labelarr(jj),110) prefix(:ii), i, j
                    else
                        www = int(log10(dble(k))) + 1
                        write(labelarr(jj),111) prefix(:ii), i, j, k
                    endif
                    jj = jj+1
                end do
            end do
        end do
    elseif ( vech(1:1) == 'l' .or. vech(1:1) == 'L' ) then
        ! lower triangle
        jj = 1
        do k = 1, d3
            do j = 1, n
                do i = j, m
                    w = int(log10(dble(i))) + 1
                    ww = int(log10(dble(j))) + 1
                    if ( d3 == 1 ) then
                        write(labelarr(jj),110) prefix(:ii), i, j
                    else
                        www = int(log10(dble(k))) + 1
                        write(labelarr(jj),111) prefix(:ii), i, j, k
                    endif
                    jj = jj+1
                end do
            end do
        end do
    else
        stop 'vech must be U or L for make_array_label'
    endif
    
else    
    ! Full array    

    jj = 1
    do k = 1, d3
        do j = 1, n
            do i = 1, m
                w = int(log10(dble(i))) + 1
                ww = int(log10(dble(j))) + 1
                if ( d3 == 1 ) then
                    write(labelarr(jj),110) prefix(:ii), i, j
                else
                    www = int(log10(dble(k))) + 1
                    write(labelarr(jj),111) prefix(:ii), i, j, k
                endif
                jj = jj+1
            end do
        end do
    end do
    
endif

110 format( A, '(', I<w>, ',', I<ww>, ')' )
111 format( A, '(', I<w>, ',', I<ww>, ',', I<www>, ')' )

end subroutine

end module