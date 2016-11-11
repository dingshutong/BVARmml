module times

implicit none

contains

subroutine diff_seconds( time1, time2, timediff )

! Input is times in seconds, e.g. from cpu_time()
! output is integer vector with days, hours, minutes, seconds and hundreths of seconds

real, intent(in) :: time1, time2
integer, intent(out) :: timediff(5)

real tdiff
integer days, hours, minutes, seconds, hundreths

tdiff = time2-time1

if ( tdiff < 0 ) stop "reverse times in diff_seconds"

days = int( tdiff/8640 )
tdiff = tdiff - days*8640
hours = int( tdiff/3600 )
tdiff = tdiff - hours*3600
minutes = int( tdiff/60 )
tdiff = tdiff - minutes*60
seconds = int( tdiff )
tdiff = tdiff - seconds
hundreths = int( tdiff*100 )

timediff = (/ days, hours, minutes, seconds, hundreths /)

end subroutine

!----------------------------------------------------------

subroutine diff_seconds_as_string( time1, time2, timediff )

! Input is times in seconds, e.g. from cpu_time()
! output is string variable on the format
! "x days, y hours, z minutes, w.v seconds"

real, intent(in) :: time1, time2
character(len=*), intent(out) :: timediff

integer i, j, k, diff(5)
logical started

if ( len(timediff) < 29 ) stop "string for time to short in diff_seconds_as_string"

call diff_seconds( time1, time2, diff )

timediff = ''

started = .false.
i = 1
if ( diff(1) > 0 ) then
    j = int(log10(dble(diff(1))))+1
    k = j + i + 4 - 1
    write(timediff(i:k),100) diff(1), " d, "
    started = .true.
    i = k+1
endif
if ( diff(2) > 0 .or. started ) then
    if ( diff(2) > 0 ) then
        j = int(log10(dble(diff(2))))+1
    else
        j = 1
    endif
    k = j + i + 4 - 1
    write(timediff(i:k),100) diff(2), " h, "
    started = .true.
    i = k+1
endif
if ( diff(3) > 0 .or. started ) then
    if ( diff(3) > 0 ) then
        j = int(log10(dble(diff(3))))+1
    else
        j = 1
    endif
    k = j + i + 4 - 1
    write(timediff(i:k),100) diff(3), " m, "
    started = .true.
    i = k+1
endif

if ( diff(4) > 0 ) then
    j = int(log10(dble(diff(4))))+1
else
    j = 1
endif
k = j + i + 3 + 2 - 1
write(timediff(i:k),101) diff(4), ".", diff(5), " s"

100 format( I<j>, A )
101 format( I<j>, A, I2.2, A )

end subroutine

!----------------------------------------------------------

function date_time_now_string()

character(len=19) date_time_now_string
integer values(8)

call date_and_time( values=values )

write(date_time_now_string,100) values(1:3), values(5:7)

100 format( I4,'-',I2.2,'-',I2.2,' ',I2,':',I2.2,':',I2.2 )

end function

!----------------------------------------------------------

function time_diff( time1, time2 )

real time_diff

integer, intent(in) :: time1(8), time2(8)
! values from DATE_AND_TIME

time_diff = 86400d0*( days_since_2008( time2(1), time2(2), time2(3) ) - &
                      days_since_2008( time1(1), time1(2), time1(3) ) ) &
          + 3600d0*( time2(5) - time1(5) )                              &
          + 60d0*( time2(6) - time1(6) )                                &
          + time2(7) - time1(7) + real( time2(8) - time1(8) )/1.0e3
          
end function

!----------------------------------------------------------

function days_since_2008( year, month, day )

! number of days since 2008-01-01
! days_since_2008( 2008, 1, 1 ) = 0

integer days_since_2008
integer, intent(in) :: year, month, day

integer i

days_since_2008 = ( year - 2008 )*365
! leap years
if ( year > 2008 ) then
    do i = 2008, year-1, 4
        days_since_2008 = days_since_2008 + 1
    end do
else if ( year < 2008 ) then
    do i = 2004, year, -4
        days_since_2008 = days_since_2008 - 1
    end do
endif

do i = 1, month-1
    select case (i)
    case( 1, 3, 5, 7, 8, 10, 12 )
        days_since_2008 = days_since_2008 + 31
    case( 2 )
        days_since_2008 = days_since_2008 + 28
        if ( mod(year,4) == 0 ) days_since_2008 = days_since_2008 + 1
    case default
        days_since_2008 = days_since_2008 + 30
    end select
end do

days_since_2008 = days_since_2008 + day - 1

end function

end module