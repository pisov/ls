program harmonics
implicit none

integer :: i, n, info
double precision :: x
double precision, dimension (:), allocatable :: d, e, b,sol
integer, dimension (:), allocatable :: ipiv
real :: start_time, end_time

! Input the problem size n
do
  write(0,*) 'Please enter problem size n ='
  read(*,*)n
  if ( n > 2 ) exit
end do

! Allocate work arrays
allocate(d(n))
allocate(e(n-1))
allocate(b(n))
allocate(sol(0:n+1))
allocate(ipiv(n))

! Initialize the interaction matrix a
d(:) =  2.d0
e(:) = -1.d0
b(:) = -1.d0 / ((n+2)*(n+1))

! Call LAPACK Linear Solver DGESV
call cpu_time(start_time)
call DPTSV(n, 1, d, e, b, n, info)
call cpu_time(end_time)

write(0,'(A18,F7.2)')'Execution time = ',end_time-start_time

!print out the result
sol(:) = 0.d0
sol(1:n) = b(1:n)
do i = 0, n+1
 x = i * 1.d0 / (n+1)
 write(*,'(2E15.5)'),x,sol(i)
end do 


! Deallocate work arrays
deallocate(d)
deallocate(e)
deallocate(b)
deallocate(ipiv)
deallocate(sol)


end program
