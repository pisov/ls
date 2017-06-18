program harmonics
implicit none

integer :: i, n, info
double precision :: x
double precision, dimension (:,:), allocatable :: a
double precision, dimension (:), allocatable :: b,sol
integer, dimension (:), allocatable :: ipiv
real :: start_time, end_time

! Input the problem size n
do
  write(0,*) 'Please enter problem size n ='
  read(*,*)n
  if ( n > 2 ) exit
end do

! Allocate work arrays
allocate(a(n,n))
allocate(b(n))
allocate(sol(0:n+1))
allocate(ipiv(n))

! Initialize the interaction matrix a
a(:,:) = 0.d0

do i = 1, n
  a(i,i) = 2.d0
  if ( i > 1) then
    a(i,i-1) = -1.d0
  end if
  if (i < n) then
    a(i,i+1) = -1.d0
  end if
  b(i) = -1.d0 / ((n+2)*(n+1))
end do

! Call LAPACK Linear Solver DGESV
call cpu_time(start_time)
call dgesv(n,1,a,n,ipiv,b,n,info)
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
deallocate(a)
deallocate(b)
deallocate(ipiv)
deallocate(sol)


end program
