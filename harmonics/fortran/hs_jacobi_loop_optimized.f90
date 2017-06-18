program harmonics
implicit none

integer :: i, j, n, info
double precision :: xi, eps, diff
double precision, dimension (:,:), allocatable :: a
double precision, dimension (:), allocatable :: b, x, xnew, sol
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
allocate(x(0:n+1))
allocate(xnew(0:n+1))

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

! Solution
eps    = 1.d-7
x(:)   = b(:)
x(0)   = 0.d0
x(n+1) = 0.d0
call cpu_time(start_time)
do
  xnew(1:n) = 0.d0
  do j = 1, n
    do i = 1, n
    if (i .NE. j) xnew(i) =  xnew(i) + a(i,j) * x(j)
    end do
!    xnew(i) = (b(i) - xnew(i))/a(i,i)
  end do
  xnew(1:n) = (b(1:n) - xnew(1:n)) / 2
  diff = dsqrt(dot_product(xnew(1:n)-x(1:n),xnew(1:n)-x(1:n)))
  if ( diff .LE. eps ) exit
  x(1:n) = xnew(1:n)
end do
call cpu_time(end_time)

write(0,'(A18,F7.2)')'Execution time = ',end_time-start_time

!print out the result
do i = 0, n+1
 xi = i * 1.d0 / (n+1)
 write(*,'(2E15.5)'),xi,x(i)
end do 


! Deallocate work arrays
deallocate(a)
deallocate(b)
deallocate(x)
deallocate(xnew)


end program
