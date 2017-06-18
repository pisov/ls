program linsolv

implicit none

integer,parameter :: nmax = 4
integer,parameter :: lda = nmax
integer :: n = nmax
real(kind=8), dimension(nmax, nmax) :: a, ainv
real(kind=8), dimension(nmax) :: x, b
integer :: i, j, info, ipiv(lda)

! Initialize the matix A
a = reshape((/2.,-1.,0.,0.,-1.,2.,-1.,0.,0.,-1.,2.,-1.,0.,0.,-1.,2./),shape(a))

!Initialize right side vector b
b = (/1.,1.,1.,1./)

!write(*,'(4F10.5)')(b(i),i=1,n)

! Solve linear system A.x = b
!       -1
! A -> A   and b -> x
! 

ainv = a
x = b

! Make LAPACK call
call dgesv(n, 1, ainv, lda, ipiv, x, n, info)

! Calculate b = A.x

b = matmul(a, x)

! Print vector b
write(*,'(4F10.5)')(b(i),i=1,n)

end program
