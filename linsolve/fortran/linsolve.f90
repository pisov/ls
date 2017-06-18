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
call dgesv(n, 1, a, lda, ipiv, b, n, info)

! Print solution vector 
write(*,'(4F10.5)')(b(i),i=1,n)

end program
