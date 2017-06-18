program linsolv

implicit none

integer,parameter :: nmax = 4
integer,parameter :: lda = nmax
integer :: n = nmax
real(kind=8), dimension(nmax, nmax) :: a, ainv, r
real(kind=8), dimension(nmax) :: x, b, d, xnew
integer :: i, j, info, ipiv(lda)
real(kind=8), parameter :: eps = 1.e-6
real(kind=8) :: diff

! Initialize the matix A
a = reshape((/2.,-1.,0.,0.,-1.,2.,-1.,0.,0.,-1.,2.,-1.,0.,0.,-1.,2./),shape(a))

!Initialize right side vector b
b = (/1.,1.,1.,1./)

!Initialize diagonal vector
d(:) = 2.0

! Initialize remainder matrix R

r = reshape((/0.,-1.,0.,0.,-1.,0.,-1.,0.,0.,-1.,0.,-1.,0.,0.,-1.,0./),shape(a)) 


! Solve linear system A.x = b
do
  xnew = (b - matmul(r, x)) / d
  diff = sqrt(dot_product(xnew - x, xnew - x))
  x = xnew
  if (diff.le.eps) exit
end do


! Print solution vector 
write(*,'(4F10.5)')(x(i),i=1,n)

end program
