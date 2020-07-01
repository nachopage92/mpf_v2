!*************************************************************************
!                   PRIOR WEIGHTFUNCTION MODULE                          !
!                   ===========================                          !
!                                                                        !
! Purpose                                                                !
! =======                                                                !
! The prior (weight function) for the maxent basis function computations !
! Author: N. Sukumar, UC Davis                                           !
! Date  : September 2007 (Version 1.4)                                   !
!                                                                        !
! Reference                                                              !
! =========                                                              !
! N. Sukumar and R. W. Wright, IJNME, Vol. 70, Number 2, 181-205, 2007   !
!                                                                        !
!*************************************************************************

MODULE priorweightfunction

! module to define prior (weight function)

implicit none     ! require explicit type declaration/checking
private           ! make everything private by default

public   weightfunction, dweightfunction          ! make these public
public   ddweightfunction                         ! make these public
public   setrmaxandpriorweight                    ! make these public
public   getrmax                                  ! make these public
public   getpriortype                             ! make these public

! types, constants ....

integer, parameter:: dp=kind(0.d0)   ! double precision type

private rmax, D, prior, small, iprior
real(dp), pointer, dimension(:)     :: rmax(:) ! support size
real(dp), pointer, dimension(:,:,:) :: D(:,:,:)! metric
character(80)                       :: prior   ! string (uniform/cubic/quartic/
                                               !         gaussian/gaussian-rbf)
integer iprior
real(dp) :: small = tiny(1.d0)

! functions, subroutines ....
CONTAINS

!*************************************************************************
! FUNCTION weightfunction(index,vec,rvec,q)
! Purpose
! =======
! Define the weight function
!*************************************************************************
!
FUNCTION weightfunction(index,vec,rvec,q)

integer , intent(in) :: index
real(dp), intent(in) :: q
real(dp), dimension(:), intent(in) :: vec, rvec
real(dp) :: beta
real(dp) :: weightfunction

  if (iprior .eq. 0) then
    weightfunction = 1.d0
  elseif (iprior .eq. 1) then  ! C^2 weight function
    if (q <= 0.5d0) then  
      weightfunction = 2.d0/3.d0 - 4.d0*q*q + 4.d0*q*q*q
    elseif (q > 0.5d0 .and. q <= 1.d0) then
      weightfunction = 4.d0/3.d0 - 4.d0*q + 4.d0*q*q - 4.d0/3.d0*q*q*q
    else 
      weightfunction = 0.d0
    endif
  elseif (iprior .eq. 2) then  ! C^2 weight function
    if (q <= 1.d0) then
      weightfunction = 1.d0 - 6.d0*q*q + 8.d0*q*q*q - 3.d0*q*q*q*q
    else
      weightfunction = 0.d0
    endif
  elseif (iprior .eq. 3) then  ! C-infty weight function
    if (q <= 1.d0) then
      weightfunction = exp(1.d0/(q*q-1.d0))
    else
      weightfunction = 0.d0
    endif
  elseif (iprior .eq. 4) then  ! Gaussian RBF
    beta = 3.d0/(rmax(index)+small)**2                        ! 
    weightfunction = exp(-beta*dot_product(rvec,rvec)) 
  endif

END FUNCTION weightfunction

!*************************************************************************
! FUNCTION dweightfunction(index,vec,rvec,q)
! Purpose
! =======
! Compute the weight function derivative
!*************************************************************************
!
FUNCTION dweightfunction(index,vec,rvec,q)

integer , intent(in) :: index
real(dp), intent(in) :: q
real(dp), dimension(:), intent(in) :: vec, rvec

real(dp) ::  t, beta
real(dp), dimension(size(vec)) :: dweightfunction, dq

dq = -rvec/(max(q,1.d-30)*rmax(index)*rmax(index))

  if (iprior .eq. 0) then
    dweightfunction = 0.d0
  elseif (iprior .eq. 1) then  ! C^2 weight function
    if (q > 0.d0 .and. q <= 0.5d0) then  
      dweightfunction = (-8.d0*q + 12.d0*q*q)*dq
    elseif (q > 0.5 .and. q <= 1.0) then
      dweightfunction = (-4.d0 + 8.d0*q - 4.d0*q*q)*dq
    else 
      dweightfunction = 0.d0
    endif
  elseif (iprior .eq. 2) then  ! C^2 weight function
    if (q > 0.d0 .and. q <= 1.d0) then
      dweightfunction = (-12.d0*q + 24.d0*q*q - 12.d0*q*q*q)*dq
    else
      dweightfunction = 0.d0
    endif
  elseif (iprior .eq. 3) then  ! C-infty weight function
    if (q > 0.d0 .and. q < 1.0) then
      t = q*q - 1.d0
      dweightfunction = (-exp(1.d0/t)*2.d0*q/(t*t))*dq
    else
      dweightfunction = 0.d0
    endif
  elseif (iprior .eq. 4) then  ! Gaussian RBF
    beta = 3.d0/(rmax(index)+small)**2 
    t = exp(-beta*dot_product(rvec,rvec))
    dweightfunction = -2.d0*beta*rvec*t 
  endif

END FUNCTION dweightfunction

!*************************************************************************
! FUNCTION ddweightfunction(index,vec,rvecin,qin)
! Purpose
! =======
! Compute the second derivatives of the weight function 
!*************************************************************************
!
FUNCTION ddweightfunction(index,vec,rvecin,qin)

integer , intent(in) :: index
real(dp), intent(in) :: qin
real(dp), dimension(:), intent(in) :: vec, rvecin

real(dp) :: q, t, beta
real(dp), dimension(size(vec)) :: dq, rvec
real(dp), dimension(size(vec),size(vec)) :: ddweightfunction
real(dp), dimension(size(vec),size(vec)) :: dq2, id, mat

q = qin
rvec = rvecin

if (q .eq. 0.d0) then                      ! find gradient close to q = 0
  rvec = matmul(D(:,:,index),(vec+1.d-30*rmax(index)))
  q = sqrt(dot_product((vec+1.d-30*rmax(index)),rvec))/rmax(index) 
endif

dq = -rvec/(q*rmax(index)*rmax(index))
id = D(:,:,index)/(rmax(index)*rmax(index))

dq2 = outerproduct(dq,dq)
  
mat = (id - dq2)/q

if (iprior .eq. 0) then
  ddweightfunction = 0.d0
elseif (iprior .eq. 1) then  ! C^2 weight function
  if (q > 0.d0 .and. q <= 0.5) then  
    ddweightfunction = (-8.d0*q + 12.d0*q*q)*mat + (-8.d0 + 24.d0*q)*dq2
  elseif (q > 0.5 .and. q <= 1.0) then
    ddweightfunction = (-4.d0 + 8.d0*q - 4.d0*q*q)*mat + (8.d0 - 8.d0*q)*dq2
  else 
    ddweightfunction = 0.d0
  endif
elseif (iprior .eq. 2) then  ! C^2 weight function
  if (q > 0.d0 .and. q <= 1.d0) then
    ddweightfunction = (-12.d0 + 24.d0*q - 12.d0*q*q)*(id-dq2) + & 
                       (-12.d0 + 48.d0*q - 36.d0*q*q)*dq2
  else
    ddweightfunction = 0.d0
  endif
elseif (iprior .eq. 3) then  ! C-infty weight function
  if (q > 0.d0 .and. q < 1.d0) then
    t = q*q - 1.d0
    ddweightfunction = -exp(1.d0/t)*2.d0*q/(t*t)*dq2 + &
                       exp(1.d0/t)*(4.d0*q*q/(t*t*t*t))*mat + &
                       exp(1.d0/t)*(2.d0*(3.d0*q*q+1.d0)/(t*t*t))*mat
  else
    ddweightfunction = 0.d0
  endif
elseif (iprior .eq. 4) then  ! Gaussian RBF
  beta = 3.d0/(rmax(index)+small)**2
  t = exp(-beta*dot_product(rvec,rvec))
  ddweightfunction = -2.d0*t*beta*identity(size(rvec)) + &
                      4.d0*t*beta*beta*outerproduct(rvec,rvec)
endif

END FUNCTION ddweightfunction

!*************************************************************************
! SUBROUTINE setrmaxandpriorweight(rsupport,dmetric,priorweight)
! Purpose
! =======
! Define the support size and the prior weight function
!*************************************************************************
!
SUBROUTINE setrmaxandpriorweight(rsupport,dmetric,priorweight)

real(dp), dimension(:), target, intent(in)     :: rsupport
real(dp), dimension(:,:,:), target, intent(in) :: dmetric
character(*), intent(in)                       :: priorweight

rmax => rsupport
D => dmetric
prior = priorweight

select case (prior)
  case ('uniform')
    iprior = 0
  case ('cubic')         ! C^2 weight function
    iprior = 1
  case ('quartic')       ! C^2 weight function
    iprior = 2
  case ('gaussian')      ! C-infty weight function
    iprior = 3
  case ('gaussian-rbf')  ! Gaussian RBF
    iprior = 4
  case default
    write(*,*)"Weight not yet coded for prior = ", prior
    stop
end select

END SUBROUTINE setrmaxandpriorweight

!*************************************************************************
! FUNCTION getrmax()
! Purpose
! =======
! Return the support sizes of the nodal weight function
!*************************************************************************
!
FUNCTION getrmax()

real(dp), dimension(size(rmax)) :: getrmax

getrmax = rmax

END FUNCTION getrmax

!*************************************************************************
! FUNCTION getpriortype()
! Purpose
! =======
! Return the prior type
!*************************************************************************
!
FUNCTION getpriortype()

character(80) :: getpriortype

getpriortype = prior

END FUNCTION getpriortype

!*************************************************************************
! FUNCTION outerproduct(a,b)
! Purpose
! =======
! Defines the outer product of two vectors
!*************************************************************************
!

FUNCTION outerproduct(a,b)

real(dp), dimension(:), intent(in) :: a,b
real(dp), dimension(size(a),size(b)) :: outerproduct

outerproduct = spread(a,dim=2,ncopies=size(b)) * &  ! continuation
               spread(b,dim=1,ncopies=size(a))

END FUNCTION outerproduct

!*************************************************************************
! FUNCTION identity(m)
! Purpose
! =======
! Return an identity matrix
!*************************************************************************
!

FUNCTION identity(m)

integer, intent(in) :: m
real(dp), dimension(m,m) :: identity
integer :: i

identity = 0.d0
do i = 1,m
  identity(i,i) = 1.d0
enddo

END FUNCTION identity

END MODULE priorweightfunction


!*************************************************************************
!                         MAXENT MODULE                                  !
!                         =============                                  !
!                                                                        !
! Purpose                                                                !
! =======                                                                !
! Module for maximum-entropy basis function computations                 !
! Author: N. Sukumar, UC Davis                                           !
! Date  : June 2008 (Version 1.4)                                        !
! References                                                             !
! ==========                                                             !
! 1. N. Sukumar, IJNME, Vol. 61, Number 12, 2159-2181, 2004              !
! 2. M. Arroyo and M. Ortiz, IJNME, Vol. 65, Number 13, 2167-2202, 2006  !
! 3. N. Sukumar, AIP Conference Proceedings, Vol. 803, 337-344, 2005     !
! 4. N. Sukumar and R. W. Wright, IJNME, Vol. 70, Number 2, 181-205, 2007!
!                                                                        !
! Dependencies                                                           !
! ============                                                           !
! o PRIORWEIGHTFUNCTION module, which is provided with this distribution !
! o L-BFGS (public-domain Fortran77 code of Liu and Nocedal), which is   !
!   included with this distribution. The original L-BFGS Fortran source  !
!   is available at http://www.netlib.org/opt/lbfgs_um.shar and also at  !
!   http://www.ece.northwestern.edu/~nocedal/lbfgs.html                  !
! o Lapack library (dpotrf and dpotri are used to compute the inverse    !
!   of a matrix). These are needed only for problems in R^d, d > 3       !
!                                                                        !
! Manual                                                                 !
! ======                                                                 !
! Please read the manual (included in the tar bundle) for installation   !
! and execution of the program (standalone or as a library).             !
!                                                                        !
!*************************************************************************

MODULE maxent

implicit none     ! require explicit type declaration/checking
private           ! make everything private by default

public drivermaxent

! types, constants ....

integer, parameter:: dp=kind(0.d0)   ! double precision type

private  n, nsd, point, method, lambda, coord, A, maxiter, epsilon
private  charlen
integer  :: n, nsd                              ! no of nodes, space-dimension
character(80)                     :: method     ! string (descent/newton/lbfgs)
real(dp), pointer, dimension(:)   :: point      ! point in nsd-dimensions
real(dp), pointer, dimension(:,:) :: coord(:,:) ! nodal coordinates (matrix)
real(dp), allocatable :: A(:,:)                 ! x^i - x (matrix)
real(dp), allocatable :: DA(:,:)                ! D*(x^i - x) (matrix)
real(dp), allocatable :: qi(:)                  ! sqrt((x^i - x)^t D (x^i - x))/rmax^i
real(dp), allocatable :: lambda(:)              ! vector (Lagrange multipliers)
integer  :: maxiter                             ! maximum number of iterations
real(dp) :: epsilon                             ! convergence tolerance
real(dp) :: charlen                             ! characteristic distance
real(dp) :: det_tol = 1.d-16                    ! tolerance on determinant
real(dp) :: small = tiny(1.d0)                  ! tolerance to avoid overflow
integer  :: error_flag = 0                      ! check if determinant of hessian is small OR
                                                ! if nonconvergence results in `maxiter'
                                                ! iterations
logical :: consist_err                          ! pass (.false.) or no-pass (.true.) 
                                                ! consistency condition 

! functions, subroutines ....
CONTAINS

!*************************************************************************
! FUNCTION phimaxent()
! Purpose
! =======
! Compute the entropy basis functions
!*************************************************************************
!
FUNCTION phimaxent()

real(dp), dimension(n) :: phimaxent  ! vector of basis function values
real(dp) :: Z
real(dp), dimension(n) :: Zcomp

Zcomp = partitionfunctioncomponents()
Z = sum(Zcomp)
phimaxent = Zcomp/Z

END FUNCTION phimaxent

!*************************************************************************
! FUNCTION dphimaxent()
! Purpose
! =======
! Compute the gradient of the entropy basis functions
!*************************************************************************
!
FUNCTION dphimaxent()

USE priorweightfunction

real(dp), dimension(nsd,n) :: dphimaxent
real(dp), dimension(n) :: phi
real(dp), dimension(nsd) :: rv1, rv2, rvec, sumw
real(dp), dimension(nsd,n) :: priormat
real(dp), dimension(nsd,nsd) :: hess, mat, res
real(dp), dimension(nsd,nsd) :: invhess, onematrix
integer :: i

onematrix = identity(nsd)
phi = phimaxent()
hess = hessian(.true.)
invhess = invmatrix(hess)

sumw = 0.d0
do i = 1,n
! adding `small' so that divide by zero never occurs
  priormat(:,i) = dweightfunction(i,A(:,i),DA(:,i),qi(i))/(weightfunction(i,A(:,i),DA(:,i),qi(i))+small)
  sumw = sumw + phi(i)*priormat(:,i)
enddo

mat = 0.d0
do i = 1,n
  rv1 = A(:,i)
  rv2 = priormat(:,i)
  mat = mat + phi(i)*outerproduct(rv1,rv2)
enddo
res = invhess - matmul(invhess,mat)

do i = 1,n
  rv1 = A(:,i)
  dphimaxent(:,i) = phi(i)*(matmul(rv1,res) + priormat(:,i) - sumw)
enddo

END FUNCTION dphimaxent

!*************************************************************************
! FUNCTION ddphimaxent()
! Purpose
! =======
! Compute the second derivatives of the entropy basis functions
!*************************************************************************
!
FUNCTION ddphimaxent()

USE priorweightfunction

real(dp), dimension(nsd,nsd,n) :: ddphimaxent
real(dp), dimension(nsd,n) :: dphi
real(dp), dimension(n) :: phi
real(dp), dimension(nsd) :: rv1, rv2, rv3, rv4, rvec, v, w
real(dp), dimension(nsd,nsd) :: gradsumw
real(dp), dimension(nsd,nsd) :: dv2
real(dp), dimension(nsd,n) :: priormat
real(dp), dimension(nsd,nsd,n) :: gradpriormat
real(dp), dimension(nsd,nsd) :: mat, res
real(dp), dimension(nsd) :: mat2
real(dp), dimension(nsd,nsd,nsd) :: gradmat, ghinv, mat1, mat3, C
real(dp), dimension(nsd,nsd) :: invhess, onematrix, ddw
real(dp)  :: wt
real(dp), dimension(nsd) :: dwt
real(dp), dimension(nsd,nsd) :: ddwt
integer :: i

onematrix = identity(nsd)
phi = phimaxent()
dphi = dphimaxent()
invhess = invhessian(.true.)

gradsumw = 0.d0
do i = 1,n
  wt   = ( weightfunction(i,A(:,i),DA(:,i),qi(i)) + 1.d-30) ! adding `small = tiny(1.d0)'
  dwt  =  dweightfunction(i,A(:,i),DA(:,i),qi(i))
  ddwt = ddweightfunction(i,A(:,i),DA(:,i),qi(i))
  priormat(:,i) = dwt/wt
  gradpriormat(:,:,i) = ddwt/wt - outerproduct(priormat(:,i),priormat(:,i))
  gradsumw = gradsumw + outerproduct(priormat(:,i),dphi(:,i)) &
             + phi(i)*gradpriormat(:,:,i)            
enddo

mat = 0.d0
do i = 1,n
  rv1 = A(:,i)
  rv2 = priormat(:,i)
  mat = mat + phi(i)*outerproduct(rv1,rv2)
enddo

mat1 = 0.d0
mat2 = 0.d0
mat3 = 0.d0
C = 0.d0
do i = 1,n
  rv1 = A(:,i)
  v = matmul(rv1,invhess)
  w = matmul(v,mat)
  rv2 = dphi(:,i)
  C = C + outerproduct3(v,w,rv2)
  rv3 = priormat(:,i)
  mat1 = mat1 + outerproduct3(rv1,rv3,rv2)
  mat2 = mat2 + phi(i)*rv3
  mat3 = mat3 + phi(i)*outerproductvecmat(rv1,gradpriormat(:,:,i))
enddo

res = transpose(matmul(invhess,mat)) - invhess
ghinv = gradhessianinv()

!
! NOTE: For more speed-up, can just do upper-half, since the matrix is symmetric; ditto in
!       functions such as vectensor3mul
!
do i = 1,n
  rv1 = A(:,i)
  v = matmul(rv1,invhess)
  ddphimaxent(:,:,i) = outerproduct(dphi(:,i),dphi(:,i))/(phi(i)+small) &
                       + phi(i)*( res &
                       + gradpriormat(:,:,i) &
                       + vectensor3mul(rv1,ghinv+C) &  
                       - vectensor3mul(v,mat1+mat3) &
                       + outerproduct(mat2,v) &
                       - gradsumw)
enddo

END FUNCTION ddphimaxent

!*************************************************************************
! FUNCTION partitionfunction()
! Purpose
! =======
! Compute the partition function
!*************************************************************************
!
FUNCTION partitionfunction()

real(dp) :: partitionfunction ! scalar output
real(dp), dimension(n) :: Zwithid

Zwithid = partitionfunctioncomponents()
partitionfunction = sum(Zwithid)

END FUNCTION partitionfunction

!*************************************************************************
! FUNCTION partitionfunctioncomponents()
! Purpose
! =======
! Compute the partition function components (vector)
!*************************************************************************
!
FUNCTION partitionfunctioncomponents()

USE priorweightfunction

real(dp), dimension(n) :: partitionfunctioncomponents ! components of Z
real(dp) w                                            ! prior weight
real(dp), dimension(n) :: lambdaAid
real(dp) rvec(nsd), r
integer :: i

lambdaAid = matmul(lambda,A)
do i = 1,n
  w = weightfunction(i,A(:,i),DA(:,i),qi(i))
  partitionfunctioncomponents(i) = w*exp(-lambdaAid(i))
enddo

END FUNCTION partitionfunctioncomponents

!*************************************************************************
! FUNCTION entropy()
! Purpose
! =======
! Compute the entropy of the basis functions
!*************************************************************************
!
FUNCTION entropy()

real(dp), dimension(n) :: phi
real(dp) :: entropy
real(dp), parameter :: TOL = 1.d-16
integer :: i

phi = phimaxent() 
entropy = 0.0
do i = 1,n
   if (phi(i) > TOL) then
      entropy = entropy - phi(i)*log(phi(i))
   endif
enddo

END FUNCTION entropy

!*************************************************************************
! FUNCTION func()
! Purpose
! =======
! Compute the objective function (inf f is solved)
!*************************************************************************
!
FUNCTION func()

real(dp) :: func !  objective f
real(dp) :: Z

Z = partitionfunction()
func = log(Z) 

END FUNCTION func

!*************************************************************************
! FUNCTION dfunc()
! Purpose
! =======
! Compute the gradient of the objective function
!*************************************************************************
!
FUNCTION dfunc()

real(dp), dimension(nsd) :: dfunc ! grad of f 
real(dp), dimension(n) :: phi

phi   = phimaxent()
dfunc = -matmul(A,phi) ! note the negative sign (see definition of Z)

END FUNCTION dfunc

!*************************************************************************
! FUNCTION hessian(flag)
! Purpose
! =======
! Compute the Hessian 
!*************************************************************************
!
FUNCTION hessian(flag)

logical, optional, intent(in) :: flag
real(dp), dimension(nsd) :: g
real(dp), dimension(nsd,nsd) :: hessian, p2
real(dp), dimension(nsd,nsd) :: mat3d
real(dp), dimension(n) :: phi
integer :: i, j, k

hessian = 0.d0
phi = phimaxent()

! ... compute upper diagonal of hessian
do i = 1,nsd
 do j = i,nsd
  do k = 1,n
    hessian(i,j) = hessian(i,j) + phi(k)*A(i,k)*A(j,k)
	enddo
 enddo
enddo
! ... compute lower diagonal
do i = 1,nsd
 do j = i+1,nsd
  hessian(j,i) = hessian(i,j)
 enddo
enddo

if (.not. present(flag)) then
  g = dfunc()
  p2 = outerproduct(g,g)
  hessian = hessian - p2
endif

END FUNCTION hessian

!*************************************************************************
! FUNCTION hessiancomp(comp)
! Purpose
! =======
! Compute the nodal contribution to the converged Hessian matrix
!*************************************************************************
!
FUNCTION hessiancomp(comp)

integer, intent(in) :: comp
real(dp), dimension(nsd,nsd) :: hessiancomp
real(dp), dimension(nsd) :: rvec
real(dp), dimension(n) :: phi
integer :: i, j

if (comp < 0 .or. comp > n) then
  write(6,*)"Invalid comp value . . aborting"
  stop
endif
  
phi = phimaxent()

rvec = A(:,comp)
hessiancomp = phi(comp)*outerproduct(rvec,rvec)

END FUNCTION hessiancomp

!*************************************************************************
! FUNCTION gradhessian()
! Purpose
! =======
! Compute the gradient of the Hessian 
!*************************************************************************

FUNCTION gradhessian()

real(dp), dimension(nsd,nsd,nsd) :: gradhessian
real(dp), dimension(nsd) :: vec1, vec2
real(dp), dimension(nsd,n) :: dphi
integer :: i

dphi = dphimaxent()
gradhessian = 0.d0
do i = 1,n
  vec1 = dphi(:,i)
  vec2 = A(:,i)
  gradhessian = gradhessian + outerproduct3(vec2,vec2,vec1)
enddo

END FUNCTION gradhessian

!*************************************************************************
! FUNCTION gradhessianinv()
! Purpose
! =======
! Compute the gradient of the inverse of the Hessian 
!*************************************************************************
!
FUNCTION gradhessianinv()

real(dp), dimension(nsd,nsd,nsd) :: gradhessianinv
real(dp), dimension(nsd,nsd) :: invhess, mat
real(dp), dimension(nsd) :: vec1, vec2
real(dp), dimension(nsd,n) :: dphi
integer :: i

dphi = dphimaxent()
invhess = invhessian(.true.)
gradhessianinv = 0.d0
do i = 1,n
  vec2 = matmul(invhess,A(:,i))
  mat = outerproduct(vec2,vec2)
  vec1 = dphi(:,i)
  gradhessianinv = gradhessianinv - outerproductmatvec(mat,vec1)
enddo

END FUNCTION gradhessianinv

!*************************************************************************
! FUNCTION invhessian(flag)
! Purpose
! =======
! Compute the inverse of the Hessian 
!*************************************************************************
!
FUNCTION invhessian(flag)

logical, optional, intent(in) :: flag
real(dp), dimension(nsd,nsd) :: mat
real(dp), dimension(nsd,nsd) :: invhessian
real(dp) :: det

if (.not. present(flag)) then
  mat = hessian()
else
  mat = hessian(flag)
endif
invhessian = invmatrix(mat)

END FUNCTION invhessian

!*************************************************************************
! FUNCTION invmatrix(mat)
! Purpose
! =======
! Compute the inverse of a matrix
!*************************************************************************
!
FUNCTION invmatrix(mat)

real(dp), dimension(:,:), intent(in) :: mat
real(dp), dimension(size(mat,DIM=1),size(mat,DIM=1)) :: invmatrix
real(dp) :: det
integer :: i, j, info, m
integer :: nn

m = size(mat,DIM=1)
INV: if (m == 1) then
  invmatrix(1,1) = 1./mat(1,1)
  det = invmatrix(1,1)
elseif (m == 2) then
  det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
!
  if (abs(det) .ge. charlen**(2*nsd)*det_tol) then
    invmatrix(1,1) = mat(2,2)/det
    invmatrix(1,2) = -mat(1,2)/det
    invmatrix(2,1) = -mat(2,1)/det
    invmatrix(2,2) = mat(1,1)/det
  endif
elseif (m == 3) then
  det = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
        - mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
        + mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))

  if (abs(det) .ge. charlen**(2*nsd)*det_tol) then
!
    invmatrix(1,1) = (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2))/det
    invmatrix(1,2) = (mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3))/det
    invmatrix(1,3) = (mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2))/det
! 
    invmatrix(2,1) = (mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3))/det
    invmatrix(2,2) = (mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1))/det
    invmatrix(2,3) = (mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3))/det
! 
    invmatrix(3,1) = (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))/det
    invmatrix(3,2) = (mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2))/det
    invmatrix(3,3) = (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1))/det

  endif
else
#ifdef USELAPBLAS
  call dpotrf('u',nsd,mat,nsd,info)
  if (info .ne. 0) then
    write(6,*)'Error in dpotrf(): Flag is ', info
    error_flag = 1
  endif
  call dpotri('u',nsd,mat,nsd,info)
  if (info .ne. 0) then
    write(6,*)'Error in dpotri(): Flag is ', info
    error_flag = 1
  endif
  invmatrix = mat
  do i = 1,nsd
    do j = 1,i-1
      invmatrix(i,j) = invmatrix(j,i)
    enddo
  enddo
#endif
endif INV

if (m .lt. 4) then
  if (abs(det) < charlen**(2*nsd)*det_tol) then
    write(6,*)"WARNING: Matrix is nearly singular"
    write(6,*)"Lambda = ", lambda
    write(6,*)"Hessian = ", mat
    write(6,*)"Det(H) ", det
    error_flag = 1
  endif
endif

END FUNCTION invmatrix

!*************************************************************************
! FUNCTION newtondirection(schemedir,iter)
! Purpose
! =======
! Compute the direction vector for Newton method or steepest descent 
!*************************************************************************
!
FUNCTION newtondirection(schemedir,iter)

character(80), intent(in) :: schemedir      ! scheme for direction
integer, intent(in) :: iter                 ! iteration number
real(dp), dimension(nsd) :: g               ! nonlinear eqs [g_i = 0]
real(dp), dimension(nsd,nsd) :: invh        ! inverse of the hessian
real(dp), dimension(nsd) :: newtondirection ! direction vector


select case (schemedir)
  case ('newton')
    g = dfunc() 
    invh = invhessian()
    if (error_flag .eq. 1) return      ! bad det in inverse
    newtondirection = -matmul(invh,g)  ! note the negative sign
  case ('descent')
    g = dfunc()  
    newtondirection = -g/sqrt(dot_product(g,g)) ! steepest descent
  case default
    write(*,*)"Scheme = ",schemedir, " not valid in newtondirection()"
    stop
end select

END FUNCTION newtondirection

!*************************************************************************
! SUBROUTINE drivermaxent_descent(maxentprint,phi)
! Purpose
! =======
! Driver subroutine for maxent basis functions
!*************************************************************************
!
SUBROUTINE drivermaxent_descent(maxentprint,phi)

logical, intent(in) :: maxentprint
real(dp), intent(out) :: phi(n)

real(dp), dimension(nsd) :: dlambda, lambda_t 
real(dp) :: errornorm
real(dp) :: alpha
integer :: iter, i

iter   = 0
lambda = 0.0                        ! initial guess
dlambda = dfunc()
errornorm = sqrt(Dot_Product(dlambda,dlambda))

do while (errornorm/charlen > epsilon .and. iter < maxiter)

  if (maxentprint) write(*,*)"ITER and ERROR: ", iter, errornorm/charlen

  lambda_t = lambda   ! save the starting lambda

  dlambda = -dlambda/errornorm   ! steepest descent direction

  alpha = stepsize(dlambda)
  lambda = lambda_t + alpha*dlambda
  dlambda = dfunc()
  errornorm = sqrt(Dot_Product(dlambda,dlambda))
  iter = iter + 1
enddo

if (iter .ge. maxiter) then
  error_flag = 2
  return
endif

phi = phimaxent()

if (errornorm/charlen <= epsilon) then
  if (maxentprint) then
    write(*,*)
    write(*,*)"*****************************************************"
    write(*,*)"*************** STEEPEST DESCENT METHOD *************"
    write(*,*)"*****************************************************"
    write(*,*)"POINT = ", point
    write(*,*)"CONVERGENCE ATTAINED IN ITERATIONS = ", iter
    write(*,*)"ASKING TOLERANCE = ", epsilon
    write(*,*)"ERROR = ",errornorm/charlen
    write(*,*)"LAGRANGE MULTIPLIERS"
    write(*,'(i5,1x,1pe20.13)')(i,lambda(i),i=1,size(lambda))
    write(*,*)"BASISFUNCTIONS"
    write(*,'(i5,1x,1pe20.13)')(i,phi(i),i=1,size(phi))
  endif
else
  if (maxentprint) then
    write(*,*)"CONVERGENCE NOT ATTAINED IN ITERATIONS = ",maxiter
    write(*,'("POINT = ",3(1pe20.13))') point
    write(*,*)"ASKING TOLERANCE = ", epsilon
    write(*,*)"ERROR = ",errornorm/charlen
    call checkconsistency
    stop
  endif
endif

END SUBROUTINE drivermaxent_descent

!*************************************************************************
! SUBROUTINE drivermaxent_newton(maxentprint,phi)
! Purpose
! =======
! Driver subroutine for maxent basis functions using Newton method
!*************************************************************************
!
SUBROUTINE drivermaxent_newton(maxentprint,phi)

logical, intent(in) :: maxentprint
real(dp), intent(out) :: phi(n)

real(dp), dimension(nsd) :: dlambda, lambda_t, g
real(dp) :: alpha
real(dp) :: errornorm
character(80) :: schemedir
integer :: iter, i

iter   = 0
lambda = 0.d0                             ! initial guess
g = dfunc()
errornorm = sqrt(Dot_Product(g,g))
schemedir = 'newton'

do while (errornorm/charlen > epsilon .and. iter < maxiter)

  if (maxentprint) write(*,*)"ITER and ERROR: ", iter, errornorm/charlen

  lambda_t = lambda
  dlambda = newtondirection(schemedir,iter)

  if (error_flag .eq. 1) return  ! bad det in inverse

  if (errornorm/charlen < 1.d-4) then
    alpha = 1.d0
  else
    alpha = stepsize(dlambda)
  endif
  lambda = lambda_t + alpha*dlambda   

  g = dfunc()
  errornorm = sqrt(Dot_Product(g,g))
  iter = iter + 1

enddo

if (iter .ge. maxiter) then
  error_flag = 2
  return
endif

phi = phimaxent()

if (errornorm/charlen <= epsilon) then
  if (maxentprint) then
    write(*,*)
    write(*,*)"*****************************************************"
    write(*,*)"***************** NEWTON METHOD *********************"
    write(*,*)"*****************************************************"
    write(*,*)"POINT = ", point
    write(*,*)"CONVERGENCE ATTAINED IN ITERATIONS = ", iter
    write(*,*)"ASKING TOLERANCE = ", epsilon
    write(*,*)"ERROR = ",errornorm/charlen
    write(*,*)"LAGRANGE MULTIPLIERS"
    write(*,'(i5,1x,1pe20.13)')(i,lambda(i),i=1,size(lambda))
    write(*,*)"BASISFUNCTIONS"
    write(*,'(i5,1x,1pe20.13)')(i,phi(i),i=1,size(phi))
  endif
else
  if (maxentprint) then
    write(*,*)"CONVERGENCE NOT ATTAINED IN ITERATIONS = ",maxiter
    write(*,'("POINT = ",3(1pe20.13))') point
    write(*,*)"ASKING TOLERANCE = ", epsilon
    write(*,*)"ERROR = ",errornorm/charlen
    call checkconsistency
    stop
  endif
endif

END SUBROUTINE drivermaxent_newton

!*************************************************************************
! FUNCTION stepsize(dlambda)
! Purpose
! =======
! Compute the step size for gradient descent or guarded Newton method
! Reference: Burden and Faires, Numerical Analysis, 8th Edition, Thomson
! Brooks/Cole, 2005.
!*************************************************************************
!
FUNCTION stepsize(dlambda)

real(dp), intent(in) :: dlambda(nsd)
real(dp) :: stepsize

real(dp), dimension(nsd) :: lambda_t 
real(dp) :: errornorm
real(dp) :: alpha, alpha1, alpha2, alpha3, g0, g1, g2, g3, h1, h2, h3
integer  :: i


alpha3 = 1.0
lambda_t = lambda
g1 = func()
lambda = lambda_t + alpha3*dlambda
g3 = func()

do while (g3 > g1 .and. alpha3 > 0.5*epsilon) 
  alpha3 = 0.5*alpha3
  lambda = lambda_t + alpha3*dlambda 
  g3 = func()
enddo

if (alpha3 > 0.5*epsilon) then
  alpha2 = 0.5*alpha3
  lambda = lambda_t + alpha2*dlambda 
  g2 = func()
  h1 = (g2 - g1)/alpha2
  h2 = (g3 - g2)/(alpha3 - alpha2)
  h3 = (h2 - h1)/alpha3
  if (h3 /= 0.) then
    alpha1 = 0.5*(alpha2 - h1/h3)
  else
    alpha1 = alpha3;
  endif
  lambda = lambda_t + alpha1*dlambda 
  g0 = func()
  if (g0 < g3) then
    alpha = alpha1;
  else
    alpha = alpha3;
  endif
else
  alpha = alpha3
endif

stepsize = alpha

END FUNCTION stepsize

!*************************************************************************
! SUBROUTINE drivermaxent_lbfgs(maxentprint,phi)
! Purpose
! =======
! Driver subroutine for maxent basis functions using LBFGS
! Reference: Liu and Nocedal, On the limited memory BFGS method for
! large scale optimization, Mathematical Programming, 45, 1989,
! pp. 503-528
!*************************************************************************
!
SUBROUTINE drivermaxent_lbfgs(maxentprint,phi)

logical, intent(in) :: maxentprint
real(dp), dimension(n), intent(out) :: phi

integer, parameter :: MSAVE = 7
integer :: NWORK
real(dp), parameter :: MACH_XTOL = 1.e-16  ! machine precision
integer :: m   ! should be between 3 and 7
real(dp) :: f
real(dp), dimension(nsd) :: g, diag
real(dp), allocatable :: w(:)
logical  :: diagco
integer, dimension(2) :: iprint
integer :: iflag, icall
integer :: iter
real(dp) :: errornorm
!
! common block
!
integer :: mp, lp
real(dp) :: gtol, stpmin, stpmax
external lb2
common /lb3/ mp, lp, gtol, stpmin, stpmax
integer i

NWORK = nsd*(2*MSAVE+1)+2*MSAVE
allocate(w(NWORK))

m = 5            ! limited memory algorithm kicks in after m iterations

iprint(1) = -1   ! negative values => do not print out information
iprint(2) = 0    ! specific information is output if iprint(1) >= 0

diagco = .false. ! no inverse hessian is supplied 
                 ! if diagco = .true., diag must contain the diagonals of the
                 ! inverse hessian

lambda = 0.0 ! initial guess
f = func()
g = dfunc()
errornorm = sqrt(dot_product(g,g))
iter = 0
iflag = 0  ! iflag = 0 => convergence, < 0 => error, > 0 => continue looping


do while (errornorm/charlen > epsilon)
  call lbfgs(nsd,m,lambda,f,g,diagco,diag,iprint,epsilon,MACH_XTOL,w,iflag)
  if (maxentprint) write(*,*)"ITER ERROR IFLAG: ", iter, errornorm, iflag
  iter = iter + 1
  if (iflag == 1) then
    f = func()
    g = dfunc()
  endif
  errornorm = sqrt(dot_product(dfunc(),dfunc()))
  if (iflag == 0 .or. errornorm/charlen <= epsilon .or. iter >= maxiter) exit
enddo

if (iflag < 0) then
  error_flag = 1
  return
endif
if (iter >= maxiter) then
  error_flag = 2
  return
endif

phi = phimaxent()
deallocate(w)

if (iflag == 0 .or. errornorm/charlen <= epsilon) then
  if (maxentprint) then
    write(*,*)
    write(*,*)"*****************************************************"
    write(*,*)"******************** L-BFGS METHOD ******************"
    write(*,*)"*****************************************************"
    write(*,*)"POINT = ", point
    write(*,*)"CONVERGENCE ATTAINED IN ITERATIONS = ", iter
    write(*,*)"ASKING TOLERANCE = ", epsilon
    write(*,*)"ERROR = ",errornorm/charlen
    write(*,*)"LAGRANGE MULTIPLIERS"
    write(*,'(i5,1x,1pe20.13)')(i,lambda(i),i=1,size(lambda))
    write(*,*)"BASISFUNCTIONS"
    write(*,'(i5,1x,1pe20.13)')(i,phi(i),i=1,size(phi))
  endif
else
  if (maxentprint) then
    write(*,*)"CONVERGENCE NOT ATTAINED IN ITERATIONS = ",maxiter
    write(*,'("POINT = ",3(1pe20.13))') point
    do i = 1,n
      write(*,*)'Node:',i,' Coord = ',coord(i,:)
    enddo
    write(*,*)"ASKING TOLERANCE = ", epsilon
    write(*,*)"ERROR = ",errornorm/charlen
    write(*,*)"IFLAG in LBFGS = ", iflag
    stop
  endif
endif

END SUBROUTINE drivermaxent_lbfgs

!*************************************************************************
! SUBROUTINE drivermaxent(ndim,nsddim,scheme,priorwt,xyz,p,dmax,dmetric,&
!                         maxit,eps,maxentprint,ierror,charlengthscale,&
!                         idet_tol,phi,dphi,ddphi)
! Purpose
! =======
! Driver subroutine for maxent basis functions, which serves as a general
! purpose interface to the MAXENT library
! 
!*************************************************************************
!
SUBROUTINE drivermaxent(ndim,nsddim,scheme,priorwt,xyz,p,dmax,dmetric,&
                        maxit,eps,maxentprint,ierror,cerror,charlengthscale,&
                        dettol,phi,dphi,ddphi)

USE priorweightfunction

integer, intent(in) :: ndim, nsddim
integer :: ierror
character(80), intent(in) :: scheme, priorwt
real(dp), intent(in) :: xyz(ndim,nsddim), p(nsddim)
real(dp), intent(in) :: dmax(ndim)
real(dp), intent(in) :: dmetric(nsddim,nsddim,ndim)
integer, intent(in)  :: maxit
real(dp), intent(in) :: eps
logical, intent(in)  :: maxentprint
logical, intent(out) :: cerror
real(dp), optional, intent(in)    :: charlengthscale, dettol
real(dp), optional, intent(out) :: phi(ndim)
real(dp), optional, intent(out) :: dphi(nsddim,ndim)
real(dp), optional, intent(out) :: ddphi(nsddim,nsddim,ndim)

integer :: i, j

!
! set `rmax', 'metric', and `prior' in the priorweightfunction module
!
call setrmaxandpriorweight(dmax,dmetric,priorwt)

n = ndim
nsd = nsddim
method = scheme
maxiter = maxit
epsilon = eps
ierror = 0
cerror = .false.
!
! set point, coord and allocate A matrix
!
call allocateandsetarrays(p,xyz,dmetric,dmax)

!
! determine characteristic length, which is needed to establish convergence criterion that is
! independent of specimen dimensions
!
if (present(charlengthscale)) then
  charlen = charlengthscale
else
  charlen = sum(getrmax())/ndim ! use average of the support sizes to make 
                                ! results independent of units
endif
if (present(dettol)) det_tol = dettol

if (maxentprint) then
  if (nsd == 1) write(*,'("POINT = ",(1pe13.6),/)')point
  if (nsd == 2) write(*,'("POINT = ",2(1pe13.6),/)')point
  if (nsd == 3) write(*,'("POINT = ",3(1pe13.6),/)')point
  if (nsd > 3) write(*,*)"POINT = ",point
  do i = 1,n
    if (nsd == 1) then
      write(*,'("NODE ",i5," : COORD = ",1(1pe13.6)," RMAX ",1pe13.6)')i,coord(i,:),dmax(i)
    else if (nsd == 2) then
      write(*,'("NODE ",i5," : COORD = ",2(1pe13.6,1x)," RMAX ",1pe13.6)')i,coord(i,:),dmax(i)
    else if (nsd == 3) then
      write(*,'("NODE ",i5," : COORD = ",3(1pe13.6,1x)," RMAX ",1pe13.6)')i,coord(i,:),dmax(i)
    else
      write(*,'("NODE ",i5," : COORD = ",4(1pe13.6,1x)," RMAX ",1pe13.6)')i,coord(i,:),dmax(i)
    endif
  enddo
  write(*,'(/,"PRIOR WEIGHTFUNCTION TYPE = ",a20,/)')getpriortype()
endif

select case (method)
  case ('descent')
    call drivermaxent_descent(maxentprint,phi)
  case ('lbfgs')
    call drivermaxent_lbfgs(maxentprint,phi)
  case ('newton')
    call drivermaxent_newton(maxentprint,phi)
  case default
    write(*,*)"Solution scheme not yet coded for scheme = ", method
    stop
end select

!
! determinant in inv hessian is too small or nonconvergence in `maxiter' iterations
!
if (error_flag .ne. 0) then
  ierror = error_flag
  call deallocatemaxentarrays
  return
endif

if (present(dphi)) then
  if (.not. present(phi)) then
    write(6,*)'Error: Must include phi in the calling subroutine too'
    stop
  endif
  dphi = dphimaxent()
  if (maxentprint) then
    call checkdconsistency()
    write(*,*)"D-BASISFUNCTIONS"
    do i = 1,n
      if (nsd == 1) write(*,'(i5,1x,1(1pe20.13,1x))')i,dphi(1:nsd,i)
      if (nsd == 2) write(*,'(i5,1x,2(1pe20.13,1x))')i,dphi(1:nsd,i)
      if (nsd == 3) write(*,'(i5,1x,3(1pe20.13,1x))')i,dphi(1:nsd,i)
      if (nsd > 3) write(*,'(i5,1x,4(1pe20.13,1x))')i,dphi(1:nsd,i)
    enddo
  endif
endif

if (present(ddphi)) then
  if (.not. present(phi)) then
    write(6,*)'Error: Must include phi in the calling subroutine too'
    stop
  endif
  ddphi = ddphimaxent()
  if (maxentprint) then
    call checkddconsistency()
    write(*,*)"DD-BASISFUNCTIONS"
    do i = 1,n
      if (nsd == 1) then
        write(*,'(i5,1x,1(1pe20.13))')i,ddphi(1,1:nsd,i)
      elseif (nsd == 2) then
        write(*,'(i5,1x,2(1pe20.13,1x))')i,ddphi(1,1:nsd,i)
        write(*,'(6x,2(1pe20.13,1x))')ddphi(2,1:nsd,i)
      elseif (nsd == 3) then
        write(*,'(i5,1x,3(1pe20.13,1x))')i,ddphi(1,1:nsd,i)
        write(*,'(6x,3(1pe20.13,1x))')ddphi(2,1:nsd,i)
        write(*,'(6x,3(1pe20.13,1x))')ddphi(3,1:nsd,i)
      else
        write(*,'(i5,1x,4(1pe20.13,1x))')i,ddphi(1,1:nsd,i)
        write(*,'(6x,4(1pe20.13,1x))')ddphi(2,1:nsd,i)
        write(*,'(6x,4(1pe20.13,1x))')ddphi(3,1:nsd,i)
        write(*,'(6x,4(1pe20.13,1x))')ddphi(4,1:nsd,i)
      endif
    enddo
  endif
endif

call deallocatemaxentarrays

END SUBROUTINE drivermaxent

!*************************************************************************
! SUBROUTINE allocateandsetarrays(p,xyz,D,rmax)
! Purpose
! =======
! Allocate and set private variables
!
!*************************************************************************
!
SUBROUTINE allocateandsetarrays(p,xyz,D,rmax)

real(dp), target, intent(in) :: p(nsd), xyz(n,nsd),D(nsd,nsd,n), rmax(n)

integer :: i

point => p
coord => xyz

allocate(A(nsd,n))
allocate(DA(nsd,n))
allocate(qi(n))
allocate(lambda(nsd))

do i = 1,n
  A(:,i)  = coord(i,:) - point(:)
  DA(:,i) = matmul(D(:,:,i),A(:,i))
  qi(i) = sqrt(dot_product(A(:,i),DA(:,i)))/rmax(i)
enddo

END SUBROUTINE allocateandsetarrays

!*************************************************************************
! SUBROUTINE deallocatemaxentarrays()
! Purpose
! =======
! Deallocate memory of private variables
!
!*************************************************************************
!
SUBROUTINE deallocatemaxentarrays()

deallocate(A)
deallocate(DA)
deallocate(qi)
deallocate(lambda)

END SUBROUTINE deallocatemaxentarrays

!*************************************************************************
! SUBROUTINE checkconsistency()
! Purpose
! =======
! Check the consistency of the basis functions
!*************************************************************************
!
SUBROUTINE checkconsistency()

real(dp), parameter :: TOL = 1.d-10
real(dp), dimension(n) :: phi
real(dp), dimension(nsd) :: sums
real(dp) :: sumphi

integer :: i

phi = phimaxent()
sumphi = sum(phi)
sums = matmul(A,phi)

if (dabs(sumphi-1.0) > TOL) then
   write(6,*)'WARNING: PU consistency error:', point
   write(6,*)"SUMPHI - 1 = ",sumphi-1.0
   write(6,*)"PHI = ", phi
   consist_err = .true.
endif

do i = 1,nsd
  if (dabs(sums(i)) > TOL) then
     write(6,*)'WARNING: Linear consistency error:', point
     write(6,*)" i = ",i," and sums(i) = ",sums(i)
     write(6,*)"PHI = ", phi
     consist_err = .true.
  endif
enddo

END SUBROUTINE checkconsistency

!*************************************************************************
! SUBROUTINE checkdconsistency()
! Purpose
! =======
! Check the consistency of the basis function derivatives
!*************************************************************************
!
SUBROUTINE checkdconsistency()

real(dp), parameter :: TOL = 1.d-8
real(dp), dimension(nsd,n) :: dphi
real(dp), dimension(nsd) :: sumders
real(dp), dimension(nsd,nsd) :: mat
real(dp) :: val

integer :: i, j

dphi = dphimaxent()
sumders = 0.d0
do i = 1,n
  sumders = sumders + dphi(:,i)
enddo

val = sum(dabs(sumders))
if (val > TOL) then
  write(6,*)'WARNING: PU deriv. const. consistency error:', point
  write(6,*)"sum of dabs(sumders) = ", val
  write(6,*)"dphi = ", dphi
  consist_err = .true.
endif

mat = 0.d0
do i = 1,n
  mat = mat + outerproduct(dphi(:,i),A(:,i))
enddo
mat = mat - identity(nsd)
val = sum(dabs(mat))

if (val > TOL) then
  write(6,*)'WARNING: PU deriv linear consistency error:', point
  write(6,*)"sum of dabs(dphix) = ", val
  write(6,*)"dphi = ", dphi
  consist_err = .true.
endif

END SUBROUTINE checkdconsistency

!*************************************************************************
! SUBROUTINE checkddconsistency()
! Purpose
! =======
! Check the consistency of the Hessian of the basis functions
!*************************************************************************
!
SUBROUTINE checkddconsistency()

real(dp), parameter :: TOL = 1.d-6
real(dp), dimension(nsd,nsd,n) :: ddphi
real(dp), dimension(nsd,nsd) :: sumhess
real(dp), dimension(nsd,nsd,nsd) :: sumhessx
real(dp) :: val

integer :: i, j

ddphi = ddphimaxent()
sumhess = 0.d0
do i = 1,n
  sumhess = sumhess + ddphi(:,:,i)
enddo
val = sum(dabs(sumhess))

if (val > TOL) then
  write(6,*)'WARNING: PU hessian const. consistency error:', point
  write(6,*)"sum of dabs(sumhess) = ", val
  write(6,*)"ddphi = ", ddphi
  consist_err = .true.
endif

sumhessx = 0.d0
do i = 1,n
   sumhessx = sumhessx + outerproductvecmat(A(:,i),ddphi(:,:,i))
enddo
val = sum(dabs(sumhessx))

if (val > TOL) then
  write(6,*)'WARNING: PU hessian linear consistency error:', point
  write(6,*)"sum of dabs(sumhessx) = ", val
  write(6,*)"ddphi = ", ddphi
  consist_err = .true.
endif

END SUBROUTINE checkddconsistency

!*************************************************************************
! FUNCTION outerproduct(a,b)
! Purpose
! =======
! Defines the outer product of two vectors
!*************************************************************************
!

FUNCTION outerproduct(a,b)

real(dp), dimension(:), intent(in) :: a, b
real(dp), dimension(size(a),size(b)) :: outerproduct

outerproduct = spread(a,dim=2,ncopies=size(b)) * &  ! continuation
               spread(b,dim=1,ncopies=size(a))

END FUNCTION outerproduct

!*************************************************************************
! FUNCTION outerproduct3(a,b,c)
! Purpose
! =======
! Defines the outer product of two vectors
!*************************************************************************
!

FUNCTION outerproduct3(a,b,c)

real(dp), dimension(:), intent(in) :: a, b, c
real(dp), dimension(size(a),size(b),size(c)) :: outerproduct3
integer :: i, j, k

do i = 1,size(a)
  do j = 1,size(b)
    do k = 1,size(c)
      outerproduct3(i,j,k) = a(i)*b(j)*c(k)
    Enddo
  enddo
enddo

END FUNCTION outerproduct3

!*************************************************************************
! FUNCTION outerproductvecmat(vec,mat,flip)
! Purpose
! =======
! Defines the outer product of a vector and a matrix (2-tensor)
!*************************************************************************
!

FUNCTION outerproductvecmat(vec,mat,flip)

logical, optional, intent(in) :: flip
real(dp), dimension(:), intent(in) :: vec
real(dp), dimension(:,:), intent(in) :: mat
real(dp), dimension(size(vec),size(vec),size(vec)) :: outerproductvecmat
integer :: i, j, k

do i = 1,size(vec)
  do j = 1,size(vec)
    do k = 1,size(vec)
      if (present(flip)) then
        outerproductvecmat(i,j,k) = vec(j)*mat(i,k)
      else
        outerproductvecmat(i,j,k) = vec(i)*mat(j,k)
      endif
    enddo
  enddo
enddo

END FUNCTION outerproductvecmat

!*************************************************************************
! FUNCTION outerproductmatvec(mat,vec,flip)
! Purpose
! =======
! Defines the outer product of a matrix (2-tensor) and a vector
!*************************************************************************
!

FUNCTION outerproductmatvec(mat,vec,flip)

logical, optional, intent(in) :: flip
real(dp), dimension(:,:), intent(in) :: mat
real(dp), dimension(:), intent(in) :: vec
real(dp), dimension(size(vec),size(vec),size(vec)) :: outerproductmatvec
integer :: i, j, k

do i = 1,size(vec)
  do j = 1,size(vec)
    do k = 1,size(vec)
      if (present(flip)) then
        outerproductmatvec(i,j,k) = mat(i,k)*vec(j)
      else
        outerproductmatvec(i,j,k) = mat(i,j)*vec(k)
      endif
    enddo
  enddo
enddo

END FUNCTION outerproductmatvec

!*************************************************************************
! FUNCTION tensor3vecmul(mat,vec)
! Purpose
! =======
! Defines the dot product of a 3-tensor and a vector
!*************************************************************************
!

FUNCTION tensor3vecmul(mat,vec)

real(dp), dimension(:,:,:), intent(in) :: mat
real(dp), dimension(:), intent(in) :: vec
real(dp), dimension(size(vec),size(vec)) :: tensor3vecmul
integer :: i, j, k

tensor3vecmul = 0.d0
do i = 1,size(vec)
  do j = 1,size(vec)
    do k = 1,size(vec)
      tensor3vecmul(i,j) = tensor3vecmul(i,j) + mat(i,j,k)*vec(k)
    enddo
  enddo
enddo

END FUNCTION tensor3vecmul

!*************************************************************************
! FUNCTION tensor3matmul(mat3,mat2)
! Purpose
! =======
! Defines the dot product of a 3-tensor and a 2-tensor
!*************************************************************************
!

FUNCTION tensor3matmul(mat3,mat2)

real(dp), dimension(:,:,:), intent(in) :: mat3
real(dp), dimension(:,:), intent(in) :: mat2
real(dp), dimension(nsd,nsd,nsd) :: tensor3matmul
integer :: i, j, k, m

tensor3matmul = 0.d0
do i = 1,nsd
  do j = 1,nsd
    do k = 1,nsd
      do m = 1,nsd
        tensor3matmul(i,j,k) = tensor3matmul(i,j,k) + mat3(i,j,m)*mat2(m,k)
      enddo
    enddo
  enddo
enddo

END FUNCTION tensor3matmul

!*************************************************************************
! FUNCTION mattensor3mul(mat2,mat3)
! Purpose
! =======
! Defines the dot product of a 2-tensor and a 3-tensor
!*************************************************************************
!

FUNCTION mattensor3mul(mat2,mat3)

real(dp), dimension(:,:), intent(in) :: mat2
real(dp), dimension(:,:,:), intent(in) :: mat3
real(dp), dimension(nsd,nsd,nsd) :: mattensor3mul
integer :: i, j, k, m

mattensor3mul = 0.d0
do i = 1,nsd
  do j = 1,nsd
    do k = 1,nsd
      do m = 1,nsd
        mattensor3mul(i,j,k) = mattensor3mul(i,j,k) + mat2(i,m)*mat3(m,j,k)
      enddo
    enddo
  enddo
enddo

END FUNCTION mattensor3mul

!*************************************************************************
! FUNCTION identity(m)
! Purpose
! =======
! Return an identity matrix
!*************************************************************************
!

FUNCTION identity(m)

integer, intent(in) :: m
real(dp), dimension(m,m) :: identity
integer :: i

identity = 0.d0
do i = 1,m
  identity(i,i) = 1.d0
enddo

END FUNCTION identity

!*************************************************************************
! FUNCTION vectensor3mul(vec,mat)
! Purpose
! =======
! Defines the dot product of a vector and a 3-tensor
!*************************************************************************
!

FUNCTION vectensor3mul(vec,mat)

real(dp), dimension(:), intent(in) :: vec
real(dp), dimension(:,:,:), intent(in) :: mat
real(dp), dimension(size(vec),size(vec)) :: vectensor3mul
integer :: i, j, k

vectensor3mul = 0.d0
do j = 1,size(vec)
  do k = 1,size(vec)
    do i = 1,size(vec)
      vectensor3mul(j,k) = vectensor3mul(j,k) + vec(i)*mat(i,j,k)
    enddo
  enddo
enddo

END FUNCTION vectensor3mul

END MODULE maxent
