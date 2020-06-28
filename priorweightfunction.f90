!*************************************************************************
!                   PRIOR WEIGHTFUNCTION MODULE                          !
!                   ===========================                          !
!                                                                        !
! Purpose                                                                !
! =======                                                                !
! The prior (weight function) for the maxent basis function computations !
! Author: N. Sukumar, UC Davis                                           !
! Date  : September 2007 (Version 1.4)                                        !
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
    beta = 3.d0/(rmax(index)+small)**2                        ! 
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
