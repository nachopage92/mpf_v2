module lme_simplificado

contains

subroutine maxent2d_shapefcn_v2(np,x_a,x,beta,maxiter,init_lam,tol,p_a,dp_a,ddp_a)
    implicit none

    integer,intent(in) :: np,maxiter
    real*8,intent(in)  :: x_a(2,np),x(2),init_lam(2),tol,beta(np)
    real*8,intent(out) :: p_a(np),dp_a(2,np),ddp_a(2,2,np)

    integer i,cont,ierr
    real*8 lam(2),dlam(2),r_norm,x_dist(2,np),J(2,2),r(2)

    !algunas definiciones
    x_dist(1,:) = x(1)-x_a(1,1:np)
    x_dist(2,:) = x(2)-x_a(2,1:np)
    
    !iniciacion de las variables para cada iteracion
    ! r   := nabla_{lambda}(logZ) = sum_a( p_a (x_a-x) )
    ! lam := mult. de Lagrange (consistencia lineal) 
    lam = init_lam
    r_norm   = huge(1d0) ! norma vectorial de r
    cont     = 0         ! contador iter
    
    !ALGORITMO DE BUSQUEDA DE RAICES (NwR)
    do while (r_norm .gt. tol)

        ! funcion auxiliar
        call gamma2D_v2

        ! invertir matriz J
        call invert(2,ierr,J) 
        if ( ierr .eq. 1 ) then 
            print*, 'error: maxent2d_shapefcn_v2, matriz j mal condicionada'
            stop
        end if

        dlam = matmul(J,r)
        lam=lam-dlam
        cont=cont+1

        ! verifica numero de iteraciones
        if( cont .gt. maxiter )then
            print*, 'ERROR: MAXENT2D_SHAPEFCN_V2, supero n. de iteraciones'
            STOP
        end if

        ! calculo norma vector r
        r_norm = sqrt(dot_product(r,r))

    end do
    
    !Gradiente espacial, primera derivada  
    dp_a = -matmul(J,x_dist) ! J ya esta invertido
    dp_a(1,1:np) = dp_a(1,1:np)*p_a
    dp_a(2,1:np) = dp_a(2,1:np)*p_a
    
    !Segunda derivada
    call hessian2D_v2
    
    return
    contains

    subroutine gamma2D_v2
        implicit none
        integer id
        real*8 temp(np),Z,gam
        do id=1,np !numerador funcion de forma
            temp(id) = exp(-beta(id) * dot_product(x_dist(:,id),x_dist(:,id)) + dot_product(lam(:),x_dist(:,id)))
        end do
        Z = sum( temp )        !funcion de particion
        p_a(:) = temp(:) / Z   !funcion de forma
        gam = log(Z)           !funcion dual de lagrange
        do id=1,2 !calculo de r (gradiente de log(z))
            r(id) = sum( p_a * x_dist(id,:)  )
        end do
        J = 0d0 !calculo de J (gradiente de r=log(z))
        do id=1,np
            J = J + p_a(id) * dyadic(x_dist(:,id),x_dist(:,id))
        end do
        J = J - dyadic(r,r)
    end subroutine

    subroutine hessian2D_v2
        implicit none
        integer id, jd
        real*8 DJ(2,2,2),aux(2,2,np),Ident(2,2)
        !calculo del tensor derivada J  
        DJ = 0d0
        do id = 1,2
            do jd = 1, np
                DJ(id,:,:) = DJ(id,:,:) + dp_a(id,jd) * dyadic(x_dist(:,jd),x_dist(:,jd))
            end do
        end do
        !calculo de producto tensorial (contraccion) entre el gradiente de J y el gradiente de la funcion de forma
        aux = 0d0
        do id = 1,2
            do jd = 1,np
                aux(:,:,jd) = aux(:,:,jd) + DJ(:,:,id) * dp_a(id,jd) 
            end do
        end do  
        !calculo de la segunda derivada de la funcion de forma
        Ident(1,1:2) = (/ 1d0 , 0d0 /)
        Ident(2,1:2) = (/ 0d0 , 1d0 /)
        do jd = 1,np
            ddp_a(:,:,jd) = -aux(:,:,jd) - dyadic( dp_a(:,jd), x_dist(:,jd) ) - p_a(jd)*Ident
            ddp_a(:,:,jd) = matmul( J, ddp_a(:,:,jd) ) ! J ya está invertido
        end do
        return
    end subroutine

end subroutine


!-----------------------------------------------------!
!    invierte una matriz de (nconsxncons) y la        ! 
!     devuelve en la misma matriz de entrada          ! 
!-----------------------------------------------------!
    subroutine invert(ncons,ierr,d) 
    
        integer, intent(in)   :: ncons
        integer, intent(out)  :: ierr
        real*8, intent(inout) :: d(ncons,ncons)
        
        integer i,j,k
        real*8 rm_,rmax,rlim
        real*8 d1(ncons,2*ncons) 
        
        rmax=0.d0
        do i=1,ncons
            d1(i,i+ncons)=1.d0 
            if (rmax.lt.dabs(d(i,i))) rmax=dabs(d(i,i)) 
            do j=1,ncons 
                d1(i,j)=d(i,j) 
                if (i.ne.j) d1(i,j+ncons)=0.d0            
            end do 
        end do 
        
        rlim=1.d-11*rmax
        
        do i=1,ncons 
            if (dabs(d1(i,i)).lt.rlim) then 
                ierr=1 
                return    
            end if 
            
            do j=i+1,ncons    
            
                rm_=d1(i,j)/d1(i,i) 
                do k=j,2*ncons
                    d1(j,k)=d1(j,k)-rm_*d1(i,k)                
                end do                    
                 
            end do 
            
        end do 
        
        do i=ncons,1,-1            
            do k=1,ncons 
                do j=i+1,ncons 
                    d1(i,k+ncons)=d1(i,k+ncons)-d1(i,j)*d1(j,k+ncons) 
                end do 
                d1(i,k+ncons)=d1(i,k+ncons)/d1(i,i) 
            end do 
        end do 
               
        ierr=0 
        
        do i=1,ncons
            do j=1,ncons
                d(i,j)=d1(i,j+ncons)
            end do 
        end do  
        
        return 
    end subroutine invert

!---------------------------------------------------------!
!    producto exterior (outer product) entre dos vectores !
!     Ignacio Apablaza. iapablaza@protonmail.com          ! 
!---------------------------------------------------------!
    pure function dyadic(a,b) result(c)
        implicit none
        real*8,intent(in)  :: a(:),b(:)
        real*8,allocatable :: c(:,:)
        integer::i,j,k,tmp(2)
        tmp=(/size(a),size(b)/)
        allocate(c(tmp(1),tmp(2)))
        do i=1,tmp(1)
            do j=1,tmp(2)
                c(i,j) = a(i)*b(j)
            end do
        end do
    return
    end function

end module lme_simplificado
