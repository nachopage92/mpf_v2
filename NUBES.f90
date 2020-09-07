module variables

    implicit none

    integer, parameter :: max_bcond   = 10

    ! derivate-type definitions

    type par_type
      character(len=100) :: test,str_npts,OS
      character(len=500) :: filename
      integer :: npt(2),mmcloud(2)
      logical :: show,cld_check
    end type par_type

    type geometria_type
      integer :: npoin,nline,ntria,mm_con(2)
      real*8, allocatable :: coord(:,:)
      integer, allocatable :: line_m(:,:),tria_m(:,:),vrtx_m(:)
      integer, allocatable :: esup1(:),esup2(:) ! linked list -> elements
      integer, allocatable :: psup1(:),psup2(:) ! linked list -> points
    end type geometria_type

    type clouds_type    ! is it better to use allocatable space???
      ! almacena multiples funciones de forma para comparar
      !real*8 :: shapes(npcloud_max,max_shsize) 
      real*8,allocatable :: pts_dist(:),inv_pts_dist(:)
      real*8 :: hcar(2)       ! hcar = ( hx , hy )
      integer,allocatable :: pts_list(:)
      integer :: np_cloud,np_first
    end type clouds_type

    type bcond_type
      integer :: id,btype
      real*8 :: value
    end type bcond_type
  
    type boundaries_type
      integer :: nbou
      type(bcond_type) :: bcond(max_bcond)
      integer :: pos_line_wall(max_bcond+1) !elemento 'linea' de contorno 
      integer :: pos_pts_wall(max_bcond+1)  !elemento 'punto' en contorno
      integer,allocatable :: pts_wall(:),line_wall(:)  !l. encadenadas
      real*8,allocatable  :: pts_norm(:,:)
      integer,allocatable :: perm(:)
    end type boundaries_type
  
    type boundary_face_type
      real*8, allocatable :: normal(:,:)
      real*8, allocatable :: center(:,:)
      real*8, allocatable :: hsize(:)
      real*8, allocatable :: length(:)
      integer, allocatable :: itype(:)
    end type boundary_face_type

    ! variables definition
    type(par_type) :: par
    type(geometria_type)          :: geometry
    type(boundary_face_type)      :: bface_data
    type(boundaries_type),target  :: boundaries
    type(clouds_type),allocatable :: clouds(:)
    character(len=100) :: msh_dir,cld_dir,sol_dir
    
    contains

    subroutine LECTURA_PARAMETROS
        implicit none
        character(len=50) :: aux
        
        ! nombre del test 
        call get_command_argument(1,par%test)

        ! numero de puntos por lado
        call get_command_argument(2,par%str_npts)

        ! nombre del archivo a leer
        par%filename=trim(par%test)//'_'//trim(par%str_npts)
        call string_to_integer_retrieve_npts_per_side
        
        ! mostrar resultados en terminal
        call get_command_argument(3,aux)
        if ( trim(aux) .eq. 'True' ) then
            par%show = .true.
        else ! diga "False" u otra cosa
            par%show = .false.
        end if
        
        ! numero minimo de puntos por nube (mmclouds(1))
        call get_command_argument(4,aux)
        read(aux,*) par%mmcloud(1)

        ! numero maximo de puntos por nube (mmclouds(2))
        call get_command_argument(5,aux)
        read(aux,*) par%mmcloud(2) 
        
        ! chequear nube calculada
        call get_command_argument(6,aux)
        if ( trim(aux) .eq. 'True' ) then
            par%cld_check = .true.
        else ! diga "False" u otra cosa
            par%cld_check = .false.
        end if
        
        ! sistema operativo
        call get_command_argument(7,par%OS)
        call directory_OS_selector

        return
        
        contains
        
        subroutine string_to_integer_retrieve_npts_per_side
            implicit none
            integer sep,siz
            ! para npts de la forma 12_34-56_78
            sep = index(par%str_npts,'-')
            siz = len(trim(par%str_npts))
            aux = par%str_npts(1:sep-1)
            read(aux,*) par%npt(1)
            aux = par%str_npts(sep+1:siz)
            read(aux,*) par%npt(2)
            return
        end subroutine string_to_integer_retrieve_npts_per_side

        subroutine directory_OS_selector
            implicit none
            if ( trim(par%OS) .eq. 'Windows' ) then ! backslash
              msh_dir = '.\DATOS\'
              cld_dir = '.\DATOS\'
              sol_dir = '.\DATOS\'
            else if ( trim(par%OS) .eq. 'Linux' ) then ! slash
              msh_dir = './DATOS/'
              cld_dir = './DATOS/'
              sol_dir = './DATOS/'
            else
              print*, 'ERROR: par%OS = ',par%OS
              print*, 'Debe indicarse SO "Windows"/"Linux". Programa detenido'
              stop
            end if
            return
        end subroutine
   
    end subroutine LECTURA_PARAMETROS


    Subroutine IMPORTAR_GEOMETRIA
        ! .............................................
        ! 1. NODOS    D.PRESCRITOS    F.PRESCRITAS
        !    n_total  nnod_d          nnod_f
        ! 2. PUNTOS     COORDENADAS
        ! do i = 1, n_total
        !   idpt   ( x_coord  y_coord )
        ! end do
        ! 3. DESPLAZAMIENTO PRESCRITO
        ! do i = 1, nnod_d
        !   idpt  ( if_xdir  if_ydir )  ( u_displ v_displ )
        ! end do
        ! 4. FUERZAS PREESCRITAS
        !   ¿? ------------------------ PENDIENTE!!!!
        ! 5. NORMALES
        ! nnod_contorno = nnod_d+nnod_f
        ! do i = 1, nnod_contorno
        !   id_contorno(i)  ( nx  ny )
        ! end do
        ! END
        ! .............................................
        Implicit None
        integer nd,nf,i,j,input,n_bound
        character (50) :: control_line,aux
        character (100) :: file_name
        logical ffail
      
        !   Lectura nombre del archivo por terminal
        file_name = trim(msh_dir)//trim(par%filename)//'.dat'
        print *,'# Reading input data...' 
    
        !   Read input geometry and boundary conditions    
        call EXISTFILE (file_name, ffail) ; if ( ffail ) stop
        
        call get_unit(input)
    
        open(input,file=file_name,action='read')
        
        read(input,*) ! NODOS  D.PRESCRITOS F.PRESCRITAS 
        read(input,*) geometry%npoin, nd, nf
        
        allocate( geometry%coord(2,geometry%npoin)  ) ! coordenadas pts totales
        
        read(input,*)! PUNTOS     COORDENADAS
        do i = 1,geometry%npoin
          read(input,*) j, geometry%coord(1:2,i)   ! total points
        end do
    
        read(input,*)! DESPLAZAMIENTOS PRESCRITOS
        do i = 1,nd
          read(input,*) ! dato no necesario
        end do

        read(input,*)! FUERZAS PRESCRITAS
        do i = 1,nf
          read(input,*) ! dato no necesario
        end do
        
        read(input,*)! NORMALES
        read(input,*) n_bound
        
        allocate(boundaries%pts_norm(2,n_bound),boundaries%perm(n_bound))

        do i = 1,n_bound
            read(input,*) boundaries%perm(i) , boundaries%pts_norm(1:2,i)
        end do

        read(input,'(a11)') control_line
        close(input)
        
        if ( control_line.eq.'END' ) then
          print *,'  input file: ',trim(file_name),' read Ok!' ; print*
        else
          print *,'  input file: ',trim(file_name)
          print *,'  !!! end of line not found in geodata file...' ; stop
        endif

        ! /////////////////////////////////////////////////////

        !   Lectura nombre del archivo por terminal
        file_name = trim(msh_dir)//trim(par%filename)//'.elem'
        print *,'# Reading input data...' 
    
        !   Read input geometry and boundary conditions    
        call EXISTFILE (file_name, ffail) ; if ( ffail ) stop
        
        call get_unit(input)
    
        open(input,file=file_name,action='read')

        read(input,*) ! 'N. Triangulos'
        read(input,*) geometry%ntria 
        allocate( geometry%tria_m(3,geometry%ntria) )
        do i = 1, geometry%ntria 
            read(input,*) geometry%tria_m(1:3,i)
        end do

        read(input,*) ! 'N. Lineas'
        read(input,*) geometry%nline
        allocate( geometry%line_m(2,geometry%ntria) )
        do i = 1,geometry%nline
            read(input,*) geometry%line_m(1:2,i) 
        end do

        read(input,'(a11)') control_line
        close(input)
        
        if ( control_line.eq.'END' ) then
          print *,'  input file: ',trim(file_name),' read Ok!' ; print*
        else
          print *,'  input file: ',trim(file_name)
          print *,'  !!! end of line not found in geodata file...' ; stop
        endif
        
        return
    end subroutine IMPORTAR_GEOMETRIA


!    ! ------------------------------------------------------------------- !
!    !                                                                     !
!    ! Subroutine GENERAR_NUBES_V1                                         !
!    !                                                                     ! 
!    ! utiliza algoritmo 'fuerza bruta' utilizando criterio por CERCANIA   !
!    ! calcula la distancia entre el punto estrella y el resto del dominio !
!    !   utiliza dos criterios:                                            ! 
!    !       1. por cercania                                               ! 
!    !       2. por numero de puntos                                       ! 
!    !                                                                     !
!    ! el orden ascendente de distancia determina utilizando la rutina     !
!    ! m_mrgrnk (calcula el ranking de una lista de numeros reales)        !
!    !                                                                     ! 
!    ! ------------------------------------------------------------------- !
!    subroutine GENERAR_NUBES_V1   ! utiliza algoritmo 'fuerza bruta'
!        use m_mrgrnk
!        implicit none
!        real*8,allocatable  :: dist(:,:),xyz(:)
!        integer,allocatable :: list(:,:)
!        real*8 diff(2)
!        integer npoin,n,i,j,npmax,npmin,np
!        ! 1. numero total de puntos
!        n = geometry%npoin
!        ! 2. calcular_distancia_entre_cada_punto
!        allocate(dist(n,n)) ; dist=0
!        do i = 1,n            ! loop sobre cada punto de colocacion
!          do j = 1,n               ! loop sobre cada punto fuente
!            diff = geometry%coord(1:2,i) - geometry%coord(1:2,j)
!            dist(i,j) = sqrt(dot_product(diff,diff))
!          end do
!        end do
!        ! 3. ordenar_segun_cercania
!        allocate(list(n,n))
!        do i = 1,n
!          call mrgrnk(dist(i,:),list(i,:))
!        end do
!        ! 4. allocar memoria
!        allocate(clouds(n))
!        ! 5. consideracion especial si el numero de punts de la nube es 
!        !       mayor que el numero total de nodos.
!        if (par%mmcloud(2).gt.geometry%npoin) par%mmcloud(2)=geometry%npoin
!        do i=1,n
!            allocate(clouds(i)%pts_dist(par%mmcloud(2)))
!            allocate(clouds(i)%inv_pts_dist(par%mmcloud(2)))
!            allocate(clouds(i)%pts_list(par%mmcloud(2)))
!        end do
!        ! 6. crear_nube_segun_criterio
!        do i=1,n
!
!          ! numero de puntos por nube: min & max
!            ! ---------------> OJO ACA
!          np = par%mmcloud(2) !numero maximo de puntos por nube
!          clouds(i)%np_cloud = np
!
!          clouds(i)%pts_list(:) = 0
!          clouds(i)%pts_list(1:np) = list(i,1:np)
!          clouds(i)%pts_dist(:) = 0d0
!          clouds(i)%pts_dist(1:np) = dist(i,list(i,1:np))
!          clouds(i)%hcar = dist(i,3)
!        end do
!
!        return
!    end subroutine GENERAR_NUBES_V1 


    ! ------------------------------------------------------------------- !
    !                                                                     !
    ! Subroutine GENERAR_NUBES_V1                                         !
    !                                                                     ! 
    ! utiliza algoritmo 'fuerza bruta' utilizando criterio por CERCANIA   !
    ! calcula la distancia entre el punto estrella y el resto del dominio !
    !   utiliza dos criterios:                                            ! 
    !       1. por cercania                                               ! 
    !       2. por numero de puntos                                       ! 
    !                                                                     ! 
    !        Criterios de creacion de nube:                               ! 
    !    1. numero de puntos, por cercania simple                         ! 
    !    3. por cercania, forma circular                                  ! 
    !    4. por cercania, forma rectangular (consideracion especial       !
    !                                        en los contornos             ! 
    !                                                                     !
    ! el orden ascendente de distancia determina utilizando la rutina     !
    ! m_mrgrnk (calcula el ranking de una lista de numeros reales)        !
    !                                                                     !
    ! En caso de priorizar tiempo en calcular, utilizar la rutina         !
    ! GENERAR_NUBES_V2, que utilizar metodos optimizados para             !
    ! generar las nubes (metodo BIN/CUADRICULAS, Rainald, L.)             !
    !                                                                     ! 
    ! ------------------------------------------------------------------- !
    subroutine GENERAR_NUBES_V1   ! utiliza algoritmo 'fuerza bruta'
        use m_mrgrnk
        implicit none
        real*8 xdis(2),diff(2),normal(2),right_lim,left_lim,up_lim,&
                &down_lim,right_offset,left_offset,up_offset,down_offset,&
                &add_left,add_right,add_up,add_down
        real*8,allocatable  :: dist(:,:),xyz(:),normaux(:,:)
        integer,allocatable :: list(:,:),list_cld(:),dist_cld(:)
        integer npoin,n,i,j,k,npmax,npmin,np,criterio,cont
        character(len=100) :: selector
        logical :: isotropic

        ! preambulo
        n = geometry%npoin
        if ( allocated(clouds) ) deallocate(clouds) ! despejar memoria
        allocate( clouds(n) )
        if (par%mmcloud(2).gt.geometry%npoin) par%mmcloud(2)=geometry%npoin ! <-- Tener en consideracion
        ! reconocer direccion del vector normal
        allocate(normaux(n,2));normaux=0
        do j=1,size(boundaries%perm)
            normaux(boundaries%perm(j),1:2) = boundaries%pts_norm(1:2,j)
        end do


        criterio = 3
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ! ALGORITMO FUERZA BRUTA

        allocate(dist(n,n)) ; dist=0
        allocate(list(n,n)) ; list=0
        allocate(list_cld(n),dist_cld(n))
        
        ! 1. NP FIJO, SEGUN CERCANIA
        IF ( CRITERIO .EQ. 1 ) THEN

            do i=1,n

                ! 1.1. algoritmo de fuerza bruta: calcula la distancia de cada
                ! punto del dominio respecto a cada punto estrella. Utilizando
                ! la subrutina MRGRNK orden la lista de puntos segun distancia
                ! ascendente
                do j = 1,n               ! loop sobre cada punto fuente
                  diff(1) = geometry%coord(1,i) - geometry%coord(1,j)
                  diff(2) = geometry%coord(2,i) - geometry%coord(2,j)
                  dist(i,j) = sqrt(dot_product(diff(1:2),diff(1:2)))
                end do

                ! ordenar datos de menor a mayor distancia
                call mrgrnk(dist(i,1:n),list(i,1:n))

                ! distancia internodal caracteristica
                clouds(i)%hcar(1) = dist(i,2) ! primer punto mas cercano respecto al n.estrella
                clouds(i)%hcar(2) = dist(i,2) ! primer punto mas cercano respecto al n.estrella

                ! 1.1. numero de puntos por nube
                np = par%mmcloud(2) !   OBSERVACION : 
                                    ! par%mmcloud(2) = max npt 
                                    !   IDEA        :
                                    ! la nube se recorta segun conveniencia
                ! 1.2. asignar memoria
                allocate(clouds(i)%pts_dist(np))
                allocate(clouds(i)%inv_pts_dist(np))
                allocate(clouds(i)%pts_list(np))

                ! 1.3. definir variables globales
                clouds(i)%np_cloud = np
                clouds(i)%pts_list(:)=0   ; clouds(i)%pts_list(1:np) = list(i,1:np)
                clouds(i)%pts_dist(:)=0d0 ; clouds(i)%pts_dist(1:np) = dist(i,list(i,1:np))

            end do


        !   2. FORMA DESCRITA: CIRCULO
        !      (DENTRO DE UN RANGO, NPT VARIABLE)
        ELSE IF ( CRITERIO .EQ. 2 ) THEN
            print*, 'criterio 2 pendiente'

        !   3. FORMA DESCRITA: RECTANGULO
        ELSE IF ( CRITERIO .EQ. 3 ) THEN

            ! numero de puntos maximos (npx*npy)
            npmax = np

            right_offset = 0.05d0 ! 5%
            left_offset  = 0.05d0
            up_offset    = 0.05d0
            down_offset  = 0.05d0

            ! creacion de nubes rectangulares
            do i = 1,n ! punto de colocacion

                ! 3.1. Definir posicion del nodo: interior, contorno, vertice
                if ( .not. (any(i.eq.geometry%line_m)) ) then
                    selector = 'interior'
                else
                    selector = 'contorno'
                end if

                do j = 1,n
                  diff(1) = geometry%coord(1,i) - geometry%coord(1,j) !pt-estrella-pt.cercano
                  diff(2) = geometry%coord(2,i) - geometry%coord(2,j)
                  dist(i,j) = sqrt(dot_product(diff(1:2),diff(1:2)))
                end do

                ! ordenar datos de menor a mayor distancia
                call mrgrnk(dist(i,1:n),list(i,1:n))

                ! distancia internodal
                clouds(i)%hcar(1) = dist(i,list(i,2))
                clouds(i)%hcar(2) = dist(i,list(i,2))

                ! 3.2. offset, limites de la nube
                right_lim = geometry%coord(1,i)+clouds(i)%hcar(1)*(1.0d0+right_offset)
                left_lim  = geometry%coord(1,i)-clouds(i)%hcar(1)*(1.0d0+left_offset)
                up_lim    = geometry%coord(2,i)+clouds(i)%hcar(2)*(1.0d0+up_offset)
                down_lim  = geometry%coord(2,i)-clouds(i)%hcar(2)*(1.0d0+down_offset)

                ! ::::::::::::::: puntos interiores :::::::::::::::::::::
                if ( trim(selector) .eq. 'interior' ) then
                    np = 0 ; list_cld = 0
                    do j = 1,n ! punto funetes
                        diff(1) = geometry%coord(1,i) - geometry%coord(1,j) ! coord x 
                        diff(2) = geometry%coord(2,i) - geometry%coord(2,j) ! coord y
                        if (  geometry%coord(1,j) .ge. left_lim  .and.&
                            & geometry%coord(1,j) .le. right_lim .and.&
                            & geometry%coord(2,j) .ge. down_lim  .and.&
                            & geometry%coord(2,j) .le. up_lim )  then
                            np = np + 1
                            dist_cld(np) = j
                        end if
                    end do

                    ! ordena de menor a mayor
                    call mrgrnk(dist(i,dist_cld(1:np)),list_cld(1:np))

                    ! 1.2. asignar memoria
                    allocate(clouds(i)%pts_dist(np))
                    allocate(clouds(i)%inv_pts_dist(np))
                    allocate(clouds(i)%pts_list(np))

                    ! 1.3. definir variables globales
                    clouds(i)%np_cloud = np
                    clouds(i)%pts_list(:)=0   ; clouds(i)%pts_list(1:np) = dist_cld(list_cld(1:np))
                    clouds(i)%pts_dist(:)=0d0 ; clouds(i)%pts_dist(1:np) = dist(i,dist_cld(list_cld(1:np)))


                ! :::::::::::::::::::::: contorno ::::::::::::::::::::::::::::
                else if ( trim(selector) .eq. 'contorno' ) then

                    normal(1:2) = normaux(i,1:2)
                    
                    ! extender limites
                    !-------------------------------
                    !   cos(l) > 0  :  +1.0 right 
                    !   cos(l) < 0  :  -1.0 left
                    if ( normal(1).gt.0 ) then
                        add_right = 0.0d0 ; add_left = -1.0d0 
                    else if ( normal(1).lt.0 ) then
                        add_right = 1.0d0 ; add_left = 0.0d0 
                    else ! normal(1)=0
                        add_right = 0.0d0 ; add_left = 0.0d0 
                    endif
                    !-------------------------------
                    !   sen(l) > 0  :  +1.0 up
                    !   sen(l) < 0  :  -1.0 down
                    !-------------------------------
                    if ( normal(2).gt.0 ) then
                        add_up = 0.0d0 ; add_down = -1.0d0 
                    else if ( normal(2).lt.0 ) then
                        add_up = 1.0d0 ; add_down = 0.0d0 
                    else ! normal(2)=0
                        add_up = 0.0d0 ; add_down = 0.0d0 
                    endif
                    add_up = add_up*clouds(i)%hcar(2)
                    add_down = add_down*clouds(i)%hcar(2)
                    add_right = add_right*clouds(i)%hcar(1)
                    add_left = add_left*clouds(i)%hcar(1)

                    np = 0 ; list_cld = 0
                    do j = 1,n ! punto funetes
                        diff(1) = geometry%coord(1,i) - geometry%coord(1,j) ! coord x 
                        diff(2) = geometry%coord(2,i) - geometry%coord(2,j) ! coord y
                        if (  geometry%coord(1,j) .ge. left_lim +add_left  .and.&
                            & geometry%coord(1,j) .le. right_lim+add_right .and.&
                            & geometry%coord(2,j) .ge. down_lim +add_down  .and.&
                            & geometry%coord(2,j) .le. up_lim   +add_up)    then
                            np = np + 1
                            dist_cld(np) = j
                        end if
                    end do

                    ! ordena de menor a mayor
                    call mrgrnk(dist(i,dist_cld(1:np)),list_cld(1:np))

                    ! 1.2. asignar memoria
                    allocate(clouds(i)%pts_dist(np))
                    allocate(clouds(i)%inv_pts_dist(np))
                    allocate(clouds(i)%pts_list(np))

                    ! 1.3. definir variables globales
                    clouds(i)%np_cloud = np
                    clouds(i)%pts_list(:)=0   ; clouds(i)%pts_list(1:np) = dist_cld(list_cld(1:np))
                    clouds(i)%pts_dist(:)=0d0 ; clouds(i)%pts_dist(1:np) = dist(i,dist_cld(list_cld(1:np)))

                else
                    print*, 'selector no ingresado'
                    print*, 'programa detenido'
                    stop

                end if

            end do

        ELSE
            print*,' INTEGER(4) :: criterio    '
            print*,' Indicar criterio 1, 2 o 3 '
            print*,' ....Programa detenido.... '
            stop

        END IF

        ! despejar memoria
        deallocate(dist,list)

        return
    end subroutine GENERAR_NUBES_V1 

!    ! -------------------------------------------------------------------!
!    !                                                                    !
!    ! Subroutine GENERAR_NUBES_V2                                        !
!    !                                                                    !
!    ! utiliza algoritmo 'fuerza bruta' enmarcado dentro de un cuadradado !
!    ! de la forma ( xc-hx*ax:xc+hx*ax , yc-hy*ay:yc+hy*ay ) donde,       !
!    !   xc,yc coordenadas nodos estrellas                                !
!    !   hx,hy distancia internodal en direccion x e y                    !
!    !   ax,ay factor de amplificacion, define el tamaño y forma          !
!    !                                                                    !
!    ! -------------------------------------------------------------------!
!    subroutine GENERAR_NUBES_V2   
!        use m_mrgrnk              
!        implicit none             
!        real*8,allocatable  :: dist(:,:),xyz(:)
!        integer,allocatable :: list(:,:)
!        real*8 diff(2)
!        integer npoin,n,i,j,npmax,npmin,np
!
!        ! numero total de puntos
!        n = geometry%npoin
!
!        ! calcular_distancia_entre_cada_punto
!        ! dist(i,j) ===== > i : nodo estrella (colocacion)
!        !                   j : nodos dominio (fuente)
!        allocate(dist(n,n)) ; dist=0
!        do i = 1,n
!          do j = 1,n
!            diff = geometry%coord(1:2,i) - geometry%coord(1:2,j)
!            dist(i,j) = sqrt(dot_product(diff,diff))
!          end do
!        end do
!
!        ! calculo de la distancia internodal
!        if ( dist_type .eq. 'regular' ) then
!            h
!        else if ( dist_type .eq. 'irregular' ) then
!            print*, 'pendiente!!'
!            stop
!        end if
!            
!        
!
!        ! ordenar_segun_cercania
!        allocate(list(n,n))
!        do i = 1,n
!          call mrgrnk(dist(i,:),list(i,:))
!        end do
!        ! allocar memoria
!        allocate(clouds(n))
!        ! consideracion especial si el numero de punts de la nube es 
!        ! mayor que el numero total de nodos.
!        if (par%mmcloud(2).gt.geometry%npoin) par%mmcloud(2)=geometry%npoin
!        do i=1,n
!            allocate(clouds(i)%pts_dist(par%mmcloud(2)))
!            allocate(clouds(i)%inv_pts_dist(par%mmcloud(2)))
!            allocate(clouds(i)%pts_list(par%mmcloud(2)))
!        end do
!        ! crear_nube_segun_criterio
!        do i=1,n
!          ! numero de puntos por nube: min & max
!
!            ! OJO ACA
!          np = par%mmcloud(2) !numero maximo de puntos por nube
!          clouds(i)%np_cloud = np
!
!          clouds(i)%pts_list(:) = 0
!          clouds(i)%pts_list(1:np) = list(i,1:np)
!          clouds(i)%pts_dist(:) = 0d0
!          clouds(i)%pts_dist(1:np) = dist(i,list(i,1:np))
!          clouds(i)%hcar = dist(i,3)
!        end do
!
!        return
!    end subroutine GENERAR_NUBES_V2 

!    subroutine GENERAR_NUBES_V3
!      implicit none
!        
!      !PARAMETER (NODOS=30000,IFIJOS=10000,NELEMENT=50000)
!      !CHARACTER FILE1*80,FILE2*80
!      !INTEGER N(3,NELEMENT),LOGI(NODOS)
!      !INTEGER CONT_NODE(2,IFIJOS),IDES(3,IFIJOS),INORM_NODE(IFIJOS)
!      !INTEGER IPRES(IFIJOS)
!      !REAL*8 X(NODOS),Y(NODOS),RDES(2,IFIJOS),RPRES(2,IFIJOS)
!      !REAL*8 CL(IFIJOS),CM(IFIJOS) 
!      integer i,npoin,ilong,isbdry(geometry%npoin)
!      real*8 x(geometry%npoin),y(geometry%npoin)
!
!      npoin = geometry%npoin
!
!      isbdry = 0
!      DO i=1,geometry%nline
!         isbdry(geometry%line_m(1,i))=1
!         isbdry(geometry%line_m(2,i))=1
!      END DO
!
!      ilong = LONG_FILE(par%filename)
!      x = geometry%coord(1,:)
!      y = geometry%coord(2,:)
!
!      call BUSCA(geometry%npoin,geometry%ntria,geometry%tria_m,isbdry,&
!           &ilong,par%filename,x,y)
!
!      return
!      contains
!
!      ! -----------------------------------------------------------
!      !  SUBROUTINE BUSCA():
!      ! determina la nube de proximidad  (Autor: Franco Perazzo M.)
!      ! output      - ICON(NNOD) : nube de proximidad
!      !             - IC(NNOD)   : numero de puntos
!      SUBROUTINE BUSCA(NNOD,NELEM,N,LOGI,ILONG,FILE,X,Y)
!
!      INTEGER I,J,K,ICOUNT,IELEM,IN,INOD,JNOD,ICOUNT1,&
!                &IP1,IPAS,NELEMC,NELEM,NNOD,ILONG
!      INTEGER N(3,NELEM),ICON(NNOD,50),IC(NNOD)
!      INTEGER LOGI(NNOD),IC2(NNOD),NE(10)
!	  REAL*8 X(NNOD),Y(NNOD),D1,DIST
!
!      CHARACTER*80 FILE
!
!      DO INOD=1,NNOD
!        IC(INOD)=1
!        IC2(INOD)=0
!      END DO
!      
!      DO IELEM=1,NELEM
!        DO I=1,3
!          IC2(N(I,IELEM))=IC2(N(I,IELEM))+1
!          DO J=1,3
!            IF (J.NE.I) THEN
!              IC(N(I,IELEM))=IC(N(I,IELEM))+1
!              ICON(N(I,IELEM),IC(N(I,IELEM)))=N(J,IELEM)
!            END IF
!          END DO
!        END DO
!      END DO
!
!      
!      DO INOD=1,NNOD
!        ICON(INOD,1)=INOD
!        DO I=IC(INOD),2,-1
!          DO J=IC(INOD),1,-1
!            IF (ICON(INOD,I).EQ.ICON(INOD,J).AND.I.NE.J) THEN
!              ICON(INOD,J)=0
!              GOTO 10
!            END IF
!          END DO
!   10     CONTINUE
!        END DO
!      END DO
!
!      DO INOD=1,NNOD
!        ICOUNT=0
!        DO I=1,IC(INOD)
!          IF (ICON(INOD,I).NE.0) THEN
!            ICOUNT=ICOUNT+1
!            ICON(INOD,ICOUNT)=ICON(INOD,I)
!          END IF          
!        END DO
!        IC(INOD)=ICOUNT
!      END DO
!      
!!CCCC----------> BUSQUEDA EN LOS CONTORNOS      
!
!      DO INOD=1,NNOD
!        
!        IF (LOGI(INOD).EQ.1) THEN
!
!          NELEMC=0
!          DO IELEM=1,NELEM
!            DO J=1,3
!              IF (N(J,IELEM).EQ.INOD) THEN
!                NELEMC=NELEMC+1
!                NE(NELEMC)=IELEM
!              END IF
!            END DO
!          END DO
!                    
!          IF (IC2(INOD).EQ.1) THEN
!
!!CCCC-----------> NO HACE NADA
!                 
!          ELSE IF (IC2(INOD).EQ.2) THEN
!
!            ICOUNT=IC(INOD)
!            DO IN=2,ICOUNT
!              IF (LOGI(ICON(INOD,IN)).EQ.0) THEN
!                DO J=1,IC(ICON(INOD,IN))
!                  IC(INOD)=IC(INOD)+1
!                  ICON(INOD,IC(INOD))=ICON(ICON(INOD,IN),J)
!                END DO
!                DO J=1,IC(ICON(INOD,IN))
!                  IC(INOD)=IC(INOD)+1
!                  ICON(INOD,IC(INOD))=ICON(ICON(INOD,IN),J)
!                END DO
!              END IF
!            END DO
!            
!          ELSE IF (IC2(INOD).GE.3) THEN
!
!            IP1=0
!            DO IN=2,IC(INOD)
!              IF (LOGI(ICON(INOD,IN)).EQ.0) IP1=IP1+1
!            END DO
!            
!            IF (IP1.EQ.1) THEN
!
!              ICOUNT=IC(INOD)
!              DO J=1,ICOUNT
!                IC(INOD)=IC(INOD)+1
!                ICON(INOD,IC(INOD))=ICON(INOD,J)
!              END DO
!              
!              DO IN=2,ICOUNT
!                IF (LOGI(ICON(INOD,IN)).EQ.0) THEN
!                  DO J=1,IC(ICON(INOD,IN))
!                    IC(INOD)=IC(INOD)+1
!                    ICON(INOD,IC(INOD))=ICON(ICON(INOD,IN),J)
!                  END DO
!                  DO J=1,IC(ICON(INOD,IN))
!                    IC(INOD)=IC(INOD)+1
!                    ICON(INOD,IC(INOD))=ICON(ICON(INOD,IN),J)
!                  END DO
!                END IF
!              END DO
!                             
!            ELSE
!              
!              ICOUNT=IC(INOD)
!              DO IN=2,ICOUNT
!                IF (LOGI(ICON(INOD,IN)).EQ.0) THEN
!                  DO J=1,IC(ICON(INOD,IN))
!                    IC(INOD)=IC(INOD)+1
!                    ICON(INOD,IC(INOD))=ICON(ICON(INOD,IN),J)
!                  END DO
!
!                  IF (IC2(INOD).EQ.30) THEN
!                    DO J=1,IC(ICON(INOD,IN))
!                      IC(INOD)=IC(INOD)+1
!                      ICON(INOD,IC(INOD))=ICON(ICON(INOD,IN),J)
!                    END DO
!                  END IF
!                  
!                END IF              
!              END DO
!              
!            END IF
!
!          END IF          
!     
!          DO J=1,IC(INOD)
!            IPAS=0
!            DO K=J+1,IC(INOD)
!              IF(ICON(INOD,J).EQ.ICON(INOD,K).AND. ICON(INOD,J).NE.0) THEN                
!                ICON(INOD,K)=0
!                IPAS=1
!              END IF
!            END DO
!            IF (IPAS.EQ.0) ICON(INOD,J)=0
!          END DO
!          
!          ICOUNT1=0
!          DO J=1,IC(INOD)
!            IF (ICON(INOD,J).NE.0) THEN
!              ICOUNT1=ICOUNT1+1
!              ICON(INOD,ICOUNT1)=ICON(INOD,J)
!            END IF
!          END DO
!          IC(INOD)=ICOUNT1
!          
!        END IF
!
!      END DO
!
!
!      DO INOD=1,NNOD
!!C          IF (IC(INOD).LT.6) WRITE(*,'(A,I8)') ' NUBE INCOMPLETA:',INOD
!
!          D1=1.0/(NNOD)**0.5
!
!1000      IF (IC(INOD).LT.6) THEN
!              IC(INOD)=1
!              ICON(INOD,1)=INOD
!              DO JNOD=1,NNOD
!                  DIST=DSQRT((X(INOD)-X(JNOD))**2.+(Y(INOD)-Y(JNOD))**2.)
!                  IF (DIST.LT.D1.AND.INOD.NE.JNOD) THEN
!                      IC(INOD)=IC(INOD)+1
!                      ICON(INOD,IC(INOD))=JNOD
!                      IF (IC(INOD).EQ.30) THEN
!                          IC(INOD)=0
!                          D1=D1*(1-.025)
!                          GOTO 1000
!                      END IF
!                  END IF
!              END DO
!
!              IF (IC(INOD).LT.9) THEN
!                  IC(INOD)=0
!                  D1=D1*(1+.05)
!                  GO TO 1000
!              END IF
!
!          END IF
!      END DO
!
!      !! exporta         
!      !OPEN(8,FILE=FILE(1:ILONG)//'.CLD',STATUS='UNKNOWN')
!      !DO INOD=1,NNOD
!      !  WRITE(8,'(15I6)') IC(INOD),(ICON(INOD,I),I=1,IC(INOD))
!      !END DO
!      !CLOSE(8)
!
!      if (allocated(clouds)) deallocate(clouds)
!      allocate(clouds(nnod))
!
!      do inod = 1,nnod
!        clouds(inod)%np_cloud = ic(inod)
!        allocate(clouds(inod)%pts_list(ic(inod)))
!        clouds(inod)%pts_list(:) = icon(inod,1:ic(inod))
!      end do
!
!
!      RETURN
!      END
!
!    end subroutine GENERAR_NUBES_V3


!    subroutine balancear_nube
!        implicit none
!        real*8 xhat_np(2),xhat_min(2),dx(2),dx_sum1(2),dx_sum2(2),exc1,exc2
!        integer np_min,np_max,np,i,j,pt
!        !calcular centro de gravedad
!        np_min = par%mmcloud(1)
!        np_max = par%mmcloud(2)
!        do i=1,geometry%npoin
!
!            ! centro de gravedad de la cantidad minima de puntos
!            xhat_min = 0d0
!            do j = 1 , np_min
!                pt = clouds(i)%pts_list(j)
!                dx = geometry%coord(1:2,pt) - geometry%coord(1:2,i)
!                xhat_min(1:2) = xhat_min(1:2) + dx
!            end do
!            xhat_min = xhat_min/dfloat(np_min)
!            np = np_min
!
!            ! incrementar puntos en la nube y verificar balanceo
!            xhat_np = xhat_min
!            dx_sum1 = xhat_np*np
!            do j = np_min+1 , np_max
!                pt = clouds(i)%pts_list(j)
!                dx = geometry%coord(1:2,pt) - geometry%coord(1:2,i)
!                
!                ! calcular distancia del vector xhat
!                dx_sum2 = dx_sum1 + dx
!                xhat_np = dx_sum2/j
!                exc1 = dot_product(xhat_min,xhat_min)
!                exc2 = dot_product(xhat_np,xhat_np)
!                if (exc2.lt.exc1) then
!                    np = j  ! incrementar tamaño de la nube
!                    xhat_min = xhat_np
!                end if
!
!                ! incrementa en cada paso
!                dx_sum1 = dx_sum2
!
!            end do
!
!            ! nuevo tamaño de nube
!            clouds(i)%np_cloud = np
!
!        end do
!
!        return
!    end subroutine balancear_nube


    subroutine EXPORTAR_NUBES
        implicit none
        integer i,k,nnod,np,input
        character(len=500) :: new_filename
        logical existence
        
        call get_unit(input)
        
        nnod = geometry%npoin
        
        new_filename = trim(par%test)//'_'//trim(par%str_npts)
        new_filename=trim(cld_dir)//trim(new_filename)//'.CLD'
        
        ! si no existe la carpeta, la crea
        Inquire( FILE = trim(cld_dir)//'README.txt', EXIST = existence )
        if ( .not. existence ) then
          call system('mkdir '//trim(cld_dir) ) 
          !call system('touch '//trim(cld_dir)//'README.txt') 
          call system('echo > '//trim(cld_dir)//'README.txt') 
          ! inquire no reconoce las carpetas, por lo que se utiliza un
          ! archivo README para verificar si existe o no la carpeta.
          ! PENDIENTE: escribir breve leyenda en archivo README
        end if
        
        ! escribir archivo
        open(unit=input,file=new_filename,action='write')
        do i=1,nnod
          write(input,*) clouds(i)%np_cloud,clouds(i)%pts_list(1:clouds(i)%np_cloud)
        end do
        close(unit=input) ! fin escritura
    
    end subroutine EXPORTAR_NUBES


    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !::::::::::::::::::::: MISC SUBROUTINES ::::::::::::::::::::::::
    !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    Subroutine STORE_LINKED_LIST ( n_ent, n_flags, tot_flags, ent_list, id_list, entities, positions )
        Implicit None
        integer, intent(in) :: n_ent,n_flags,tot_flags(n_flags),ent_list(n_ent),id_list(n_ent)
        integer, intent(out) :: entities(n_ent),positions(n_flags+1)
        ! Local variables
        integer :: i,iflag
        
        positions(1) = 0
        
        do i = 2,n_flags+1  ! flag positions in the entities list
          positions(i) = positions(i-1) + tot_flags(i-1)
        end do
        
        do i = 1,n_ent
          iflag = id_list(i)
          positions(iflag) = positions(iflag) + 1
          entities(positions(iflag)) = ent_list(i)
        end do
        
        do i = n_flags,2,-1
         positions(i) = positions(i-1)
        end do
        
        positions(1) = 0
        
    End Subroutine STORE_LINKED_LIST


    Subroutine EXISTFILE ( fname, ffail )
    ! Checks whether a file exists before opening
        character (len = *) , intent(in) :: fname
        logical, intent(out) :: ffail
        inquire (file=fname, exist=ffail) ! check
        if ( ffail) then        
          ffail = .false.
        else
          print 10, fname
          ffail = .true.
        endif
 10     format('   !!! Cannot find file ',a) 
    End Subroutine EXISTFILE


    Subroutine get_unit ( iunit )
    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 September 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
        implicit none
        
        integer ( kind = 4 ) i
        integer ( kind = 4 ) ios
        integer ( kind = 4 ) iunit
        logical lopen
        
        iunit = 0
        
        do i = 1, 99
        
          if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
        
            inquire ( unit = i, opened = lopen, iostat = ios )
        
            if ( ios == 0 ) then
              if ( .not. lopen ) then
                iunit = i
                return
              end if
            end if
        
          end if
        
        end do
        
        return
    end subroutine


    !-----------------------------------------------------------------!
    ! CALCULA LA CANTIDAD DE LETRAS DEL NOMBRE DEL ARCHIVO DE ENTRADA !
    !-----------------------------------------------------------------!
    FUNCTION LONG_FILE(FILE)
    CHARACTER FILE*80
    INTEGER LONG_FILE,I
    DO I=1,80
      IF (FILE(I:I).EQ.' ') THEN
        LONG_FILE=I-1
        RETURN
      END IF
    END DO
    WRITE(10,'(A,A)') 'ERROR....( hay un error con el nombre del file)'
    STOP
    END

    subroutine print_header
        implicit none
        print*, '                                                 '
        print*, ':::::::::::::::::::::::::::::::::::::::::::::::::'
        print*, '       Programa Fortran 90                       '
        print*, '               Generador de nubes 2D             '
        print*, '                                                 '
    end subroutine print_header


    subroutine print_end_of_program
        implicit none 
        print*, '                                                 '
        print*, '       Fin Programa Fortran 90                   '
        print*, '               Generador de nubes 2D             '
        print*, ':::::::::::::::::::::::::::::::::::::::::::::::::'
        print*, '                                                 '
    end subroutine print_end_of_program

end module variables


! ::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::
! ::::::::::::::::::::::::::::::::::::::::::::::


!   GENERADOR DE NUBES 2D

! genera nube de aproximacion a partir de la
! siguiente informacion de entrada: 


PROGRAM GENERADOR_NUBES_2D 

    USE VARIABLES
    
    implicit none
    integer k
    character(len=20) :: generador_nube
    
    call print_header
    
 !.................. informacion ingresada ..........................!
 !  - par%test      (character(100)) : nombre del test implementado  !
 !  - par%npt       (integer(2))     : numero de puntos              !
 !  - par%mmclouds  (integer(2))     : numero min. & max. pt. x nube !
 !  - par%cld_check (logical)        : revisar calidad de la nube    !
 !  - par%show      (logical)        : mostrar info en terminal      !
 !  - par%OS        (character(100)) : Sistema operativo utilizado   !
    CALL LECTURA_PARAMETROS 
    
    CALL IMPORTAR_GEOMETRIA  
                             
    CALL GENERAR_NUBES_V1   ! utiliza algoritmo 'fuerza bruta'
    !CALL GENERAR_NUBES_V2   ! utiliza informacion de la triangulacion
    
    CALL EXPORTAR_NUBES
    
    call print_end_of_program


END PROGRAM GENERADOR_NUBES_2D
