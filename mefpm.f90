program Maximum_Entropy_Finite_Point_Method_Solver

    IMPLICIT NONE

    INTEGER,PARAMETER :: MAXCLD = 16  ! numero maximo permitido por nubes

    CHARACTER(LEN=100) :: FILENAME
    
    INTEGER :: IWEIGT          ! selector de 
    INTEGER :: NNOD,NFIXP,NFIXF,NNORM

    INTEGER,ALLOCATABLE :: IFIXP_NODE(:),IFIXF_NODE(:)
    INTEGER,ALLOCATABLE :: INORM_NODE(:),BAN(:,:)

    REAL*8 BET                 ! parametro de solver 
    REAL*8 YOUNG,POISSON       ! parametros del material (elasticidad lineal)

    REAL*8,ALLOCATABLE :: X(:),Y(:)
    REAL*8,ALLOCATABLE :: RFIXP_VALUE(:),RFIXV_VALUE(:)
    REAL*8,ALLOCATABLE :: RFIXF_VALUEX(:),RFIXF_VALUEY(:)
    REAL*8,ALLOCATABLE :: RNORM_VALUEX(:),RNORM_VALUEY(:)

    INTEGER,ALLOCATABLE :: CL_CONECT(:,:),NPUNTOS(:)

    REAL*8,ALLOCATABLE :: RM(:,:,:),DMAX(:)


!-------------> LEE EL NOMBRE DEL PROBLEMA
      
    OPEN(1,FILE='NUBESPUNT.DAT',STATUS='OLD')
    READ(1,'(A)') FILENAME
    CLOSE(1)

!-------------> LECTURA PARAMETROS GENERALES, GEOMETRIA, Y NUBES

    CALL READ_GENERAL_DATA
    
    CALL READ_GEOM_DATA

    CALL READ_CLOUDS

!-------------> CALCULO DE LA FUNCION DE FORMA

    CALL N_DN_D2N

    CONTAINS

!----------------------------------------------------!
!     SUBRUTINA DE LECTURA DE DATOS DEL PROBLEMA     !
!----------------------------------------------------!     
    SUBROUTINE READ_GENERAL_DATA
    
    !---------> LEE DATOS GENERALES DEL PROBLEMA
          
    OPEN(1,FILE=TRIM(FILENAME)//'.GEN',STATUS='OLD')
    
    READ(1,*)
    !---------> TIPO DE FUNCION DE PONDERACION            
    READ(1,*) IWEIGT
    READ(1,*)
    !---------> PARAMETRO DEL SOLVER
    READ(1,*) BET  
    READ(1,*)
    !---------> PROPIEDADES DEL MATERIAL
    READ(1,*)
    READ(1,*) YOUNG
    READ(1,*) 
    READ(1,*) POISSON    
    
    CLOSE(1)
    
    RETURN
    END

!----------------------------------------------------
!C    SUBRUTINA DE LECTURA DE DATOS DEL PROBLEMA    C
!----------------------------------------------------      
    SUBROUTINE READ_GEOM_DATA
    
    INTEGER I,NOD

    !-------------> LEE DATOS DE LA GEOMETRIA
    
    OPEN(1,FILE=trim(FILENAME)//'.dat',STATUS='OLD')
    
    !--------> DATOS DE CONTROL
    
    READ(1,*)
    READ(1,*) NNOD,NFIXP,NFIXF

    !--------> ASIGNAR MEMORIA
    ALLOCATE( X(NNOD),Y(NNOD) )
    ALLOCATE( IFIXP_NODE(NFIXP),RFIXP_VALUE(NFIXP),RFIXV_VALUE(NFIXP),BAN(6,NFIXP) )
    ALLOCATE( IFIXF_NODE(NFIXF),RFIXF_VALUEX(NFIXF),RFIXF_VALUEY(NFIXF) )
    
    !---------> COORDENADAS DE LOS PUNTOS
    READ(1,*)
    I=0
    DO WHILE (I.LT.NNOD) 
      I=I+1 
      READ(1,*) NOD,X(NOD),Y(NOD)
    END DO
    
    !---------> CONDICIONES DE CONTORNO
    !---------> PUNTO CON VALOR DE DESPLAZAMIENTO PRESCRITO
    I=0
    READ(1,*)
    DO WHILE (I.LT.NFIXP) 
      I=I+1 
      READ(1,*) IFIXP_NODE(I),BAN(1,I),BAN(2,I),RFIXP_VALUE(I),RFIXV_VALUE(I)
    END DO
          
    !-------> PUNTO CON VALOR DE FUERZA PRESCRITA
    I=0
    READ(1,*)
    DO WHILE (I.LT.NFIXF) 
      I=I+1 
      READ(1,*) IFIXF_NODE(I),RFIXF_VALUEX(I),RFIXF_VALUEY(I)
    END DO
    
    !-------> NODO CON VALOR NORMAL PRESCRITO
    I=0
    READ(1,*)
    READ(1,*)NNORM
    !--------> ASIGNAR MEMORIA VARIABLES DE NORMALES
    ALLOCATE( INORM_NODE(NNORM),RNORM_VALUEX(NNORM),RNORM_VALUEY(NNORM) )
    DO WHILE (I.LT.NNORM) 
      I=I+1 
      READ(1,*) INORM_NODE(I),RNORM_VALUEX(I),RNORM_VALUEY(I)
    END DO
    
    CLOSE(1)
    
    WRITE(*,'(A10,I5)') 'NNOD  :', NNOD 
    WRITE(*,'(A10,I5)') 'NFIXP :', NFIXP
    WRITE(*,'(A10,I5)') 'NFIXF :', NFIXF 
    WRITE(*,'(A10,I5)') 'NNORM :', NNORM
    
    CLOSE(1)
    
    RETURN
    END

      
!------------------------------------------------------------------!
!    LEE LAS NUBES DE PUNTOS DE UN ARCHIVO CON EXTENSION ".CLD"    !
!     CON EL SIGUIENTE FORMATO:  NPTOS , PT_1 , PT_2 , PT_3 ,      !
!     ... , PT_NPTOS                                               !
!------------------------------------------------------------------!
    SUBROUTINE READ_CLOUDS
    
    INTEGER INOD,J
    
    ALLOCATE( CL_CONECT(MAXCLD,NNOD),NPUNTOS(NNOD) )

    OPEN(1,FILE=TRIM(FILENAME)//'.CLD',STATUS='OLD')
    
    DO INOD=1,NNOD
      READ(1,*) NPUNTOS(INOD),(CL_CONECT(J,INOD),J=1,NPUNTOS(INOD))
    END DO

    CLOSE(1)
    
    RETURN
    END


!-----------------------------------------------------------
!C    CALCULO DE LAS FUNCIONES DE FORMA Y SUS DERIVADAS    C
!-----------------------------------------------------------
    SUBROUTINE N_DN_D2N
    
    INTEGER,PARAMETER :: NCONS = 6   ! [ 1 ; x ; y ; x^2 ; xy ; y^2 ]
    INTEGER INOD,I,J,IERR,J1,J2,J3
    
    REAL*8 DISMAX,DX,DY,DIS,BETA,ALFA,PED,PEX,SUMA,XSTAR,YSTAR
    REAL*8 X(NNOD),Y(NNOD)
    REAL*8 C(MAXCLD,NCONS),D(NCONS,NCONS),RM_P(NCONS,MAXCLD),W(MAXCLD)

    ALLOCATE( RM(NCONS,MAXCLD,NNOD) , DMAX(NNOD) )
    
    DO INOD=1,NNOD 
    
        XSTAR=X(CL_CONECT(1,INOD))
        YSTAR=Y(CL_CONECT(1,INOD))
        
        !--------------> CALCULO DE LA DISTANCIA AL NODO MAS LEJANO DE LA NUBE   
        DMAX(INOD)=0.D0
        DO J=1,NPUNTOS(INOD)  
            
            DX=X(CL_CONECT(J,INOD))-XSTAR 
            DY=Y(CL_CONECT(J,INOD))-YSTAR 
            DIS=SQRT(DX**2+DY**2) 
            IF (DMAX(INOD).LT.DIS) DMAX(INOD)=DIS 
            
        END DO   
        
        !--------------> CALCULO DE LOS PARAMETROS DE LA FUNCION DE PESO
        
        BETA=dmax(inod)*1.1
        ALFA=BETA/2.0
        PED=DEXP(-(BETA/ALFA)**2) 
        
        DO J=1,NPUNTOS(INOD)  
          
            DX=X(CL_CONECT(J,INOD))-XSTAR 
            DY=Y(CL_CONECT(J,INOD))-YSTAR 
            
            
            IF (IWEIGT.EQ.0) THEN
              
                DIS=SQRT(DX**2+DY**2)           
                PEX=DEXP(-(DIS/ALFA)**2)             
                W(J)=(PEX-PED)/(1.D0-PED) 
                
            ELSE IF (IWEIGT.EQ.1) THEN
              
                W(J)=1.D0
                !W(J)=.020D0
                !IF (J.EQ.1) W(J)=1.D0
            
            ELSE IF (IWEIGT.EQ.2) THEN
            
                DIS=SQRT(DX**2+DY**2)
                W(J)=1.D0-DIS/BETA
              
            ELSE IF (IWEIGT.EQ.3) THEN  
            
                DIS=SQRT(DX**2+DY**2)
                W(J)=1.D0-(DIS/BETA)**2
            
            ELSE IF (IWEIGT.EQ.4) THEN  
            
                DIS=SQRT(DX**2+DY**2)
                W(J)=.5D0*(1+DCOS(3.141592654*DIS/BETA)) 
            
            ELSE IF (IWEIGT.EQ.5) THEN  
            
                DIS=SQRT(DX**2+DY**2)
                W(J)=1.D0-6.D0*(DIS/BETA)**2+8.D0*(DIS/BETA)**3-3.D0*(DIS/BETA)**4
              
            END IF
            
            dx=dx/dmax(inod)
            dy=dy/dmax(inod)
            C(J,1)=1.D0 
            C(J,2)=DX
            C(J,3)=DY
            C(J,4)=(DX)**2.D0 
            C(J,5)=DX*DY 
            C(J,6)=(DY)**2.D0
            !c          C(J,7)=DX**3.D0
            !c          C(J,8)=(DX**2.D0)*DY
            !c          C(J,9)=(DY**2.D0)*DX
            !c          C(J,10)=DY**3.D0
        
        END DO 
        
        !--------------> CALCULO DE LA MATRIZ  "CT * C" 
        
        DO J1=1,NCONS
            DO J2=1,NCONS 
                SUMA=0d0
                DO J3=1,NPUNTOS(INOD) 
                    SUMA=SUMA+C(J3,J2)*C(J3,J1)*W(J3) 
                END DO 
                D(J1,J2)=SUMA 
            END DO 
        END DO
        
        
        !--------------> INVIERTE LA MATRIZ  "D"  
        
        CALL INVERT(NCONS,IERR,D)
        IF (IERR.EQ.1) THEN          
            STOP 'PROBLEMAS EN LAS NUBES...'
        !c         write(*,*) 'PROBLEMAS EN LAS NUBES...',inod
        END IF 
        
        !--------------> CALCULO DE  "D(-1) * CT" 
        
        DO J1=1,NCONS
            DO J2=1,NPUNTOS(INOD)
                SUMA=0d0
                DO J3=1,NCONS 
                    SUMA=SUMA+D(J1,J3)*C(J2,J3) 
                END DO 
                RM_P(J1,J2)=SUMA*W(J2)
            END DO
        END DO
        
        !--------------> LIMPIA LAS FUNCIONES DE FORMA
        
        DO I=1,NPUNTOS(INOD)
            DO J=1,10
                RM(J,I,INOD)=0.D0
            END DO                    
        END DO
        
        !--------------> GUARDA LOS VALORES NECESARIOS
        
        DO I=1,NPUNTOS(INOD)
        !c          write(*,*) rm(1,i,inod)
        RM(1,I,INOD)=RM_P(1,I)
        RM(2,I,INOD)=RM_P(2,I)/dmax(inod)
        RM(3,I,INOD)=RM_P(3,I)/dmax(inod)
        RM(4,I,INOD)=RM_P(4,I)*2/dmax(inod)**2
        RM(5,I,INOD)=RM_P(5,I)/dmax(inod)**2
        RM(6,I,INOD)=RM_P(6,I)*2/dmax(inod)**2
        END DO
    
    !c        write(*,*) rm(1,i,inod)
    END DO
    
    RETURN
    END



!-----------------------------------------------------!
!    INVIERTE UNA MATRIZ DE (NCONSxNCONS) Y LA        ! 
!     DEVUELVE EN LA MISMA MATRIZ DE ENTRADA          ! 
!-----------------------------------------------------!
    SUBROUTINE INVERT(NCONS,IERR,D) 
    
    INTEGER, INTENT(IN)   :: NCONS
    INTEGER, INTENT(OUT)  :: IERR
    REAL*8, INTENT(INOUT) :: D(NCONS,NCONS)
    
    INTEGER I,J,K
    REAL*8 RM,RMAX,RLIM
    REAL*8 D1(NCONS,2*NCONS) 
    
    RMAX=0.D0
    DO I=1,NCONS
        D1(I,I+NCONS)=1.D0 
        IF (RMAX.LT.DABS(D(I,I))) RMAX=DABS(D(I,I)) 
        DO J=1,NCONS 
            D1(I,J)=D(I,J) 
            IF (I.NE.J) D1(I,J+NCONS)=0.D0            
        END DO 
    END DO 
    
    RLIM=1.D-11*RMAX
    
    DO I=1,NCONS 
        IF (DABS(D1(I,I)).LT.RLIM) THEN 
            IERR=1 
            RETURN    
        END IF 
        
        DO J=I+1,NCONS    
        
            RM=D1(I,J)/D1(I,I) 
            DO K=J,2*NCONS
                D1(J,K)=D1(J,K)-RM*D1(I,K)                
            END DO                    
             
        END DO 
        
    END DO 
    
    DO I=NCONS,1,-1            
        DO K=1,NCONS 
            DO J=I+1,NCONS 
                D1(I,K+NCONS)=D1(I,K+NCONS)-D1(I,J)*D1(J,K+NCONS) 
            END DO 
            D1(I,K+NCONS)=D1(I,K+NCONS)/D1(I,I) 
        END DO 
    END DO 
           
    IERR=0 
    
    DO I=1,NCONS
        DO J=1,NCONS
            D(I,J)=D1(I,J+NCONS)
            write(*,*) d(i,j)
        END DO 
    END DO  
    
    RETURN 
    END 
      

end program Maximum_Entropy_Finite_Point_Method_Solver



