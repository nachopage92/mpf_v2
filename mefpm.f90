program Finite_Point_Method_Solver

    IMPLICIT NONE

    INTEGER,PARAMETER :: MAXCLD = 16  ! numero maximo permitido por nubes
    INTEGER,PARAMETER :: NCONS  = 6   ! [ 1 ; x ; y ; x^2 ; xy ; y^2 ]
                                      ! phi,dphix,dphiy,ddphixx,ddphixy,ddphiyy

    CHARACTER(LEN=100) :: FILENAME
    
    INTEGER :: IWEIGT          ! selector de funcion de peso
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

    INTEGER,ALLOCATABLE :: CL_CONECTN(:,:),IFIXP_NODEN(:),NPUNTOSN(:)
    REAL*8,ALLOCATABLE :: RFIXP_VALUEN(:)

    REAL*8,ALLOCATABLE :: S(:,:),S1(:,:),DMAX(:)
    REAL*8,ALLOCATABLE :: RM_FWLS(:,:,:),RM_MaxEnt(:,:,:),RM(:,:,:)
    REAL*8,ALLOCATABLE :: F(:),PRESS(:),DU(:),DV(:),U(:),V(:)



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

    CALL SOLVE_POT

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
    END SUBROUTINE READ_GENERAL_DATA

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
    END SUBROUTINE READ_GEOM_DATA

      
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
    END SUBROUTINE READ_CLOUDS

      

!:::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::

    SUBROUTINE SOLVE_POT
    
    REAL*8 RES(NNOD),RESB(NNOD),PK(NNOD),PKB(NNOD),APK(NNOD),APKB(NNOD)
    REAL*8 ESFX(NNOD),ESFY(NNOD),TAU(NNOD)
    REAL*8 DESPU,DESPV,DS1,DS2,DS3,DUDX,DVDX,DUDY,DVDY
    REAL*8 ESFXEXA(NNOD),ESFYEXA(NNOD),TAUEXA(NNOD)
    INTEGER I,INOD

!--------------> ASIGNAR MEMORIA
      
    ALLOCATE( S(MAXCLD*2,NNOD*2),S1(MAXCLD*2,NNOD*2) )
    ALLOCATE( F(NNOD*2),PRESS(NNOD*2),DU(NNOD),DV(NNOD),U(NNOD),V(NNOD) )

!--------------> ENSAMBLA LA MATRIZ DE RIGIDEZ
      
    CALL ENSAMB_LAPLACE
      
!--------------> IMPONE EL TERMINO "F" DEL SISTEMA DE ECUACIONES
    
    DO INOD=1,NNOD*2
        F(INOD)=0.D0      
    END DO

!--------------> IMPONE EL VALOR DE LAS FUERZAS PRESCRITAS

    DO I=1,NFIXF
        INOD=IFIXF_NODE(I)

!--------------------------------------------------------------------CCC
!C               ENSAMBLANDO H CARACTERISTICO                          C
!--------------------------------------------------------------------CCC

!cccc-----> ec.contorno-h/2*(ec.dominio)=0
         
        F(INOD*2-1)=RFIXF_VALUEX(I)*((1.d0-(POISSON**2.d0))/YOUNG)
        F(INOD*2)=RFIXF_VALUEY(I)*((1.d0-(POISSON**2.d0))/YOUNG)
         
!c         write(*,*) inod,rfixf_valuex(i),rfixf_valuey(i)       !f(inod*2-1),f(inod*2)

!--------------------------------------------------------------------CCC
!C             SIN ENSAMBLAR H CARACTERISTICO                          C
!--------------------------------------------------------------------CCC

!c         F(INOD*2-1)=RFIXF_VALUEX(I)*((1-(POISSON**2))/YOUNG)
!c         F(INOD*2)=RFIXF_VALUEY(I)*((1-(POISSON**2))/YOUNG)
         
    END DO

!!--------------> IMPONE EL VALOR DE LOS DESPLAZAMIENTOS PRESCRITOS
!
    DO I=1,NFIXP
        INOD=IFIXP_NODE(I)
        IF (BAN(1,I).EQ.1) THEN
            F(INOD*2-1)=RFIXP_VALUE(I)
        END IF
        IF (BAN(2,I).EQ.1) THEN
            F(INOD*2)=RFIXV_VALUE(I)
        END IF
        ! write(*,*) f(inod*2-1),f(inod*2)
    END DO

!!--------------> PRECONDICIONA LA MATRIZ DE RIGIDEZ Y
!!-------------->      EL TERMINO INDEPENDIENTE
!
    CALL PRECF
      
!--------------> REALIZA CAMBIO DE VARIABLES
      
    CALL CAMBIOVAR
      
!--------------> RESUELVE EL SISTEMA DE ECUACIONES SIN ENSAMBLAR
!--------------> LOS NODOS PRESCRITOS
      
!       CALL GRADCONJ(S,PRESS,F,NNOD*2,CL_CONECTN,NPUNTOSN, &
!     &  IFIXP_NODEN,NFIXP*2,RFIXP_VALUEN)

!--------------> RESUELVE EL SISTEMA DE ECUACIONES ENSAMBLANDO
!--------------> LOS NODOS PRESCRITOS

!    CALL GRADCONJ(S,PRESS,F,NNOD*2,MAXCLD*2,CL_CONECTN &
!    &  ,NPUNTOSN,IFIXP_NODEN,0,RFIXP_VALUEN)

    CALL BCG2()     !   BICGSTAB_v2

    !CALL UMF()      !   Unsymmetric-pattern MultiFrontal method 
      
    OPEN(1,FILE=TRIM(FILENAME)//'2d.res',STATUS='UNKNOWN')
    
    DO INOD=1,NNOD
        DU(INOD)=PRESS(INOD*2-1)
        DV(INOD)=PRESS(INOD*2)
        ! write(*,*) inod,DU(INOD),DV(INOD)
    END DO
    
!--------------> IMPRESION DE LOS DESPLAZAMIENTOS
    WRITE(1,'(A15,I4,I6,3I4)') 'DESPLAZAMIENTOS',2,0,2,1,0
    
    DO INOD=1,NNOD
    
       CALL EVAVECT(INOD,DESPU,DESPV,1)
       WRITE(1,'(I6,2E25.13)') INOD,DESPU,DESPV
    END DO

!--------------> IMPRESION DE LOS ESFUERZOS
      
    WRITE(1,'(A15,I4,I6,3I4)') 'ESFUERZOS      ',2,0,3,1,0
    
    DS1 = (1.D0-POISSON)/2.d0
    DS3 = YOUNG/(1.d0-(POISSON**2.d0))
    
    DO INOD=1,NNOD
        
        CALL EVAVECT(INOD,DUDX,DVDX,2)
        CALL EVAVECT(INOD,DUDY,DVDY,3)
        
        ESFX(INOD)=DS3*(DUDX+POISSON*DVDY)
        ESFY(INOD)=DS3*(POISSON*DUDX+DVDY)
        TAU(INOD)=DS3*(DS1*(DVDX+DUDY))
        
        WRITE(1,'(I6,3E25.13)') INOD,ESFX(INOD),ESFY(INOD),TAU(INOD)
        
    END DO
    
    CLOSE(1)

!!----C---> TENSIONES (exactas,numericas) en un plano a una distancia 
!!----C     de 12 unidades del borde superior(libre) de la placa
!
!    OPEN(1,FILE=TRIM(FILENAME)//'2d.res-1',STATUS='UNKNOWN')
!    
!    
!    DO INOD=1,NNOD
!        if (x(inod).eq.12) then
!        
!        ESFXEXA(INOD)=-20/3.14159265/12*(cos(atan(y(inod)/12)))**4
!        TAUEXA(INOD)=-20/3.14159265/12*sin(atan(y(inod)/12))*(cos(atan(y(inod)/12)))**3
!        
!        WRITE(1,'(5E16.6)')y(inod),esfxexa(inod),esfx(inod),tauexa(inod),tau(inod)
!        end if
!         
!    END DO
!    
!    CLOSE(1)
    
    RETURN
    
    END



    SUBROUTINE ENSAMB_LAPLACE
    
    INTEGER I,INOD,J
    
    REAL*8 HCAR,DS1,AF

!--------------> LIMPIA LA MATRIZ DE RIGIDEZ
      
    DO INOD=1,NNOD
      DO J=1,NPUNTOS(INOD)*2
         S(J,INOD*2-1)=0.D0
         S(J,INOD*2)=0.D0
         S1(J,INOD*2-1)=0.D0
         S1(J,INOD*2)=0.D0
      END DO
    END DO
    
    DS1=(1.D0-POISSON)/2.d0
      
!--------------> ENSAMBLA LAS ECUACIONES EN EL DOMINIO
      
    DO INOD=1,NNOD
    
      ! CALL H_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCAR,HCARX,HCARY,X,Y)
      ! write(*,*) HCAR,HCARX,HCARY
       
      DO J=1,NPUNTOS(INOD)
         S(J*2-1,INOD*2-1)=(RM(4,J,INOD)+(DS1*RM(6,J,INOD)))
         S(J*2,INOD*2-1)=(POISSON*RM(5,J,INOD)+(DS1*RM(5,J,INOD)))
         S(J*2-1,INOD*2)=S(J*2,INOD*2-1)
         S(J*2,INOD*2)=(RM(6,J,INOD)+(DS1*RM(4,J,INOD)))
      END DO
    
    END DO
      
!----C-------> ENSAMBLA LAS ECUACIONES EN EL CONTORNO

    DO I=1,NNORM
    
        INOD=INORM_NODE(I)
        
        CALL H_CARACT(INOD,HCAR)
        !CALL HM_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCARM,X,Y)
        
        DO J=1,NPUNTOS(INOD)

!--------------------------------------------------------CCC
!C            ENSAMBLANDO H CARACTERISTICO                 C
!--------------------------------------------------------CCC            

!cccc-----> ec.contorno-h/2*(ec.dominio)=0

            AF=0.0

            S(J*2-1,INOD*2-1)=-((RM(4,J,INOD)+(DS1*RM(6,J,INOD)))*(AF*HCAR/2.D0))+ &
               & ((RM(2,J,INOD)*RNORM_VALUEX(I))+(DS1*RM(3,J,INOD)*RNORM_VALUEY(I)))

            S(J*2,INOD*2-1)=-((POISSON*RM(5,J,INOD)+(DS1*RM(5,J,INOD)))*(AF*HCAR/2.D0))+ &
               & ((POISSON*RM(3,J,INOD)*RNORM_VALUEX(I))+(DS1*RM(2,J,INOD)*RNORM_VALUEY(I)))

            S(J*2-1,INOD*2)=-((POISSON*RM(5,J,INOD)+(DS1*RM(5,J,INOD)))*(AF*HCAR/2.D0))+ &
               & ((DS1*RM(3,J,INOD)*RNORM_VALUEX(I))+(POISSON*RM(2,J,INOD)*RNORM_VALUEY(I)))
            
            S(J*2,INOD*2)=-((RM(6,J,INOD)+(DS1*RM(4,J,INOD)))*(AF*HCAR/2.D0))+ &
               & ((DS1*RM(2,J,INOD)*RNORM_VALUEX(I))+(RM(3,J,INOD)*RNORM_VALUEY(I)))


!--------------------------------------------------------CC
!C              SIN ENSAMBLAR H CARACTERISTICO            C
!--------------------------------------------------------CC

!c           S(J*2-1,INOD*2-1)=RM(2,J,INOD)*RNORM_VALUEX(I)+ &
!c     &     DS1*RM(3,J,INOD)*RNORM_VALUEY(I)

!c           S(J*2,INOD*2-1)=POISSON*RM(3,J,INOD)*RNORM_VALUEX(I)+ &
!c     &     DS1*RM(2,J,INOD)*RNORM_VALUEY(I)

!c           S(J*2-1,INOD*2)=DS1*RM(3,J,INOD)*RNORM_VALUEX(I)+ &
!c     &     POISSON*RM(2,J,INOD)*RNORM_VALUEY(I)
           
!c           S(J*2,INOD*2)=DS1*RM(2,J,INOD)*RNORM_VALUEX(I)+ &
!c     &     RM(3,J,INOD)*RNORM_VALUEY(I)
           
        END DO
    END DO
      
           
!------------> ENSAMBLA LAS ECUACIONES EN LOS NODOS
!------------> CON DESPLAZAMIENTOS PRESCRITOS

    DO I=1,NFIXP

        INOD=IFIXP_NODE(I)

        IF (BAN(1,I).EQ.1) THEN
            DO J=1,NPUNTOS(INOD)
                S(J*2-1,INOD*2-1)=RM(1,J,INOD)
                S(J*2,INOD*2-1)=0.0
            END DO
        END IF
           
        IF (BAN(2,I).EQ.1) THEN
            DO J=1,NPUNTOS(INOD)
                S(J*2-1,INOD*2)=0.0
                S(J*2,INOD*2)=RM(1,J,INOD)                
            END DO
        END IF

    END DO
    
    !open (unit=10,file="matriz de rigidez",iostat=k,action="write")
    !do inod=1,nnod*2
    !      write(10,'(18e17.7)') (s(j,inod),j=1,2*nnod)
    !end do
    !close (unit=10)      
      

    RETURN
    END SUBROUTINE ENSAMB_LAPLACE

!-----------------------------------------------!
!       PRECONDICIONA LA MATRIZ DE RIGIDEZ      !
!          Y EL TERMINO INDEPENDIENTE           !
!-----------------------------------------------!

    SUBROUTINE PRECF
    
    INTEGER INOD,J
    REAL*8 FS1,FS2
    
    DO INOD=1,NNOD
        FS1=S(1,INOD*2-1)
        FS2=S(2,INOD*2)         
        F(INOD*2-1)=F(INOD*2-1)/FS1
        F(INOD*2)=F(INOD*2)/FS2
        DO J=1,NPUNTOS(INOD)*2
            S(J,INOD*2-1)=S(J,INOD*2-1)/FS1
            S(J,INOD*2)=S(J,INOD*2)/FS2
        END DO
    END DO
    
    RETURN
    END SUBROUTINE PRECF


!-----------------------------------------------!
!         REALIZA CAMBIO DE VARIABLES           !
!-----------------------------------------------!
      SUBROUTINE CAMBIOVAR

      INTEGER INOD,I,J

      ALLOCATE(IFIXP_NODEN(NFIXP*2),RFIXP_VALUEN(NFIXP*2))
      ALLOCATE(CL_CONECTN(MAXCLD*2,NNOD*2),NPUNTOSN(NNOD*2))
      
      DO INOD=1,NFIXP
        IFIXP_NODEN(INOD*2-1)=IFIXP_NODE(INOD)*2-1
        IFIXP_NODEN(INOD*2)=IFIXP_NODE(INOD)*2
        RFIXP_VALUEN(INOD*2-1)=RFIXP_VALUE(INOD)
        RFIXP_VALUEN(INOD*2)=RFIXV_VALUE(INOD)
      END DO
      
      DO J=1,NNOD
        NPUNTOSN(J*2-1)=NPUNTOS(J)*2
        NPUNTOSN(J*2)=NPUNTOS(J)*2
         DO I=1,NPUNTOS(J)
            CL_CONECTN(I*2-1,J*2-1)=CL_CONECT(I,J)*2-1
            CL_CONECTN(I*2,J*2-1)=CL_CONECT(I,J)*2
            CL_CONECTN(I*2-1,J*2)=CL_CONECT(I,J)*2-1
            CL_CONECTN(I*2,J*2)=CL_CONECT(I,J)*2
         END DO
      END DO

      RETURN
      END SUBROUTINE CAMBIOVAR


!-------------------------------------------!
!     CALCULO DEL H CARACTERISTICO MINIMO   !
!-------------------------------------------!
    SUBROUTINE H_CARACT(INOD,HCAR)
    
    INTEGER INOD,J
    REAL*8 H,HCAR,HCARX,HCARY,HX,HY
            
    HCAR=1.D20        
    DO J=2,NPUNTOS(INOD)
        H=((((X(CL_CONECT(1,INOD))-X(CL_CONECT(J,INOD)))**2.+ &
          & (Y(CL_CONECT(1,INOD))-Y(CL_CONECT(J,INOD)))**2.)**.5))
        IF (H.LT.HCAR) HCAR=H
    END DO

    !HCARX=HCAR
    !HCARY=HCAR
    !
    !DO J=2,NPUNTOS(INOD)
    !   HX=DABS (X(CL_CONECT(J,INOD))-X(CL_CONECT(1,INOD)))
    !   HY=DABS (Y(CL_CONECT(J,INOD))-Y(CL_CONECT(1,INOD)))
    !   H=(HX**2+HY**2)**0.5
    !   IF (H.LT.HCAR) HCAR=H
    !   IF (HX.LT.HCARX) HCARX=HX
    !   IF (HY.LT.HCARY) HCARY=HY
    !   IF ((HX.GT.(HCAR/2)).AND.(HX.LT.HCARX)) HCARX=HX
    !   IF ((HY.GT.(HCAR/2)).AND.(HY.LT.HCARY)) HCARY=HY
    !END DO
        
    RETURN
    END SUBROUTINE H_CARACT

!----------------------------------------------!
!     SOLVER. GRADIENTES BICONJUGADOS PARA     !
!     MATRICES NO SIMETRICAS Y NO DEFINIDAS    !
!     POSITIVAS                                !
!                                              ! 
!     MATRIZ ALMACENADA EN FILAS               !
!                                    15/01/98  !
!----------------------------------------------!
    SUBROUTINE GRADCONJ(S,PRESS,B,NNOD,MAXCLD,CL_CONECT,NPUNTOS, &
     &  IFIX_NODE,NFIX,RFIX_VALUE)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: NNOD,NFIX,MAXCLD
    INTEGER,INTENT(IN) :: CL_CONECT(MAXCLD,NNOD),IFIX_NODE(NFIX),NPUNTOS(NNOD)
    REAL*8,INTENT(IN)  :: RFIX_VALUE(NFIX),S(MAXCLD,NNOD),B(NNOD)
    REAL*8,INTENT(OUT) :: PRESS(NNOD)
    
    REAL*8 RR_1,RR_2,CONJERR,ALF,BET,PKBAPK,FRR_1
    REAL*8 RES(NNOD),RESB(NNOD),PK(NNOD),PKB(NNOD),APK(NNOD),APKB(NNOD)

    INTEGER IFIX,IK,INOD,K
         
    CONJERR=1.D-14
    
    DO IFIX=1,NFIX
        PRESS(IFIX_NODE(IFIX))=RFIX_VALUE(IFIX)
    END DO
    
    CALL RESIDUO(RES,S,PRESS,CL_CONECT,NPUNTOS,NNOD,MAXCLD &
     &  ,IFIX_NODE,NFIX,RFIX_VALUE,0)
    
    DO INOD=1,NNOD 
        RES(INOD)=B(INOD)-RES(INOD) 
    END DO
    
    DO IFIX=1,NFIX 
        RES(IFIX_NODE(IFIX))=0.D0
    END DO
    
    DO INOD=1,NNOD
        RESB(INOD)=RES(INOD)
    END DO
    
    RR_1=0.D0
    DO INOD=1,NNOD
        RR_1=RR_1+RES(INOD)*RESB(INOD)
    END DO
    K=0
    DO WHILE (DABS(RR_1).GT.CONJERR.AND.K.LT.30000)
    
        K=K+1
                
        RR_1=0.D0
        DO INOD=1,NNOD
            RR_1=RR_1+RES(INOD)*RESB(INOD)
        END DO
        
        IF (K.EQ.1) THEN
            
            DO IK=1,NNOD
              PK(IK)=RES(IK)
              PKB(IK)=RESB(IK)
            END DO
            
        ELSE
        
            BET=RR_1/RR_2
            
            DO IK=1,NNOD
                PK(IK)=RES(IK)+BET*PK(IK)
                PKB(IK)=RESB(IK)+BET*PKB(IK)            
            END DO
            
        END IF
        
        CALL RESIDUO(APK,S,PK,CL_CONECT,NPUNTOS,NNOD,MAXCLD &
          &    ,IFIX_NODE,NFIX,RFIX_VALUE,0)        
        CALL RESIDUO(APKB,S,PKB,CL_CONECT,NPUNTOS,NNOD,MAXCLD &
          &    ,IFIX_NODE,NFIX,RFIX_VALUE,1)
        
        CALL SUBPKBAPK(PKBAPK,APK,PKB,NNOD)
        
        ALF=RR_1/PKBAPK
        
        DO IK=1,NNOD
          
            PRESS(IK)=PRESS(IK)+ALF*PK(IK)
            RES(IK)=RES(IK)-ALF*APK(IK)
            RESB(IK)=RESB(IK)-ALF*APKB(IK)
          
        END DO
                
        RR_2=RR_1
        write(20,*) k,rr_1
        
    END DO
    
    WRITE(*,'(A,I5)')'  ITERACIONES DE GBC....',K
    
    RETURN
    END

!------------------------------------!
!     FUNCION AUXILIAR DEL SOLVER    !
!     CALCULA EL RESIDUO             !      
!------------------------------------!      
    PURE SUBROUTINE RESIDUO(RES,S,PRESS,CL_CONECT,NPUNTOS,NNOD,MAXCLD &
    &  ,IFIX_NODE,NFIX,RFIX_VALUE,IFLAG)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: NNOD,NFIX,IFLAG,MAXCLD
    INTEGER,INTENT(IN) :: CL_CONECT(MAXCLD,NNOD),IFIX_NODE(NFIX),NPUNTOS(NNOD)
    REAL*8,INTENT(IN)  :: S(MAXCLD,NNOD),PRESS(NNOD),RFIX_VALUE(NFIX)
    REAL*8,INTENT(OUT) :: RES(NNOD)
    
    INTEGER IFIX,INOD,J
    
    DO INOD=1,NNOD
        RES(INOD)=0.D0
    END DO
    
    IF (IFLAG.EQ.0) THEN
    
        DO INOD=1,NNOD
            DO J=1,NPUNTOS(INOD)
                RES(INOD)=RES(INOD)+S(J,INOD)*PRESS(CL_CONECT(J,INOD))
            END DO
        END DO
        
    ELSE IF (IFLAG.EQ.1) THEN
    
        DO INOD=1,NNOD
          
            DO J=1,NPUNTOS(INOD)
                RES(CL_CONECT(J,INOD))=RES(CL_CONECT(J,INOD)) &
                 &+S(J,INOD)*PRESS(INOD)
            END DO
            
        END DO
        
    END IF
    
    DO IFIX=1,NFIX
        RES(IFIX_NODE(IFIX))=0.D0
    END DO      
    
    RETURN
    END

!----------------------------------------!
!     FUNCION AUXILIAR [2] DEL SOLVER    ! 
!----------------------------------------!
    PURE SUBROUTINE SUBPKBAPK(PKBAPK,APK,PK,NNOD)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN) :: NNOD
    REAL*8,INTENT(OUT) :: PKBAPK
    REAL*8,INTENT(IN)  :: APK(NNOD),PK(NNOD)
    
    INTEGER INOD

    PKBAPK=0.D0
    DO INOD=1,NNOD
        PKBAPK=PKBAPK+PK(INOD)*APK(INOD)
    END DO
    
    RETURN
    END


!:::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine BCG2

        integer, parameter :: l=1
        logical okprint,nonzero
        integer mxmv,info,ldw
        real*8 tol,rhs(NNOD*2),work(NNOD*2,2*l+5)
        character(3) typestop

        okprint = .true.
        tol = 1.d-14
        PRESS = 0.d0  ! initial guess
        nonzero = .false.
        typestop = 'max'
        mxmv = 100000
        ldw = NNOD*2*(2*l+5) ! o mas grande

! okprint == (input) LOGICAL. If okprint=.true. residual norm
!            will be printed to *
! l       == (input) INTEGER BiCGstab's dimension <= 2
!            Set l=2 for highly nonsymmetric problems
! n       == (input) INTEGER size of the system to solve 
! x       == (input/output) DOUBLE PRECISION array dimension n
!            initial guess on input, solution on output
! rhs     == (input) DOUBLE PRECISION array dimension n
!            right-hand side (rhs) vector
! matvec  == (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
! nonzero == (input) LOGICAL tells
!            BiCGstab(\ell) if initial guess in x is zero or not. 
!            If nonzero is .FALSE., one MATVEC call is saved.
! tol     == (input/output) DOUBLE PRECISION tolerance for all possible
!            stopping criteria (see the 'typestop' parameter)
!            On output, if info=0 or 1, tol is actually achieved
!            residual reduction or whatever (see the 'typestop' parameter)
! typestop== (input) CHARACTER*3 stopping criterion (||.|| denotes 
!            the 2-norm):
!            typestop='rel' -- relative stopping crit.: ||res|| < tol*||res0||
!            typestop='abs' -- absolute stopping crit.: ||res||<tol
!            typestop='max' -- maximum  stopping crit.: max(abs(res))<tol
! NOTE(for typestop='rel' and 'abs'): To save computational work, the value of 
!            residual norm used to check the convergence inside the main iterative 
!            loop is computed from 
!            projections, i.e. it can be smaller than the true residual norm
!            (it may happen when e.g. the 'matrix-free' approach is used).
!            Thus, it is possible that the true residual does NOT satisfy
!            the stopping criterion ('rel' or 'abs').
!            The true residual norm (or residual reduction) is reported on 
!            output in parameter TOL -- this can be changed to save 1 MATVEC
!            (see comments at the end of the subroutine)
! mxmv   ==  (input/output) INTEGER.  On input: maximum number of matrix 
!            vector multiplications allowed to be done.  On output: 
!            actual number of matrix vector multiplications done
! work   ==  (workspace) DOUBLE PRECISION array dimension (n,2*l+5))
! ldw    ==  (input) INTEGER size of work, i.e. ldw >= n*(2*l+5)
! info   ==  (output) INTEGER.  info = 0 in case of normal computations
!            and 
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (try to enlarge
!            parameter l=\ell to get rid of this)
! ----------------------------------------------------------
        call bicgstab2 (okprint,l,NNOD*2,PRESS,F,matvec,nonzero,tol,&
                        &typestop,mxmv,work,ldw,info)

        return
    end subroutine BCG2

    subroutine matvec(n,x,y)
    !using as EXTERNAL matrix vector subroutine
    !to deliver y:=A*x by CALL matvec(n,x,y)
        integer,intent(in) :: n
        real*8,intent(in)  :: x(n)
        real*8,intent(out) :: y(n)
        integer i,j
        y = 0.d0
        do i = 1, nnod*2
            do j=1,npuntosn(i)
                y(i) = y(i) + S(j,i)*x(CL_CONECTN(j,i))
            end do
        end do
        return
    end subroutine matvec


!----------------------------------------------------!
!     EVALUA UNA FUNCION VECTORIAL Y SUS DERIVADA    !
!----------------------------------------------------!
    PURE SUBROUTINE EVAVECT(INOD,VALX,VALY,ITIPE)
    
    INTEGER,INTENT(IN) :: INOD,ITIPE
    REAL*8,INTENT(OUT) :: VALX,VALY
    INTEGER J

!--------------> ITIPE = 1  FUNCION DE FORMA
!--------------> ITIPE = 2  DERIVADA RESPECTO A X   d /dx
!--------------> ITIPE = 3  DERIVADA RESPECTO A Y   d /dy
!--------------> ITIPE = 4  DERIVADA SEGUNDA RESPECTO A  X   d2 /dx2
!--------------> ITIPE = 5  DERIVADA SEGUNDA RESPECTO A X,Y  d2 /dxdy
!--------------> ITIPE = 6  DERIVADA SEGUNDA RESPECTO A  Y   d2 /dy2

    VALX=0.D0
    VALY=0.D0
    DO J=1,NPUNTOS(INOD)
        VALX=VALX+RM(ITIPE,J,INOD)*DU(CL_CONECT(J,INOD))
        VALY=VALY+RM(ITIPE,J,INOD)*DV(CL_CONECT(J,INOD))
    END DO
    
    RETURN
    END


!:::::::::::::::::::::::::::::::::::::::::::::::::::
!:::::::::::::::::::::::::::::::::::::::::::::::::::

!----------------------------------------------------------!
!     CALCULO DE LAS FUNCIONES DE FORMA Y SUS DERIVADAS    !
!----------------------------------------------------------!

    subroutine N_DN_D2N
        integer i
        CALL N_DN_D2N_FWLS
        CALL N_DN_D2N_MaxEnt
        allocate( RM(NCONS,MAXCLD,NNOD) ) ; RM=0.d0
        do i = 1,nnod
            if ( any( i .eq. inorm_node) ) then ! pertenece al contorno : fwls
               RM(1:NCONS,1:MAXCLD,i) = RM_FWLS(1:NCONS,1:MAXCLD,i) 
            else ! pertenece al interior
               RM(1:NCONS,1:MAXCLD,i) = RM_MaxEnt(1:NCONS,1:MAXCLD,i) 
               !RM(1:NCONS,1:MAXCLD,i) = RM_FWLS(1:NCONS,1:MAXCLD,i) 
            end if
            !print*, 'inod (nodo estrella) = ',i  ;   print*,'phi:'
            !print*, RM(1,:,i)                    ;   print*, ''
        end do

        ! despejar espacio
        deallocate(RM_FWLS,RM_MaxEnt)
        
        return
    end subroutine N_DN_D2N


    SUBROUTINE N_DN_D2N_MaxEnt
        use maxent , only : drivermaxent
        integer maxit,ierror,inod,j,np
        character(80) scheme,priorwt
        logical cerror,maxentprint
        real*8 charlengthscale,dettol,eps,p(2),hcar
        real*8,allocatable :: phi(:),dphi(:,:),ddphi(:,:,:),&
            &dmetric(:,:,:),xyz(:,:),dmax_(:)

        ALLOCATE( RM_MaxEnt(NCONS,MAXCLD,NNOD) ) ; RM_MaxEnt = 0.d0
        ! NOTA : DMAX calculado anteriormente en N_DN_D2N

                                    ! OPCIONES DISPONIBLES:
        scheme = 'newton'           ! descent,lbfgs,newton
        priorwt = 'gaussian-rbf'    ! uniform,cubi,quartic,gaussian,gaussian-rbf

        eps    = 1.d-14      ! convergence tolerance 
        dettol = 1.d-16      ! determinant tolerance with default of 10d−16
        maxit  = 10000       ! maximum number of iteration
        maxentprint = .true. ! print convergence information

        do inod = 1,nnod ! recorre cada nodo estrella del dominio 

            np = npuntos(inod)
        
            ! asignar memoria
            allocate(phi(np),dphi(2,np),ddphi(2,2,np),&
                &dmetric(2,2,np),xyz(np,2),dmax_(np))
        
            xyz=0.d0 ; phi=0.d0 ; dphi=0.d0 ; ddphi=0.d0 ! init var

            ! length scale for domain so that convergence criterion 
            ! is independent of the magnitude of nodal coordinates.
            call H_CARACT(INOD,HCAR)
            charlengthscale = hcar
            
            p = (/ x(cl_conect(1,inod)) , y(cl_conect(1,inod)) /) ! punto estrella
            do j = 1,np
                dmetric(1:2,1:2,j)=reshape((/1.d0,0.d0,0.d0,1.d0/),(/2,2/))    ! isotropic
                xyz(j,1:2) = (/ x(cl_conect(j,inod)) , y(cl_conect(j,inod)) /) ! coord cld

                dmax_(j) = dmax(cl_conect(j,inod))
            end do

            if ( any( inod .eq. inorm_node ) ) then ! is boundary?
                ! WARNING : Matrix is nearly singular
                !call  drivermaxent(np,2,scheme,priorwt,xyz,p,dmax_,dmetric,&
                !                    & maxit,eps,maxentprint,ierror,cerror,&
                !                    & charlengthscale=charlengthscale,&
                !                    & dettol=dettol,phi=phi)
                phi = 0.d0 ; dphi = 0.d0 ; ddphi = 0.d0
            else 
                call  drivermaxent(np,2,scheme,priorwt,xyz,p,dmax_,dmetric,&
                                    & maxit,eps,maxentprint,ierror,cerror,&
                                    & charlengthscale=charlengthscale,&
                                    & dettol=dettol,phi=phi,dphi=dphi,ddphi=ddphi)
            end if

            RM_MaxEnt(1,1:np,inod) = phi(1:np)
            RM_MaxEnt(2,1:np,inod) = dphi(1,1:np)
            RM_MaxEnt(3,1:np,inod) = dphi(2,1:np)
            RM_MaxEnt(4,1:np,inod) = ddphi(1,1,1:np)
            RM_MaxEnt(5,1:np,inod) = ddphi(2,1,1:np) ! o ddphi(1,2,1:np)
            RM_MaxEnt(6,1:np,inod) = ddphi(2,2,1:np)

            ! despejar memoria
            deallocate(phi,dphi,ddphi,dmetric,xyz,dmax_)

        end do

        return
    END SUBROUTINE N_DN_D2N_MaxEnt


    SUBROUTINE N_DN_D2N_FWLS
    
        INTEGER INOD,I,J,IERR,J1,J2,J3
        
        REAL*8 DISMAX,DX,DY,DIS,BETA,ALFA,PED,PEX,SUMA,XSTAR,YSTAR
        REAL*8 C(MAXCLD,NCONS),D(NCONS,NCONS),RM_P(NCONS,MAXCLD),W(MAXCLD)

        ALLOCATE( RM_FWLS(NCONS,MAXCLD,NNOD) , DMAX(NNOD) )

        !--------------> LIMPIA LAS FUNCIONES DE FORMA
        RM_FWLS=0.D0
        
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
                
                DX=DX/DMAX(INOD)
                DY=DY/DMAX(INOD)
                C(J,1)=1.D0 
                C(J,2)=DX
                C(J,3)=DY
                C(J,4)=(DX)**2.D0 
                C(J,5)=DX*DY 
                C(J,6)=(DY)**2.D0
                !C(J,7)=DX**3.D0
                !C(J,8)=(DX**2.D0)*DY
                !C(J,9)=(DY**2.D0)*DX
                !C(J,10)=DY**3.D0
            
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
                PRINT*, 'PROBLEMAS EN LAS NUBES... inod = ',inod
                STOP 
            
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
            
            !--------------> GUARDA LOS VALORES NECESARIOS
            
            DO I=1,NPUNTOS(INOD)
                RM_FWLS(1,I,INOD)=RM_P(1,I)
                RM_FWLS(2,I,INOD)=RM_P(2,I)/dmax(inod)
                RM_FWLS(3,I,INOD)=RM_P(3,I)/dmax(inod)
                RM_FWLS(4,I,INOD)=RM_P(4,I)*2/dmax(inod)**2
                RM_FWLS(5,I,INOD)=RM_P(5,I)/dmax(inod)**2
                RM_FWLS(6,I,INOD)=RM_P(6,I)*2/dmax(inod)**2
            END DO
        
            !write(*,*) 'funcion de forma, inod =',inod
            !write(*,*) rm(1,1:npuntos(inod),inod) ; write(*,*)
        END DO
        
        RETURN
    END SUBROUTINE N_DN_D2N_FWLS

!-----------------------------------------------------!
!    INVIERTE UNA MATRIZ DE (NCONSxNCONS) Y LA        ! 
!     DEVUELVE EN LA MISMA MATRIZ DE ENTRADA          ! 
!-----------------------------------------------------!
    SUBROUTINE INVERT(NCONS,IERR,D) 
    
        INTEGER, INTENT(IN)   :: NCONS
        INTEGER, INTENT(OUT)  :: IERR
        REAL*8, INTENT(INOUT) :: D(NCONS,NCONS)
        
        INTEGER I,J,K
        REAL*8 RM_,RMAX,RLIM
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
            
                RM_=D1(I,J)/D1(I,I) 
                DO K=J,2*NCONS
                    D1(J,K)=D1(J,K)-RM_*D1(I,K)                
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
            END DO 
        END DO  
        
        RETURN 
    END SUBROUTINE INVERT

end program Finite_Point_Method_Solver



