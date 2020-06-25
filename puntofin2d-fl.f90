     
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                   C
!C     METODO NUMERICO SIN MALLA DE MINIMOS CUADRADOS PONDERADOS     C
!C             FUNCION BASE DE INTERPOLACION CUADRATICA              C
!C                                                                   C
!C                           2-DIMENSIONES                           C
!C                           TENSION PLANA                           C
!C                                                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (NODOS=20000)          !NUMERO MAXIMO DE PUNTOS
      PARAMETER (IFIJOS=20000)         !NUMERO MAXIMO DE PUNTOS PRESCRITOS
      PARAMETER (ITERAC=20000)         !NUMERO MAXIMO DE ITERACIONES

      CHARACTER FILE*80

      INTEGER NNOD ,NFIXP,NFIXF,NNORM,ILONG,IWEIGT
      INTEGER CL_CONECT(16,NODOS),NPUNTOS(NODOS)
      INTEGER IFIXP_NODE(IFIJOS),IFIXF_NODE(IFIJOS)
      INTEGER INORM_NODE(IFIJOS),BAN(6,IFIJOS)

      REAL*8 BET,YOUNG,POISSON,dmax(nodos)      
      REAL*8 X(NODOS),Y(NODOS),RM(10,16,NODOS),RMA(1,16,NODOS)
      REAL*8 PRESS(NODOS),U(NODOS),V(NODOS)
      REAL*8 RFIXP_VALUE(IFIJOS),RFIXV_VALUE(IFIJOS),RES(ITERAC)
      REAL*8 RNORM_VALUEY(IFIJOS),RNORM_VALUEX(IFIJOS)
      REAL*8 RFIXF_VALUEX(IFIJOS),RFIXF_VALUEY(IFIJOS)

!CCCC----------> LEE EL NOMBRE DEL PROBLEMA
      
      OPEN(1,FILE='NUBESPUNT.DAT',STATUS='OLD')
      READ(1,'(A)') FILE
      CLOSE(1)
   
!CCCC----------> LECTURA DE DATOS

!CCCC-----> LEE DATOS GENERALES      

      CALL READ_GENERAL_DATA(ILONG,IWEIGT,BET,YOUNG,POISSON,FILE)
            
!CCCC-----> LEE DATOS DE GEOMETRIA Y COND. DE CONTORNO
      
      CALL READ_GEOM_DATA(NODOS,IFIJOS,NNOD,NFIXP,NFIXF,NNORM,ILONG &
     &   ,IFIXP_NODE,IFIXF_NODE,INORM_NODE,X,Y,BAN,RFIXP_VALUE, &
     &   RFIXV_VALUE,RNORM_VALUEX,RNORM_VALUEY,RFIXF_VALUEX, &
     &          RFIXF_VALUEY,FILE)

!CCCC----------> LECTURA DE NUBES PARA CALCULO DE LA FUNCION DE FORMA PARA ALISADO

      CALL CLOUDS_DATA(ILONG,NNOD,CL_CONECT,NPUNTOS,X,Y,FILE,0)
!
!!CCCC----------> CALCULO DE LA FUNCION DE FORMA PARA ALISADO
!      
!      NCONS=6
!      iweigt=0
!      CALL N_DN_D2N(NCONS,NNOD,IWEIGT,NPUNTOS,CL_CONECT,X,Y,dmax,RM)
!      DO INOD=1,NNOD
!        DO I=1,NPUNTOS(INOD)
!           RMA(1,I,INOD)=RM(1,I,INOD)
!!c           write(*,*) rma(1,i,inod)
!        END DO
!      END DO

!CCCC----------> LECTURA DE NUBES PARA CALCULO DE LA FUNCION DE FORMA Y SUS DERIVADAS
      
      CALL CLOUDS_DATA(ILONG,NNOD,CL_CONECT,NPUNTOS,X,Y,FILE,0)
      
!CCCC----------> CALCULO DE LAS FUNCIONES DE FORMA Y SUS DERIVADAS

      NCONS=6
      iweigt=0
      CALL N_DN_D2N(NCONS,NNOD,IWEIGT,NPUNTOS,CL_CONECT,X,Y,dmax,RM)
      
!CCCC----------> CALCULO DEL FLUJO POTENCIAL
      
      CALL SOLVE_POT(ILONG,NNOD,NNORM,NFIXF,NFIXP,NPUNTOS,CL_CONECT &
     &  ,INORM_NODE,IFIXP_NODE,IFIXF_NODE,U,V,PRESS,X,Y,dmax,RM,RMA &
     &  ,RNORM_VALUEX,RNORM_VALUEY,RFIXF_VALUEX,RFIXF_VALUEY,BAN,RFIXP_VALUE &
     &  ,RFIXV_VALUE,POISSON,YOUNG,FILE)
      
      END

!CCCC----------------------------------------------------------------CCCC
!C    CALCULA LA CANTIDAD DE LETRAS DEL NOMBRE DEL ARCHIVO DE ENTRADA   C
!CCCC----------------------------------------------------------------CCCC
      FUNCTION LONG_FILE(FILE)

      CHARACTER FILE*80

      INTEGER LONG_FILE
      
      DO I=1,80
        IF (FILE(I:I).EQ.' ') THEN
          LONG_FILE=I-1
          RETURN
        END IF
      END DO

      WRITE(10,'(A,A)') 'ERROR...( hay un error con el nombre del file)'
      STOP
      
      END

!CCCC--------------------------------------------CCCC
!C    SUBRUTINA DE LECTURA DE DATOS DEL PROBLEMA    C
!CCCC--------------------------------------------CCCC      
      SUBROUTINE READ_GENERAL_DATA(ILONG,IWEIGT,BET,YOUNG,POISSON,FILE)
      
      CHARACTER FILE*80

      INTEGER ILONG,LONG_FILE,IWEIGT

      REAL*8 BET,YOUNG,POISSON
      
      ILONG=LONG_FILE(FILE)

!CCCC----------> LEE DATOS GENERALES DEL PROBLEMA
      
      OPEN(1,FILE=FILE(1:ILONG)//'.GEN',STATUS='OLD')

      READ(1,*)
!CCCC-----> TIPO DE FUNCION DE PONDERACION            
      READ(1,*) IWEIGT
      READ(1,*)
!CCCC-----> PARAMETRO DEL SOLVER
      READ(1,*) BET  
      READ(1,*)
!CCCC-----> PROPIEDADES DEL MATERIAL
      READ(1,*)
      READ(1,*) YOUNG
      READ(1,*) 
      READ(1,*) POISSON    

      CLOSE(1)

      RETURN
      END

!CCCC--------------------------------------------CCCC
!C    SUBRUTINA DE LECTURA DE DATOS DEL PROBLEMA    C
!CCCC--------------------------------------------CCCC      
      SUBROUTINE READ_GEOM_DATA(NODOS,IFIJOS,NNOD,NFIXP,NFIXF &
     &  ,NNORM,ILONG,IFIXP_NODE,IFIXF_NODE,INORM_NODE,X,Y &
     &  ,BAN,RFIXP_VALUE,RFIXV_VALUE,RNORM_VALUEX,RNORM_VALUEY &
     &  ,RFIXF_VALUEX,RFIXF_VALUEY,FILE)

      CHARACTER FILE*80
      
      INTEGER NNOD,NFIX,ILONG,NFIXP,NFIXV,NNORM
      INTEGER IFIXP_NODE(IFIJOS),IFIXF_NODE(IFIJOS)
      INTEGER INORM_NODE(IFIJOS),BAN(6,IFIJOS)
      
      REAL*8 X(NODOS),Y(NODOS)
      REAL*8 RFIXP_VALUE(IFIJOS),RFIXV_VALUE(IFIJOS)
      REAL*8 RFIXF_VALUEX(IFIJOS),RFIXF_VALUEY(IFIJOS)
      REAL*8 RNORM_VALUEX(IFIJOS),RNORM_VALUEY(IFIJOS)
      
!CCCC----------> LEE DATOS DE LA GEOMETRIA

      OPEN(1,FILE=FILE(1:ILONG)//'.dat',STATUS='OLD')

!CCCC-----> DATOS DE CONTROL

	READ(1,*)
	READ(1,*)NNOD,NFIXP,NFIXF

!CCCC-----> COORDENADAS DE LOS PUNTOS
      READ(1,*)
      I=0
      DO WHILE (I.LT.NNOD) 
        I=I+1 
        READ(1,*) NOD,X(NOD),Y(NOD)
      END DO

!CCCC-----> CONDICIONES DE CONTORNO
!CCCC-----> PUNTO CON VALOR DE DESPLAZAMIENTO PRESCRITO
      I=0
      READ(1,*)
      DO WHILE (I.LT.NFIXP) 
        I=I+1 
        READ(1,*) IFIXP_NODE(I),BAN(1,I),BAN(2,I),RFIXP_VALUE(I), &
     &  RFIXV_VALUE(I)
      END DO
      
!CCCC---> PUNTO CON VALOR DE FUERZA PRESCRITA
      I=0
      READ(1,*)
      DO WHILE (I.LT.NFIXF) 
        I=I+1 
        READ(1,*) IFIXF_NODE(I),RFIXF_VALUEX(I),RFIXF_VALUEY(I)
      END DO

!CCCC---> NODO CON VALOR NORMAL PRESCRITO
      I=0
      READ(1,*)
	READ(1,*)NNORM
      DO WHILE (I.LT.NNORM) 
        I=I+1 
        READ(1,*) INORM_NODE(I),RNORM_VALUEX(I),RNORM_VALUEY(I)
      END DO

      CLOSE(1)

	WRITE(*,*)NNOD,NFIXP,NFIXF,NNORM
!c	STOP

      CLOSE(1)

      RETURN
      END
      
!CCCC-------------------------------CCCC
!C    LECTURA Y GENERACION DE NUBES    C
!CCCC-------------------------------CCCC      
      SUBROUTINE CLOUDS_DATA(ILONG,NNOD,CL_CONECT &
     &  ,NPUNTOS,X,Y,FILE,NPN)

      CHARACTER FILE*80

      INTEGER ILONG,NNOD,NPN
      INTEGER CL_CONECT(16,NNOD),NPUNTOS(NNOD)

      REAL*8 X(NNOD),Y(NNOD)

      CALL READ_CLOUDS(ILONG,NNOD,CL_CONECT,NPUNTOS,FILE,NPN)

      RETURN
      END

      
!CCCC------------------------------------------------------------CCCC
!C    LEE LAS NUBES DE PUNTOS DE UN ARCHIVO CON EXTENSION ".CLD"    C
!C     CON EL SIGUIENTE FORMATO:  NPTOS , PT_1 , PT_2 , PT_3 ,      C
!C     ... , PT_NPTOS                                               C
!CCCC------------------------------------------------------------CCCC
      SUBROUTINE READ_CLOUDS(ILONG,NNOD,CL_CONECT,NPUNTOS,FILE,NPN)

      CHARACTER FILE*80
      
      INTEGER ILONG,NNOD,NPN
      INTEGER CL_CONECT(16,NNOD),NPUNTOS(NNOD)

      IF (NPN.EQ.0) THEN
         
         OPEN(1,FILE=FILE(1:ILONG)//'.CLD',STATUS='OLD')

      ELSE IF (NPN.EQ.1) THEN

         OPEN(1,FILE=FILE(1:ILONG)//'.CLD1',STATUS='OLD')

      END IF
      
      DO INOD=1,NNOD

        READ(1,*) NPUNTOS(INOD),(CL_CONECT(J,INOD),J=1,NPUNTOS(INOD))

      END DO
      CLOSE(1)

      RETURN
      END

!CCCC---------------------------------------------------CCCC
!C    CALCULO DE LAS FUNCIONES DE FORMA Y SUS DERIVADAS    C
!CCCC---------------------------------------------------CCCC
      SUBROUTINE N_DN_D2N(NCONS,NNOD,IWEIGT,NPUNTOS,CL_CONECT,X,Y,dmax,RM)

      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD),IWEIGT,NCONS

      REAL*8 DISMAX,DX,DY,DIS,BETA,ALFA,PED,PEX
      REAL*8 X(NNOD),Y(NNOD),RM(10,16,NNOD),dmax(nnod)
      REAL*8 C(16,10),D(10,10),RM_P(10,16),W(16)
      
      DO INOD=1,NNOD 

        XSTAR=X(CL_CONECT(1,INOD))
        YSTAR=Y(CL_CONECT(1,INOD))
        
!CCCC----------> CALCULO DE LA DISTANCIA AL NODO MAS LEJANO DE LA NUBE   
          
	  dmax(inod)=0.d0
        DO J=1,NPUNTOS(INOD)  
          
          DX=X(CL_CONECT(J,INOD))-XSTAR 
          DY=Y(CL_CONECT(J,INOD))-YSTAR 
          DIS=SQRT(DX**2+DY**2) 
          IF (dmax(inod).LT.DIS) dmax(inod)=DIS 
          
        END DO   
        
!CCCC----------> CALCULO DE LOS PARAMETROS DE LA FUNCION DE PESO
        
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
!c            W(J)=.020D0
!c            IF (J.EQ.1) W(J)=1.D0

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
        
!CCCC----------> CALCULO DE LA MATRIZ  "CT * C" 
        
        DO J1=1,NCONS
          DO J2=1,NCONS 
            SUM=0 
            DO J3=1,NPUNTOS(INOD) 
              SUM=SUM+C(J3,J2)*C(J3,J1)*W(J3) 
            END DO 
            D(J1,J2)=SUM 
          END DO 
        END DO

        
!CCCC----------> INVIERTE LA MATRIZ  "D"  
        
        CALL INVERT(NCONS,IERR,D)
        IF (IERR.EQ.1) THEN          
          STOP 'PROBLEMAS EN LAS NUBES...'
!c         write(*,*) 'PROBLEMAS EN LAS NUBES...',inod
        END IF 
        
!CCCC----------> CALCULO DE  "D(-1) * CT" 
        
        DO J1=1,NCONS
          DO J2=1,NPUNTOS(INOD)
            SUM=0 
            DO J3=1,NCONS 
              SUM=SUM+D(J1,J3)*C(J2,J3) 
            END DO 
            RM_P(J1,J2)=SUM*W(J2)
          END DO
        END DO
        
!CCCC----------> LIMPIA LAS FUNCIONES DE FORMA

        DO I=1,NPUNTOS(INOD)
          DO J=1,10
            RM(J,I,INOD)=0.D0
          END DO                    
        END DO
        
!CCCC----------> GUARDA LOS VALORES NECESARIOS
        
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


!CCCC------------------------------------------------CCCC
!C    INVIERTE UNA MATRIZ DE (NCONSxNCONS) Y LA         C
!C     DEVUELVE EN LA MISMA MATRIZ DE ENTRADA           C
!CCCC------------------------------------------------CCCC      
      SUBROUTINE INVERT(NCONS,IERR,D) 

      REAL*8 RMAX,RLIM
      REAL*8 D(10,10),D1(10,20) 
 
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
!c            write(*,*) d(i,j)
         END DO 
      END DO  
    
      RETURN 
      END 
      
!CCCC---------------------------------CCCCCCC
!C    CALCULO DEL H CARACTERISTICO MINIMO   C
!CCCC---------------------------------CCCCCCC
      SUBROUTINE H_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCAR,X,Y)

      INTEGER INOD,NNOD,J
      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD)
      
      REAL*8 H,HCAR,HCARX,HCARY,HX,HY
      REAL*8 X(NNOD),Y(NNOD)
              
      HCAR=1.D20        
      DO J=2,NPUNTOS(INOD)
        H=((((X(CL_CONECT(1,INOD))-X(CL_CONECT(J,INOD)))**2.+ &
     &  (Y(CL_CONECT(1,INOD))-Y(CL_CONECT(J,INOD)))**2.)**.5))
        IF (H.LT.HCAR) HCAR=H
      END DO

!c      HCARX=HCAR
!c      HCARY=HCAR
      
!C      DO J=2,NPUNTOS(INOD)
!C         HX=DABS (X(CL_CONECT(J,INOD))-X(CL_CONECT(1,INOD)))
!C         HY=DABS (Y(CL_CONECT(J,INOD))-Y(CL_CONECT(1,INOD)))
!C         H=(HX**2+HY**2)**0.5
!C         IF (H.LT.HCAR) HCAR=H
!C         IF (HX.LT.HCARX) HCARX=HX
!C         IF (HY.LT.HCARY) HCARY=HY
!C         IF ((HX.GT.(HCAR/2)).AND.(HX.LT.HCARX)) HCARX=HX
!C         IF ((HY.GT.(HCAR/2)).AND.(HY.LT.HCARY)) HCARY=HY
!C      END DO
        
      RETURN
      END

!CCCC-----------------------------------CCCC
!C    CALCULO DEL H CARACTERISTICO MAXIMO  C
!CCCC-----------------------------------CCCC
      SUBROUTINE HM_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCARM,X,Y)

      INTEGER INOD,NNOD,J
      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD)
      
      REAL*8 H,HCARM,HX,HY
      REAL*8 X(NNOD),Y(NNOD)
              
      HCARM=1.D-20        
      DO J=2,NPUNTOS(INOD)
        H=((((X(CL_CONECT(1,INOD))-X(CL_CONECT(J,INOD)))**2.+ &
     &  (Y(CL_CONECT(1,INOD))-Y(CL_CONECT(J,INOD)))**2.)**.5))
        IF (H.GT.HCARM) HCARM=H
      END DO
        
      RETURN
      END

!CCCC-------------------------------------CCCC
!C    ENSAMBLA LAS ECUACIONES                C
!C    EN EL DOMINIO Y EN LOS CONTORNOS       C
!CCCC-------------------------------------CCCC
      SUBROUTINE ENSAMB_LAPLACE(NNOD,NFIXP,NNORM,INORM_NODE,NPUNTOS &
     & ,IFIXP_NODE,CL_CONECT,BAN,RNORM_VALUEX,RNORM_VALUEY,RM,S,S1,X,Y,dmax,YOUNG &
     & ,POISSON,DS1,DS2)
      
      INTEGER NNOD,NNORM,NFIXP
      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD),INORM_NODE(NNORM)
      INTEGER IFIXP_NODE(NFIXP),BAN(6,NFIXP)
      
      REAL*8 HCAR,DS1,DS2,YOUNG,POISSON,AF
      REAL*8 S(32,NNOD*2),S1(32,NNOD*2),RM(10,16,NNOD),X(NNOD),Y(NNOD),dmax(nnod)
      REAL*8 RNORM_VALUEX(NNORM),RNORM_VALUEY(NNORM)

!CCCC----------> LIMPIA LA MATRIZ DE RIGIDEZ
      
      DO INOD=1,NNOD
        DO J=1,NPUNTOS(INOD)*2
           S(J,INOD*2-1)=0.D0
           S(J,INOD*2)=0.D0
           S1(J,INOD*2-1)=0.D0
           S1(J,INOD*2)=0.D0
        END DO
      END DO

      DS1=(1-POISSON)/2
      
!CCCC----------> ENSAMBLA LAS ECUACIONES EN EL DOMINIO
      
      DO INOD=1,NNOD

!c        CALL H_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCAR,HCARX,HCARY,X,Y)

!c        write(*,*) HCAR,HCARX,HCARY
         
         
        DO J=1,NPUNTOS(INOD)
        
!c        k=cl_conect(j,inod)
        
           S(J*2-1,INOD*2-1)=(RM(4,J,INOD)+(DS1*RM(6,J,INOD)))
           S(J*2,INOD*2-1)=(POISSON*RM(5,J,INOD)+(DS1*RM(5,J,INOD)))
           S(J*2-1,INOD*2)=S(J*2,INOD*2-1)
           S(J*2,INOD*2)=(RM(6,J,INOD)+(DS1*RM(4,J,INOD)))

!c           write(*,*) s(j,inod),rm(1,j,inod)
        END DO
      END DO
      
!CCCCC-------> ENSAMBLA LAS ECUACIONES EN EL CONTORNO

      DO I=1,NNORM
         INOD=INORM_NODE(I)

         CALL H_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCAR,X,Y)
!C         CALL HM_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCARM,X,Y)
         

         DO J=1,NPUNTOS(INOD)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C            ENSAMBLANDO H CARACTERISTICO                 C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            

!cccc-----> ec.contorno-h/2*(ec.dominio)=0

           AF=0.0
		            
           S(J*2-1,INOD*2-1)=-((RM(4,J,INOD)+(DS1*RM(6,J,INOD)))*(AF*HCAR/2.D0))+ &
     &               ((RM(2,J,INOD)*RNORM_VALUEX(I))+ &
     &               (DS1*RM(3,J,INOD)*RNORM_VALUEY(I)))

           S(J*2,INOD*2-1)=-((POISSON*RM(5,J,INOD)+(DS1*RM(5,J,INOD)))*(AF*HCAR/2.D0))+ &
     &               ((POISSON*RM(3,J,INOD)*RNORM_VALUEX(I))+ &
     &               (DS1*RM(2,J,INOD)*RNORM_VALUEY(I)))

           S(J*2-1,INOD*2)=-((POISSON*RM(5,J,INOD)+(DS1*RM(5,J,INOD)))*(AF*HCAR/2.D0))+ &
     &               ((DS1*RM(3,J,INOD)*RNORM_VALUEX(I))+ &
     &               (POISSON*RM(2,J,INOD)*RNORM_VALUEY(I)))
           
           S(J*2,INOD*2)=-((RM(6,J,INOD)+(DS1*RM(4,J,INOD)))*(AF*HCAR/2.D0))+ &
     &               ((DS1*RM(2,J,INOD)*RNORM_VALUEX(I))+ &
     &               (RM(3,J,INOD)*RNORM_VALUEY(I)))


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C              SIN ENSAMBLAR H CARACTERISTICO            C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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
      
           
!CCCCC-------> ENSAMBLA LAS ECUACIONES EN LOS NODOS
!CCCCC-------> CON DESPLAZAMIENTOS PRESCRITOS

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
      END

      
!CCCC---------------------------------------CCCC
!C    SOLVER. GRADIENTES BICONJUGADOS PARA     C
!C    MATRICES NO SIMETRICAS Y NO DEFINIDAS    C
!C    POSITIVAS                                C
!C                                             C 
!C    MATRIZ ALMACENADA EN FILAS               C
!C                                   15/01/98  C
!CCCC---------------------------------------CCCC
      SUBROUTINE GRADCONJ(S,PRESS,B,NNOD,CL_CONECT,NPUNTOS, &
     &  IFIX_NODE,NFIX,RFIX_VALUE,PK,PKB,APK,APKB,RES,RESB)
      
      INTEGER NNOD,NFIX
      INTEGER CL_CONECT(32,NNOD),IFIX_NODE(NFIX),NPUNTOS(NNOD*2)

      REAL*8 RR_1,RR_2,CONJERR,ALF,BET,PKBAPK,FRR_1
      REAL*8 RES(NNOD),RESB(NNOD),PK(NNOD),PKB(NNOD),APK(NNOD) &
     & ,APKB(NNOD)
      REAL*8 B(NNOD),PRESS(NNOD),S(32,NNOD),RFIX_VALUE(NFIX)
           
      !CONJERR=1.D-6
      !CONJERR=1.D-14
      CONJERR=1.D-18

      DO IFIX=1,NFIX
        PRESS(IFIX_NODE(IFIX))=RFIX_VALUE(IFIX)
      END DO
      
      CALL RESIDUO(RES,S,PRESS,CL_CONECT,NPUNTOS,NNOD &
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
      
      RR_1=FRR_1(RES,RESB,NNOD)
      K=0
      DO WHILE (DABS(RR_1).GT.CONJERR.AND.K.LT.30000)

        K=K+1
                
        RR_1=FRR_1(RES,RESB,NNOD)
        
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
        
        CALL RESIDUO(APK,S,PK,CL_CONECT,NPUNTOS,NNOD &
     &    ,IFIX_NODE,NFIX,RFIX_VALUE,0)        
        CALL RESIDUO(APKB,S,PKB,CL_CONECT,NPUNTOS,NNOD &
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

      
!CCCC-----------------------------CCCC
!C    FUNCION AUXILIAR DEL SOLVER    C
!C    CALCULA EL RESIDUO             C      
!CCCC-----------------------------CCCC      
      SUBROUTINE RESIDUO(RES,S,PRESS,CL_CONECT,NPUNTOS,NNOD &
     &  ,IFIX_NODE,NFIX,RFIX_VALUE,IFLAG)

      INTEGER NNOD,NFIX
      INTEGER CL_CONECT(32,NNOD),IFIX_NODE(NFIX),NPUNTOS(NNOD*2)
      
      REAL*8 S(32,NNOD),RES(NNOD),PRESS(NNOD),RFIX_VALUE(NFIX)
      
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
     &        +S(J,INOD)*PRESS(INOD)
          END DO
          
        END DO
        
      END IF
      
      DO IFIX=1,NFIX
        RES(IFIX_NODE(IFIX))=0.D0
      END DO      
      
      RETURN
      END

!CCCC---------------------------------CCCC
!C    FUNCION AUXILIAR [1] DEL SOLVER    C 
!CCCC---------------------------------CCCC
      FUNCTION FRR_1(RES,RESB,NNOD)

      INTEGER NNOD
      
      REAL*8 FRR_1
      REAL*8 RES(NNOD),RESB(NNOD)

      FRR_1=0.D0
      DO INOD=1,NNOD
        FRR_1=FRR_1+RES(INOD)*RESB(INOD)
      END DO

      RETURN
      END
        
!CCCC---------------------------------CCCC
!C    FUNCION AUXILIAR [2] DEL SOLVER    C 
!CCCC---------------------------------CCCC
      SUBROUTINE SUBPKBAPK(PKBAPK,APK,PK,NNOD)

      INTEGER NNOD

      REAL*8 PKBAPK
      REAL*8 APK(NNOD),PK(NNOD)
      
      PKBAPK=0.D0
      DO INOD=1,NNOD
        PKBAPK=PKBAPK+PK(INOD)*APK(INOD)
      END DO

      RETURN
      END

!CCCC----------------------------------------CCCC
!C      PRECONDICIONA LA MATRIZ DE RIGIDEZ      C
!C         Y EL TERMINO INDEPENDIENTE           C
!CCCC----------------------------------------CCCC

      SUBROUTINE PRECF(NNOD,NPUNTOS,S,F)

      INTEGER NNOD,NPUNTOS(NNOD)
      REAL*8 FS1,FS2
      REAL*8 S(32,NNOD*2),F(NNOD*2)

      
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
      END
      
!CCCC----------------------------------------CCCC
!C        REALIZA CAMBIO DE VARIABLES           C 
!CCCC----------------------------------------CCCC

      SUBROUTINE CAMBIOVAR(NFIXP,IFIXP_NODE,IFIXP_NODEN,RFIXP_VALUE, &
     & RFIXV_VALUE,RFIXP_VALUEN,NNOD,NPUNTOS,NPUNTOSN,CL_CONECT, &
     & CL_CONECTN)

      INTEGER NFIXP,NNOD
      INTEGER NPUNTOS(NNOD),NPUNTOSN(NNOD*2)
      INTEGER IFIXP_NODE(NFIXP),IFIXP_NODEN(NFIXP*2)
      INTEGER CL_CONECT(16,NNOD),CL_CONECTN(32,NNOD*2)

      REAL*8 RFIXP_VALUE(NFIXP),RFIXP_VALUEN(NFIXP*2)
      REAL*8 RFIXV_VALUE(NFIXP)

      
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
      END
      
!CCCC---------------------------------------------CCCC
!C    EVALUA UNA FUNCION VECTORIAL Y SUS DERIVADA    C
!CCCC---------------------------------------------CCCC
      SUBROUTINE EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,VALX,VALY,RM,U,V &
     & ,ITIPE,dmax)
      
      INTEGER NNOD,INOD,ITIPE
      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD)

      REAL*8 VALX,VALY
      REAL*8 RM(10,16,NNOD),U(NNOD),V(NNOD),dmax(nnod)

!CCCC----------> ITIPE = 1  FUNCION DE FORMA
!CCCC----------> ITIPE = 2  DERIVADA RESPECTO A X   d /dx
!CCCC----------> ITIPE = 3  DERIVADA RESPECTO A Y   d /dy
!CCCC----------> ITIPE = 4  DERIVADA SEGUNDA RESPECTO A  X   d2 /dx2
!CCCC----------> ITIPE = 5  DERIVADA SEGUNDA RESPECTO A X,Y  d2 /dxdy
!CCCC----------> ITIPE = 6  DERIVADA SEGUNDA RESPECTO A  Y   d2 /dy2

      VALX=0.D0
      VALY=0.D0
      DO J=1,NPUNTOS(INOD)
        VALX=VALX+RM(ITIPE,J,INOD)*U(CL_CONECT(J,INOD))
        VALY=VALY+RM(ITIPE,J,INOD)*V(CL_CONECT(J,INOD))
      END DO

      RETURN
      END

      
!CCCC-------------------------------------------CCCC
!C    EVALUA UNA FUNCION ESCALAR Y SUS DERIVADA    C
!CCCC-------------------------------------------CCCC
      SUBROUTINE EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,VAL,RMA,U,ITIPE)
      
      INTEGER NNOD,INOD,ITIPE
      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD)

      REAL*8 VAL
      REAL*8 U(NNOD),RMA(1,16,NNOD)

!CCCC----------> ITIPE = 1  FUNCION DE FORMA
!CCCC----------> ITIPE = 2  DERIVADA RESPECTO A X   d /dx
!CCCC----------> ITIPE = 3  DERIVADA RESPECTO A Y   d /dy
!CCCC----------> ITIPE = 4  DERIVADA SEGUNDA RESPECTO A  X   d2 /dx2
!CCCC----------> ITIPE = 5  DERIVADA SEGUNDA RESPECTO A X,Y  d2 /dxdy
!CCCC----------> ITIPE = 6  DERIVADA SEGUNDA RESPECTO A  Y   d2 /dy2

      VAL=0.D0
      DO J=1,NPUNTOS(INOD)
        VAL=VAL+RMA(ITIPE,J,INOD)*U(CL_CONECT(J,INOD))
      END DO

      RETURN
      END

!CCCC-------------------------------------------CCCC
!C    RESUELVE UN SISTEMA DE ECUACIONES CON        C
!C    GAUSS SEIDEL PARA ENCONTRAR LOS PARAMETROS   C
!CCCC-------------------------------------------CCCC
      SUBROUTINE TRANSFVECT(NNOD,NPUNTOS,CL_CONECT,RM &
     &  ,UB,VB,WB,UPARAUX,VPARAUX,WPARAUX)

      INTEGER NNOD
      INTEGER KITER,INOD,INORMV,IFLAG
      INTEGER NPUNTOS(NNOD),CL_CONECT(16,NNOD)
      REAL*8 RESDEN,RESNUM,RES,UN,UT1,UT2
      REAL*8 RINODU,RINODV,RINODW,VP1,VP2,VP3
      REAL*8 UB(NNOD),VB(NNOD),WB(NNOD)
      REAL*8 UPARAUX(NNOD),VPARAUX(NNOD),WPARAUX(NNOD)
      REAL*8 RM(10,16,NNOD)
           
      RESDEN=1.D10
      DO INOD=1,NNOD
        UPARAUX(INOD)=UB(INOD)
        VPARAUX(INOD)=VB(INOD)
        WPARAUX(INOD)=WB(INOD)
        RES=UB(INOD)**2+VB(INOD)**2+WB(INOD)**2
        IF (RESDEN.GT.RES) RESDEN=RES
      END DO

      KITER=0
      RES=1.D-9
      DO WHILE(KITER.LT.30.AND.RES.GT.1.D-20)
        
        RES=0.D0
        KITER=KITER+1
        DO INOD=1,NNOD
          
          RINODU=0.D0
          RINODV=0.D0
          RINODW=0.D0
          
          DO J=2,NPUNTOS(INOD)
            RINODU=RINODU+RM(1,J,INOD)*UPARAUX(CL_CONECT(J,INOD))
            RINODV=RINODV+RM(1,J,INOD)*VPARAUX(CL_CONECT(J,INOD))
            RINODW=RINODW+RM(1,J,INOD)*WPARAUX(CL_CONECT(J,INOD))
          END DO
          
          VP1=(UB(INOD)-RINODU)/RM(1,1,INOD)
          VP2=(VB(INOD)-RINODV)/RM(1,1,INOD)
          VP3=(WB(INOD)-RINODW)/RM(1,1,INOD)
          
          RESNUM=((VP1-UPARAUX(INOD))**2.D0 &
     &      +(VP2-VPARAUX(INOD))**2.D0 &
     &	    +(VP3-WPARAUX(INOD))**2.D0)
          
          IF (RES.LT.RESNUM) RES=RESNUM
          
          UPARAUX(INOD)=VP1
          VPARAUX(INOD)=VP2
          WPARAUX(INOD)=VP3  
          
          
        END DO
        
      END DO
      
      
      RETURN
      END      

      
      SUBROUTINE SOLVE_POT(ILONG,NNOD,NNORM,NFIXF,NFIXP,NPUNTOS,CL_CONECT &
     &  ,INORM_NODE,IFIXP_NODE,IFIXF_NODE,U,V,PRESS,X,Y,dmax,RM,RMA,RNORM_VALUEX &
     &  ,RNORM_VALUEY,RFIXF_VALUEX,RFIXF_VALUEY,BAN,RFIXP_VALUE,RFIXV_VALUE &
     &  ,POISSON,YOUNG,FILE)

      PARAMETER (NODOS=20000)          ! MAXIMA DIMENSION DE LA MAT
      PARAMETER (IFIJOS=20000)
      PARAMETER (ITERAC=20000)         !NUMERO MAXIMO DE ITERACIONES

      CHARACTER FILE*80
      
      INTEGER NNOD,NNORM,NFIXF,NFIXP,ILONG
      INTEGER NPUNTOS(NNOD),NPUNTOSN(NODOS)
      INTEGER CL_CONECT(16,NNOD),CL_CONECTN(32,NODOS)
      INTEGER IFIXP_NODE(NFIXP),IFIXF_NODE(NFIXF),INORM_NODE(NNORM)
      INTEGER IFIXP_NODEN(IFIJOS),BAN(6,IFIJOS)
      
      REAL*8 BET,YOUNG,POISSON,dmax(nnod)
      REAL*8 DESPU,DESPV,DS1,DS2,DS3,DUDX,DVDX,DUDY,DVDY
      REAL*8 DU_DX(NODOS),DV_DX(NODOS),DU_DY(NODOS),DV_DY(NODOS)
      REAL*8 DU_DXAU(NODOS),DV_DXAU(NODOS),DU_DYAU(NODOS),DV_DYAU(NODOS)
      REAL*8 DUDXAU,DVDXAU,DUDYAU,DVDYAU
      REAL*8 D2UDX2,D2VDX2,D2UDXDY,D2VDXDY,D2UDY2,D2VDY2
      REAL*8 D3UDX2DY,D3VDX2DY,D3UDY2DX,D3VDY2DX
      REAL*8 D2ESFX,D2ESFY,D2TAU
      REAL*8 DESFXDX,DESFYDY,DTAUDX,DTAUDY
      REAL*8 ESFXA,ESFYA,TAUA
      REAL*8 ESFX(NODOS),ESFY(NODOS),TAU(NODOS)
      REAL*8 ESFXEXA(NODOS),ESFYEXA(NODOS),TAUEXA(NODOS)
      REAL*8 ESFXAU(NODOS),ESFYAU(NODOS),TAUAU(NODOS)
      REAL*8 ESFXPA(NODOS),ESFYPA(NODOS),TAUPA(NODOS)
      REAL*8 S(32,NODOS*2),S1(32,NODOS*2),RM(10,16,NNOD),RMA(1,16,NNOD)
      REAL*8 X(NNOD),Y(NNOD),A(NODOS),B(NODOS)
      REAL*8 RNORM_VALUEX(NNORM),RNORM_VALUEY(NNORM)
      REAL*8 RFIXF_VALUEX(NFIXF),RFIXF_VALUEY(NFIXF)
      REAL*8 RFIXP_VALUE(NFIXP),RFIXV_VALUE(NFIXP)
      REAL*8 RFIXP_VALUEN(IFIJOS)
      REAL*8 F(NODOS*2),PRESS(NNOD*2),DU(NODOS),DV(NODOS),U(NNOD),V(NNOD)
      REAL*8 DESNUMX(NODOS),DESNUMY(NODOS),DESEXAX(NODOS),DESEXAY(NODOS)
      REAL*8 DELTAX(NODOS),DELTAY(NODOS),DELTAXE(NODOS),DELTAYE(NODOS),DELTAXYE(NODOS)
      REAL*8 ERRX(NODOS,2),ERRY(NODOS,2),ERRGL(NODOS,5)
      REAL*8 ERRXAU(NODOS,2),ERRYAU(NODOS,2),ERRGLAUX(NODOS,2),ERRGLAUY(NODOS,2)
      REAL*8 RES(NODOS),RESB(NODOS),PK(NODOS),PKB(NODOS)
      REAL*8 APK(NODOS),APKB(NODOS),PREC(NODOS)
     
!CCCC----------> ENSAMBLA LA MATRIZ DE RIGIDEZ
      
      CALL ENSAMB_LAPLACE(NNOD,NFIXP,NNORM,INORM_NODE,NPUNTOS &
     & ,IFIXP_NODE,CL_CONECT,BAN,RNORM_VALUEX,RNORM_VALUEY,RM,S,S1,X,Y,dmax, &
     & YOUNG,POISSON,DS1,DS2)
      
!CCCC----------> IMPONE EL TERMINO "F" DEL SISTEMA DE ECUACIONES
    
      DO INOD=1,NNOD*2
        F(INOD)=0.D0      
      END DO

!CCCC----------> IMPONE EL VALOR DE LAS FUERZAS PRESCRITAS

      DO I=1,NFIXF
         INOD=IFIXF_NODE(I)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C               ENSAMBLANDO H CARACTERISTICO                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!cccc-----> ec.contorno-h/2*(ec.dominio)=0
         
         F(INOD*2-1)=RFIXF_VALUEX(I)*((1-(POISSON**2))/YOUNG)
         F(INOD*2)=RFIXF_VALUEY(I)*((1-(POISSON**2))/YOUNG)
         
!c         write(*,*) inod,rfixf_valuex(i),rfixf_valuey(i)       !f(inod*2-1),f(inod*2)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C             SIN ENSAMBLAR H CARACTERISTICO                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!c         F(INOD*2-1)=RFIXF_VALUEX(I)*((1-(POISSON**2))/YOUNG)
!c         F(INOD*2)=RFIXF_VALUEY(I)*((1-(POISSON**2))/YOUNG)
         
      END DO

!CCCC----------> IMPONE EL VALOR DE LOS DESPLAZAMIENTOS PRESCRITOS

      DO I=1,NFIXP
         INOD=IFIXP_NODE(I)
         IF (BAN(1,I).EQ.1) THEN
            F(INOD*2-1)=RFIXP_VALUE(I)
         END IF
         IF (BAN(2,I).EQ.1) THEN
            F(INOD*2)=RFIXV_VALUE(I)
         END IF
!c         write(*,*) f(inod*2-1),f(inod*2)
      END DO

!CCCC----------> PRECONDICIONA LA MATRIZ DE RIGIDEZ Y
!CCCC---------->      EL TERMINO INDEPENDIENTE

      CALL PRECF(NNOD,NPUNTOS,S,F)
      
!CCCC----------> REALIZA CAMBIO DE VARIABLES
      
      CALL CAMBIOVAR(NFIXP,IFIXP_NODE,IFIXP_NODEN,RFIXP_VALUE, &
     &  RFIXV_VALUE,RFIXP_VALUEN,NNOD,NPUNTOS,NPUNTOSN,CL_CONECT, &
     &  CL_CONECTN)
      
!CCCC----------> RESUELVE EL SISTEMA DE ECUACIONES SIN ENSAMBLAR
!CCCC----------> LOS NODOS PRESCRITOS
      
       CALL GRADCONJ(S,PRESS,F,NNOD*2,CL_CONECTN,NPUNTOSN, &
     &  IFIXP_NODEN,NFIXP*2,RFIXP_VALUEN,PK,PKB,APK, &
     &  APKB,RES,RESB)

!CCCC----------> RESUELVE EL SISTEMA DE ECUACIONES ENSAMBLANDO
!CCCC----------> LOS NODOS PRESCRITOS

!      CALL GRADCONJ(S,PRESS,F,NNOD*2,CL_CONECTN,NPUNTOSN &
!     &  ,IFIXP_NODEN,0,RFIXP_VALUEN,PK,PKB,APK &
!     &  ,APKB,RES,RESB)

       
      OPEN(1,FILE=FILE(1:ILONG)//'2d.res',STATUS='UNKNOWN')

      DO INOD=1,NNOD
        DU(INOD)=PRESS(INOD*2-1)
        DV(INOD)=PRESS(INOD*2)
!c          write(*,*) inod,DU(INOD),DV(INOD)
      END DO
    
!CCCC----------> IMPRESION DE LOS DESPLAZAMIENTOS
      WRITE(1,'(A15,I4,I6,3I4)') 'DESPLAZAMIENTOS',2,0,2,1,0
      
      DO INOD=1,NNOD

         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DESPU,DESPV,RM,DU,DV,1,dmax)
         WRITE(1,'(I6,2E25.13)') INOD,DESPU,DESPV
      END DO

!CCCC----------> IMPRESION DE LOS ESFUERZOS
      
      WRITE(1,'(A15,I4,I6,3I4)') 'ESFUERZOS      ',2,0,3,1,0

      DS3=YOUNG/(1-(POISSON**2))
      
      DO INOD=1,NNOD
        
        CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDX,DVDX,RM,DU,DV,2,dmax)
        CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDY,DVDY,RM,DU,DV,3,dmax)

        ESFX(INOD)=DS3*(DUDX+POISSON*DVDY)
        ESFY(INOD)=DS3*(POISSON*DUDX+DVDY)
        TAU(INOD)=DS3*(DS1*(DVDX+DUDY))
        
        WRITE(1,'(I6,3E25.13)') INOD,ESFX(INOD),ESFY(INOD),TAU(INOD)
        
      END DO
      
      CLOSE(1)

!CCCCC---> TENSIONES (exactas,numericas) en un plano a una distancia 
!CCCCC     de 12 unidades del borde superior(libre) de la placa

      OPEN(1,FILE=FILE(1:ILONG)//'2d.res-1',STATUS='UNKNOWN')


      DO INOD=1,NNOD
		if (x(inod).eq.12) then

		ESFXEXA(INOD)=-20/3.14159265/12*(cos(atan(y(inod)/12)))**4
		TAUEXA(INOD)=-20/3.14159265/12*sin(atan(y(inod)/12))*(cos(atan(y(inod)/12)))**3
        
		WRITE(1,'(5E16.6)')y(inod),esfxexa(inod),esfx(inod),tauexa(inod),tau(inod)
		end if
         
      END DO


      CLOSE(1)


!CCCCC-------->  CALCULO DE TENSIONES ALISADAS
      
!c      OPEN(1,FILE=FILE(1:ILONG)//'2d.res-2',STATUS='UNKNOWN')

!c      WRITE(1,'(A15,I4,I6,3I)') 'residuo-contorno',2,0,3,1,0

!CCCCC-------->  SUAVIZADO DE TENSIONES
!c      SUM=0.0
!c	SUM1=0.0
!c      sum2=0.0

!c      DO INOD=1,NNOD

!c        CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDX,DVDX,RM,DU,DV,2)
!c        CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDY,DVDY,RM,DU,DV,3)

!c        ESFX(INOD)=DS3*(DUDX+POISSON*DVDY)
!c        ESFY(INOD)=DS3*(POISSON*DUDX+DVDY)
!c        TAU(INOD)=DS3*DS1*(DVDX+DUDY)

!c        CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,ESFXA,RMA,ESFX,1)
!c        CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,ESFYA,RMA,ESFY,1)
!c        CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,TAUA,RMA,TAU,1)

!c        errgl(inod,1)=ESFX(INOD)
!c        errgl(inod,2)=ESFXA
!c        errgl(inod,3)=(ESFXA-ESFX(INOD))                               !(TAUA-TAU(INOD))
!c        SUM=SUM+ERRGL(INOD,1)
!c        sum2=sum2+dabs(ESFX(INOD))

!c        errgl(inod,1)=((ESFXA-ESFX(INOD))**2+(ESFYA-ESFY(INOD))**2+(2*(TAUA-TAU(INOD))**2))**0.5

!c        WRITE(1,'(I6,3E16.6)') INOD,errgl(inod,1),errgl(inod,2),errgl(inod,3)

!c        sum2=sum2+(esfx(inod)**2+esfy(inod)**2+tau(inod)**2)

!c      END DO
      
!c      sum2=sum2**0.5
!c      write(*,*) sum2

!c      DO INOD=1,NNOD
!c         ERRGL(INOD,2)=ERRGL(INOD,1)/sum2
!c         WRITE(1,'(I6,3E16.6)') INOD,errgl(inod,1),errgl(inod,2)
!c      END DO

!c      CLOSE(1)

      
!CCCCC-------> Residuo contorno

!c      DO INOD=1,NNOD
      
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDX,DVDX,RM,DU,DV,2,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDY,DVDY,RM,DU,DV,3,dmax)

!c         DU_DX(INOD)=DUDX
!c         DV_DX(INOD)=DVDX
!c         DU_DY(INOD)=DUDY
!c         DV_DY(INOD)=DVDY

!c         CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,DUDXAU,RMA,DU_DX,1)
!c         CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,DVDXAU,RMA,DV_DX,1)
!c         CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,DUDYAU,RMA,DU_DY,1)
!c         CALL EVAESC(INOD,NNOD,NPUNTOS,CL_CONECT,DVDYAU,RMA,DV_DY,1)

!c         ESFXAU(INOD)=DS3*(DUDXAU+POISSON*DVDYAU)
!c         ESFYAU(INOD)=DS3*(POISSON*DUDXAU+DVDYAU)
!c         TAUAU(INOD)=DS3*DS1*(DVDXAU+DUDYAU)
        
!c         WRITE(1,'(I6,3E16.6)') INOD,ESFXAU(INOD),ESFX(INOD)

!c      END DO

      
!c      DO I=1,NNORM
!c         INOD=INORM_NODE(I)

!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDXAU,DVDXAU,RM,DU_DX,DV_DX,1,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDYAU,DVDYAU,RM,DU_DY,DV_DY,1,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDX2,D2VDX2,RM,DU_DXAU,DV_DXAU,2,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDXDY,D2VDXDY,RM,DU_DYAU,DV_DYAU,2,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDY2,D2VDY2,RM,DU_DYAU,DV_DYAU,3,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDX2,D2VDX2,RM,DU,DV,4,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDXDY,D2VDXDY,RM,DU,DV,5,dmax)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDY2,D2VDY2,RM,DU,DV,6,dmax)

!c         CALL H_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCAR,X,Y)
!c	   AF=1.0


!cccc-----> residuo contorno=h/2*(ec.dominio)
         
!c         ERRXAU(INOD,1)= ((DUDXAU*RNORM_VALUEX(I))+(DS1*DUDYAU*RNORM_VALUEY(I))+ &
!c     &             (POISSON*DVDYAU*RNORM_VALUEX(I))+(DS1*DVDXAU*RNORM_VALUEY(I))) 
!c         ERRXAU(INOD,2)= (DABS(D2UDX2+(DS1*D2UDY2)+(POISSON*D2VDXDY)+(DS1*D2VDXDY)))
!c	    ERRXAU(INOD,2)= ESFX(INOD)*RNORM_VALUEX(I)+TAU(INOD)*RNORM_VALUEY(I)


!c         ERRYAU(INOD,1)= ((POISSON*DUDXAU*RNORM_VALUEY(I))+(DS1*DUDYAU*RNORM_VALUEX(I))+ &
!c     &             (DVDYAU*RNORM_VALUEY(I))+(DS1*DVDXAU*RNORM_VALUEX(I)))
         
!c         ERRYAU(INOD,2)= (DABS((POISSON*D2UDXDY)+(DS1*D2UDXDY)+D2VDY2+(DS1*D2VDX2)))
!C	    ERRYAU(INOD,2)= ESFY(INOD)*RNORM_VALUEY(I)+TAU(INOD)*RNORM_VALUEX(I)


!c         SUM=SUM+(ERRXAU(INOD,2))**2.0+(ERRYAU(INOD,2))**2.0
!C	   SUM1=SUM1+ERRYAU(INOD,2)
!C	    IF (X(INOD).EQ.8.0) THEN
!c         WRITE(1,'(I6,3E16.6)') INOD,ERRXAU(INOD,2),ERRYAU(INOD,2)
!C         WRITE(1,'(I6,3E16.6)') INOD,((ERRXAU(INOD,2))**2+(ERRYAU(INOD,2))**2)**0.5

!C         END IF
!c      END DO
!c
!c	WRITE(*,*) SUM**0.5
!C	WRITE(*,*) SUM1
         
!c      DO I=1,NNORM
!c         INOD=INORM_NODE(I)

!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDX,DVDX,RM,DU,DV,2)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,DUDY,DVDY,RM,DU,DV,3)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDX2,D2VDX2,RM,DU,DV,4)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDXDY,D2VDXDY,RM,DU,DV,5)
!c         CALL EVAVECT(INOD,NNOD,NPUNTOS,CL_CONECT,D2UDY2,D2VDY2,RM,DU,DV,6)


!c         ERRX(INOD,1)= (DUDX*RNORM_VALUEX(I)+DS1*DUDY*RNORM_VALUEY(I)+ &
!c     &             POISSON*DVDY*RNORM_VALUEX(I)+DS1*DVDX*RNORM_VALUEY(I))

!C         ERRX(INOD,2)=(D2UDX2+DS1*D2UDY2+POISSON*D2VDXDY+DS1*D2VDXDY)
         
!c         ERRY(INOD,1)= (POISSON*DUDX*RNORM_VALUEY(I)+DS1*DUDY*RNORM_VALUEX(I)+ &
!c     &             DVDY*RNORM_VALUEY(I)+DS1*DVDX*RNORM_VALUEX(I))

!C         ERRY(INOD,2)=(POISSON*D2UDXDY+DS1*D2UDXDY+D2VDY2+DS1*D2VDX2)

!c      END DO

!c      DO I=1,NFIXF
!c         INOD=IFIXF_NODE(I)

!c         ERRX(INOD,1)=ERRX(INOD,1)-(RFIXF_VALUEX(I)/DS3)
!c         ERRY(INOD,1)=ERRY(INOD,1)-(RFIXF_VALUEY(I)/DS3)

!c      END DO

!C      DO I=1,NFIXF
!C         INOD=IFIXF_NODE(I)

!C         ERRXAU(INOD,1)=ERRXAU(INOD,1)-(RFIXF_VALUEX(I)/DS3)
!C         ERRYAU(INOD,1)=ERRYAU(INOD,1)-(RFIXF_VALUEY(I)/DS3)

!C      END DO

      
!c      DO I=1,NNORM
!c         INOD=INORM_NODE(I)
!c         CALL H_CARACT(INOD,NNOD,NPUNTOS,CL_CONECT,HCAR,HCARX,HCARY,X,Y)
!c         AF=1.0

!c         ERRGLAUX(INOD,1)= ERRX(INOD,1)                      !-(AF*HCAR/2*ERRX(INOD,2))
!c         ERRGLAUY(INOD,1)= ERRY(INOD,1)                      !-(AF*HCAR/2*ERRY(INOD,2))

!C         ERRGLAUX(INOD,2)= ERRXAU(INOD,1)-(AF*HCAR/2*ERRXAU(INOD,2))
!C         ERRGLAUY(INOD,2)= ERRYAU(INOD,1)-(AF*HCAR/2*ERRYAU(INOD,2))


!C         ERRGLAUX= ((ERRXAU(INOD,1)-ERRX(INOD,1))/(ERRXAU(INOD,2)-ERRX(INOD,2)))
!C         IF (ERRGLAUX.GT.1.0) THEN
!C            ERRGLAUX=1.0
!C         END IF
!C         IF (ERRGLAUX.LT.0.0) THEN
!C            ERRGLAUX=0.0
!C         END IF
         
         
!C         ERRGLAUY= ((ERRYAU(INOD,1)-ERRY(INOD,1))/(ERRYAU(INOD,2)-ERRY(INOD,2)))
!C         IF (ERRGLAUY.GT.1.0) THEN
!C            ERRGLAUY=1.0
!C         END IF
!C         IF (ERRGLAUY.LT.0.0) THEN
!C            ERRGLAUY=0.0
!C         END IF
         
!C         ERRGL(INOD,1)=  ERRGLAUX(INOD,1)-ERRGLAUX(INOD,2)
!C         ERRGL(INOD,2)=  ERRGLAUY(INOD,1)-ERRGLAUY(INOD,2)

         
!c         WRITE(1,'(I6,3E16.6)') INOD,ERRGLAUX(INOD,1),ERRGLAUY(INOD,1)
         
!c      END DO

!C      DO I=1,NNORM
!C         INOD=INORM_NODE(I)
!C         SUM=SUM+ERRGL(INOD)
!C      END DO
      
!C      SUM=SUM/NNORM
!C      write(*,*) SUM
      
      CLOSE(1)

      RETURN
      END
