!##################################################################################################
!esta subrotina faz resolve as integrais de contorno para um ponto de gauss
!##################################################################################################
    
SUBROUTINE INTEGR_CONTORNO(JJ,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
    PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
IMPLICIT NONE

INTEGER::JJ,NPT,PSUPG,MON,TMON,J,K,i,L,SHAPE
INTEGER,DIMENSION(3)::CPRES
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::AUX1,AUX2,XT,YT,DiST,RAIOi,W,JAC,N1,N2,ALFA,ALFA2,q,ALFAC,dc
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(3)::VSUPPRES,MATRIZR2,INT43,MATRIZD5,INT53,INT83,INTPEN13,INTPEN23,INTPEN33
REAL*8,DIMENSION(3,2)::MATRIZD3,MATRIZD4
REAL*8,DIMENSION(3,3)::MATRIZN3,AUX33
REAL*8,DIMENSION(2,2)::MATRIZN4,MATRIZN5
REAL*8,DIMENSION(1,3*PSUPG)::MATRIZFi5,AUX13P1
REAL*8,DIMENSION(2,3*PSUPG)::MATRIZFi3,MATRIZFi4,MATRIZN4FI3,MATRIZN5Fi4,AUX23P,MATRIZN4Fi4,MATRIZN5Fi3,MATRIZFi6
REAL*8,DIMENSION(3,3*PSUPG)::AUX33P1,AUX33P2,INT333N,AUX33P3,INT733N,MATRIZFi1PEN,INTPEN133N,INTPEN233N,INTPEN333N
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY,DDMX,DDMY,DDMXY
REAL*8,DIMENSION(TMON,PSUPG)::AiB,dAiXB,dAiYB,AidBX,AidBY,Sb
REAL*8,DIMENSION(PSUPG,PSUPG)::Sa,DIST3
REAL*8,DIMENSION(PSUPG)::Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY,DIST2,R,dRx,dRy,ddRx,ddRy,ddRxy
REAL*8,DIMENSION(1)::XB,YB
REAL*8,DIMENSION(3*NPT,3*NPT)::KGLOBAL
REAL*8,DIMENSION(3*NPT)::FGLOBAL
REAL*8,DIMENSION(3,4)::MATRIZD7

q=1.03
ALFAC=1.42

!calcular a distancia do ponto de gauss ao ponto base
AUX1=(X(JJ)-XT)**2.D0
AUX2=(Y(JJ)-YT)**2.D0
Dist=DSQRT(AUX1+AUX2)
!avaliar a função peso(spline de 4 ordem) no ponto de gauss (função peso centrada no ponto base)
W=1.D0-6.D0*((Dist/RAIOi)**2.D0)+8.D0*((Dist/RAIOi)**3.D0)-3.D0*((Dist/RAIOi)**4.D0)

XB(1)=XT
YB(1)=YT
CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XB,YB,M,DMX,DMY,DDMX,DDMY,DDMXY)

SELECT CASE(SHAPE)
    CASE(1)
        Fi=0.D0
        dFiX=0.D0
        dFiY=0.D0
        DO K=1,PSUPG
            DO L=1,TMON
                Fi(K)=M(1,L)*AiB(L,K)+Fi(K)
                DFiX(K)=DMX(1,L)*AiB(L,K)+DFiX(K)
                DFiY(K)=DMY(1,L)*AiB(L,K)+DFiY(K)
                DFiX(K)=M(1,L)*dAiXB(L,K)+DFiX(K)
                DFiY(K)=M(1,L)*dAiyB(L,K)+DFiY(K)
                DFiX(K)=M(1,L)*AidBx(L,K)+DFiX(K)
                DFiY(K)=M(1,L)*AidBy(L,K)+DFiY(K)
            ENDDO
        ENDDO
    CASE(2)

        Fi=0.D0
        dFiX=0.D0
        dFiY=0.D0
        ddFix=0.D0
        ddFiY=0.D0
        ddFixy=0.D0

        DO K=1,PSUPG
            DO L=1,TMON
                Fi(K)=M(1,L)*Sb(L,K)+Fi(K)
            ENDDO
        ENDDO        
 
    DO K=1,PSUPG
        DIST2(K)=((X(FG(K))-XT)**2+(Y(FG(K))-YT)**2)
    ENDDO

    DO K=1,PSUPG
        DO L=1,PSUPG
            DIST3(K,L)=((X(FG(L))-X(FG(K)))**2+(Y(FG(L))-Y(FG(K)))**2)
        ENDDO
    ENDDO
    
    dc=DSQRT(dist3(1,3))
    
    DO K=1,PSUPG
        R(K)=((DIST2(K))+(ALFAC*dc)**2)**q
    ENDDO

    DO K=1,PSUPG
        DO L=1,PSUPG
            Fi(K)=R(L)*Sa(L,K)+Fi(K)
        ENDDO
    ENDDO    


    DO K=1,PSUPG
        !R(1,J)=((DIST(1,J))+(ALFAC)**2)**q
        dRx(K)=2.d0*(-X(FG(K))+XT)*q*((DIST2(K))+(ALFAC*dc)**2)**(q-1.)
        dRy(K)=2.d0*(-Y(FG(K))+YT)*q*((DIST2(K))+(ALFAC*dc)**2)**(q-1.)
        ddRx(K)=(2.d0*q*(DIST2(K)+(ALFAC*dc)**2)**(q-1.))+(4.D0*q*(q-1)*(DIST2(K)+(ALFAC*dc)**2)**(q-2.))*(-X(FG(K))+XT)**2
        ddRy(K)=(2.d0*q*(DIST2(K)+(ALFAC*dc)**2)**(q-1.))+(4.D0*q*(q-1)*(DIST2(K)+(ALFAC*dc)**2)**(q-2.))*(-Y(FG(K))+YT)**2
        ddRxy(K)=4.D0*q*(q-1.D0)*(DIST2(K)+(ALFAC*dc)**2)**(q-2.)*(-X(FG(K))+XT)*(-Y(FG(K))+YT)
    ENDDO        
            
    DO K=1,PSUPG
        DO L=1,PSUPG
            dFix(K)=dRx(L)*Sa(L,K)+dFix(K)
            dFiy(K)=dRy(L)*Sa(L,K)+dFiy(K)    
    
            ddFix(K)=ddRx(L)*Sa(L,K)+ddFix(K)
            ddFiy(K)=ddRY(L)*Sa(L,K)+ddFiy(K)    
            ddFixy(K)=ddRXY(L)*Sa(L,K)+ddFiXY(K)    
        ENDDO
    ENDDO        
            
    DO K=1,PSUPG
        DO L=1,TMON
            dFix(K)=dMx(1,L)*Sb(L,K)+dFix(K)
            dFiy(K)=dMy(1,L)*Sb(L,K)+dFiy(K)
    
            ddFix(K)=ddMX(1,L)*Sb(L,K)+ddFix(K)
            ddFiy(K)=ddMY(1,L)*Sb(L,K)+ddFiy(K)
            ddFixy(K)=ddMXY(1,L)*Sb(L,K)+ddFixy(K)
        ENDDO
    ENDDO        
        
ENDSELECT
IF(CPRES(1).EQ.1)THEN
    !resolve aqui as integrais com rotacao normal prescrita

!montar a matriz Fi3, matriz Fi4
    MATRIZFi3=0.D0
    MATRIZFi4=0.D0
    MATRIZFi1PEN=0.D0
    
    DO K=1,3*PSUPG,3                               
	    MATRIZFI3(1,K)=N1*DFiX((k+2)/3)
		MATRIZFI3(1,K+1)=N2*DFiX((k+2)/3)
		MATRIZFI3(2,K)=N1*DFiY((k+2)/3)
		MATRIZFI3(2,K+1)=N2*DFiY((k+2)/3)

        MATRIZFI4(1,K)=-N2*DFiX((K+2)/3)           
		MATRIZFI4(1,K+1)=N1*DFiX((K+2)/3)
		MATRIZFI4(2,K)=-N2*DFiY((K+2)/3)
		MATRIZFI4(2,K+1)=N1*DFiY((K+2)/3)
        
        MATRIZFi1PEN(1,K)=Fi((K+2)/3)   
        MATRIZFi1PEN(2,K+1)=Fi((K+2)/3)   
        MATRIZFi1PEN(3,K+2)=Fi((K+2)/3)   
    ENDDO    

    !iniciando solução da terceira integral da equação 20 para rotacao normal prescrita
    MATRIZN4Fi3=0.D0
    MATRIZN5Fi4=0.D0

    !multiplicando a matriz N4 pela matriz Fi3 e matriz N5 pela matriz Fi4
    MATRIZN4Fi3=MATMUL(MATRIZN4,MATRIZFi3)
    MATRIZN5Fi4=MATMUL(MATRIZN5,MATRIZFi4)
    !somando MatrizN4*matrizFi3+MatrizN5*matrizFi4
    AUX23P=MATRIZN4Fi3+MATRIZN5Fi4

    AUX33P1=0.D0
    !calculando MatrizD3*(MatrizN4*matrizFi3+MatrizN5*matrizFi4)
    AUX33P1=MATMUL(MATRIZD3,AUX23P)
    AUX33P2=0.D0
    !calculando MatrizN3*MatrizD3*(MatrizN4*matrizFi3+MatrizN5*matrizFi4)
    AUX33P2=MATMUL(MATRIZN3,AUX33P1)

    INT333N=0.D0
    !calculando o resultado da integral para o ponto de gauss
    INT333N=-W*JAC*AUX33P2

    !computar a terceira integral na matriz KGLOBAL
    DO K=1,PSUPG
        DO J=1,3
  		    KGLOBAL(JJ*3-2,FG(K)*3-3+J)=INT333N(1,K*3-3+J)+KGLOBAL(JJ*3-2,FG(K)*3-3+J)	
		    KGLOBAL(JJ*3-1,FG(K)*3-3+J)=INT333N(2,K*3-3+J)+KGLOBAL(JJ*3-1,FG(K)*3-3+J)
		    KGLOBAL(JJ*3,FG(K)*3-3+J)=INT333N(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J)
	    ENDDO														  		 	   
    ENDDO  
   
    !iniciando a solução da quarta integral da equação 20
    INT43=0.D0
    !produto de matrizN3*matrizR2
    INT43=MATMUL(MATRIZN3,MATRIZR2)
    !produto de W*matrizN3*matrizR2*JAC
    INT43=W*JAC*INT43

    !computando a quarta integral na matriz FGLOBAL
    FGLOBAL(JJ*3-2)=INT43(1)+FGLOBAL(JJ*3-2)
    FGLOBAL(JJ*3-1)=INT43(2)+FGLOBAL(JJ*3-1)
    FGLOBAL(JJ*3)=INT43(3)+FGLOBAL(JJ*3)
   
    !iniciando a solução da integral da penalidade (matriz K)
    IF (ABS(ALFA).GE.ALFA2) THEN
        
        AUX33=0.D0
        DO i=1,3
            DO J=1,2
                DO K=1,2 !estou reaproveitando matrizes... por isso as dimensoes parecem erradas
                    AUX33(i,J)=MATRIZN3(i,K)*MATRIZN4(K,J)+AUX33(I,J)
                ENDDO
            ENDDO  
        ENDDO
        
        INTPEN133N=0.D0
        DO i=1,3
            DO J=1,3*PSUPG
                DO K=1,3
                    INTPEN133N(i,J)=AUX33(i,K)*MATRIZFi1PEN(K,J)+INTPEN133N(i,J)
                ENDDO
                INTPEN133N(i,J)=INTPEN133N(i,J)*W*ALFA*JAC
            ENDDO
        ENDDO
        
        !computar a integral da penalidade referente a rotacao normal na matriz KGLOBAL
        DO K=1,PSUPG
            DO J=1,3
      		    KGLOBAL(JJ*3-2,FG(K)*3-3+J)=INTPEN133N(1,K*3-3+J)+KGLOBAL(JJ*3-2,FG(K)*3-3+J)	
    		    KGLOBAL(JJ*3-1,FG(K)*3-3+J)=INTPEN133N(2,K*3-3+J)+KGLOBAL(JJ*3-1,FG(K)*3-3+J)
    		    KGLOBAL(JJ*3,FG(K)*3-3+J)=INTPEN133N(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J)
    	    ENDDO														  		 	   
        ENDDO  
    
        !iniciando a solução da integral da penalidade (vetor F)
        INTPEN13(1)=ALFA*W*N1*VSUPPRES(1)*JAC
        INTPEN13(2)=ALFA*W*N2*VSUPPRES(1)*JAC
        INTPEN13(3)=0.D0
        !computando a solucao da integral da penalidade na matriz FGLOBAL
        FGLOBAL(JJ*3-2)=INTPEN13(1)+FGLOBAL(JJ*3-2)
        FGLOBAL(JJ*3-1)=INTPEN13(2)+FGLOBAL(JJ*3-1)
        !(O terceiro termo sempre é nulo em intpen13
    ENDIF
    
ELSEIF(CPRES(1).EQ.0)THEN
    !solução da quinta integral para momento normal prescrito
    AUX33=0.D0
    DO i=1,3
        DO j=1,3
            AUX33(I,j)=W*MATRIZN3(i,J)
        ENDDO
    ENDDO
    INT53=0.D0
    DO i=1,3
        INT53(i)=AUX33(i,1)*VSUPPRES(1)*JAC
    ENDDO
    !computando a quinta integral na matriz FGLOBAL
    FGLOBAL(JJ*3-2)=INT53(1)+FGLOBAL(JJ*3-2)
    FGLOBAL(JJ*3-1)=INT53(2)+FGLOBAL(JJ*3-1)
    !int53(3) sempre é nulo
    FGLOBAL(JJ*3)=INT53(3)+FGLOBAL(JJ*3)    
    
ELSE
    write(*,*) 'algum erro'
ENDIF

IF(CPRES(2).EQ.1)THEN
    !resolve aqui as integrais com rotacao tangencial prescrita
!montar a matriz Fi3, matriz Fi4
    MATRIZFi3=0.D0
    MATRIZFi4=0.D0
    MATRIZFi1PEN=0.D0
    DO K=1,3*PSUPG,3
	    MATRIZFI3(1,K)=N1*DFiX((k+2)/3)
		MATRIZFI3(1,K+1)=N2*DFiX((k+2)/3)
		MATRIZFI3(2,K)=N1*DFiY((k+2)/3)
		MATRIZFI3(2,K+1)=N2*DFiY((k+2)/3)

        MATRIZFI4(1,K)=-N2*DFiX((K+2)/3)
		MATRIZFI4(1,K+1)=N1*DFiX((K+2)/3)
		MATRIZFI4(2,K)=-N2*DFiY((K+2)/3)
		MATRIZFI4(2,K+1)=N1*DFiY((K+2)/3)
        
        MATRIZFi1PEN(1,K)=Fi((K+2)/3)   
        MATRIZFi1PEN(2,K+1)=Fi((K+2)/3)   
        MATRIZFi1PEN(3,K+2)=Fi((K+2)/3)           
    ENDDO    
    
    !iniciando solução da terceira integral da equação 20 para rotacao tangencial prescrita
    MATRIZN4Fi4=0.D0
    MATRIZN5Fi3=0.D0
    DO i=1,2
        DO J=1,3*PSUPG
            DO K=1,2
                !multiplicando a matriz N4 pela matriz Fi4 e matriz N5 pela matriz Fi3
                MATRIZN4Fi4(i,J)=MATRIZN4(i,K)*MATRIZFi4(K,J)+MATRIZN4Fi4(i,J)
                MATRIZN5Fi3(i,J)=MATRIZN5(i,K)*MATRIZFi3(K,J)+MATRIZN5Fi3(i,J)
            ENDDO
            !somando MatrizN4*matrizFi3+MatrizN5*matrizFi4
            AUX23P(i,J)=MATRIZN4Fi4(i,J)+MATRIZN5Fi3(i,J)
        ENDDO
    ENDDO
    AUX33P1=0.D0
    DO i=1,3
        DO J=1,3*PSUPG
            DO K=1,2
                !calculando MatrizD4*(MatrizN4*matrizFi3+MatrizN5*matrizFi4)
                AUX33P1(i,J)=MATRIZD4(i,K)*AUX23P(K,J)+AUX33P1(i,J)
            ENDDO
        ENDDO
    ENDDO
    AUX33P2=0.D0
    DO i=1,3
        DO J=1,3*PSUPG
            DO K=1,3
                !calculando MatrizN3*MatrizD4*(MatrizN4*matrizFi3+MatrizN5*matrizFi4)
                AUX33P2(i,J)=MATRIZN3(i,K)*AUX33P1(K,J)+AUX33P2(i,J)
            ENDDO
        ENDDO
    ENDDO
    INT333N=0.D0
    DO i=1,3
        DO J=1,3*PSUPG
            !calculando o resultado da integral para o ponto de gauss
            INT333N(i,J)=-W*AUX33P2(i,J)*JAC
        ENDDO
    ENDDO
    !computar a terceira integral na matriz KGLOBAL
    DO K=1,PSUPG
        DO J=1,3
  		    KGLOBAL(JJ*3-2,FG(K)*3-3+J)=INT333N(1,K*3-3+J)+KGLOBAL(JJ*3-2,FG(K)*3-3+J)	
		    KGLOBAL(JJ*3-1,FG(K)*3-3+J)=INT333N(2,K*3-3+J)+KGLOBAL(JJ*3-1,FG(K)*3-3+J)
		    KGLOBAL(JJ*3,FG(K)*3-3+J)=INT333N(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J)
	    ENDDO														  		 	   
    ENDDO
    
    IF (ABS(ALFA).GE.ALFA2) THEN
        !iniciando a solução da integral da penalidade (matriz K) para rotação tangencial prescrita
        AUX33=0.D0
        DO i=1,3
            DO J=1,2
                DO K=1,2 !estou reaproveitando matrizes... por isso as dimensoes parecem erradas
                    AUX33(i,J)=MATRIZN3(i,K)*MATRIZN5(K,J)+AUX33(I,J)
                ENDDO
            ENDDO  
        ENDDO
        
        INTPEN233N=0.D0
        DO i=1,3
            DO J=1,3*PSUPG
                DO K=1,3
                    INTPEN233N(i,J)=AUX33(i,K)*MATRIZFi1PEN(K,J)+INTPEN233N(i,J)
                ENDDO
                INTPEN233N(i,J)=INTPEN233N(i,J)*W*ALFA*JAC
            ENDDO
        ENDDO
        
        !computar a integral da penalidade referente a rotacao tangencial na matriz KGLOBAL
        DO K=1,PSUPG
            DO J=1,3
  	    	    KGLOBAL(JJ*3-2,FG(K)*3-3+J)=INTPEN233N(1,K*3-3+J)+KGLOBAL(JJ*3-2,FG(K)*3-3+J)	
	    	    KGLOBAL(JJ*3-1,FG(K)*3-3+J)=INTPEN233N(2,K*3-3+J)+KGLOBAL(JJ*3-1,FG(K)*3-3+J)
	    	    KGLOBAL(JJ*3,FG(K)*3-3+J)=INTPEN233N(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J)
	        ENDDO														  		 	   
        ENDDO  

        !iniciando a solução da integral da penalidade (vetor F)
        INTPEN23(1)=-ALFA*W*N2*VSUPPRES(2)*JAC
        INTPEN23(2)=ALFA*W*N1*VSUPPRES(2)*JAC
        INTPEN23(3)=0.D0
        !computando a solucao da integral da penalidade na matriz FGLOBAL
        FGLOBAL(JJ*3-2)=INTPEN23(1)+FGLOBAL(JJ*3-2)
        FGLOBAL(JJ*3-1)=INTPEN23(2)+FGLOBAL(JJ*3-1)
        !(O terceiro termo sempre é nulo em intpen23
    ENDIF
    
ELSEIF(CPRES(2).EQ.0)THEN
    !solução da quinta integral para momento tangencial prescrito
    AUX33=0.D0
    DO i=1,3
        DO j=1,3
            AUX33(I,j)=W*MATRIZN3(i,J)
        ENDDO
    ENDDO
    INT53=0.D0
    DO i=1,3
        INT53(i)=AUX33(i,2)*VSUPPRES(2)*JAC
    ENDDO
    !computando a quinta integral na matriz FGLOBAL
    FGLOBAL(JJ*3-2)=INT53(1)+FGLOBAL(JJ*3-2)
    FGLOBAL(JJ*3-1)=INT53(2)+FGLOBAL(JJ*3-1)
    !int53(3) sempre é nulo
    FGLOBAL(JJ*3)=INT53(3)+FGLOBAL(JJ*3)    
    
ELSE
    write(*,*) 'algum erro'
    !resolve aqui as integrais com momento tangencial prescrito
ENDIF

IF(CPRES(3).EQ.1)THEN
    !resolve aqui as integrais com deslocamento tranversal prescrito
    !montar matriz Fi5 e Fi6
    MATRIZFI5=0.D0
    MATRIZFI6=0.D0
    MATRIZFi1PEN=0.D0
    DO K=1,3*PSUPG,3
		MATRIZFI5(1,K)=N1*Fi((K+2)/3)
		MATRIZFI5(1,K+1)=N2*Fi((K+2)/3)
        
		MATRIZFI6(1,K+2)=DFiX((K+2)/3)
		MATRIZFI6(2,K+2)=DFiY((K+2)/3)
        
        MATRIZFi1PEN(1,K)=Fi((K+2)/3)   
        MATRIZFi1PEN(2,K+1)=Fi((K+2)/3)   
        MATRIZFi1PEN(3,K+2)=Fi((K+2)/3)   
    ENDDO
    AUX13P1=0.D0
    !produto da primeira linha da matriz N4 pela matriz Fi6
    DO i=1,3*PSUPG
        DO J=1,2
            AUX13P1(1,i)=MATRIZN4(1,J)*MATRIZFi6(J,i)+AUX13P1(1,i)
        ENDDO
    ENDDO
    !somar o produto da primeira linha da matriz N4 pela matriz Fi6 com a matriz Fi5
    DO i=1,3*PSUPG
        AUX13P1(1,i)=AUX13P1(1,i)+MATRIZFi5(1,i)
    ENDDO
    AUX33P3=0.D0
    DO i=1,3*PSUPG
        AUX33P3(3,i)=MATRIZD5(3)*AUX13P1(1,i)
    ENDDO
    INT733N=0.D0
    DO i=1,3*PSUPG
        INT733N(3,i)=-W*AUX33P3(3,i)*JAC
    ENDDO
    !computar o resultado da sétima integral na matriz KGLOBAL
    DO K=1,PSUPG
        DO J=1,3
            !o resultado das duas primeiras linhas da sétima integral é nulo
  		    KGLOBAL(JJ*3,FG(K)*3-3+J)=INT733N(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J)
	    ENDDO														  		 	   
    ENDDO     
    
    !iniciando a solução da integral da penalidade (matriz K) para deslocamento transversal prescrito
    IF (ABS(ALFA).GE.ALFA2) THEN
        INTPEN333N=0.D0
            DO J=1,PSUPG
               INTPEN333N(3,J*3)=MATRIZFi1PEN(3,J*3)*W*ALFA*JAC
            ENDDO
        
        !computar a integral da penalidade referente ao deslocamento transversal na matriz KGLOBAL
        DO K=1,PSUPG
            DO J=1,3
!  	    	    KGLOBAL(JJ*3-2,FG(K)*3-3+J)=INTPEN233N(1,K*3-3+J)+KGLOBAL(JJ*3-2,FG(K)*3-3+J)	
!	    	    KGLOBAL(JJ*3-1,FG(K)*3-3+J)=INTPEN233N(2,K*3-3+J)+KGLOBAL(JJ*3-1,FG(K)*3-3+J)
	    	    KGLOBAL(JJ*3,FG(K)*3-3+J)=INTPEN333N(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J)
	        ENDDO														  		 	   
        ENDDO  
 
        !iniciando a solução da integral da penalidade (vetor F)
        INTPEN33(1)=0.D0
        INTPEN33(2)=0.D0
        INTPEN33(3)=ALFA*W*VSUPPRES(3)*JAC
        !computando a solucao da integral da penalidade na matriz FGLOBAL
        FGLOBAL(JJ*3)=INTPEN33(3)+FGLOBAL(JJ*3)
        !(o primeiro e o segundo termo sempre são nulos em intpen33    
    ENDIF
    
ELSEIF(CPRES(3).EQ.0)THEN
    !resolve aqui as integrais com força transversal prescrita
    !solução da oitava integral 
    INT83=0.D0
    INT83(3)=W*VSUPPRES(3)*JAC
    
    !computar a oitava integral na matriz FGLOBAL
    FGLOBAL(JJ*3)=INT83(3)+FGLOBAL(JJ*3)  
ELSE
    write(*,*) 'algum erro'
    
ENDIF


RETURN
ENDSUBROUTINE INTEGR_CONTORNO
