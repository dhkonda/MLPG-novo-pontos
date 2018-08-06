!##################################################################################################
!esta subrotina faz resolve as integrais de contorno para um ponto de gauss
!##################################################################################################
    !sub rotina nova
    
SUBROUTINE INTEGR_CONTORNO_2(JJ,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
    PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
IMPLICIT NONE

INTEGER::JJ,NPT,PSUPG,MON,TMON,J,K,i,L,SHAPE
INTEGER,DIMENSION(3)::CPRES
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::AUX1,AUX2,XT,YT,DiST,RAIOi,W,JAC,N1,N2,ALFA,ALFA2,dc,q,ALFAC
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(3)::VSUPPRES,MATRIZR2,INT43,MATRIZD5,INT53,INT83,INTPEN13,INTPEN23,INTPEN33
REAL*8,DIMENSION(3,2)::MATRIZD3,MATRIZD4
REAL*8,DIMENSION(3,3)::MATRIZN3,AUX33
REAL*8,DIMENSION(2,2)::MATRIZN4,MATRIZN5
REAL*8,DIMENSION(1,3*PSUPG)::MATRIZFi5,AUX13P1
REAL*8,DIMENSION(2,3*PSUPG)::MATRIZFi3,MATRIZFi4,MATRIZN4FI3,MATRIZN5Fi4,AUX23P,MATRIZN4Fi4,MATRIZN5Fi3,MATRIZFi6
REAL*8,DIMENSION(3,3*PSUPG)::AUX33P1,AUX33P2,INT333N,AUX33P3,INT733N,MATRIZFi1PEN,INTPEN133N,INTPEN233N,INTPEN333N,AUXDESLP
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY,DDMX,DDMY,DDMXY
REAL*8,DIMENSION(TMON,PSUPG)::AiB,dAiXB,dAiYB,AidBX,AidBY,Sb
REAL*8,DIMENSION(PSUPG,PSUPG)::Sa,DIST3
REAL*8,DIMENSION(PSUPG)::Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY,DIST2,R,dRx,dRy,ddRx,ddRy,ddRxy
REAL*8,DIMENSION(1)::XB,YB
REAL*8,DIMENSION(3*NPT,3*NPT)::KGLOBAL
REAL*8,DIMENSION(3*NPT)::FGLOBAL
REAL*8,DIMENSION(3,4)::MATRIZD7
REAL*8,DIMENSION(4,3*PSUPG)::MATRIZFi7

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
        
        
MATRIZFi7=0.D0
DO K=1,3*PSUPG,3
    MATRIZFI7(1,K)=CPRES(1)*((N1**2)*dFiX((K+2)/3)+(N1*N2*dFiY((K+2)/3)))
    MATRIZFI7(1,K+1)=CPRES(1)*((N1*N2)*dFiX((K+2)/3)+((N2**2)*dFiY((K+2)/3)))
    MATRIZFI7(2,K)=CPRES(2)*((N2**2)*dFiX((K+2)/3)-(N1*N2*dFiY((K+2)/3)))   
    MATRIZFI7(2,K+1)=CPRES(2)*(-(N1*N2)*dFiX((K+2)/3)+((N1**2)*dFiY((K+2)/3)))
    MATRIZFI7(3,K)=CPRES(1)*(-(N1*N2)*dFiX((K+2)/3)+((N1**2)*dFiY((K+2)/3)))+CPRES(2)*(-(N1*N2)*dFiX((K+2)/3)-((N2**2)*dFiY((K+2)/3)))
    MATRIZFI7(3,K+1)=CPRES(1)*(-(N2**2)*dFiX((K+2)/3)+((N1*N2)*dFiY((K+2)/3)))+CPRES(2)*((N1**2)*dFiX((K+2)/3)+((N1*N2)*dFiY((K+2)/3)))
    MATRIZFI7(4,K)=CPRES(1)*N1*FI((K+2)/3)
    MATRIZFI7(4,K+1)=CPRES(1)*N2*FI((K+2)/3)
    MATRIZFI7(4,K+2)=CPRES(3)*(N1*DFiX((K+2)/3)+N2*DFiY((K+2)/3))
ENDDO
AUXDESLP=0.D0

AUXDESLP=MATMUL(MATRIZD7,MATRIZFI7)
AUXDESLP=MATMUL(MATRIZN3,AUXDESLP)
AUXDESLP=-W*JAC*AUXDESLP          !RESULTADO DA TERCEIRA E DA SETIMA INTEGRAL  

    !computar a terceira e semetima integral na matriz KGLOBAL    !AINDA FALTA COMPUTAR UMA INTEGRAL... O TERMO Q VAI NO VETOR LIVRE Q VEM DE Mn
DO K=1,PSUPG
    DO J=1,3
	    KGLOBAL(JJ*3-2,FG(K)*3-3+J)=AUXDESLP(1,K*3-3+J)+KGLOBAL(JJ*3-2,FG(K)*3-3+J)	
	    KGLOBAL(JJ*3-1,FG(K)*3-3+J)=AUXDESLP(2,K*3-3+J)+KGLOBAL(JJ*3-1,FG(K)*3-3+J)
	    KGLOBAL(JJ*3,FG(K)*3-3+J)=AUXDESLP(3,K*3-3+J)+KGLOBAL(JJ*3,FG(K)*3-3+J) !É  A LINHA RESPONSAVEL PELA 7 INTEGRAL
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

    
IF(CPRES(1).EQ.0)THEN
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

ENDIF

IF(CPRES(2).EQ.0)THEN
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
    
ENDIF
    
IF(CPRES(3).EQ.0)THEN
    !resolve aqui as integrais com força transversal prescrita
    !solução da oitava integral 
    INT83=0.D0
    INT83(3)=W*VSUPPRES(3)*JAC
    
    !computar a oitava integral na matriz FGLOBAL
    FGLOBAL(JJ*3)=INT83(3)+FGLOBAL(JJ*3)  
    
ENDIF

ENDSUBROUTINE INTEGR_CONTORNO_2
