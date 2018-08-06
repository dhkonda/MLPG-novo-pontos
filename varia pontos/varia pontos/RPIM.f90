!O PONTO BASE PRECISA FAZER PARTE DO SUPORTE NESTE PROCESSO DE APROXIMAÇÃO
!pois a matriz R deve ser simétrica (deve ter algum jeito de fazer excluindo o pontos base, mas depois eu vejo isso)    
    !ESTA SUBROTINA NAO PODE MAIS SER USADA NESTE PROGRAMA (USAR A SUBROTINA QUE EVITA A INVERSAO DE MATRIZES)
!    CALL RPIM(i,XB,YB,NPT,NPSUP,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
    
SUBROUTINE RPIM2(JJ,XB,YB,NPT,NPSUP,X,Y,FG,MON,TMON,Fi,dFiX,dFiY,ddFiX,ddFiY,ddFiXY)
    IMPLICIT NONE                           
    
INTEGER::i,J,JJ,NPT,NPSUP,MON,TMON
INTEGER,DIMENSION(NPSUP)::FG
REAL*8::q,ALFAC,dc,UM,ZERO,ZEROUM,ZERO2,ZERO3,ZERO4,XB,YB
REAL*8,DIMENSION(1)::XBM,YBM
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(NPSUP,NPSUP)::R,DIST
REAL*8,DIMENSION(1,NPSUP)::dRX,dRY,ddRx,ddRy,ddRxy
REAL*8,DIMENSION(NPSUP,TMON)::P
REAL*8,DIMENSION(1,TMON)::dPx,dPy
REAL*8,DIMENSION(NPSUP)::Fi,dFiX,dFiY,ddFiX,ddFiY,ddFiXY
REAL*8,DIMENSION(NPSUP+TMON)::FiRPIM,dFiRPIMX,dFiRPIMY,ddFiRPIMX,ddFiRPIMY,ddFiRPIMxy
REAL*8,DIMENSION(TMON+NPSUP,TMON+NPSUP)::G,Gi,IDENT !,dGx,dGy,dGix,dGiy
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY,DDMX,DDMY,DDMXY    
q=1.03
ALFAC=1.42

!calculo da Matriz R
DO i=1,NPSUP
    DO J=1,NPSUP
        DIST(i,J)=((X(FG(J))-X(FG(i)))**2+(Y(FG(J))-Y(FG(i)))**2)
    ENDDO
ENDDO

dc=DSQRT(dist(1,2))

DO i=1,NPSUP
    DO J=1,NPSUP
        R(i,J)=((DIST(i,J))+(ALFAC*dc)**2)**q
    ENDDO
ENDDO

!!calculo da derivada da matriz R
!!derivada em X
!DO i=1,NPSUP
!    DO J=1,NPSUP
!        dRx(i,J)=2.d0*(X(FG(J))-X(FG(i)))*q*((DIST(i,J))+(ALFAC*DIST(1,2))**2)**(q-1.)
!    ENDDO
!ENDDO
!!derivada em y
!DO i=1,NPSUP
!    DO J=1,NPSUP
!        dRy(i,J)=2.d0*(Y(FG(J))-Y(FG(i)))*q*((DIST(i,J))+(ALFAC*DIST(1,2))**2)**(q-1.)
!    ENDDO
!ENDDO

!calculo da matriz P e suas derivadas
DO i=1,NPSUP
    CALL MONTA_BASE_MONOMIAL(MON,TMON,1,X(FG(i)),Y(FG(i)),M,DMX,DMY,DDMX,DDMY,DDMXY)     
    DO J=1,TMON
        P(i,J)=M(1,J)
        !dPx(i,J)=dMX(1,J)
        !dPy(i,J)=dMY(1,J)
    ENDDO    
ENDDO

!montagem da matriz G
G=0.D0
DO I=1,NPSUP
    DO J=1,NPSUP
        G(i,J)=R(i,J)
    ENDDO
ENDDO
DO i=1,NPSUP
    DO J=NPSUP+1,NPSUP+TMON
        G(i,J)=P(I,J-NPSUP)
    ENDDO
ENDDO
DO i=NPSUP+1,NPSUP+TMON
    DO J=1,NPSUP
        G(i,J)=P(J,i-NPSUP)
    ENDDO
ENDDO

!!montagem da matriz G derivada
!dGx=0.d0
!dGy=0.d0
!DO I=1,NPSUP
!    DO J=1,NPSUP
!        dGx(i,J)=dRx(i,J)
!        dGy(i,J)=dRy(i,J)
!    ENDDO
!ENDDO
!DO i=1,NPSUP
!    DO J=NPSUP+1,NPSUP+TMON
!        dGx(i,J)=dPx(I,J-NPSUP)
!        dGy(i,J)=dPy(I,J-NPSUP)
!    ENDDO
!ENDDO
!DO i=NPSUP+1,NPSUP+TMON
!    DO J=1,NPSUP
!        dGx(i,J)=dPx(J,i-NPSUP)
!        dGy(i,J)=dPy(J,i-NPSUP)
!    ENDDO
!ENDDO

!calculo de G inversa
DO i=1,TMON+NPSUP
    IDENT=0.d0
    IDENT(i,i)=1.d0
    CAll r8mat_fs(TMON+NPSUP,G,IDENT(:,i),Gi(:,i),JJ)
ENDDO

!!calculo da derivada da inversa de G
!dGiX=0.d0
!dGiY=0.d0
!dGiX=-MATMUL(Gi,MATMUL(dGX,Gi))
!dGiY=-MATMUL(Gi,MATMUL(dGY,Gi))

!Calculo de fi
FiRPIM=0.D0
DO J=1,NPSUP
        DIST(1,J)=((X(FG(J))-XB)**2+(Y(FG(J))-YB)**2)
ENDDO
DO J=1,NPSUP
    R(1,J)=((DIST(1,J))+(ALFAC*dc)**2)**q
ENDDO
DO i=1,NPSUP+TMON
    DO J=1,NPSUP
        FiRPIM(i)=R(1,J)*Gi(J,i)+FiRPIM(i)    
    ENDDO
ENDDO

XBM(1)=XB
YBM(1)=YB

CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XBM,YBM,M,DMX,DMY,DDMX,DDMY,DDMXY)
DO i=1,NPSUP+TMON
    DO J=NPSUP+1,NPSUP+TMON
        FiRPIM(i)=M(1,J-NPSUP)*Gi(J,i)+FiRPIM(i)
    ENDDO
ENDDO

!calculo da derivada de Fi
!composto por duas parcelas
dFiRPIMX=0.D0
dFiRPIMY=0.D0
ddFiRPIMX=0.D0
ddFiRPIMY=0.D0
ddFiRPIMXY=0.D0

DO J=1,NPSUP
        DIST(1,J)=((X(FG(J))-XB)**2+(Y(FG(J))-YB)**2)
ENDDO

DO J=1,NPSUP
    !R(1,J)=((DIST(1,J))+(ALFAC)**2)**q
    dRx(1,J)=2.d0*(-X(FG(J))+XB)*q*((DIST(1,J))+(ALFAC*dc)**2)**(q-1.)
    dRy(1,J)=2.d0*(-Y(FG(J))+YB)*q*((DIST(1,J))+(ALFAC*dc)**2)**(q-1.)
    ddRx(1,J)=(2.d0*q*(DIST(1,J)+(ALFAC*dc)**2)**(q-1.))+(4.D0*q*(q-1)*(DIST(1,J)+(ALFAC*dc)**2)**(q-2.))*(-X(FG(J))+XB)**2
    ddRy(1,J)=(2.d0*q*(DIST(1,J)+(ALFAC*dc)**2)**(q-1.))+(4.D0*q*(q-1)*(DIST(1,J)+(ALFAC*dc)**2)**(q-2.))*(-Y(FG(J))+YB)**2
    ddRxy(1,J)=4.D0*q*(q-1.D0)*(DIST(1,J)+(ALFAC*dc)**2)**(q-2.)*(-X(FG(J))+XB)*(-Y(FG(J))+YB)
ENDDO

DO i=1,NPSUP+TMON
    DO J=1,NPSUP
        dFiRPIMx(i)=dRx(1,J)*Gi(J,i)+dFiRPIMx(i)
        dFiRPIMy(i)=dRy(1,J)*Gi(J,i)+dFiRPIMy(i)    

        ddFiRPIMx(i)=ddRx(1,J)*Gi(J,i)+ddFiRPIMx(i)
        ddFiRPIMy(i)=ddRY(1,J)*Gi(J,i)+ddFiRPIMy(i)    
        ddFiRPIMxy(i)=ddRXY(1,J)*Gi(J,i)+ddFiRPIMXY(i)    
    ENDDO
ENDDO

XBM(1)=XB
YBM(1)=YB

CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XBM,YBM,M,DMX,DMY,DDMX,DDMY,DDMXY)
DO i=1,NPSUP+TMON
    DO J=NPSUP+1,NPSUP+TMON
        dFiRPIMx(i)=dMx(1,J-NPSUP)*Gi(J,i)+dFiRPIMx(i)
        dFiRPIMy(i)=dMy(1,J-NPSUP)*Gi(J,i)+dFiRPIMy(i)

        ddFiRPIMx(i)=ddMX(1,J-NPSUP)*Gi(J,i)+ddFiRPIMx(i)
        ddFiRPIMy(i)=ddMY(1,J-NPSUP)*Gi(J,i)+ddFiRPIMy(i)
        ddFiRPIMxy(i)=ddMXY(1,J-NPSUP)*Gi(J,i)+ddFiRPIMxy(i)
    ENDDO
ENDDO


!TESTE DE FiRPIM
!UM=0.D0
!ZERO=0.D0
!ZEROUM=0.D0
!ZERO2=0.D0
!ZERO3=0.D0
!ZERO4=0.D0
!DO i=1,NPSUP !+TMON
!    UM=FiRPIM(i)+UM
!    ZERO=dFiRPIMX(i)+ZERO
!    ZEROUM=dFiRPIMY(i)+ZEROUM
!    ZERO2=ddFiRPIMX(i)+ZERO2
!    ZERO3=ddFiRPIMy(i)+ZERO3
!    ZERO4=ddFiRPIMXY(i)+ZERO4
!ENDDO
!UM=1.D0-UM
!WRITE(*,*) JJ,UM,'Fi'
!WRITE(*,*) JJ,ZERO,'dFiX'
!WRITE(*,*) JJ,ZEROUM,'dFiY'
!WRITE(*,*) JJ,ZERO2,'ddFiX'
!WRITE(*,*) JJ,ZERO3,'ddFiY'
!WRITE(*,*) JJ,ZERO4,'ddFixY'

DO i=1,NPSUP
    Fi(i)=FiRPIM(i)
    dFiX(i)=dFiRPIMX(i)
    dFiY(i)=dFiRPIMY(i)
    ddFiX(i)=ddFiRPIMX(i)
    ddFiY(i)=ddFiRPIMY(i)
    ddFiXY(i)=ddFiRPIMXY(i)
ENDDO


END SUBROUTINE RPIM2