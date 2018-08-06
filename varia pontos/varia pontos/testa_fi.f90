SUBROUTINE TESTA_FI(NPT,NPSUP,X,Y,MON,TMON,SHAPE)
    IMPLICIT NONE
    
INTEGER:: NPT,NPSUP,MON,TMON,i,J,K,SHAPE
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(NPSUP)::FG
REAL*8::UM,ZERO,ZEROUM,AC,max,min,maxdx,mindx,maxddx,minddx,maxdy,mindy,maxddy,minddy,minddxy,maxddxy
REAL*8,DIMENSION(NPT)::X,Y,Z,Zap,dX,dY,dXap,dYap,ddxap,ddyap,ddx,ddy,ddxy,ddxyap
REAL*8,DIMENSION(NPSUP)::DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY

!DESL=0.D0
Zap=0.D0
dXap=0.D0
dYap=0.D0
ddxap=0.d0
ddyap=0.d0
ddxyap=0.d0

!!CALCULO DO VALOR ANALITICO --- FUNCAO=XCOS(X)YSEN(y)
!DO i=1,NPT
!    Z(i)=X(i)*COS(X(i))*Y(i)*SIN(Y(i))
!    dX(i)=COS(X(i))*Y(i)*SIN(Y(i))-X(i)*SIN(X(i))*Y(i)*SIN(Y(i))
!    dY(i)=X(i)*COS(X(i))*SIN(Y(i))+X(i)*COS(X(i))*Y(i)*COS(Y(i))
!ENDDO

!!CALCULO DO VALOR ANALITICO --- FUNCAO=X^2+y^2
!DO i=1,NPT
!    Z(i)=X(i)**3.+Y(i)**3.
!    dX(i)=3.*(x(i))**2.
!    dY(i)=3.*y(i)**2.
!ENDDO

shape=2

!CALCULO DO VALOR ANALITICO --- FUNCAO=X^2+y^2
DO i=1,NPT
    Z(i)=1.d0+2.d0*x(i)**2+2.d0*y(i)**2 
    dX(i)=4.d0*x(i)**1
    dY(i)=4.d0*y(i)**1
    ddx(i)=4.d0 !12.d0*x(i)**1
    ddy(i)=4.d0 !12.d0*y(i)**1
    ddxy(i)=0.d0
ENDDO
max=maxval(z)
min=minval(z)
maxdx=maxval(dx)
mindx=minval(dx)
maxddx=maxval(ddx)
minddx=minval(ddx)
maxdy=maxval(dy)
mindy=minval(dy)
maxddy=maxval(ddy)
minddy=minval(ddy)
maxddxy=maxval(ddxy)
minddxy=minval(ddxy)

DO i=1,NPT
     
    CALL MONTA_SUPORTE(X(i),Y(i),X,Y,NPT,NPSUP,NO,FG,DG)                
    
    !aplicar o MMQM no nós do contorno para  obter as funções de forma para impor as condicoes de contorno
if(shape.eq.1)then                    
            CALL MINIMOS_QUADRADOS(i,X(i),Y(i),X,Y,NPSUP,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
elseif(shape.eq.2)then
            CALL RPIM(i,X(i),Y(i),NPT,NPSUP,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
endif    
    
    DO J=1,NPSUP
        Zap(i)=Fi(J)*Z(FG(J))+Zap(i)
        dXap(i)=dFix(J)*Z(FG(J))+dXAP(i)
        dYap(i)=dFiy(J)*Z(FG(J))+dYAP(i)
        ddxap(i)=ddFix(j)*Z(FG(J))+ddxAP(i)
        ddyap(i)=ddFiy(j)*Z(FG(J))+ddyAP(i)
    ENDDO    
ENDDO

!do i=1,npt
!    write(*,*) i,z(I),zap(I)
!enddo
!write(*,*) "u aaproximado"
!pause
!
!do i=1,npt
!    write(*,*) i,dx(I),dxap(I)
!enddo
!write(*,*) "dx aaproximado"
!pause
!
!do i=1,npt
!    write(*,*) i,dy(I),dyap(I)
!enddo
!write(*,*) "dy aaproximado"
!pause
!calculando o erro do valor aproximado e imprimindo em porcentagem.
AC=0.D0
DO i=1,NPT
    UM=ZAP(i)-Z(i)
!    UM=UM*100./Z(i)
!    WRITE(*,*) I,UM,z(I)
    AC=AC+ABS(UM)
ENDDO
WRITE(*,*) AC/npt,'ERRO AC DA FUNCAO APROXIMADA'
write(*,*) max,'maximo'
write(*,*) min,'minimo'
PAUSE
AC=0.D0
DO i=1,NPT
    UM=dXap(i)-dX(i)
!    UM=UM*100./dX(i)
!    WRITE(*,*) I,UM,dx(i)
    AC=AC+ABS(UM)
ENDDO
WRITE(*,*) AC/npt,'ERRO AC DA FUNCAO dX APROXIMADA'
write(*,*) maxdx,'maximo'
write(*,*) mindx,'minimo'
PAUSE
AC=0.D0
DO i=1,NPT
    UM=dYap(i)-dY(i)
!    UM=UM*100./dY(i)
!    WRITE(*,*) I,UM,dy(I)
    AC=AC+ABS(UM)
ENDDO
WRITE(*,*) AC/npt,'ERRO AC DA FUNCAO dY APROXIMADA'
write(*,*) maxdy,'maximo'
write(*,*) mindy,'minimo'
pause
AC=0.D0
DO i=1,NPT
    UM=ddxap(i)-ddx(i)
!    UM=UM*100./dY(i)
!    WRITE(*,*) I,UM,dy(I)
    AC=AC+ABS(UM)
ENDDO
WRITE(*,*) AC/npt,'ERRO AC DA FUNCAO ddx APROXIMADA'
write(*,*) maxddx,'maximo'
write(*,*) minddx,'minimo'
pause
AC=0.D0
DO i=1,NPT
    UM=ddYap(i)-ddY(i)
!    UM=UM*100./dY(i)
!    WRITE(*,*) I,UM,dy(I)
    AC=AC+ABS(UM)
ENDDO
WRITE(*,*) AC/npt,'ERRO AC DA FUNCAO ddY APROXIMADA'
write(*,*) maxddy,'maximo'
write(*,*) minddy,'minimo'
pause
AC=0.D0
DO i=1,NPT
    UM=ddxYap(i)-ddxY(i)
!    UM=UM*100./dY(i)
!    WRITE(*,*) I,UM,dy(I)
    AC=AC+ABS(UM)
ENDDO
WRITE(*,*) AC/npt,'ERRO AC DA FUNCAO ddXY APROXIMADA'
write(*,*) maxddxy,'maximo'
write(*,*) minddxy,'minimo'


ENDSUBROUTINE TESTA_FI