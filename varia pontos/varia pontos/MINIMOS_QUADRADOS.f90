!##########################################################################################
!SUBROTINA PARA MONTAR A FUNÇÃO DE FORMA PARA O PONTO BASE XT,YT
!##########################################################################################    
SUBROUTINE  MINIMOS_QUADRADOS(JJ,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)  
IMPLICIT NONE

INTEGER::i,J,K,NPT,PSUPG,MON,TMON,JJ
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::XT,YT,RSG,UM,ZERO,ZEROUM
REAL*8,DIMENSION(1)::XBM,YBM
REAL*8,DIMENSION(PSUPG)::XB,YB,DG,W,DWX,DWY,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(:,:),ALLOCATABLE::M,DMX,DMY,DDMX,DDMY,DDMXY
REAL*8,DIMENSION(TMON,PSUPG)::B,dBx,dBy,AiB,dAiXB,dAiYB,AidBX,AidBY
REAL*8,DIMENSION(TMON,TMON)::A,Ai,IDENT,dAx,dAy,dAiX,dAiY,AUX1,AUX2

!pause

DO i=1,PSUPG
    XB(i)=X(FG(i))
    YB(i)=Y(FG(i))
ENDDO

!montando a matriz P(armazenada em M) para o suporte do ponto XT,YT

ALLOCATE(M(PSUPG,TMON),DMX(PSUPG,TMON),DMY(PSUPG,TMON),DDMX(PSUPG,TMON),DDMY(PSUPG,TMON),DDMXY(PSUPG,TMON))

CALL MONTA_BASE_MONOMIAL(MON,TMON,PSUPG,XB,YB,M,DMX,DMY,DDMX,DDMY,DDMXY)

!calculo do raio da função peso para aplicação do MMQM no suporte do ponto de Gauss
RSG=1.1*DG(PSUPG)
!montando a matriz W (spline de 4 ordem) e suas derivadas
!WRITE(*,*) RSG
DO i=1,PSUPG
    W(i)=1.D0-6.D0*((DG(i)/RSG)**2)+8.D0*((DG(i)/RSG)**3)-3.D0*((DG(i)/RSG)**4)
    DWX(i)=(12.D0*(XT-X(FG(i)))/RSG**2)*(-1.D0+2.D0*((DG(i))/RSG)-((DG(i))/RSG)**2) !essas derivadas estão corretas!!! 
    DWY(i)=(12.D0*(YT-Y(FG(i)))/RSG**2)*(-1.D0+2.D0*((DG(i))/RSG)-((DG(i))/RSG)**2)
ENDDO

!construção da matriz B e suas derivadas
DO i=1,TMON
    DO J=1,PSUPG
        B(i,J)=M(J,i)*W(J)
        dBx(i,J)=M(J,i)*DWX(J)
        dBy(i,J)=M(J,i)*DWY(J)
    ENDDO
ENDDO
!construção da matriz A e suas derivadas
A=0.
dAx=0.
dAy=0.
A=MATMUL(B,M)
dAx=MATMUL(dBX,M)
dAy=MATMUL(dBy,M)

CALL TESTA_CONDICIONAMENTO(A,TMON,jj)

!determinando a matriz A inversa
DO i=1,TMON
    IDENT=0.
    IDENT(i,i)=1.
    CAll r8mat_fs(TMON,A,IDENT(:,i),Ai(:,i),JJ)
!    CALL SVDSOLVE(A,IDENT(:,i),tmon,1,Ai(:,i))    
ENDDO
!calculo das derivadas da matriz A inversa
AUX1=0.
AUX2=0.
AUX1=MATMUL(dAx,Ai)
AUX2=MATMUL(dAy,Ai)
dAiX=0.
dAiY=0.
dAiX=MATMUL(-Ai,AUX1)
dAiY=MATMUL(-Ai,AUX2)
!calculo dos parametros que serão utilizados na integração
AiB=0.
dAiXB=0.
dAiYB=0.
AidBX=0.
AidBY=0.

AiB=MATMUL(Ai,B)
dAiXB=MATMUL(dAiX,B)
dAiYB=MATMUL(dAiY,B)
AidBX=MATMUL(Ai,dBX)
AidBY=MATMUL(Ai,dbY)

DEALLOCATE(M,DMX,DMY,DDMX,DDMY,DDMXY)

!Obtenção da função de forma e suas derivadas.

!construindo a base monomial para o ponto de gauss XT,YT
XBM(1)=XT
YBM(1)=YT

ALLOCATE(M(1,TMON),DMX(1,TMON),DMY(1,TMON),DDMX(1,TMON),DDMY(1,TMON),DDMXY(1,TMON))

CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XBM,YBM,M,DMX,DMY,DDMX,DDMY,DDMXY)

dFix=0.D0
dFiy=0.D0
Fi=0.D0
UM=0.D0
ZERO=0.D0
ZEROUM=0.D0
DO J=1,PSUPG
    DO K=1,TMON
        !dFi é composto por uma soma de tres termos
        !calculo do primeiro termo
        dFiX(J)=DMX(1,K)*AiB(K,J)+dFiX(J)
        dFiY(J)=DMY(1,K)*AiB(K,J)+dFiY(J)
        !calculo do segundo termo
        dFiX(J)=M(1,K)*dAiXB(K,J)+dFiX(J)
        dFiY(J)=M(1,K)*dAiYB(K,J)+dFiY(J)
        !calculo do terceiro termo
        dFiX(J)=M(1,K)*AidBX(K,J)+dFiX(J)
        dFiY(J)=M(1,K)*AidBY(K,J)+dFiY(J)
        !calculo de Fi
        Fi(J)=M(1,K)*AiB(K,J)+Fi(J)
    ENDDO
    UM=UM+Fi(J)
    ZERO=dFiX(J)+ZERO
    ZEROUM=dFiY(J)+ZEROUM
ENDDO
UM=UM-1.D0

!
!dfix=-dfix
!dfiy=-dfiy
!IF(ABS(ZEROUM).GT.1.E-2)THEN
!    write(*,*) ZEROUM
!    WRITE(*,*) 'PROBLEMA NA dFiY DO PONTO',JJ
!    WRITE(*,*) 'AUMENTE OS PONTOS ADICIONAIS NO SUPORTE - cod xYX'
!    PAUSE
!ENDIF
!IF(ABS(ZERO).GT.1.E-2)THEN
!    write(*,*) ZERO
!    WRITE(*,*) 'PROBLEMA NA dFiX DO PONTO',I
!    WRITE(*,*) 'AUMENTE OS PONTOS ADICIONAIS NO SUPORTE - cod xY'
!    PAUSE
!ENDIF
!IF(ABS(UM).GT.0.1)THEN
!    write(*,*) um
!    WRITE(*,*) 'PROBLEMA NA FI DO PONTO',I
!    WRITE(*,*) 'AUMENTE OS PONTOS ADICIONAIS NO SUPORTE - cod x'
!    PAUSE
!ENDIF

DEALLOCATE(M,DMX,DMY,DDMX,DDMY,DDMXY)


RETURN
END SUBROUTINE MINIMOS_QUADRADOS