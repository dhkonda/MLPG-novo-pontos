!  MLPGversao20.f90 
!
!  FUNCTIONS:
!  MLPGversao20 - Entry point of console application.
!
!****************************************************************************
!
!  PROGRAM: MLPGversao20
!
!  PURPOSE:  Entry point for the console application.
!   
!****************************************************************************

SUBROUTINE MLPG_NEW(adp,DESLOCAMENTO,NNO,NDESL)
USE ALOCAVEIS
IMPLICIT NONE
!n�o esquece de corrigir o problema dos angulos

CHARACTER*20::ARQS
CHARACTER::TRANS
INTEGER::NNO,NDESL
INTEGER:: I,J,K,L,NNOC,NDOM,MAXSUP,NCANT,MON,NPG,NPRES,ADP,NPT,PSUPG,TMON,NPSUP,INFO,SHAPE!TIPOSUP,iDCMENOR,
REAL*8::  E,NI,T,CARGA,RAIOi,D !,DPLINTAUX,TETADPLAUX !,TETA3,TETA4,TETA2,TETA1,TETAi,TETAF,JAC,XT,YT
REAL*8::  UM,ALFA,FAT,teste,DESLOCAMENTO !X1,Y1,X2,Y2,X3,Y3,X4,Y4,XV1,YV1,XV2,YV2,N1,N2,
!INTEGER,DIMENSION(3)::CPRES1,CPRES2,CPRES
!INTEGER,DIMENSION(2)::IDLMENOR
INTEGER,DIMENSION(:),ALLOCATABLE::iDC,FG,IPIV,NODUPLO !C,NO,CA,NOPRES,
INTEGER,DIMENSION(:,:),ALLOCATABLE::iDCadj,Cadj !TIPOPRES,
REAL*8,DIMENSION(:),ALLOCATABLE::XC,YC,XG,WG,DG,FGLOBAL,XGLOBAL,DESL !X,Y,
REAL*8,DIMENSION(:,:),ALLOCATABLE::VCANT,KGLOBAL,M,DMX,DMY !,VALORP !VPRES,
!vari�veis calculadas na rotina minimos quadrados
REAL*8,DIMENSION(:),ALLOCATABLE::Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY
!matrizes constantes utilizadas na integra��o
REAL*8,DIMENSION(2,2)::MATRIZD2 !,MATRIZN4,MATRIZN5
REAL*8,DIMENSION(3,3)::MATRIZD1 !,MATRIZN3
REAL*8,DIMENSION(3)::MATRIZR1,MATRIZR2,MATRIZD5 !,VSUPPRES,VSUPPRES1,VSUPPRES2
REAL*8,DIMENSION(2)::DPLINT !,TETADPL
REAL*8,DIMENSION(3,2)::MATRIZD3,MATRIZD4
REAL*8,DIMENSION(3,4)::MATRIZD7
REAL*8,PARAMETER:: ALFA2=1.E-5 !determina o tipo de imposicao de condicao de contorno para deslocamento

SHAPE=2 !SHAPE=1 MLS ---- SHAPE=2 RPIM

write(*,*) 'carregando dados de entrada'
CALL INPUT(NNOC,NDOM,MAXSUP,NCANT,MON,NPG,E,NI,T,CARGA,NPRES,ADP,NPT,PSUPG,TMON,D,ARQS,NPSUP,ALFA,FAT) !,X,Y,C,NO,CA,NOPRES,TIPOPRES,VPRES,

PSUPG=NPSUP

!subrotina input verificada - ok!
write(*,*) 'abrindo arquivo de saida'
CALL OPEN_OUTPUT(ARQS,MON,ALFA,ALFA2,SHAPE,E,NI,T,CARGA)

ALLOCATE(iDC(NCANT),XC(NCANT),YC(NCANT),VCANT(NCANT,2),iDCadj(NCANT,2),Cadj(NCANT,2),NODUPLO(NNOC))

write(*,*) 'calculando cantos'
CALL CANTOS(NCANT,NNOC,C,NO,NPT,XC,YC,X,Y,VCANT,iDC,CA,iDCadj,Cadj,NODUPLO)
!vcant,iDCadj,Cadj ok

ALLOCATE(XG(NPG),WG(NPG),KGLOBAL(3*NPT,3*NPT),FGLOBAL(3*NPT))

write(*,*) 'carregando pontos de integracao'    
CALL Calcula_Pontos_Gauss(NPG,XG,WG)

write(*,*) 'calculando matrizes nao variaveis'
CALL MATRIZES_FIXAS(D,Ni,MATRIZD1,T,CARGA,MATRIZR1,MATRIZR2,MATRIZD2,MATRIZD3,MATRIZD4,MATRIZD5,MATRIZD7)


!iniciando a matriz KGLOBAL e FGLOBAL
KGLOBAL=0.D0
FGLOBAL=0.D0

write(*,*) 'montando matriz de rigidez e vetor de forcas'
CALL INTEGRACAO(PSUPG,NPSUP,NPT,FAT,NCANT,XC,YC,VCANT,iDCadj,NNOC,iDC,NPRES,NPG,XG,WG,SHAPE,MON,TMON,&
    FGLOBAL,KGLOBAL,CARGA,MATRIZR1,MATRIZR2,MATRIZD5,MATRIZD2,MATRIZD1,MATRIZD3,MATRIZD4,MATRIZD7,ALFA2,ALFA,NODUPLO)

    
write(*,*) 'testando condicionamento da matriz de rigidez'
CALL TESTA_CONDICIONAMENTO(KGLOBAL,3*NPT,0)

ALLOCATE(XGLOBAL(3*NPT),IPIV(3*NPT))
XGLOBAL=0.D0

write(*,*) 'resolvendo o sistema de equacoes'
J=2 !J=1 usa o solver do lapack... j=2 utiliza o meu solver ou o svd
IF(J.EQ.1)THEN
    CALL DGETRF(3*NPT,3*NPT,KGLOBAL,3*NPT,IPIV,INFO)
    IF(INFO.EQ.0)THEN
        CONTINUE
    ELSE
        WRITE(*,*) 'PROBLEMA EM DGETRF'
        PAUSE
    ENDIF
    TRANS='N'
    XGLOBAL=FGLOBAL
    CALL DGETRS(TRANS,3*NPT,1,KGLOBAL,3*NPT,IPIV,XGLOBAL,3*NPT,INFO)
    IF(INFO.EQ.0)THEN
        CONTINUE
    ELSE
        WRITE(*,*) 'PROBLEMA EM DGETRS'
        PAUSE
    ENDIF
ELSEIF(J.EQ.2)THEN
    i=999 
    CALL r8mat_fs (3*NPT,KGLOBAL,FGLOBAL,XGLOBAL,i)
!    i=998
    IF(i.EQ.998)THEN
       
        CALL SVDSOLVE(kglobal,fglobal,3*npt,1,XGLOBAL)

    ENDIF
ENDIF

!interpolar os valores ficticios para obter os valores reais
ALLOCATE(DESL(3*NPT))

write(*,*) 'interpolando valores finais'
CALL INTERPOLA_VALORES_FINAIS(NPT,NPSUP,X,Y,MON,TMON,XGLOBAL,DESL,SHAPE)

write(*,*) 'gravando resultados'
CALL OUTPUT(ARQS,NCANT,iDC,XC,YC,Cadj,NNOC,NDOM,NPG,ADP,NPT,DESL,NO,PSUPG)

write(*,*) 'chegou ao fim'

!testar a fun��o fi e derivadas

!CALL TESTA_FI(NPT,NPSUP,X,Y,MON,TMON,SHAPE)

DESLOCAMENTO=DESL(NNO*NDESL)

DEALLOCATE(C,CA,CADJ,DESL,FGLOBAL,IDC,IDCADJ,IPIV,KGLOBAL)
DEALLOCATE(NO,NODUPLO,VCANT,WG,X,XC,XG,XGLOBAL,Y,YC)
DEALLOCATE(TIPOPRES,NOPRES,VPRES)


!READ(*,*)

ENDSUBROUTINE MLPG_NEW

    

