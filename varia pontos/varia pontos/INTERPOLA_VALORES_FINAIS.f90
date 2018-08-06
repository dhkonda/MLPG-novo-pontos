SUBROUTINE INTERPOLA_VALORES_FINAIS(NPT,NPSUP,X,Y,MON,TMON,XGLOBAL,DESL,shape)
    IMPLICIT NONE

INTEGER:: NPT,NPSUP,MON,TMON,i,J,K,SHAPE
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(NPSUP)::FG
REAL*8,DIMENSION(3*NPT)::DESL,XGLOBAL
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(NPSUP)::DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY,DDMX,DDMY,DDMXY
REAL*8,DIMENSION(TMON,NPSUP)::AiB,dAiXB,dAiYB,AidBX,AidBY,Sb
REAL*8,DIMENSION(NPSUP,NPSUP)::Sa

DESL=0.D0

DO i=1,NPT
     
    CALL MONTA_SUPORTE(X(i),Y(i),X,Y,NPT,NPSUP,NO,FG,DG)                
    
    !aplicar o MMQM no nós do contorno para  obter as funções de forma para impor as condicoes de contorno

if(shape.eq.1)then
    CALL MINIMOS_QUADRADOS(i,X(i),Y(i),X,Y,NPSUP,FG,NPT,MON,TMON,DG,AiB,dAixB,dAiyB,AidBx,AidBy,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
elseif(shape.eq.2)then
            CALL RPIM(i,X(i),Y(i),NPT,NPSUP,X,Y,FG,MON,TMON,Sa,Sb,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
endif    
    DO J=1,NPSUP
        DESL(i*3-2)=Fi(J)*XGLOBAL(FG(J)*3-2)+DESL(i*3-2)
        DESL(i*3-1)=Fi(J)*XGLOBAL(FG(J)*3-1)+DESL(i*3-1)
        DESL(i*3)=Fi(J)*XGLOBAL(FG(J)*3)+DESL(i*3)
!        write(*,*) fi(j),xglobal(fg(j)*3-2)
    ENDDO
ENDDO
    
    
    
    
END SUBROUTINE INTERPOLA_VALORES_FINAIS