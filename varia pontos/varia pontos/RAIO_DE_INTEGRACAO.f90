!################################################################################
!SUBROTINA QUE CALCULA O TAMANHO DO RAIO DE INTEGRAÇÃO DO PONTO i
!################################################################################
SUBROUTINE RAIO_DE_INTEGRACAO(i,NPT,X,Y,NNOC,NODUPLO,RAIOi)
IMPLICIT NONE

INTEGER::i,J,NPT,NNOC
INTEGER,DIMENSION(NNOC)::NODUPLO
REAL*8::AUX1,AUX2,RAIOi,A1
REAL*8,DIMENSION(NPT)::X,Y,DiST

RAIOi=1.d308
DO J=1,NPT
    AUX1=(X(i)-X(J))**2
    AUX2=(Y(i)-Y(J))**2
    DiST(J)=DSQRT(AUX1+AUX2)
    IF(DiST(J).GT.0)THEN
        IF(RAIOi.GT.DiST(J))THEN
            RAIOi=DiST(J)
        ENDIF
    ENDIF
ENDDO

!calculo do raio de integração para os nós duplos
IF(i.LT.NNOC)THEN
    IF(NODUPLO(i).EQ.1)THEN
        RAIOi=1.d308
        DO J=1,NPT
            AUX1=(X(i)-X(J))**2
            AUX2=(Y(i)-Y(J))**2
            DiST(J)=DSQRT(AUX1+AUX2)
            IF(DiST(J).GT.0)THEN
                
                IF(RAIOi.GT.DiST(J))THEN
                    A1=RAIOi
                    RAIOi=DiST(J)
                ELSEIF(A1.GT.DIST(J))THEN
                    A1=DIST(J)
                ENDIF
            ENDIF        
        ENDDO
        RAIOI=A1
    ENDIF
ENDIF


ENDSUBROUTINE RAIO_DE_INTEGRACAO
