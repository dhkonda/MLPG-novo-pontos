    
!##################################################################
!subrotina para imprimir o andamento do programa
!##################################################################
SUBROUTINE EVO(i,NPT,PORC,PULA2)
IMPLICIT NONE

INTEGER::I,NPT,PORC,PORC2,PULA,PULA2
REAL::J
J=i*100./NPT

PORC2=INT(J)

IF(PORC.NE.PORC2)THEN
    WRITE(*,11) PORC2,'%'
    PORC=PORC2
ENDIF
    
PULA=MOD(PORC2,10)
IF(PULA.EQ.0.AND.PULA2.NE.PORC2)THEN
    PULA2=PORC2
    WRITE(*,13)
13 FORMAT (/)
ENDIF
    
10 FORMAT (F7.3,A1,\)
11 FORMAT (i4,A1,\)
   


RETURN
END SUBROUTINE EVO