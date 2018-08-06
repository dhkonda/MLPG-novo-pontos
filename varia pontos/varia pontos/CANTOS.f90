!###################################################################
!ESTA SUBROTINA CALCULA OS CANTOS ADJACENTES A UM CANTO E TAMBÉM
!OS VERSORES DOS LADOS DA PLACA (SENTIDO ANTI-HORÁRIO POSITIVO)
!###################################################################

SUBROUTINE CANTOS(NCANT,NNOC,C,NO,NPT,XC,YC,X,Y,VCANT,iDC,CA,iDCadj,Cadj,NODUPLO)
USE TIPO
IMPLICIT NONE

INTEGER::i,J,NCANT,NNOC,NPT
INTEGER,DIMENSION(NPT)::NO
INTEGER,DIMENSION(NNOC)::C,CA,NODUPLO
INTEGER,DIMENSION(NCANT)::iDC
INTEGER,DIMENSION(NCANT,2)::iDCadj
REAL*8::AUX1,AUX2,AUX3
REAL*8,DIMENSION(NCANT)::XC,YC,DD
REAL*8,DIMENSION(NPT)::X,Y
INTEGER,DIMENSION(NCANT,2)::Cadj
REAL*8, DIMENSION(NCANT,2)::Xadj,Yadj
REAL*8,DIMENSION(NCANT,2)::VCANT
REAL*8,DIMENSION(NNOC)::XD,YD
REAL*8,DIMENSION(NCANT,NPT)::DIST
TYPE(GROUP),DIMENSION(NPT)::DORD
REAL*8::AJUSTEDD

AJUSTEDD=0.1D0
NODUPLO=0

J=0
DO i=1,NNOC
	IF (C(I).EQ.1) THEN
		!contador para identificar a quantidade de cantos
        J=J+1
		IF(J.GT.NCANT)THEN
			WRITE(*,*) 'ERRO NO LANCAMETO DOS CANTOS DA PLACA'
			READ(*,*)
        ENDIF
        !identificando se o Nó(i) é nó de canto
		iDC(J)=NO(I)
   	    !salvando as coordenadas dos nós de canto
        XC(J)=X(i)
		YC(J)=Y(i)
        !Salvando o canto anterior como canto adjacente 1.
		Cadj(J,1)=CA(i)
	ENDIF
ENDDO

DO i=1,NCANT		 !encontrando o canto posterior a um canto
	DO J=1,NCANT
		IF(Cadj(i,1).EQ.iDC(J)) THEN
            !salvando o canto posterior como canto adjacente 2
			Cadj(J,2)=iDC(i)
			!salvando a identidade dos cantos adjacentes
            iDCadj(i,1)=J
			iDCadj(J,2)=i
		ENDIF
	ENDDO																  
ENDDO
DO i=1,NCANT		 
		DO J=1,NCANT
			IF(iDC(i).EQ.Cadj(J,1)) THEN
				!identificando e salvando as coordenadas dos cantos adjacentes 1
                Xadj(J,1)=XC(i)
				Yadj(J,1)=YC(i)
            ELSE IF (iDC(i).EQ.Cadj(J,2)) THEN
                !identificando e salvando as coordenadas dos cantos adjacentes 2
				Xadj(J,2)=XC(i)
				Yadj(J,2)=YC(i)
			ENDIF			
        ENDDO
ENDDO

DO i=1,NCANT						!impressao da conectividade dos cantos da placa e calculo dos versores dos lados da placa.
    AUX1=(Xadj(i,2)-XC(i))**2
    AUX2=((Yadj(i,2))-YC(i))**2
    AUX3=DSQRT(AUX1+AUX2)
    !calculando os versores dos lados da placa. Orientação positiva anti-horária
	VCANT(i,1)=(Xadj(i,2)-XC(i))/AUX3
	VCANT(i,2)=	((Yadj(i,2))-YC(i))/AUX3
ENDDO

!definindo as novas coordenadas dos nós duplos
XD=0.D0
YD=0.D0
    !encontrar o ponto mais proximo aos cantos
DO i=1,NCANT
    DO J=1,NPT
        DIST(i,j)=DSQRT((XC(i)-X(J))**2+(YC(i)-Y(J))**2)
        DORD(J)%VALOR=DIST(i,J)
        DORD(J)%ORDEM=NO(J)
        DORD(J)%ORDEM2=J
    ENDDO
    CALL QuickSort(DORD,NPT) !me interessa só a 3 menos distancia, já que as duas primeiras são do nó duplo
    DD(i)=DORD(3)%VALOR
ENDDO
    
!DO i=1,NNOC
!    IF(C(i).EQ.1)THEN !é canto...
!        !encontrar o versor da direcao deste canto
!        DO J=1,NCANT
!            IF(iDC(J).EQ.NO(i))THEN
!                !CALCULO A NOVA COORDENADA DO PONTO AQUI
!                XD(i)=X(i)+DD(J)*VCANT(J,1)*AJUSTEDD
!                YD(i)=Y(i)+DD(J)*VCANT(J,2)*AJUSTEDD
!            ENDIF
!        ENDDO
!    ENDIF
!ENDDO

!ajuste dos nós de canto sem a flag C(i)=1
DO i=1,NNOC
    IF(C(i).NE.1)THEN
        DO J=1,NCANT
            IF(XC(J).EQ.X(i).AND.YC(J).EQ.Y(i))THEN
                XD(i)=X(i)-DD(J)*VCANT(iDCadj(J,1),1)*AJUSTEDD
                YD(i)=Y(i)-DD(J)*VCANT(iDCadj(J,1),2)*AJUSTEDD
            ENDIF
        ENDDO    
    ENDIF
ENDDO

DO i=1,NNOC
    IF(XD(i).NE.0.OR.YD(i).NE.0.)THEN
        X(i)=XD(i)
        Y(i)=YD(i)
        NODUPLO(i)=1
    ENDIF    
ENDDO


RETURN
ENDSUBROUTINE CANTOS    
