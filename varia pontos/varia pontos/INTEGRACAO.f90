SUBROUTINE INTEGRACAO(PSUPG,NPSUP,NPT,FAT,NCANT,XC,YC,VCANT,iDCadj,NNOC,iDC,NPRES,NPG,XG,WG,SHAPE,MON,TMON,&    
        FGLOBAL,KGLOBAL,CARGA,MATRIZR1,MATRIZR2,MATRIZD5,MATRIZD2,MATRIZD1,MATRIZD3,MATRIZD4,MATRIZD7,ALFA2,ALFA,NODUPLO)
        
    USE ALOCAVEIS
    IMPLICIT NONE

INTEGER::i,J,K,L,PSUPG,NPT,NCANT,NNOC,TIPOSUP,iDCMENOR,NPRES,NPG,INFO,SHAPE,MON,TMON,NPSUP,PORC,PULA2
INTEGER,DIMENSION(2)::IDLMENOR
INTEGER,DIMENSION(3)::CPRES1,CPRES2,CPRES
INTEGER,DIMENSION(:),ALLOCATABLE::FG
INTEGER,DIMENSION(NNOC)::NODUPLO
INTEGER,DIMENSION(NCANT)::iDC
INTEGER,DIMENSION(NCANT,2)::iDCadj
REAL*8::RAIOi,FAT,TETA1,TETA2,TETA3,TETA4,X1,Y1,X2,Y2,X3,Y3,X4,Y4,TETAi,TETAF,JAC,XT,YT,CARGA,XV1,YV1,XV2,YV2,N1,N2,DPLINTAUX,TETADPLAUX
REAL*8,DIMENSION(2)::DPLINT,TETADPL
REAL*8,DIMENSION(NPG)::XG,WG
REAL*8,DIMENSION(NCANT)::XC,YC
REAL*8,DIMENSION(NCANT,2)::VCANT
REAL*8,DIMENSION(TMON,NPSUP)::AiB,dAixB,dAiyB,AidBx,AidBy,Sb
REAL*8,DIMENSION(NPSUP,NPSUP)::Sa
REAL*8,DIMENSION(:),ALLOCATABLE::DG,Fi,dFiX,dFiY,ddFiX,ddFiY,ddFiXY
REAL*8,DIMENSION(NNOC,3)::VALORP
REAL*8,DIMENSION(3)::VSUPPRES,VSUPPRES1,VSUPPRES2
REAL*8,DIMENSION(3)::MATRIZR1,MATRIZR2,MATRIZD5
REAL*8,DIMENSION(3*NPT)::FGLOBAL
REAL*8,DIMENSION(3*NPT,3*NPT)::KGLOBAL
REAL*8,DIMENSION(2,2)::MATRIZD2,MATRIZN4,MATRIZN5
REAL*8,DIMENSION(3,2)::MATRIZD3,MATRIZD4
REAL*8,DIMENSION(3,3)::MATRIZD1,MATRIZN3
REAL*8,DIMENSION(3,4)::MATRIZD7
REAL*8,PARAMETER::Pi=4.d0*datan(1.d0) !3.14159265358979323846d0
REAL*8:: ALFA2,ALFA
CHARACTER*20::ARQS

WRITE(*,23) 
PORC=0
PULA2=0
DO i=1,NPT

    ALLOCATE(FG(PSUPG),DG(PSUPG),Fi(PSUPG),dFix(PSUPG),dFiy(PSUPG),ddFiX(PSUPG),ddFiY(PSUPG),ddFiXY(PSUPG))
    
    CALL RAIO_DE_INTEGRACAO(i,NPT,X,Y,NNOC,NODUPLO,RAIOi)

    raioi=FAT*raioi

    CALL CLASSIFICA_PONTO(i,NCANT,NPT,X,Y,XC,YC,VCANT,RAIOi,iDCadj,NNOC,C,TIPOSUP,TETA1,TETA2,TETA3,TETA4,DPLINT,TETADPL,&
        X1,Y1,X2,Y2,X3,Y3,X4,Y4,iDCMENOR,iDLMENOR)
 
    CALL TIPO_CONTORNO(TIPOSUP,NNOC,NPT,NO,iDLMENOR,NCANT,NPRES,NOPRES,VPRES,TIPOPRES,iDC,iDCadj,iDCMENOR,&
    CPRES1,CPRES2,VSUPPRES1,VSUPPRES2,VALORP)

    CALL MONTA_SUPORTE(X(i),Y(i),X,Y,NPT,PSUPG,NO,FG,DG)

if(shape.eq.1)then
            CALL MINIMOS_QUADRADOS(i,X(i),Y(i),X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY) !Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY
elseif(shape.eq.2)then
            CALL RPIM(i,X(i),Y(i),NPT,PSUPg,X,Y,FG,MON,TMON,Sa,Sb,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
endif

    
    
    DO J=1,NPG
        DO K=1,NPG
            
            !a integral de domínio de teta2 à teta1 é feita para todos os TIPOSUP
            INFO=1
            CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)    
            
            CALL T_COORD_CIRCULO(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,XT,YT,RAIOi,WG(J),WG(K),JAC)

!            CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
!
!if(shape.eq.1)then
!            CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY) !Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

            CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)
      
            IF(TIPOSUP.EQ.2)THEN

                !integrar o triangulo
                !para tiposup2 a integral é feita de teta1 à teta2

                INFO=2
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
 
                DPLINTAUX=DPLINT(1)
                TETADPLAUX=TETADPL(1)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

!                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
! 
!if(shape.eq.1)then 
!                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                    FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)
                
                !integrais de contorno apenas com o laço J 
                IF(K.EQ.1)THEN
                
                    XV1=X1
                    YV1=Y1
                    XV2=X2
                    YV2=Y2
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
!                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
!
!if(shape.eq.1)then
!                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif                    
                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES1(L)
			            VSUPPRES(L)=VSUPPRES1(L)
		            ENDDO
                    
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
                    
                ENDIF
                
            ELSEIF(TIPOSUP.EQ.3) THEN

                !integrar o primeiro triangulo
                !Para tiposup3 o primeiro triangulo é integrado de teta1 à teta3

                INFO=3
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi) 
                
                DPLINTAUX=DPLINT(1)
                TETADPLAUX=TETADPL(1)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

!                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
!
!if(shape.eq.1)then                
!                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                    FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)                

                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X1
                    YV1=Y1
                    XV2=X3
                    YV2=Y3
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
!                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
!
!if(shape.eq.1)then
!                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES1(L)
			            VSUPPRES(L)=VSUPPRES1(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
                
                ENDIF

                !integrar o segundo triangulo
                !Para tiposup3 o segundo triangulo é integrado de teta3 à teta2                

                INFO=31
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)
                
                DPLINTAUX=DPLINT(2)
                TETADPLAUX=TETADPL(2)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

!                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
!
!if(shape.eq.1)then                
!                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                    FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)                
           
                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X3
                    YV1=Y3
                    XV2=X2
                    YV2=Y2
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
!                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
!
!if(shape.eq.1)then                    
!                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES2(L)
			            VSUPPRES(L)=VSUPPRES2(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
                                
                ENDIF
        
            ELSEIF(TIPOSUP.EQ.4)THEN
                
                !integrar o setor circular
                !a integral é feita de teta3 à teta4
                INFO=4
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)                
                
                CALL T_COORD_CIRCULO(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,XT,YT,RAIOi,WG(J),WG(K),JAC)
                
!                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
!
!if(shape.eq.1)then                
!                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                    FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)
                             
                !integrar o primeiro triangulo
                !Para tiposup4 o primeiro triangulo é integrado de teta1 à teta3                
                INFO=41
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)             
                
                DPLINTAUX=DPLINT(1)
                TETADPLAUX=TETADPL(1)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

!                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
!   
!if(shape.eq.1)then                
!                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                    FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)                

                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X1
                    YV1=Y1
                    XV2=X3
                    YV2=Y3
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
!                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
!          
!if(shape.eq.1)then                    
!                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES1(L)
			            VSUPPRES(L)=VSUPPRES1(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
                                          
                ENDIF

                !integrar o segundo triangulo
                !Para tiposup4 o segundo triangulo é integrado de teta4 à teta2                
                INFO=42
                CALL CORRIGE_TETA(INFO,TETA1,TETA2,TETA3,TETA4,TETAi,TETAF,Pi)  
    
                DPLINTAUX=DPLINT(2)
                TETADPLAUX=TETADPL(2)
                
                CALL T_COORD_TRI(X(i),Y(i),XG(J),XG(K),TETAi,TETAF,WG(J),WG(K),DPLINTAUX,TETADPLAUX,XT,YT,JAC)

!                CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)
!
!if(shape.eq.1)then                
!                CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                CALL INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,MON,&
                    FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)    

                IF(K.EQ.1)THEN
                !INTEGRAR O CONTORNO UTILIZANDO O LAÇO J
                   
                    XV1=X4
                    YV1=Y4
                    XV2=X2
                    YV2=Y2
                    CALL T_COORD_CONTORNO(XG(J),XV1,YV1,XV2,YV2,WG(J),XT,YT,JAC)
                                    
!                    CALL MONTA_SUPORTE(XT,YT,X,Y,NPT,PSUPG,NO,FG,DG)                
!                
!if(shape.eq.1)then                    
!                    CALL MINIMOS_QUADRADOS(i,XT,YT,X,Y,PSUPG,FG,NPT,MON,TMON,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!elseif(shape.eq.2)then
!            CALL RPIM(i,XT,YT,NPT,PSUPg,X,Y,FG,MON,TMON,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
!endif

                    CALL MATRIZES_NORMAIS(XV1,YV1,XV2,YV2,MATRIZN3,MATRIZN4,MATRIZN5,N1,N2)
                    
                    DO L=1,3
			            CPRES(L)=CPRES2(L)
			            VSUPPRES(L)=VSUPPRES2(L)
		            ENDDO
                    !AGORA É SÓ CHAMAR UMA ROTINA QUE RESOLVE AS INTEGRAIS DE CONTORNO
                    CALL INTEGR_CONTORNO(i,NPT,X,Y,XT,YT,RAIOi,JAC,CPRES,VSUPPRES,MATRIZN3,MATRIZN4,MATRIZN5,Sa,Sb,SHAPE,&
                        PSUPG,FG,N1,N2,MON,TMON,AiB,dAixB,dAiyB,AidBx,AidBy,MATRIZD3,KGLOBAL,FGLOBAL,MATRIZR2,MATRIZD4,MATRIZD5,ALFA,ALFA2,MATRIZD7)
                                          
                ENDIF

            ENDIF
 
        ENDDO
    ENDDO    
    
    !aplicar o MMQM no nós do contorno para  obter as funções de forma para impor as condicoes de contorno
    
    DEALLOCATE(FG,DG,Fi,dFix,dFiy,ddFiX,ddFiY,ddFiXY)
    
    !aplicar condições de contorno (somente para nós do contorno)
    CALL CONDICAO_DE_CONTORNO(ALFA,ALFA2,i,NNOC,C,NCANT,iDLMENOR,NO,NPT,iDC,X,Y,NPSUP,MON,TMON,KGLOBAL,FGLOBAL,VALORP,TIPOPRES,VCANT,SHAPE)

CALL EVO(i,NPT,PORC,PULA2)

ENDDO

WRITE(*,23)
23 FORMAT (/)
    
END SUBROUTINE INTEGRACAO