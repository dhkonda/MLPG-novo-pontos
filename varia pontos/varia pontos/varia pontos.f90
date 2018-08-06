!  variapontos.f90 
!
!  FUNCTIONS:
!  variapontos - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: variapontos
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program variapontos

    implicit none
    INTEGER ADP,NNO,NDESL,STOPADP,I
    REAL*8::DESLOCAMENTO
    CHARACTER*20:: ARQTESTE
    ! Variables

    ! Body of variapontos
    NNO=181
    NDESL=3
    
    ARQTESTE='VARSUP.TXT'
    OPEN (3,FILE=ARQTESTE,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
    
    WRITE(*,*) 'COMEÇAR O SUPORTE COM QUANTOS PONTOS ADICIONAIS?'
    READ(*,*) ADP
    WRITE(*,*) 'QUANDO PARAR?'
    READ(*,*) STOPADP
    DO I=ADP,STOPADP
        CALL MLPG_NEW(i,DESLOCAMENTO,NNO,NDESL)    
        WRITE(3,10) I,DESLOCAMENTO
10 FORMAT (I3,5X,E12.5)        
        CLOSE ( 1, STATUS='KEEP') 
        CLOSE ( 2, STATUS='KEEP') 
        
    ENDDO
    PAUSE

    end program variapontos

