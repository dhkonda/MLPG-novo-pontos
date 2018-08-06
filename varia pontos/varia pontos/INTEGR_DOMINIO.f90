     
!#############################################################################################
!Esta subrotina faz a integração de domínio para um ponto de gauss
!#############################################################################################
SUBROUTINE INTEGR_DOMINIO(i,NPT,X,Y,XT,YT,RAIOi,MATRIZD1,AiB,dAixB,dAiyB,AidBx,AidBy,Sa,Sb,TMON,PSUPG,&
    MON,FG,MATRIZR1,MATRIZD2,CARGA,JAC,KGLOBAL,FGLOBAL,SHAPE)
IMPLICIT NONE

INTEGER::i,J,K,L,NPT,TMON,PSUPG,MON,SHAPE
INTEGER,DIMENSION(PSUPG)::FG
REAL*8::RAIOi,XT,YT,Dist,W,DWX,DWY,CARGA,JAC,q,ALFAC,dc
REAL*8,DIMENSION(NPT)::X,Y
REAL*8,DIMENSION(3,3)::MATRIZW1,MATRIZD1,AUX33
REAL*8,DIMENSION(TMON,PSUPG)::AiB,dAiXB,dAiYB,AidBX,AidBY,Sb !B,dBx,dBy,
REAL*8,DIMENSION(PSUPG,PSUPG)::Sa,DIST3
REAL*8,DIMENSION(PSUPG)::ddFix,ddFiy,ddFixy,dFix,dFiy,Fi,DIST2,R,dRx,dRy,ddRx,ddRy,ddRxy
REAL*8,DIMENSION(1,TMON)::M,DMX,DMY,DDMX,DDMY,DDMXY
REAL*8,DIMENSION(2,3*PSUPG)::MATRIZFI2
REAL*8,DIMENSION(3,3*PSUPG)::MATRIZFi1,INT133N,INT633N
REAL*8,DIMENSION(3)::MATRIZR1
REAL*8,DIMENSION(3,2)::MATRIZW2,AUX32
REAL*8,DIMENSION(2,2)::MATRIZD2
REAL*8,DIMENSION(3*NPT,3*NPT)::KGLOBAL
REAL*8,DIMENSION(3*NPT)::FGLOBAL
REAL*8,DIMENSION(1)::XB,YB

q=1.03d0
ALFAC=1.42d0

!calcular a distancia do ponto de gauss ao ponto base
Dist=DSQRT((X(i)-XT)**2+(Y(i)-YT)**2)

!avaliar a função peso(spline de 4 ordem) no ponto de gauss (função peso centrada no ponto base)
W=1.-6.*((Dist/RAIOi)**2)+8.*((Dist/RAIOi)**3)-3.*((Dist/RAIOi)**4)
DWX=-(12.D0*(X(i)-XT)/RAIOi**2)*(-1.D0+2.D0*((Dist)/RAIOi)-((Dist)/RAIOi)**2)
DWY=-(12.D0*(Y(i)-YT)/RAIOi**2)*(-1.D0+2.D0*((Dist)/RAIOi)-((Dist)/RAIOi)**2)

!montar a matriz w1(3x3) para o ponto de gauss
MATRIZW1=0.
MATRIZW1(1,1)=DWX
MATRIZW1(1,3)=DWY
MATRIZW1(2,2)=DWY
MATRIZW1(2,3)=DWX

XB(1)=XT
YB(1)=YT
CALL MONTA_BASE_MONOMIAL(MON,TMON,1,XB,YB,M,DMX,DMY,DDMX,DDMY,DDMXY)

SELECT CASE(SHAPE)
    CASE(1)

        Fi=0.D0
        dFiX=0.D0
        dFiY=0.D0
        DO K=1,PSUPG
            DO L=1,TMON
                Fi(K)=M(1,L)*AiB(L,K)+Fi(K)
                DFiX(K)=DMX(1,L)*AiB(L,K)+DFiX(K)
                DFiY(K)=DMY(1,L)*AiB(L,K)+DFiY(K)
                DFiX(K)=M(1,L)*dAiXB(L,K)+DFiX(K)
                DFiY(K)=M(1,L)*dAiyB(L,K)+DFiY(K)
                DFiX(K)=M(1,L)*AidBx(L,K)+DFiX(K)
                DFiY(K)=M(1,L)*AidBy(L,K)+DFiY(K)
            ENDDO
        ENDDO
CASE(2)

        Fi=0.D0
        dFiX=0.D0
        dFiY=0.D0
        ddFix=0.D0
        ddFiY=0.D0
        ddFixy=0.D0

        DO K=1,PSUPG
            DO L=1,TMON
                Fi(K)=M(1,L)*Sb(L,K)+Fi(K)
            ENDDO
        ENDDO        
 
    DO K=1,PSUPG
        DIST2(K)=((X(FG(K))-XT)**2+(Y(FG(K))-YT)**2)
    ENDDO
    
    DO K=1,PSUPG
        DO L=1,PSUPG
            DIST3(K,L)=((X(FG(L))-X(FG(K)))**2+(Y(FG(L))-Y(FG(K)))**2)
        ENDDO
    ENDDO
    
    dc=DSQRT(dist3(1,3))
    DO K=1,PSUPG
        R(K)=((DIST2(K))+(ALFAC*dc)**2)**q
    ENDDO

    DO K=1,PSUPG
        DO L=1,PSUPG
            Fi(K)=R(L)*Sa(L,K)+Fi(K)
        ENDDO
    ENDDO    

    DO K=1,PSUPG
        !R(1,J)=((DIST(1,J))+(ALFAC)**2)**q
        dRx(K)=2.d0*(-X(FG(K))+XT)*q*((DIST2(K))+(ALFAC*dc)**2)**(q-1.)
        dRy(K)=2.d0*(-Y(FG(K))+YT)*q*((DIST2(K))+(ALFAC*dc)**2)**(q-1.)
        ddRx(K)=(2.d0*q*(DIST2(K)+(ALFAC*dc)**2)**(q-1.))+(4.D0*q*(q-1)*(DIST2(K)+(ALFAC*dc)**2)**(q-2.))*(-X(FG(K))+XT)**2
        ddRy(K)=(2.d0*q*(DIST2(K)+(ALFAC*dc)**2)**(q-1.))+(4.D0*q*(q-1)*(DIST2(K)+(ALFAC*dc)**2)**(q-2.))*(-Y(FG(K))+YT)**2
        ddRxy(K)=4.D0*q*(q-1.D0)*(DIST2(K)+(ALFAC*dc)**2)**(q-2.)*(-X(FG(K))+XT)*(-Y(FG(K))+YT)
    ENDDO        
            
    DO K=1,PSUPG
        DO L=1,PSUPG
            dFix(K)=dRx(L)*Sa(L,K)+dFix(K)
            dFiy(K)=dRy(L)*Sa(L,K)+dFiy(K)    
    
            ddFix(K)=ddRx(L)*Sa(L,K)+ddFix(K)
            ddFiy(K)=ddRY(L)*Sa(L,K)+ddFiy(K)    
            ddFixy(K)=ddRXY(L)*Sa(L,K)+ddFiXY(K)    
        ENDDO
    ENDDO        
            
    DO K=1,PSUPG
        DO L=1,TMON
            dFix(K)=dMx(1,L)*Sb(L,K)+dFix(K)
            dFiy(K)=dMy(1,L)*Sb(L,K)+dFiy(K)
    
            ddFix(K)=ddMX(1,L)*Sb(L,K)+ddFix(K)
            ddFiy(K)=ddMY(1,L)*Sb(L,K)+ddFiy(K)
            ddFixy(K)=ddMXY(1,L)*Sb(L,K)+ddFixy(K)
        ENDDO
    ENDDO        

ENDSELECT        
        
!montar a matriz Fi1 (3x3*psupg) 
MATRIZFI1=0.D0
DO J=1,3*PSUPG,3                   
    MATRIZFI1(1,J)=dFiX((J+2)/3)
    MATRIZFI1(2,J+1)=dFiY((J+2)/3)
    MATRIZFI1(3,J)=dFiY((J+2)/3)
    MATRIZFI1(3,J+1)=dFiX((J+2)/3)
ENDDO
!a primeira e a segunda integral da eq.20 já podem ser resolvidas aqui

!montagem da matriz W2
MATRIZW2=0.D0
MATRIZW2(1,1)=W
MATRIZW2(2,2)=W
MATRIZW2(3,1)=DWX
MATRIZW2(3,2)=DWY

!montar a matriz Fi2 (2x3*psupg)
MATRIZFI2=0.D0
DO J=1,3*PSUPG,3                
    MATRIZFI2(1,J)=Fi((J+2)/3)
    MATRIZFI2(2,J+1)=Fi((J+2)/3)
    MATRIZFI2(1,J+2)=dFiX((J+2)/3)
    MATRIZFI2(2,J+2)=dFiY((J+2)/3)
ENDDO

!a sexta integral da eq.20 já pode ser resolvida aqui. A nona integral também!

!RESOLVENDO AS INTEGRAIS
! primeira integral
AUX33=0.D0
!produto da matriz W1 pela matriz D1
AUX33=MATMUL(MATRIZW1,MATRIZD1)

!produto de AUX33(MATRIZW1*MATRIZD1) pela matriz Fi1
INT133N=0.D0
INT133N=MATMUL(AUX33,MATRIZFi1)
!multiplicando o resultado da integral no ponto de gauss pelo jacobiano da transformada e pelo peso do ponto de gauss
INT133N=JAC*INT133N

!salvar o resultado da integral 1 (int133n) na matriz K global
DO J=1,PSUPG
    DO K=1,3
        !primeira linha relativa ao nó i da matriz K global
        KGLOBAL(i*3-2,FG(J)*3-3+K)=INT133N(1,J*3-3+K)+KGLOBAL(i*3-2,FG(J)*3-3+K)
        !segunda linha relativa ao nó i da matriz K global
        KGLOBAL(i*3-1,FG(J)*3-3+K)=INT133N(2,J*3-3+K)+KGLOBAL(i*3-1,FG(J)*3-3+K)        
        !terceira linha relativa ao nó i da matriz K global
        KGLOBAL(i*3,FG(J)*3-3+K)=INT133N(3,J*3-3+K)+KGLOBAL(i*3,FG(J)*3-3+K)       
    ENDDO
ENDDO

!segunda integral
FGLOBAL(3*i-2)=-MATRIZW1(1,1)*MATRIZR1(1)*JAC+FGLOBAL(3*i-2)
FGLOBAL(3*i-1)=-MATRIZW1(2,2)*MATRIZR1(2)*JAC+FGLOBAL(3*i-1)

!sexta integral
AUX32=0.D0
!produto de matriz W2 por matriz D2
AUX32=MATMUL(MATRIZW2,MATRIZD2)
!produto de AUX32(MATRIZW2*MATRIZD2) pela matriz Fi2
INT633N=0.D0
INT633N=MATMUL(AUX32,MATRIZFi2)
!multiplicando o resultado da integral no ponto de gauss pelo jacobiano da transformada e pelo peso do ponto de gauss
INT633N=JAC*INT633N
!salvar o resultado da integral 6 (int633n) na matriz K global 
DO J=1,PSUPG  !tem erro aqui
    DO K=1,3
        !primeira linha relativa ao nó i da matriz K global
        KGLOBAL(i*3-2,FG(J)*3-3+K)=INT633N(1,J*3-3+K)+KGLOBAL(i*3-2,FG(J)*3-3+K)
        !segunda linha relativa ao nó i da matriz K global
        KGLOBAL(i*3-1,FG(J)*3-3+K)=INT633N(2,J*3-3+K)+KGLOBAL(i*3-1,FG(J)*3-3+K)        
        !terceira linha relativa ao nó i da matriz K global
        KGLOBAL(i*3,FG(J)*3-3+K)=INT633N(3,J*3-3+K)+KGLOBAL(i*3,FG(J)*3-3+K)       
    ENDDO
ENDDO

!nona integral
FGLOBAL(3*i)=W*CARGA*JAC+FGLOBAL(3*i)

RETURN
END SUBROUTINE INTEGR_DOMINIO

    