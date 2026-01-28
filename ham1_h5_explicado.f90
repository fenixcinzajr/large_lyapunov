program kutta

!****************************************************************************


!################################################################
!								                              #
!	INTEGRACAO DE E.D.O. VIA RUNGE - KUTTA DE ORDEM 4 .	      #
!								                              #
!	CALCULO DA SECAO DE POINCARE INTEGRANDO NO TEMPO .	      #
!
!   Essencialmente são geradas duas saídas: t102.dat que contém
!os pontos com energia maior do que zero e t103.dat que contém
!os pontos com energia menor do que zero.
!
!     								                          #
!################################################################
!                           
      Implicit real*8(A-H,P-Z)

      EXTERNAL dx,dp
!     Runge-Kutta de ordem 6 com integração no tempo "t" e gerando espaço de fase (p,q)
      OPEN(UNIT=1,FILE='t1.dat')
      OPEN(UNIT=2,FILE='t2.dat')
      OPEN(UNIT=3,FILE='t3.dat')
      OPEN(UNIT=4,FILE='t4.dat')
      OPEN(UNIT=5,FILE='t5.dat')
      OPEN(UNIT=6,FILE='t6.dat')
      OPEN(UNIT=7,FILE='t7.dat')
      OPEN(UNIT=8,FILE='t8.dat')
      OPEN(UNIT=101,FILE='t101.dat')
      OPEN(UNIT=102,FILE='t102.dat')
      OPEN(UNIT=103,FILE='t103.dat')
      OPEN(UNIT=104,FILE='t104.dat')
      OPEN(UNIT=105,FILE='t105.dat')
      
      q=1.d0       !parâmetro qsi
      a=1.d0       !parâmetro eta
      xe=1.d0      !parâmetro xe
      e0=2.5d0     !parâmetro epsilon0
      w=sqrt(3.d0) !parâmetro da frequência ômega
      
      z1=e0
      z2=w
      z3=a
      z4=q
      z5=xe
      z6=0.d0
      
pi=4.D0*DATAN(1.D0)
doispi=pi+pi
pimenos=-pi
nh1=10000        !defino nh1, um número de ajuste arbitrário que será inversamente proporcional ao passo - se quiser diminuir o passo, aumenta-se nh1
h=2.d0*pi/(dfloat(nh1)*w) !passo h é definido tal que w.t será múltiplo inteiro de 2.pi, portanto depende da frequência w.
tfinal=1000.d0   !tempo final de integração
NITER=tfinal/h   !número de iterações necessárias para que se atinja o tempo final
div=NITER/nh1    !número de divisões contidas no tempo final (é o número de cruzamentos com a seção de Poincaré para o tempo dado)
ndiv=int(div)    !tornando ndiv um inteiro. Isso permite integrar o sistema até o último cruzamento com a seção de Poincaré 
nmax=ndiv*nh1    !número de iterações necessárias até o último cruzamento com a seção de Poincaré
titeracao=nmax*h !novo tempo de integração (menor do que o tempo final) tal que o sistema seja integrado até o último cruzamento com a s.P.
coletas=dfloat(nmax)/dfloat(nh1) !quantas vezes o sistema cruza com a seção de Poincaré no tempo dado
tempocalculado=nmax*h !novo tempo de integração
x=0.d0



        write(*,*)"passo =",h
        write(*,*)"iteracoes =",nmax
        write(*,*)"ndiv =",ndiv
        write(*,*)"xi =",xi
        write(*,*)"xf =",xf
        write(*,*)"tempocalculado =",tempocalculado
        write(*,*)"quantos threads serao usados?"
        read(*,*)nt


!---------------------
!---------------------
ix=21
xi=-0.9d0
xf=0.1d0
px=(xf-xi)/dfloat(ix)
!---------------------
ip=21
psi=0.2d0
psf=1.2d0
pp=(psf-psi)/dfloat(ip)


!$OMP Parallel Do default(private) shared(psi,ip,pp,xi,ix,px,nmax,niter,h,nh1,w,nma,nme,z1,z2,z3,z4,z5,z6)&
!$omp Num_threads(nt) schedule(dynamic)
do j=1,ix
    write(*,*)j,"de",ix
    x=xi+dfloat(j)*px
    x0=x
             
    do k=1,ip
        x=x0
        p=psi+dfloat(k)*pp
        p0=p
        t=0.d0
        do i1=1,8
            DO i=1,nmax
                call rungekutta(t,x,p,i,h,z1,z2,z3,z4,z5,z6)
                !quando i=nh1, entao t=i*(h)=nh1*(h)=nh1*(2*pi/(nh1*w))=2*pi/w Assim w.t é multiplo de 2.pi
                iresto=mod(i,nh1)
                if (iresto.eq.0.d0) then
                    !condição de escape/dissociação:
                    if (x.gt.20.d0.and.p.gt.0.00006d0) then
                        write(1,*)x0,p0
                        goto 103
                    end if
                    !coletando apenas quando i1>2 para que i1=1 seja considerado como um transiente (desprezo as coletas iniciais até "titeracao"):
                    if (i1.ge.2) then
                        ener=(exp(-2*x) - 2.d0*exp(-x))/2.d0 + p**2/2.d0
                        if (ener.gt.0.d0) then
                            !$omp critical
                            write(102,*)x,p
                            !$omp end critical
                        else
                            !$omp critical
                            write(103,*)x,p
                            !$omp end critical
                        end if
                    end if
                end if
            END do
        end do
        !$omp critical
        write(2,*)x0,p0
        write(3,*)x,p
        !$omp end critical
103 continue
    end do
end do 
!$OMP END PARALLEL DO

end program kutta

!############################################################################      


!############################################################################      
!	DERIVADA DE x : dx/dt = del(H)/del(p)

      REAL*8 FUNCTION dx(T,x,p,z1,z2,z3,z4,z5,z6)
      IMPLICIT REAL*8(A-H,P-Z)
      
      dx = p
      
      RETURN
      END
      
!############################################################################      
!	DERIVADA DE p : dp/dt = -del(H)/del(x)

      REAL*8 FUNCTION dp(T,x,p,z1,z2,z3,z4,z5,z6)
      IMPLICIT REAL*8(A-H,P-Z)
      
      q=1.d0
      a=1.d0
      xe=1.d0
      e0=2.5d0
      w=sqrt(3.d0)
      !z1=e0
      !z2=w
      !z3=a
      !z4=q
      !z5=xe
      dp=dexp(-2.d0*x)-dexp(-x)+e0*dsin(w*t)*(a*dcos(a*(x+xe))*dexp(-q*(x+xe)**4)-dsin(a*(x+xe))*4.d0*q*((x+xe)**3)*dexp(-q*(x+xe)**4))/a
      

      RETURN
      END
      

      
subroutine rungekutta(t,x,p,i,h,z1,z2,z3,z4,z5,z6)
implicit real*8 (A-H,P-Z)

    ak1=dx(T,x,p,z1,z2,z3,z4,z5,z6)
    al1=dp(T,x,p,z1,z2,z3,z4,z5,z6)
    
    aux1=t+h/2.d0
    aux2=x+(h/2.d0)*ak1
    aux3=p+(h/2.d0)*al1
    
    ak2=dx(aux1,aux2,aux3,z1,z2,z3,z4,z5,z6)
    al2=dp(aux1,aux2,aux3,z1,z2,z3,z4,z5,z6)
    
    aux4=x+(h/2.d0)*ak2
    aux5=p+(h/2.d0)*al2
    
    ak3=dx(aux1,aux4,aux5,z1,z2,z3,z4,z5,z6)
    al3=dp(aux1,aux4,aux5,z1,z2,z3,z4,z5,z6)
    
    aux6=t+h
    aux7=x+h*ak3
    aux8=p+h*al3
    
    ak4=dx(aux6,aux7,aux8,z1,z2,z3,z4,z5,z6)
    al4=dp(aux6,aux7,aux8,z1,z2,z3,z4,z5,z6)
    
    x=x+(ak1+2.d0*ak2+2.d0*ak3+ak4)*h/6.d0
    p=p+(al1+2.d0*al2+2.d0*al3+al4)*h/6.d0
    t=dfloat(i)*h

return
end
    
