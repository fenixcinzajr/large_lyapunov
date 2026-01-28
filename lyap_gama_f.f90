
!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  calcular o expoente de Lyapunov para sistemas contínuos (EDOs).
!
!  Autor: Allan Kardec de Almeida Junior - fenixcinza@yahoo.com.br
!
!****************************************************************************

    program ODE
Implicit real*8(A-H,P-Z)
    parameter n=2
    parameter nn=6
    external fcn,dx,dp
    
    dimension y(nn),znorm(n),gsc(n),cum(n),c(24),w(NN,9)
    
    
    
z1=2.5d0
z2=sqrt(3.d0)
z3=1.d0
z4=1.d0
z5=1.d0
z6=0.3d0   
    
wfreq=sqrt(3.d0)
pi=4.D0*DATAN(1.D0)
nh1=10000
h=2.d0*pi/(dfloat(nh1)*wfreq)
tfinal=1000.d0
NITER=tfinal/h
div=NITER/nh1
ndiv=int(div)
nmax=ndiv*nh1
titeracao=nmax*h
coletas=dfloat(nmax)/dfloat(nh1)
transiente=dfloat(nmax)/2.d0
ntrans=int(transiente)
tempocalculado=nmax*h

write(*,*)"passo =",h
write(*,*)"iteracoes =",nmax
write(*,*)"ndiv =",ndiv
write(*,*)"quantos threads serao usados?"
read(*,*)nt
nt=7
    
    OPEN(UNIT=100,FILE='t100.dat')
    OPEN(UNIT=101,FILE='t101.dat')
    OPEN(UNIT=102,FILE='t102.dat')
    OPEN(UNIT=103,FILE='exp_lyapunov_gama.dat')
    
iz=8
zi=0.1d0
zf=0.15d0
pz=(zf-zi)/float(iz)

!$OMP Parallel Do default(private) shared(iz,zi,pz,nmax,h,nh1,wfreq,z1,z2,z3,z4,z5,z6)&
!$omp Num_threads(nt) schedule(dynamic)
do l1001=0,iz
zz=zi+float(l1001)*pz

!transiente:
    x=-0.4d0
    p=0.7d0
    t=0.d0

do j=1,100
do i=1,nmax
    call rungekutta(t,x,p,i,h,z1,z2,z3,z4,z5,zz)
end do
end do
x0=x
p0=p


!acabou o transiente


    y(1)=x0
    y(2)=p0

    
    neq=nn
    x=0.0
    ind=1
    nw=nn
    tol=1.d-4
    nstep=nmax
    stpsze=h
    i0=nh1

    do 10 i=n+1,nn
        y(i)=0.0
10  continue
    do 20 i=1,n
        y((n+1)*i)=1.0
        cum(i)=0.0
20  continue
    
    !type*, 'tol,nstep,stpsze,i0?'
    !accept*,tol,nstep,stpsze,i0

    neq=nn
    x=0.0
    ind=1
    nw=nn
    tol=1.d-4
    nstep=nmax
    stpsze=h
    i0=nh1

    do 100 i=1,nstep
        xend=stpsze*float(i)
        
        call dverk (neq, fcn, x, y, xend, tol, ind, c, nw, w, zz)
        
        znorm(1)=0.0
        do 30 j=1,n
            znorm(1)=znorm(1)+y(n*j+1)**2
30      continue
        znorm(1)=sqrt(znorm(1))
        do 40 j=1,n
            
            y(n*j+1)=y(n*j+1)/znorm(1)
40      continue
        
        
        do 80 j=2,n
            
            do 50 k=1,(j-1)
                gsc(k)=0.0
                do 50 L=1,n
                    gsc(k)=gsc(k)+y(n*L+j)*y(n*L+k)
50          continue
            
            
            do 60 k=1,n
                do 60 L=1,(j-1)
                    y(n*k+j)=y(n*k+j)-gsc(L)*y(n*k+L)
60          continue
            
            
            znorm(j)=0.0
            do 70 k=1,n
                znorm(j)=znorm(j)+y(n*k+j)**2
70          continue
            znorm(j)=sqrt(znorm(j))
            
            
            do 80 k=1,n
                y(n*k+j)=y(n*k+j)/znorm(j)
80      continue
        
        
        do 90 k=1,n
            cum(k)=cum(k)+log(znorm(k))/log(2.)
90      continue
        
        !if ((mod(i,i0).eq.0).and.(i/nh1.gt.265)) type*,y(1),y(2),(cum(k)/x,k=1,n)
        if (mod(i,i0).eq.0) then
            ener=(exp(-2*y(1)) - 2.d0*exp(-y(1)))/2.d0 + y(2)**2/2.d0
            !$omp critical
                write(100,*),zz,ener
                write(101,*),zz,y(1)
                write(102,*),zz,y(2)
            !$omp end critical
            if (i.eq.nmax) then
                !$omp critical
                if (cum(1).ge.cum(2)) then
                    explyap=cum(1)/x
                else
                    explyap=cum(2)/x
                end if
                write(103,*),zz,explyap
                type*,zz,explyap
                !$omp end critical
            end if
        end if
        
100 continue 
end do
!$OMP END PARALLEL DO
    pause
    call exit
  
    end
    

    
    subroutine fcn (n,x,y,yprime,zz)
    IMPLICIT REAL*8(A-H,P-Z)
    dimension y(6),yprime(6)
    !equações de movimento:
    
      q=1.d0
      a=1.d0
      xe=1.d0
      e0=2.5d0
      wfreq=sqrt(3.d0)
      z6=0.3d0

    
yprime(1)=y(2)
yprime(2)=exp(-2.d0*y(1))-exp(-y(1))+e0*sin(wfreq*x)*(a*cos(a*(y(1)+xe))*exp(-q*(y(1)+xe)**4)-sin(a*(y(1)+xe))*4.d0*q*((y(1)+xe)**3)*exp(-q*(y(1)+xe)**4))/a-zz*y(2)
    
do 10 i=0,1
    !equações de movimento linearizadas:
        

yprime(3+i)=y(5+i)
        
!(zz*Sin(x*wfreq)*(a*Cos(a*(y(1) + xe)) - 4.d0*q*(y(1) + xe)**3*Sin(a*(y(1) + xe))))/(a*Exp(q*(y(1) + xe)**4)) + &
yprime(5+i)=     (-Exp(-y(1)) + (2.d0*(-1.d0 + Exp(y(1))))/Exp(2.d0*y(1)) + (e0*Sin(x*wfreq)*&
            (-4.d0*a*q*(y(1) + xe)**3*Cos(a*(y(1) + xe)) - a**2*Sin(a*(y(1) + xe)) - 12.d0*q*(y(1) + xe)**2*Sin(a*(y(1) + xe))))/(a*Exp(q*(y(1) + xe)**4)) - &
        (4.d0*e0*q*(y(1) + xe)**3*Sin(x*wfreq)*(a*Cos(a*(y(1) + xe)) - 4.d0*q*(y(1) + xe)**3*Sin(a*(y(1) + xe))))/(a*Exp(q*(y(1) + xe)**4)))*y(3+i) - zz*y(5+i)
10  continue
    
    return
    end
    


      
      
      
      
      subroutine dverk (n, fcn, x, y, xend, tol, ind, c, nw, w, zz)
      integer n, ind, nw, k
      double precision x, y(n), xend, tol, c(1), w(nw,9), temp
!c
!c***********************************************************************
!c                                                                      *
!c note added 11/14/85.                                                 *
!c                                                                      *
!c if you discover any errors in this subroutine, please contact        *
!c                                                                      *
!c        kenneth r. jackson                                            *
!c        department of computer science                                *
!c        university of toronto                                         *
!c        toronto, ontario,                                             *
!c        canada   m5s 1a4                                              *
!c                                                                      *
!c        phone: 416-978-7075                                           *
!c                                                                      *
!c        electronic mail:                                              *
!c        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
!c        csnet:  krj@toronto                                           *
!c        arpa:   krj.toronto@csnet-relay                               *
!c        bitnet: krj%toronto@csnet-relay.arpa                          *
!c                                                                      *
!c dverk is written in fortran 66.                                      *
!c                                                                      *
!c the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
!c set for a  vax  in  double  precision.  they  should  be  reset,  as *
!c described below, if this program is run on another machine.          *
!c                                                                      *
!c the c array is declared in this subroutine to have one element only, *
!c although  more  elements  are  referenced  in this subroutine.  this *
!c causes some compilers to issue warning messages.  there is,  though, *
!c no  error  provided  c is declared sufficiently large in the calling *
!c program, as described below.                                         *
!c                                                                      *
!c the following external statement  for  fcn  was  added  to  avoid  a *
!c warning  message  from  the  unix  f77 compiler.  the original dverk *
!c comments and code follow it.                                         *
!c                                                                      *
!c***********************************************************************
!c
      external fcn
!c
!c***********************************************************************
!c                                                                      *
!c     purpose - this is a runge-kutta  subroutine  based  on  verner's *
!c fifth and sixth order pair of formulas for finding approximations to *
!c the solution of  a  system  of  first  order  ordinary  differential *
!c equations  with  initial  conditions. it attempts to keep the global *
!c error proportional to  a  tolerance  specified  by  the  user.  (the *
!c proportionality  depends  on the kind of error control that is used, *
!c as well as the differential equation and the range of integration.)  *
!c                                                                      *
!c     various options are available to the user,  including  different *
!c kinds  of  error control, restrictions on step sizes, and interrupts *
!c which permit the user to examine the state of the  calculation  (and *
!c perhaps make modifications) during intermediate stages.              *
!c                                                                      *
!c     the program is efficient for non-stiff systems.  however, a good *
!c variable-order-adams  method  will probably be more efficient if the *
!c function evaluations are very costly.  such a method would  also  be *
!c more suitable if one wanted to obtain a large number of intermediate *
!c solution values by interpolation, as might be the case  for  example *
!c with graphical output.                                               *
!c                                                                      *
!c                                    hull-enright-jackson   1/10/76    *
!c                                                                      *
!c***********************************************************************
!c                                                                      *
!c     use - the user must specify each of the following                *
!c                                                                      *
!c     n  number of equations                                           *
!c                                                                      *
!c   fcn  name of subroutine for evaluating functions - the  subroutine *
!c           itself must also be provided by the user - it should be of *
!c           the following form                                         *
!c              subroutine fcn(n, x, y, yprime)                         *
!c              integer n                                               *
!c              double precision x, y(n), yprime(n)                     *
!c                      *** etc ***                                     *
!c           and it should evaluate yprime, given n, x and y            *
!c                                                                      *
!c     x  independent variable - initial value supplied by user         *
!c                                                                      *
!c     y  dependent variable - initial values of components y(1), y(2), *
!c           ..., y(n) supplied by user                                 *
!c                                                                      *
!c  xend  value of x to which integration is to be carried out - it may *
!c           be less than the initial value of x                        *
!c                                                                      *
!c   tol  tolerance - the subroutine attempts to control a norm of  the *
!c           local  error  in  such  a  way  that  the  global error is *
!c           proportional to tol. in some problems there will be enough *
!c           damping  of  errors, as well as some cancellation, so that *
!c           the global error will be less than tol. alternatively, the *
!c           control   can   be  viewed  as  attempting  to  provide  a *
!c           calculated value of y at xend which is the exact  solution *
!c           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
!c           is proportional to tol.  (the norm  is  a  max  norm  with *
!c           weights  that  depend on the error control strategy chosen *
!c           by the user.  the default weight for the k-th component is *
!c           1/max(1,abs(y(k))),  which therefore provides a mixture of *
!c           absolute and relative error control.)                      *
!c                                                                      *
!c   ind  indicator - on initial entry ind must be set equal to  either *
!c           1  or  2. if the user does not wish to use any options, he *
!c           should set ind to 1 - all that remains for the user to  do *
!c           then  is  to  declare c and w, and to specify nw. the user *
!c           may also  select  various  options  on  initial  entry  by *
!c           setting ind = 2 and initializing the first 9 components of *
!c           c as described in the next section.  he may also  re-enter *
!c           the  subroutine  with ind = 3 as mentioned again below. in *
!c           any event, the subroutine returns with ind equal to        *
!c              3 after a normal return                                 *
!c              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
!c              -1, -2, or -3 after an error condition (see below)      *
!c                                                                      *
!c     c  communications vector - the dimension must be greater than or *
!c           equal to 24, unless option c(1) = 4 or 5 is used, in which *
!c           case the dimension must be greater than or equal to n+30   *
!c                                                                      *
!c    nw  first dimension of workspace w -  must  be  greater  than  or *
!c           equal to n                                                 *
!c                                                                      *
!c     w  workspace matrix - first dimension must be nw and second must *
!c           be greater than or equal to 9                              *
!c                                                                      *
!c     the subroutine  will  normally  return  with  ind  =  3,  having *
!c replaced the initial values of x and y with, respectively, the value *
!c of xend and an approximation to y at xend.  the  subroutine  can  be *
!c called  repeatedly  with new values of xend without having to change *
!c any other argument.  however, changes in tol, or any of the  options *
!c described below, may also be made on such a re-entry if desired.     *
!c                                                                      *
!c     three error returns are also possible, in which  case  x  and  y *
!c will be the most recently accepted values -                          *
!c     with ind = -3 the subroutine was unable  to  satisfy  the  error *
!c        requirement  with a particular step-size that is less than or *
!c        equal to hmin, which may mean that tol is too small           *
!c     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
!c        probably  means  that the requested tol (which is used in the *
!c        calculation of hmin) is too small                             *
!c     with ind = -1 the allowed maximum number of fcn evaluations  has *
!c        been  exceeded,  but  this  can only occur if option c(7), as *
!c        described in the next section, has been used                  *
!c                                                                      *
!c     there are several circumstances that will cause the calculations *
!c to  be  terminated,  along with output of information that will help *
!c the user determine the cause of  the  trouble.  these  circumstances *
!c involve  entry with illegal or inconsistent values of the arguments, *
!c such as attempting a normal  re-entry  without  first  changing  the *
!c value of xend, or attempting to re-enter with ind less than zero.    *
!c                                                                      *
!c***********************************************************************
!c                                                                      *
!c     options - if the subroutine is entered with ind = 1, the first 9 *
!c components of the communications vector are initialized to zero, and *
!c the subroutine uses only default values  for  each  option.  if  the *
!c subroutine  is  entered  with ind = 2, the user must specify each of *
!c these 9 components - normally he would first set them all  to  zero, *
!c and  then  make  non-zero  those  that  correspond to the particular *
!c options he wishes to select. in any event, options may be changed on *
!c re-entry  to  the  subroutine  -  but if the user changes any of the *
!c options, or tol, in the course of a calculation he should be careful *
!c about  how  such changes affect the subroutine - it may be better to *
!c restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
!c program  -  the information is available to the user, but should not *
!c normally be changed by him.)                                         *
!c                                                                      *
!c  c(1)  error control indicator - the norm of the local error is  the *
!c           max  norm  of  the  weighted  error  estimate  vector, the *
!c           weights being determined according to the value of c(1) -  *
!c              if c(1)=1 the weights are 1 (absolute error control)    *
!c              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
!c                 control)                                             *
!c              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
!c                 (relative  error  control,  unless abs(y(k)) is less *
!c                 than the floor value, abs(c(2)) )                    *
!c              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
!c                 (here individual floor values are used)              *
!c              if c(1)=5 the weights are 1/abs(c(k+30))                *
!c              for all other values of c(1), including  c(1) = 0,  the *
!c                 default  values  of  the  weights  are  taken  to be *
!c                 1/max(1,abs(y(k))), as mentioned earlier             *
!c           (in the two cases c(1) = 4 or 5 the user must declare  the *
!c           dimension of c to be at least n+30 and must initialize the *
!c           components c(31), c(32), ..., c(n+30).)                    *
!c                                                                      *
!c  c(2)  floor value - used when the indicator c(1) has the value 3    *
!c                                                                      *
!c  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
!c           to be abs(c(3)) - otherwise it uses the default value      *
!c              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
!c           where dwarf is a very small positive  machine  number  and *
!c           rreb is the relative roundoff error bound                  *
!c                                                                      *
!c  c(4)  hstart specification - if not zero, the subroutine  will  use *
!c           an  initial  hmag equal to abs(c(4)), except of course for *
!c           the restrictions imposed by hmin and hmax  -  otherwise it *
!c           uses the default value of hmax*(tol)**(1/6)                *
!c                                                                      *
!c  c(5)  scale specification - this is intended to be a measure of the *
!c           scale of the problem - larger values of scale tend to make *
!c           the method more reliable, first  by  possibly  restricting *
!c           hmax  (as  described  below) and second, by tightening the *
!c           acceptance requirement - if c(5) is zero, a default  value *
!c           of  1  is  used.  for  linear  homogeneous  problems  with *
!c           constant coefficients, an appropriate value for scale is a *
!c           norm  of  the  associated  matrix.  for other problems, an *
!c           approximation to  an  average  value  of  a  norm  of  the *
!c           jacobian along the trajectory may be appropriate           *
!c                                                                      *
!c  c(6)  hmax specification - four cases are possible                  *
!c           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
!c              min(abs(c(6)),2/abs(c(5)))                              *
!c           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
!c           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
!c              2/abs(c(5))                                             *
!c           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
!c              of 2                                                    *
!c                                                                      *
!c  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
!c           error  return with ind = -1 will be caused when the number *
!c           of function evaluations exceeds abs(c(7))                  *
!c                                                                      *
!c  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
!c           interrupt   the  calculations  after  it  has  chosen  its *
!c           preliminary value of hmag, and just before choosing htrial *
!c           and  xtrial  in  preparation for taking a step (htrial may *
!c           differ from hmag in sign, and may  require  adjustment  if *
!c           xend  is  near) - the subroutine returns with ind = 4, and *
!c           will resume calculation at the point  of  interruption  if *
!c           re-entered with ind = 4                                    *
!c                                                                      *
!c  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
!c           interrupt   the  calculations  immediately  after  it  has *
!c           decided whether or not to accept the result  of  the  most *
!c           recent  trial step, with ind = 5 if it plans to accept, or *
!c           ind = 6 if it plans to reject -  y(*)  is  the  previously *
!c           accepted  result, while w(*,9) is the newly computed trial *
!c           value, and w(*,2) is the unweighted error estimate vector. *
!c           the  subroutine  will  resume calculations at the point of *
!c           interruption on re-entry with ind = 5 or 6. (the user  may *
!c           change ind in this case if he wishes, for example to force *
!c           acceptance of a step that would otherwise be rejected,  or *
!c           vice versa. he can also restart with ind = 1 or 2.)        *
!c                                                                      *
!c***********************************************************************
!c                                                                      *
!c  summary of the components of the communications vector              *
!c                                                                      *
!c     prescribed at the option       determined by the program         *
!c           of the user                                                *
!c                                                                      *
!c                                    c(10) rreb(rel roundoff err bnd)  *
!c     c(1) error control indicator   c(11) dwarf (very small mach no)  *
!c     c(2) floor value               c(12) weighted norm y             *
!c     c(3) hmin specification        c(13) hmin                        *
!c     c(4) hstart specification      c(14) hmag                        *
!c     c(5) scale specification       c(15) scale                       *
!c     c(6) hmax specification        c(16) hmax                        *
!c     c(7) max no of fcn evals       c(17) xtrial                      *
!c     c(8) interrupt no 1            c(18) htrial                      *
!c     c(9) interrupt no 2            c(19) est                         *
!c                                    c(20) previous xend               *
!c                                    c(21) flag for xend               *
!c                                    c(22) no of successful steps      *
!c                                    c(23) no of successive failures   *
!c                                    c(24) no of fcn evals             *
!c                                                                      *
!c  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
!c                                                                      *
!c***********************************************************************
!c                                                                      *
!c  an overview of the program                                          *
!c                                                                      *
!c     begin initialization, parameter checking, interrupt re-entries   *
!c  ......abort if ind out of range 1 to 6                              *
!c  .     cases - initial entry, normal re-entry, interrupt re-entries  *
!c  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
!c  v........abort if n.gt.nw or tol.le.0                               *
!c  .        if initial entry without options (ind .eq. 1)              *
!c  .           set c(1) to c(9) equal to zero                          *
!c  .        else initial entry with options (ind .eq. 2)               *
!c  .           make c(1) to c(9) non-negative                          *
!c  .           make floor values non-negative if they are to be used   *
!c  .        end if                                                     *
!c  .        initialize rreb, dwarf, prev xend, flag, counts            *
!c  .     case 2 - normal re-entry (ind .eq. 3)                         *
!c  .........abort if xend reached, and either x changed or xend not    *
!c  .        re-initialize flag                                         *
!c  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
!c  v        transfer control to the appropriate re-entry point.......  *
!c  .     end cases                                                  .  *
!c  .  end initialization, etc.                                      .  *
!c  .                                                                v  *
!c  .  loop through the following 4 stages, once for each trial step .  *
!c  .     stage 1 - prepare                                          .  *
!c***********error return (with ind=-1) if no of fcn evals too great .  *
!c  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
!c  .        calc hmin, scale, hmax                                  .  *
!c***********error return (with ind=-2) if hmin .gt. hmax            .  *
!c  .        calc preliminary hmag                                   .  *
!c***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
!c  .        calc hmag, xtrial and htrial                            .  *
!c  .     end stage 1                                                .  *
!c  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
!c  .     stage 3 - calc the error estimate                          .  *
!c  .     stage 4 - make decisions                                   .  *
!c  .        set ind=5 if step acceptable, else set ind=6            .  *
!c***********interrupt no 2 if requested....................re-entry.v  *
!c  .        if step accepted (ind .eq. 5)                              *
!c  .           update x, y from xtrial, ytrial                         *
!c  .           add 1 to no of successful steps                         *
!c  .           set no of successive failures to zero                   *
!c**************return(with ind=3, xend saved, flag set) if x .eq. xend *
!c  .        else step not accepted (ind .eq. 6)                        *
!c  .           add 1 to no of successive failures                      *
!c**************error return (with ind=-3) if hmag .le. hmin            *
!c  .        end if                                                     *
!c  .     end stage 4                                                   *
!c  .  end loop                                                         *
!c  .                                                                   *
!c  begin abort action                                                  *
!c     output appropriate  message  about  stopping  the  calculations, *
!c        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
!c        previous xend,  no of  successful  steps,  no  of  successive *
!c        failures, no of fcn evals, and the components of y            *
!c     stop                                                             *
!c  end abort action                                                    *
!c                                                                      *
!c***********************************************************************
!c
!c     ******************************************************************
!c     * begin initialization, parameter checking, interrupt re-entries *
!c     ******************************************************************
!c
!c  ......abort if ind out of range 1 to 6
         if (ind.lt.1 .or. ind.gt.6) go to 500
!c
!c        cases - initial entry, normal re-entry, interrupt re-entries
         go to (5, 5, 45, 1111, 2222, 2222), ind
!c        case 1 - initial entry (ind .eq. 1 or 2)
!c  .........abort if n.gt.nw or tol.le.0
    5       if (n.gt.nw .or. tol.le.0.d0) go to 500
            if (ind.eq. 2) go to 15
!c              initial entry without options (ind .eq. 1)
!c              set c(1) to c(9) equal to 0
               do 10 k = 1, 9
                  c(k) = 0.d0
   10          continue
               go to 35
   15       continue
!c              initial entry with options (ind .eq. 2)
!c              make c(1) to c(9) non-negative
               do 20 k = 1, 9
                  c(k) = dabs(c(k))
   20          continue
!c              make floor values non-negative if they are to be used
               if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
                  do 25 k = 1, n
                     c(k+30) = dabs(c(k+30))
   25             continue
   30          continue
   35       continue
!c           initialize rreb, dwarf, prev xend, flag, counts
            c(10) = 2.d0**(-56)
            c(11) = 1.d-35
!c           set previous xend initially to initial value of x
            c(20) = x
            do 40 k = 21, 24
               c(k) = 0.d0
   40       continue
            go to 50
!c        case 2 - normal re-entry (ind .eq. 3)
!c  .........abort if xend reached, and either x changed or xend not
   45       if (c(21).ne.0.d0 .and.&
                             (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!c           re-initialize flag
            c(21) = 0.d0
            go to 50
!c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!c           transfer control to the appropriate re-entry point..........
!c           this has already been handled by the computed go to        .
!c        end cases                                                     v
   50    continue
!c
!c     end initialization, etc.
!c
!c     ******************************************************************
!c     * loop through the following 4 stages, once for each trial  step *
!c     * until the occurrence of one of the following                   *
!c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!c     *        stage 4                                                 *
!c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!c     *        requested, in stage 1 or stage 4                        *
!c     ******************************************************************
!c
99999 continue
!c
!c        ***************************************************************
!c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!c        * and some parameter  checking,  and  end  up  with  suitable *
!c        * values of hmag, xtrial and htrial in preparation for taking *
!c        * an integration step.                                        *
!c        ***************************************************************
!c
!c***********error return (with ind=-1) if no of fcn evals too great
            if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
               ind = -1
               return
  100       continue
!c
!c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
            if (ind .eq. 6) go to 105
               call fcn(n, x, y, w(1,1),zz)
               c(24) = c(24) + 1.d0
  105       continue
!c
!c           calculate hmin - use default unless value prescribed
            c(13) = c(3)
            if (c(3) .ne. 0.d0) go to 165
!c              calculate default value of hmin
!c              first calculate weighted norm y - c(12) - as specified
!c              by the error control indicator c(1)
               temp = 0.d0
               if (c(1) .ne. 1.d0) go to 115
!c                 absolute error control - weights are 1
                  do 110 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  110             continue
                  c(12) = temp
                  go to 160
  115          if (c(1) .ne. 2.d0) go to 120
!c                 relative error control - weights are 1/dabs(y(k)) so
!c                 weighted norm y is 1
                  c(12) = 1.d0
                  go to 160
  120          if (c(1) .ne. 3.d0) go to 130
!c                 weights are 1/max(c(2),abs(y(k)))
                  do 125 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(2))
  125             continue
                  c(12) = dmin1(temp, 1.d0)
                  go to 160
  130          if (c(1) .ne. 4.d0) go to 140
!c                 weights are 1/max(c(k+30),abs(y(k)))
                  do 135 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  135             continue
                  c(12) = dmin1(temp, 1.d0)
                  go to 160
  140          if (c(1) .ne. 5.d0) go to 150
!c                 weights are 1/c(k+30)
                  do 145 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  145             continue
                  c(12) = temp
                  go to 160
  150          continue
!c                 default case - weights are 1/max(1,abs(y(k)))
                  do 155 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  155             continue
                  c(12) = dmin1(temp, 1.d0)
  160          continue
               c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
  165       continue
!c
!c           calculate scale - use default unless value prescribed
            c(15) = c(5)
            if (c(5) .eq. 0.d0) c(15) = 1.d0
!c
!c           calculate hmax - consider 4 cases
!c           case 1 both hmax and scale prescribed
               if (c(6).ne.0.d0 .and. c(5).ne.0.d0)&
                                         c(16) = dmin1(c(6), 2.d0/c(5))
!c           case 2 - hmax prescribed, but scale not
               if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
!c           case 3 - hmax not prescribed, but scale is
               if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
!c           case 4 - neither hmax nor scale is provided
               if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
!c
!c***********error return (with ind=-2) if hmin .gt. hmax
            if (c(13) .le. c(16)) go to 170
               ind = -2
               return
  170       continue
!c
!c           calculate preliminary hmag - consider 3 cases
            if (ind .gt. 2) go to 175
!c           case 1 - initial entry - use prescribed value of hstart, if
!c              any, else default
               c(14) = c(4)
               if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
               go to 185
  175       if (c(23) .gt. 1.d0) go to 180
!c           case 2 - after a successful step, or at most  one  failure,
!c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!c              overflow. then avoid reduction by more than half.
               temp = 2.d0*c(14)
               if (tol .lt. (2.d0/.9d0)**6*c(19))&
                                 temp = .9d0*(tol/c(19))**(1./6.)*c(14)
               c(14) = dmax1(temp, .5d0*c(14))
               go to 185
  180       continue
!c           case 3 - after two or more successive failures
               c(14) = .5d0*c(14)
  185       continue
!c
!c           check against hmax
            c(14) = dmin1(c(14), c(16))
!c
!c           check against hmin
            c(14) = dmax1(c(14), c(13))
!c
!c***********interrupt no 1 (with ind=4) if requested
            if (c(8) .eq. 0.d0) go to 1111
               ind = 4
               return
!c           resume here on re-entry with ind .eq. 4   ........re-entry..
 1111       continue
!c
!c           calculate hmag, xtrial - depending on preliminary hmag, xend
            if (c(14) .ge. dabs(xend - x)) go to 190
!c              do not step more than half way to xend
               c(14) = dmin1(c(14), .5d0*dabs(xend - x))
               c(17) = x + dsign(c(14), xend - x)
               go to 195
  190       continue
!c              hit xend exactly
               c(14) = dabs(xend - x)
               c(17) = xend
  195       continue
!c
!c           calculate htrial
            c(18) = c(17) - x
!c
!c        end stage 1
!c
!c        ***************************************************************
!c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!c        * stage 3. w(*,9) is temporary storage until finally it holds *
!c        * ytrial.                                                     *
!c        ***************************************************************
!c
            temp = c(18)/1398169080000.d0
!c
            do 200 k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
  200       continue
            call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2),zz)
!c
            do 205 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0&
                                + w(k,2)*298276070400.d0  )
  205       continue
            call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3),zz)
!c
            do 210 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0&
                                     - w(k,2)*3728450880000.d0&
                                     + w(k,3)*3495422700000.d0 )
  210       continue
            call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4),zz)
!c
            do 215 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0&
                                     + w(k,2)*12816549900000.d0&
                                     - w(k,3)*9284716546875.d0&
                                     + w(k,4)*1237962206250.d0 )
  215       continue
            call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5),zz)
!c
            do 220 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0&
                                     - w(k,2)*11185352640000.d0&
                                     + w(k,3)*9172628850000.d0&
                                     - w(k,4)*427218330000.d0&
                                     + w(k,5)*482505408000.d0  )
  220       continue
            call fcn(n, x + c(18), w(1,9), w(1,6),zz)
!c
            do 225 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0&
                                     + w(k,2)*2311639545600.d0&
                                     - w(k,3)*1322092233000.d0&
                                     - w(k,4)*453006781920.d0&
                                     + w(k,5)*326875481856.d0  )
  225       continue
            call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7),zz)
!c
            do 230 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0&
                                     - w(k,2)*9754668000000.d0&
                                     + w(k,3)*7897110375000.d0&
                                     - w(k,4)*192082660000.d0&
                                     + w(k,5)*400298976000.d0&
                                     + w(k,7)*201586000000.d0  )
  230       continue
            call fcn(n, x + c(18), w(1,9), w(1,8),zz)
!c
!c           calculate ytrial, the extrapolated approximation and store
!c              in w(*,9)
            do 235 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0&
                                     + w(k,3)*545186250000.d0&
                                     + w(k,4)*446637345000.d0&
                                     + w(k,5)*188806464000.d0&
                                     + w(k,7)*15076875000.d0&
                                     + w(k,8)*97599465000.d0   )
  235       continue
!c
!c           add 7 to the no of fcn evals
            c(24) = c(24) + 7.d0
!c
!c        end stage 2
!c
!c        ***************************************************************
!c        * stage 3 - calculate the error estimate est. first calculate *
!c        * the  unweighted  absolute  error  estimate vector (per unit *
!c        * step) for the unextrapolated approximation and store it  in *
!c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!c        * specified by the error  control  indicator  c(1).  finally, *
!c        * modify  this result to produce est, the error estimate (per *
!c        * unit step) for the extrapolated approximation ytrial.       *
!c        ***************************************************************
!c
!c           calculate the unweighted absolute error estimate vector
            do 300 k = 1, n
               w(k,2) = (   w(k,1)*8738556750.d0&
                         + w(k,3)*9735468750.d0&
                         - w(k,4)*9709507500.d0&
                         + w(k,5)*8582112000.d0&
                         + w(k,6)*95329710000.d0&
                         - w(k,7)*15076875000.d0&
                         - w(k,8)*97599465000.d0)/1398169080000.d0
  300       continue
!c
!c           calculate the weighted max norm of w(*,2) as specified by
!c           the error control indicator c(1)
            temp = 0.d0
            if (c(1) .ne. 1.d0) go to 310
!c              absolute error control
               do 305 k = 1, n
                  temp = dmax1(temp,dabs(w(k,2)))
  305          continue
               go to 360
  310       if (c(1) .ne. 2.d0) go to 320
!c              relative error control
               do 315 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/y(k)))
  315          continue
               go to 360
  320       if (c(1) .ne. 3.d0) go to 330
!c              weights are 1/max(c(2),abs(y(k)))
               do 325 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))&
                                  / dmax1(c(2), dabs(y(k))) )
  325          continue
               go to 360
  330       if (c(1) .ne. 4.d0) go to 340
!c              weights are 1/max(c(k+30),abs(y(k)))
               do 335 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))&
                                  / dmax1(c(k+30), dabs(y(k))) )
  335          continue
               go to 360
  340       if (c(1) .ne. 5.d0) go to 350
!c              weights are 1/c(k+30)
               do 345 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
  345          continue
               go to 360
  350       continue
!c              default case - weights are 1/max(1,abs(y(k)))
               do 355 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))&
                                  / dmax1(1.d0, dabs(y(k))) )
  355          continue
  360       continue
!c
!c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!c              - est is intended to be a measure of the error  per  unit
!c              step in ytrial
            c(19) = temp*c(14)*c(15)
!c
!c        end stage 3
!c
!c        ***************************************************************
!c        * stage 4 - make decisions.                                   *
!c        ***************************************************************
!c
!c           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
!c
!c***********interrupt no 2 if requested
            if (c(9) .eq. 0.d0) go to 2222
               return
!c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
!c
            if (ind .eq. 6) go to 410
!c              step accepted (ind .eq. 5), so update x, y from xtrial,
!c                 ytrial, add 1 to the no of successful steps, and set
!c                 the no of successive failures to zero
               x = c(17)
               do 400 k = 1, n
                  y(k) = w(k,9)
  400          continue
               c(22) = c(22) + 1.d0
               c(23) = 0.d0
!c**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1.d0
                  return
  405          continue
               go to 420
  410       continue
!c              step not accepted (ind .eq. 6), so add 1 to the no of
!c                 successive failures
               c(23) = c(23) + 1.d0
!c**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
!c
!c        end stage 4
!c
      go to 99999
!c     end loop
!c
!c  begin abort action
  500 continue
!c
      write(6,505) ind, tol, x, n, c(13), xend, nw, c(16), c(20),&
           c(22), c(23), c(24), (y(k), k = 1, n)
  505 format( /// 1h0, 58hcomputation stopped in dverk with the followin&
     g values -&
        / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,&
               1pd22.15&
        / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =,&
               1pd22.15&
        / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =,&
               1pd22.15&
        / 1h0, 14x, 27hno of successful steps    =, 0pf8.0&
        / 1h , 14x, 27hno of successive failures =, 0pf8.0&
        / 1h , 14x, 27hno of function evals      =, 0pf8.0&
        / 1h0, 23hthe components of y are&
        // (1h , 1p5d24.15)                                           )
!c
      stop
!c
!c  end abort action
!c
      end

      
      
      
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
      gama=0.3d0
      !z1=e0
      !z2=w
      !z3=a
      !z4=q
      !z5=xe
      !z6=gama
      dp=dexp(-2.d0*x)-dexp(-x)+e0*dsin(w*t)*(a*dcos(a*(x+xe))*dexp(-q*(x+xe)**4)-dsin(a*(x+xe))*4.d0*q*((x+xe)**3)*dexp(-q*(x+xe)**4))/a-z6*p
      
      

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