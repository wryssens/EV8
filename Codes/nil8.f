c______________________________________________________________________________
      program nil8


c-- Vendredi 7 janvier 2005


c    ***********************************************************************
c    *                                                                     *
c    * Copyright  P. Bonche, H. Flocard, P.H. Heenen                       *
c    *            in case of trouble, take it easy...                      *
c    *            and (please) let us know                                 *
c    *                                                                     *
c    ***********************************************************************
c    *                                                                     *
c    * THIS PROGRAM CAN BE OBTAINED FREELY FROM ANY OF THE ABOVE PERSONS   *
C    *   Reference to them should be made in all publications presenting   *
c    *             results obtained with this program.                     *
c    *                                                                     *
c    ***********************************************************************

c    *********************************************************************
c    *                                                                   *
c    * nilson program                                                    *
c    * skyrme interaction,time-reversal symmetry                         *
c    *   available forces:                                               *
c    *        SkM* SIII Ska SGII RATP T6 Sly4 Sly5 Sly6 Sly7 SkP         *
c    * approximate correction for the c.m. motion and coulomb exchange   *
c    * nucleus symmetric through x=0,y=0,z=0 plane symmetries            *
c    * in mome the slicing is done along the z-axis if icqx = 5 0r 6     *
c    *                                       y-axis           3 or 4     *
c    *                                       x-axis           1 or 2     *
c    *                                                                   *
c    *********************************************************************

c      ****************************************************************
c      *                                                              *
c      * this subroutine reads the data and setups the starting wave- *
c      * functions                                                    *
c      * this version can only start from nilson wave functions.      *
c      * fort.13 contains in output the final wave-functions.         *
c      *                                                              *
c      * to compile properly the program requires a param8.h file     *
c      *   containing the folowing information:                       *
c      * parameter (mx=14,my=14,mz=14,mc=8,mv=mx*my*mz,mq=4*mv,mw=40) *
c      * parameter (meven=5,modd=4)                                   *
c      *   where                                                      *
c      *   mx, my,mz are : numbers of points in each directions       *
c      *   meven, modd   : numbers of major hamonic oscillator shells *
c      *                   meven is either equal to modd or to modd+1 *
c      *                                                              *
c      ****************************************************************

c            *******************************************************
c            *                      data                           *
c            *  - head                                -  20a4   -  *
c            *  - nwaven,nwavep,npn,npp               -  4i5    -  *
c            *  - njmunu,ncm2,nmass                   -  3i5    -  *
c            *  - afor                                -    a4   -  *
c            *      if (afor.eq.'XXXX')                         -  *
c            *       - t0,x0,                         -  2e15.8 -  *
c            *       - t1,x1,t2,x2                    -  4e15.8 -  *
c            *       - t3a,x3a,yt3a                   -  3e15.8 -  *
c            *       - t3b,x3b,yt3b                   -  3e15.8 -  *
c            *       - wso,wsoq                       -  2e15.8 -  *
c            *  - nx,ny,nz                            -  3i5    -  *
c            *  - dx                                  -   e15.8 -  *
c            *  - iq1,iq2                             -  2i5    -  *
c            *  - alpha,qqq                           -  2e15.8 -  *
c            *  - rcut                                -   e15.8 -  *
c            *                                                     *
c            *******************************************************

c   **************************************************************
c   * Program logic                                              *
c   *   - nil8                                                   *
c   *     - lecfor           skyrme force                        *
c   *       - blockdata case                                     *
c   *     - nilson           initial nilson orbitals             *
c   *     - evolve           calculation of the energy           *
c   *       - inisp                                              *
c   *       - class                                              *
c   *       - gaphf                                              *
c   *       - ortho          orthogonalization of the spwf       *
c   *       - densit         density                             *
c   *       - newpot         H.F. potential                      *
c   *         - mome                                             *
c   *         - vbound       coulomb boudary conditions          *
c   *         - vcoul        coulomb potential                   *
c   *         - vcal         skyrme contributions                *
c   *       - figaro         printout                            *
c   *         - fprte                                            *
c   *         - fprtj                                            *
c   *         - fprtz                                            *
c   *         - fprtm                                            *
c   *     - wrini                                                *
c   *     - writ8            storage of the spwf on fort.13      *
c   * and also :                                                 *
c   *     - diagon                                               *
c   *     - deriv, der, lapla                                    *
c   *     - scopy, sdot, saxpy sscal                             *
c   **************************************************************

c Number of orbitals following the Nilsson filling
c    n+/n- number of positive/parity parity orbitals per (sub)shells
c    N+/N- cumulative number for the corresponding filling
c    N     total number of particles
c    ( + - ) numbers of positive/negative orbitals 
c            (2 nucleons per orbits due to time reversal invariance)

c            ______________________________________________________
c           |    orbitals        |  occ  |  total   | shell closure|
c           |                    | n+ n- |  N+   N- |   N  ( +  -) |
c           |______________________________________________________|
c           |  1 s               |  2    |   2      |   2          |
c           |           1 p      |     6 |        6 |   8  ( 1, 3) |
c           |  2 sd              | 12    |  14      |              |
c           |           2 f7/2   |     8 |       14 |  28  ( 7, 7) |
c           |           2 p,f5/2 |    12 |       26 |              |
c           |  3 g9/2            | 10    |  24      |  50  (12,13) |
c           |  3 sdg7/2          | 20    |  44      |              |
c           |           3 h11/2  |    12 |       38 |  82  (22,19) |
c           |           3 pfh9/2 |    30 |       68 |              |
c           |  4 i13/2           | 14    |  58      | 126  (29,34) |
c           |  4 g9/2i11/2       | 22    |  80      |              |
c           |           4 j15/2  |    16 |       84 | 164  (42,42) |
c           |  4 sdg7/2          | 20    | 100      |              |
c           |           4 h11/2  |    12 |       96 | 196  (50,48) |
c           |           4 pfj13/2|    34 |      130 |              |
c           |           4 h9/2   |    10 |      140 | 240  (50,70) |
c           |  5 k17/2           | 18    | 118      |              |
c           |  5 i13/2           | 14    | 132      | 272  (66,70) |
c           |______________________________________________________|



c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (t100=100.0d0,t180=180.0d0)
      parameter (thbar=6.58218d0,txmn=1.044673d0)
      character*4 afor,head

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz
      common /nxyz / dx,dv
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /wave / p1(mx,my,mz),p2(mx,my,mz),p3(mx,my,mz),p4(mx,my,mz)
     1              ,w1(mx,my,mz),w2(mx,my,mz),w3(mx,my,mz),w4(mx,my,mz)


c.... read formats
  101 format (20a4)
  102 format (5e15.8)
  103 format (6i5)
  105 format (i5,5e15.8)

c.... write formats
  200 format (/,'  __________________________________________________ ',
     1        /,' |                                                  |',
     2        /,' |  program nil8                    version 1.0.0   |',
     3        /,' |  Friday, January 7th, 2005                       |',
     4        /,' |  Copyright  P. Bonche, H. Flocard, P.H. Heenen   |',
     5        /,' |    in case of trouble, take it easy...           |',
     6        /,' |    and (please) let us know                      |',
     7        /,' |__________________________________________________|')
  201 format (//,' parameters: mx,my,mz ',3i5,
     2         /,'                   mw ',i5)
  203 format (//,' ',78('_'),/,' run information',//,20a4)
  207 format (/,'  njmunu = ',i5)
  208 format (/,'  njmunu = ',i5,' -- spin tensor is self-consistent ')
  209 format (  '  ncm2   = ',i5)
  210 format (  '  ncm2   = ',i5,' -- 2-body c.m. is self-consistent ')
  211 format (  '  nmass  = ',i5,' (m_n = m_p)')
  212 format (  '  nmass  = ',i5,' experimental nucleon masses')
  230 format (/,'  Skyrme force: ',a4,'  (new set of parameters)')
  231 format (/,'  Skyrme force: ',a4)
  232 format (/,'  ____________________________________________     ',
     1        /,' |        !!!!  ATTENTION PLEASE   !!!!       |    ',
     2        /,' |    non standard choice of njmunu :',i2,'       |',
     3        /,' |                        instead of ',i2,'       |',
     4        /,' |  The program will resume at your own risk  |    ',
     5        /,' |____________________________________________|    ',/)
  233 format (/,'  ____________________________________________     ',
     1        /,' |        !!!!  ATTENTION PLEASE   !!!!       |    ',
     2        /,' |    non standard choice of ncm2   :',i2,'       |',
     3        /,' |                        instead of ',i2,'       |',
     4        /,' |  The program will resume at your own risk  |    ',
     5        /,' |____________________________________________|    ',/)
  234 format (/,'  ____________________________________________     ',
     1        /,' |        !!!!  ATTENTION PLEASE   !!!!       |    ',
     2        /,' |    non standard choice of nmass  :',i2,'       |',
     3        /,' |                        instead of ',i2,'       |',
     4        /,' |  The program will resume at your own risk  |    ',
     5        /,' |____________________________________________|    ',/)
  236 format (/,'  ____________________________________________     ',
     1        /,' |        !!!!  ATTENTION PLEASE   !!!!       |    ',
     2        /,' |          The program was run with          |    ',
     3        /,' |     a non standard choice of njmunu        |',
     4        /,' |____________________________________________|    ',/)
  237 format (/,'  ____________________________________________     ',
     1        /,' |        !!!!  ATTENTION PLEASE   !!!!       |    ',
     2        /,' |          The program was run with          |    ',
     3        /,' |     a non standard choice of ncm2          |',
     4        /,' |____________________________________________|    ',/)
  238 format (/,'  ____________________________________________     ',
     1        /,' |        !!!!  ATTENTION PLEASE   !!!!       |    ',
     2        /,' |          The program was run with          |    ',
     3        /,' |     a non standard choice of nmass         |',
     4        /,' |____________________________________________|    ',/)
  240 format   ('  t0 =',f13.6,' x0 =',f12.6,/,
     1          '  t1 =',f13.6,' x1 =',f12.6,/,
     2          '  t2 =',f13.6,' x2 =',f12.6,/,
     3          '  t3a=',f13.6,' x3a=',f12.6,'  sigma =',f10.6,/,
     4          '  t3b=',f13.6,' x3b=',f12.6,'  sigmb =',f10.6,/,
     5          '   w =',f13.6,' wq =',f12.6)
  250 format (/,'  nx=',i5,6x,'  ny=',i5,6x,'  nz=',i5,/,
     1          '  dx=',f8.3,' Fm')
  260 format (/,'  number of wave-functions n=',i5,'  p=',i5)
  261 format   ('  mass =',i5,'  charge =',i5,'  n =',i5)
  262 format (//,' ',78('_'))
  264 format   (' iq1 iq2=',2i6)
  265 format   (' parameters al0=',f8.3,'    q=',f8.3,'  gamma=',f8.3)
  266 format   ('            alx=',f8.3,'  aly=',f8.3,'    alz=',f8.3)
  267 format (/,' order of the quadrupole moments i=',i5,/,
     1          '  i=1 x,y,z    i=2 x,z,y    i=3 y,x,z',/,
     2          '  i=4 y,z,x    i=5 z,x,y    i=6 z,y,x')
  268 format   ('  rcut = ',f8.3,/,'  acut = ',f8.3)
  280 format (/,'  characteristics of the neutron wave-functions')
  282 format (/,'  characteristics of the proton wave-functions')
  283 format   ('  blocks ((+++),(--+),(-+-),(+--)) ',
     1                  ' ((++-),(---),(-++),(+-+)',/,21x,i5,22x,i5)
  299 format (//,' ',78('_'))


c........................................................................
      read  101,head
      print 200
      print 201,mx,my,mz,mw
      print 203,head

c........................................................................
      pi   = tt4 * atan2(one,one)
      traf = pi/t180

      do i=1,mq
        w1(i,1,1) = zero
        p1(i,1,1) = zero
      enddo

c............ number of neutron and proton wave-functions:  nwaven,nwavep
c                                  number of neutrons and protons npn,npp
c                                  npn should be less than 2*nwaven
c                                  npp should be less than 2*nwavep
c   a nucleon and its time time-reversed for each wave-function
c   a wave function is made of four 3d components,
c   they are plus (plus)  the   real    part of the spin up (down),
c            plus (minus) the imaginary part of the spin up (down),
c            plus (minus) the   real    part of the spin down (up)
c            plus (plus)  the imaginary part of the spin down (up)
c            of the s.p.states (time reversed states).
c       the first component of the s.p. reference states is by convention
c                             symmetric through x=0, y=0 plane inversions
c       if kparz (+ or -1) is the parity of the first of these components
c       across the z plane, the parities of the 4 components with respect
c       to the x, y and z plane are:  (+,+, kparz)  spin  up  real
c                                     (-,-, kparz)  spin  up  imaginary
c                                     (-,+,-kparz)  spin down real
c                                     (+,-,-kparz)  spin down imaginary
      read  103,nwaven,nwavep,npn,npp
      print 260,nwaven,nwavep
      if (nwaven*nwavep.eq.0) stop 'nwave'
      if (npp.eq.0) npp=2*nwavep
      if (npn.eq.0) npn=2*nwaven
      ato = npp + npn
      print 261,npp+npn,npp,npn
      if (npp.gt.2*nwavep) stop 'npp'
      if (npn.gt.2*nwaven) stop 'npn'

c...................... t0*(1+x0*ps)*del+(t1/2)*(1+x1*ps)*(k2*del+del*k2)
c                 +t2*(1+x2*ps)*(k*del*k)+(t3/6)*(1+x3*ps)*(rho**yt3)*del
c                                      +i*wso*(sig1+sig2)*(k*cross*del*k)
      read  103,mjmunu,mcm2,mmass
      if (mjmunu.eq.0) print 207,mjmunu
      if (mjmunu.eq.1) print 208,mjmunu
      if   (mcm2.eq.0) print 209,mcm2
      if   (mcm2.eq.1) print 210,mcm2
      if  (mmass.eq.0) print 211,mmass
      if  (mmass.eq.1) print 212,mmass

      read  101,afor
      call lecfor (kfor)
      if (kfor.eq.0) then
        read  102,t0,x0
        read  102,t1,x1,t2,x2
        read  102,t3a,x3a,yt3a
        read  102,t3b,x3b,yt3b
        read  102,wso,wsoq
        njmunu = mjmunu
        ncm2   = mcm2
        nmass  = mmass
        mjmunu = 0
        mcm2   = 0
        mmass  = 0
        print 230,afor
      else
        print 231,afor
        if (mjmunu.ne.njmunu) then
          print 232,mjmunu,njmunu
          njmunu = mjmunu
          mjmunu = 1
        else
          mjmunu = 0
        endif
        if (mcm2.ne.ncm2) then
          print 233,mcm2,ncm2
          ncm2 = mcm2
          mcm2 = 1
        else
          mcm2 = 0
        endif
        if (mmass.ne.nmass) then
          print 234,mmass,nmass
          nmass = mmass
          mmass = 1
        else
          mmass = 0
        endif
      endif
      call lecmas
      print 240,t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                           ,t3b,x3b,yt3b,wso,wsoq

c.......... characteristics of the mesh, the time step and the total time
      read  103,nx,ny,nz,npx,npy,npz
      read  102,dx
      print 250,nx,ny,nz,dx
      if (nx.ne.mx)  stop
      if (ny.ne.my)  stop
      if (nz.ne.mz)  stop
      if (npx.gt.mc) stop
      if (npy.gt.mc) stop
      if (npz.gt.mc) stop
      nnx = mx + npx
      nny = my + npy
      nnz = mz + npz

c.......................................... characteristics of the output
c              npair=0 h.f. filling

c............. characteristics of the starting point    (modein=0   stop)
c       modein = 1-6 eigenvectors of the nilsson hamiltonian
c           paramaters kappa and mu in data in subroutine nilson
c               alpha must be positif = m*omega0/hbar
c                 qqq must be larger than one
c                 ho+ = ho0*qqq**(-2.*cos(xgam)/3.)
c                 hom = ho0*qqq**(-2.*cos(xgam-2*pi/3)/3.)
c                 ho- = ho0*qqq**(-2.*cos(xgam+2*pi/3)/3.)
c        thus  ho0**3 = ho+*hom*ho-
c        if xgam =  0 prolate nucleus qqq=omega perp/omega axis=ho-/ho+
c        if xgam = 60  oblate nucleus qqq=omega axis/omega perp=ho-/ho+

c        the decreasing order of the quadrupole moments was determined by
c        icqx (1=x,y,z-2=x,z,y-3=y,x,z-4=y,z,x-5=z,x,y-6=z,y,x)
c        in this program, icqx is no longer used

c   qx     = -delq2*(iq1-iq2)/2.
c   qy     = -delq2*(iq1+2*iq2)/2.
c   qz     =  delq2*(2*iq1+iq2)/2.

c   q0     =  delq2*sqrt(iq1**2+iq2**2+iq1*iq2)

c   note: qx = -q0*cos(gam+60)
c         qy = -q0*cos(gam-60)
c         qz =  q0*cos(gam)
c        gam = atan2 (qx-qy,sqrt(3.)*qz)
c            = atan2 (iq2*sqrt(3.),2*iq1+iq2)
c            = 0 at the spherical point (iq1=iq2=0)

      read  103,iq1,iq2
      read  102,alpha,qqq
      if ((iq1.eq.0).and.(iq2.eq.0)) qqq = one
      xgam = zero
      if ((iq1.ne.0).or.(iq2.ne.0)) then
        ca = iq2 * sqrt(tt3)
        sa = 2 * iq1 + iq2
        xgam = atan2(ca,sa)/traf
      endif

      print 262
      print 264,iq1,iq2
      print 265,alpha,qqq,xgam
      if (alpha.lt.zero) stop
      if   (qqq.lt.one)  stop
      xgam = xgam*traf
      alz  = alpha*qqq**(-two*cos(xgam)/tt3)
      alx  = alpha*qqq**(-two*cos(xgam-two*pi/tt3)/tt3)
      aly  = alpha*qqq**(-two*cos(xgam+two*pi/tt3)/tt3)
      print 266,alx,aly,alz

      call nilson (alx,aly,alz)

      do i=1,mq
        w1(i,1,1) = zero
      enddo

c.....................................   calculation of icqx and of iaxis
      if (iq1.lt.0) go to 36
      icqx=5
      if (iq2.ge.0) go to 42
      icqx=6
      if (iq1+iq2.lt.0) icqx=4
      go to 42
   36 icqx=3
      if (iq2.lt.0) go to 42
      icqx=1
      if (iq1+iq2.ge.0) icqx=2
   42 continue

c................... rcut= radius of the spherical cut-off function on q2
      read 102,rcut
      if (rcut.eq.zero) rcut=t100
      acut = two*dx
      print 262
      print 267,icqx
      print 268,rcut,acut

c................... determination of the vectors having a given symmetry
      do it=1,2
        if (it.eq.1) print 280
        if (it.eq.2) print 282
        print 283,(npar(i,it),i=1,2)
      enddo

c................................... standard output with skyrme energies
      call evolve

c.................. storage of the wave-functions and print of input data
      print 299

      call wrini

      open (13,form='unformatted')
      call writ8
      close (13)

      if (mjmunu.eq.1) print 236
      if   (mcm2.eq.1) print 237
      if  (mmass.eq.1) print 238

      stop
      end

c______________________________________________________________________________
      subroutine lecfor (kfor)

      implicit real*8 (a-h,o-z)
      character*4 afor

      parameter (zero=0.0d0)

      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /skm  / prm(20,16)


c..............................................................................
c           t_0      t_1       t_2       t_3   
c  1  SkM*  -2645.0    410.0    -135.0    15595.0 
c  2  SIII  -1128.75   395.0     -95.0    14000.0 
c  3  Ska   -1602.78   570.88    -67.7     8000.0 
c  4  SGII  -2645.00   340.0     -41.9    15595.0 
c  5  RATP  -2160.00   513.0     121.0    11600.0 
c  6  T6    -1794.20   294.0    -294.0    12817.0 
c  7  Sly4  -2488.913  486.818  -546.395  13777.0 
c  8  Sly5  -2483.450  484.230  -556.690  13757.0 
c  9  Sly6  -2479.500  462.180  -448.610  13673.0 
c 10  Sly7  -2480.800  461.290  -433.930  13669.0 
c 11  SkP   -2931.70   320.62   -337.41   18708.97
c 12  St3a
c 13  St3b
c 14  Smnp

c          x_0      x_1      x_2      x_3 
c  SkM*   0.090    0.000    0.000    0.000
c  SIII   0.450    0.000    0.000    1.000
c  Ska   -0.020    0.000    0.000   -0.286
c  SGII   0.090   -0.0588   1.425    0.060
c  RATP   0.418   -0.360   -2.290    0.586
c  T6     0.392   -0.500   -0.500    0.500
c  Sly4   0.834   -0.344   -1.0      1.354
c  Sly5   0.776   -0.317   -1.0      1.263
c  Sly6   0.825   -0.465   -1.0      1.355
c  Sly7   0.848   -0.492   -1.0      1.393
c  SkP    0.29215  0.65318 -0.53732  0.18103 

c          disper    W_so   W_so
c  SkM*   1.0  6.0  130.0  130.0
c  SIII   1.0  1.0  120.0  120.0
c  Ska    1.0  3.0  125.0  125.0
c  SGII   1.0  6.0  105.0  105.0
c  RATP   1.0  5.0  120.0  120.0
c  T6     1.0  3.0  107.0  107.0
c  Sly4   1.0  6.0  123.0  123.0
c  Sly5   1.0  6.0  125.0  125.0
c  Sly6   1.0  6.0  122.0  122.0
c  Sly7   1.0  6.0  125.0  125.0
c  SkP    1.0  6.0  100.0  100.0


c..............................................................................
      kfor = 0

      if (afor.eq.'Skm*') kfor =  1
      if (afor.eq.'SIII') kfor =  2
      if (afor.eq.'Ska ') kfor =  3
      if (afor.eq.'SGII') kfor =  4
      if (afor.eq.'RATP') kfor =  5
      if (afor.eq.'T6  ') kfor =  6
      if (afor.eq.'Sly4') kfor =  7
      if (afor.eq.'Sly5') kfor =  8
      if (afor.eq.'Sly6') kfor =  9
      if (afor.eq.'Sly7') kfor = 10
      if (afor.eq.'SkP ') kfor = 11
      if (afor.eq.'St3a') kfor = 12
      if (afor.eq.'St3b') kfor = 13
      if (afor.eq.'Smnp') kfor = 14

      if (kfor.eq.0) return

      t0   = prm( 1,kfor)
      t1   = prm( 2,kfor)
      t2   = prm( 3,kfor)
      t3a  = prm( 4,kfor)
      t3b  = prm( 5,kfor)
      x0   = prm( 6,kfor)
      x1   = prm( 7,kfor)
      x2   = prm( 8,kfor)
      x3a  = prm( 9,kfor)
      x3b  = prm(10,kfor)
      yt3a = prm(11,kfor)/prm(12,kfor)
      yt3b = prm(13,kfor)/prm(14,kfor)

      wso  = prm(16,kfor)
      wsoq = prm(17,kfor)

      ncm2   = 0
      njmunu = 0
      nmass  = 0
      if (prm(18,kfor).ne.zero) ncm2   = 1
      if (prm(19,kfor).ne.zero) njmunu = 1
      if (prm(20,kfor).ne.zero) nmass  = 1

      return
      end

c______________________________________________________________________________
      subroutine lecmas

      implicit real*8 (a-h,o-z)
      character*4 afor

      parameter (two=2.0d0)
      parameter (hhbar=6.58218d0,xxmn =1.044673d0)
      parameter (hc= 197.327053d0,xmasn=939.565360d0,xmasp=938.272029)
      parameter (clum=29.9792458d0)
c.. these values are taken from the Particle Physics Booklet (july 2004)

      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /skm  / prm(20,16)


c............................................................. constants
c   hbar   (MeV*10-22sec)
c          (average between neutron and proton)
c   hbm  = hbar**2/m

      if (nmass.eq.0) then
        hbar =  hhbar   ! = 6.58218d0
        xmn  =  xxmn    ! = 1.044673d0
        hbm(1)  = hbar*hbar/xmn
        hbm(2)  = hbar*hbar/xmn
        xmasrd  = two*xmasn*xmasp/(xmasn+xmasp)
        xm(1)   = xmasrd
        xm(2)   = xmasrd
      else
        hbarc  = hc             ! = 197.327 053 d0 
        hbar   = hc/clum        ! =  29.979 245 8d0
        xm(1)  = xmasn          ! = 939.565 330 d0
        xm(2)  = xmasp          ! = 938.271 998 d0 
        hbm(1) = hc**2 / xmasn
        hbm(2) = hc**2 / xmasp      
      endif
      xm(3)   = npn * xm(1) + npp * xm(2)

      return
      end

c______________________________________________________________________________
      blockdata case

      implicit real*8 (a-h,o-z)

      common /skm  / prm (20,16)

      data prm
c............................................................... SkM*  1
     1 / -2645.0d0,    410.0d0,   -135.0d0,   15595.0d0,       0.0d0,
     1       0.090d0,    0.000d0,    0.000d0,     0.000,       0.0d0,
     1       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     1     130.0d0,    130.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... SIII  2
     2   -1128.75d0,   395.0d0,    -95.0d0,   14000.0d0,       0.0d0,
     2       0.450d0,    0.000d0,    0.000d0,     1.000,       0.0d0,
     2       1.0d0,      1.0d0,      0.0d0,       1.0d0,       0.0d0,
     2     120.0d0,    120.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... Ska   3
     3   -1602.78d0,   570.88d0,   -67.7d0,    8000.0d0,       0.0d0,
     3      -0.020d0,    0.000d0,    0.000d0,    -0.286d0,     0.0d0,
     3       1.0d0,      3.0d0,      0.0d0,       1.0d0,       0.0d0,
     3     125.0d0,    125.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... SGII  4
     4   -2645.00d0,   340.0d0,    -41.9d0,   15595.0d0,       0.0d0,
     4       0.090d0,   -0.0588d0,   1.425d0,     0.0604d0,    0.0d0,
     4       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     4     105.0d0,    105.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... RATP  5
     5   -2160.00d0,   513.0d0,    121.0d0,   11600.0d0,       0.0d0,
     5       0.418d0,   -0.360d0,   -2.290d0,     0.586d0,     0.0d0,
     5       1.0d0,      5.0d0,      0.0d0,       1.0d0,       0.0d0,
     5     120.0d0,    120.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... T6    6
     6   -1794.20d0,   294.0d0,   -294.0d0,   12817.0d0,       0.0d0,
     6       0.392d0,   -0.500d0,   -0.500d0,     0.500,       0.0d0,
     6       1.0d0,      3.0d0,      0.0d0,       1.0d0,       0.0d0,
     6     107.0d0,    107.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... Sly4  7
     7   -2488.913d0,  486.818d0, -546.395d0, 13777.0d0,       0.0d0,
     7       0.834d0,   -0.344d0,   -1.0d0,       1.354,       0.0d0,
     7       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     7     123.0d0,    123.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... Sly5  8
     8   -2483.450d0,  484.230d0, -556.690d0, 13757.0d0,       0.0d0,
     8       0.776d0,   -0.317d0,   -1.0d0,       1.263d0,     0.0d0,
     8       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     8     125.0d0,    125.0d0,      0.0d0,       1.0d0,       0.0d0,
c............................................................... Sly6  9
     9   -2479.500d0,  462.180d0, -448.610d0, 13673.0d0,       0.0d0,
     9       0.825d0,   -0.465d0,   -1.0d0,       1.355d0,     0.0d0,
     9       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     9     122.0d0,    122.0d0,      1.0d0,       0.0d0,       0.0d0,
c............................................................... Sly7 10
     A   -2480.800d0,  461.290d0, -433.930d0, 13669.0d0,       0.0d0,
     A       0.848d0,   -0.492d0,   -1.0d0,       1.393d0,     0.0d0,
     A       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     A     125.0d0,    125.0d0,      1.0d0,       1.0d0,       0.0d0,
c............................................................... SkP  11
     B   -2931.70d0,   320.62d0,  -337.41d0,  18708.97d0,      0.0d0,
     B       0.29215d0,  0.65318d0, -0.53732d0,   0.18103d0,   0.0d0,
     B       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     B     100.0d0,    100.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... St3a 12
     C   -1855.90419330d0,     481.06391932d0,  -348.59605803d0,
     C   14067.9736d0,       -4276.44686d0,
     C       0.70980632d0,       0.01695512d0,    -0.85227559d0,
     C       1.54481568d0,       3.34586349d0,
     C       1.0d0,              3.0d0,
     C       2.0d0,              3.0d0,              0.0d0,
     C     120.44329826d0,     120.44329826d0,
     C       0.0d0,              1.0d0,              0.0d0,
c............................................................... St3b 13
     D   -1855.22721423d0,     481.21024666d0,  -250.92888862d0,
     D   14048.405767d0,     -4253.789521d0,
     D       0.73031876d0,      -0.13489175d0,    -0.69796468d0,
     D       1.60209914d0,       3.46279658d0,
     D       1.0d0,              3.0d0,
     D       2.0d0,              3.0d0,              0.0d0,
     D     117.70472089d0,     117.70472089d0,
     D       0.0d0,              0.0d0,              0.0d0,
c............................................................... Smnp 14
     E   -2496.4745d0, 487.5245d0,-598.5310d0,13899.3413d0,    0.0d0,
     E       0.8605d0,   0.1765d0,  -1.0d0,       1.2614d0,    0.0d0,
     E       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     E     123.0547d0, 123.0547d0,   0.0d0,       0.0d0,       1.0d0,
c............................................................... S... 15
     F       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     F       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     F       0.0d0,      1.0d0,      0.0d0,       1.0d0,       0.0d0,
     F       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c............................................................... S... 16
     G       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     G       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     G       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     G       0.0d0,      1.0d0,      0.0d0,       1.0d0,       0.0d0 /

      end

c______________________________________________________________________________
      subroutine nilson (hox,hoy,hoz)

      implicit real*8 (a-h,o-z)
      character*4 afor
      include 'param8.h'

      parameter (mblc=meven+modd+1)
      parameter (ms=(mblc*(mblc**2-1))/6)
      parameter (ndim=(mblc*(mblc-1))/2)
      parameter (mqa=(ms*(3*mblc**2-2))/10)

      parameter (zero=0.0d0,one=1.0d0,tp5=0.5d0,two=2.0d0,tt3=3.0d0)
      parameter (t50=50.0d0)

      common /big  / h(ndim,ndim),s(ndim,ndim),d(ndim),wd(ndim)
      common /noyau/ nwn,nwp,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /orbit/ e(ms),he(mblc,mz,3),nor(ms),npa(ms),ntrs(ms)
     2              ,nx(ms),ny(ms),nz(ms),nsi(mblc+1,4),ns(mblc+1)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /stor / ark(mq,mw)
      common /wave / wf(mx,my,mz*4),psi(mx,my,mz*4)

      dimension a(mqa),irep(mblc+1)
      dimension xk(4),xmu(4),cf(2)
      equivalence (a,wf)

      data ca,cb /0.986d0,0.14d0/
      data xk,xmu/0.08d0,0.08d0,0.0637d0,0.0637d0
     1           ,0.0d0 ,0.0d0 ,  0.42d0,  0.60d0   /

c......................... mz must be larger or equal than both mx and my

c     neven   nodd   nvec+   nvec-   nblc    ms   mblc   ndim     mqa
c       1       0      1       0       1      1     2       1       1
c       1       1      0       3       2      4     3       3      10
c       2       1      6       0       3     10     4       6      46
c       2       2      0      10       4     20     5      10     146
c       3       2     15       0       5     35     6      15     371
c       3       3      0      21       6     56     7      21     812
c       4       3     28       0       7     84     8      28    1592
c       4       4      0      36       8    120     9      36    2892
c       5       4     45       0       9    165    10      45    4917
c       5       5      0      55      10    220    11      55    7942
c       6       5     66       0      11    286    12      66   12298
c       6       6      0      78      12    364    13      78   18382
c       7       6     91       0      13    455    14      91   26663
c       7       7      0     105      14    560    15     105   37688


c     parameter (meven=7,modd=6,mblc=meven+modd+1)
c     parameter (ms=(mblc*(mblc**2-1))/6)
c     parameter (ndim=(mblc*(mblc-1))/2)
c     parameter (mqa=(ms*(3*mblc**2-2))/10)
c  meven = number of even harmonic oscillator shells
c  modd  = number of odd  harmonic oscillator shells
c ! Note that neven must be either equal to nodd or to nodd + 1
c     mblc = meven + modd + 1
c     ms   = (mblc*(mblc**2-1))/6)
c     ndim = (mblc*(mblc   -1))/2)
c     mqa  = (ms*(3*mblc**2-2))/10)


c............. this subroutine determines the starting point by selecting
c              the nwn(nwp) lowest eigenstates of the nilsson hamiltonian
c              h = sum(1 to 3) of hoi*(ni+1/2)-xk*ho0*(2*l*s+xmu*l**2)
c              the quantities hox(hoy,hoz) are in fact m*omegax(y,z)/hbar
c              are given in data 
c              (G.Gustafson, I.L.Lamm, B.Nilsson and S.G.Nilsson, 
c               Ark Fys 36(1966)613)
c              the operator l is the stretched angular momentum.
c              the operator l**2 is corrected so that its trace over each
c               major shell is zero.(see also copybook n0 17)


  100 format (/,' standard order nilsson hamiltonian')
  101 format (/,' neutron levels kappa=',e10.3,' mu=',e9.2,
     1          ' al0n=al0*(1+',f5.2,'*(n-z)/a)',/,
     2          ' (n0,nor,energy/(hbar*omega0),parity)',/,' ')
  102 format (/,' proton  levels kappa=',e10.3,' mu=',e9.2,
     1          ' al0p=al0*(',f6.3,'-',f5.2,'*(n-z)/a)',/,
     2          ' (n0,nor,energy/(hbar*omega0),parity)',/,' ')
  103 format (' (',2i4,f8.3,i3,') (',2i4,f8.3,i3,') (',2i4,f8.3,i3,')')

  201 format (//,'    _______________________________',
     1         /,'   |                               |',
     2         /,'   | neven = ',i2,'  nodd = ',i2,'         |',
     3         /,'   | wrong values !!               |',
     4         /,'   | neven must be either equal to |',
     5         /,'   |       nodd or to nodd + 1     |'
     6         /,'   |_______________________________|')
  202 format (//,'    ______________________________ ',
     1         /,'   |                              |',
     2         /,'   | meven = ',i2,'  modd = ',i2,'        |',
     3         /,'   | mblc  = ',i2,'  too small !!     |',
     4         /,'   | mblc must be larger or equal |',
     5         /,'   |      than meven + modd       |'
     6         /,'   |______________________________|')
  203 format (//,'    ______________________________________________ ',
     1         /,'   |                                              |',
     2         /,'   | mqa = ',i6,'  is too large                   |',
     3         /,'   | mqa must be smaller or equal than 8*mx*my*mz |',
     4         /,'   |     = ',i6,'                                 |',
     5         /,'   | recompile with a bigger mesh                 |',
     6         /,'   |  or without the equivalence statement        |',
     7         /,'   |     in the subroutine                        |'
     8         /,'   |______________________________________________|')



c................................... neven = number of even parity shells
c                                    nodd  = number of odd  parity shells
      neven = meven
      nodd  = modd
      nmax  = max(neven,nodd)
      if ((neven.ne.nodd).and.(neven.ne.nodd+1)) then
        print 201,neven,nodd
        stop " neven & nodd "
      endif
      if (mblc.le.neven+nodd) then
        print 202,neven,nodd,mblc
        stop " mblc is too small "
      endif
      if (mqa.gt.2*mq) then
        print 203,mqa,2*mq
        stop " recompile without other parameters in nilson"
      endif

c..................................................... ordering the basis
      nvec  = 0
      nblc  = 0
      ns(1) = 0

c..................................................... loop on the blocks
c                                                             even parity
      do ni=1,neven
        n    = 2*(ni-1)
        nblc = nblc + 1
        nsi(nblc,1) = ((n+2)*(n+4))/8
        nsi(nblc,2) = ((n+2)*n)/8
        nsi(nblc,3) = nsi(nblc,2)
        nsi(nblc,4) = nsi(nblc,2)
c............................................................ sub-block 1
        do i=1,ni
          nij = ni - i + 1
          do j=1,nij
            nvec = nvec + 1
            nx(nvec) = 2*(i-1)
            ny(nvec) = 2*(j-1)
            nz(nvec) = n - nx(nvec) - ny(nvec)
          enddo
        enddo
c............................................................ sub-block 2
        ni1 = ni - 1
        if (ni1.ne.0) then
          do i=1,ni1
            nij = ni - i
            do j=1,nij
              nvec = nvec + 1
              nx(nvec) = 2*(i-1) + 1
              ny(nvec) = 2*(j-1) + 1
              nz(nvec) = n - nx(nvec) - ny(nvec)
            enddo
          enddo
c............................................................ sub-block 3
          do i=1,ni1
            nij = ni - i
            do j=1,nij
              nvec = nvec + 1
              nx(nvec) = 2*(i-1) + 1
              ny(nvec) = 2*(j-1)
              nz(nvec) = n - nx(nvec) - ny(nvec)
            enddo
          enddo
c............................................................ sub-block 4
          do i=1,ni1
            nij = ni - i
             do j=1,nij
              nvec = nvec + 1
              nx(nvec) = 2*(i-1)
              ny(nvec) = 2*(j-1) + 1
              nz(nvec) = n - nx(nvec) - ny(nvec)
            enddo
          enddo
        endif
        ns(nblc+1) = ns(nblc) + ((n+1)*(n+2))/2
      enddo

c............................................................. odd parity
      do ni=1,nodd
        n    = 2*ni - 1
        nblc = nblc + 1
        nsi(nblc,1) = ((n+1)*(n+3))/8
        nsi(nblc,2) = ((n-1)*(n+1))/8
        nsi(nblc,3) = nsi(nblc,1)
        nsi(nblc,4) = nsi(nblc,1)
c............................................................ sub-block 1
        do i=1,ni
          nij = ni - i + 1
          do j=1,nij
            nvec = nvec + 1
            nx(nvec) = 2*(i-1)
            ny(nvec) = 2*(j-1)
            nz(nvec) = n - nx(nvec) - ny(nvec)
          enddo
        enddo
c............................................................ sub-block 2
        ni1 = ni - 1
        if (ni1.ne.0) then
          do i=1,ni1
            nij = ni - i
            do j=1,nij
              nvec = nvec + 1
              nx(nvec) = 2*i-1
              ny(nvec) = 2*j-1
              nz(nvec) = n - nx(nvec) - ny(nvec)
            enddo
          enddo
        endif
c............................................................ sub-block 3
        do i=1,ni
          nij = ni - i + 1
          do j=1,nij
            nvec = nvec + 1
            nx(nvec) = 2*i-1
            ny(nvec) = 2*(j-1)
            nz(nvec) = n - nx(nvec) - ny(nvec)
          enddo
        enddo
c............................................................ sub-block 4
        do i=1,ni
          nij=ni - i + 1
          do j=1,nij
            nvec = nvec + 1
            nx(nvec) = 2*(i-1)
            ny(nvec) = 2*j-1
            nz(nvec) = n - nx(nvec) - ny(nvec)
          enddo
        enddo
        ns(nblc+1) = ns(nblc) + ((n+1)*(n+2))/2
      enddo

c............................. one dimensionnal oscillator wave-functions
      xis   =     (npn-npp)
      xis   = xis/(npn+npp)
      cf(2) =  ca-cb*xis
      cf(1) = one+cb*xis
c.................................................... loop on the isospin
      print 100
      nwave = 0
      do 15 it=1,2
      nn = max(mx,my,mz)
      do ndd=1,3
        if (ndd.eq.1) xho = sqrt(hox)*dx
        if (ndd.eq.2) xho = sqrt(hoy)*dx
        if (ndd.eq.3) xho = sqrt(hoz)*dx
        xho = xho * sqrt(cf(it))
        do j=1,nmax
          nn1 = 2*j - 1
          nn2 = 2*j
          if (j.eq.1) then
            do i=1,nn
              x = xho*(i-tp5)
              he(1,i,ndd) = exp(-x*x/two)
              he(2,i,ndd) = x*he(1,i,ndd)
            enddo
          else
          do i=1,nn
              x = (xho*(i-tp5))**2
              he(nn1,i,ndd) = he(nn1-2,i,ndd)*x
              he(nn2,i,ndd) = he(nn2-2,i,ndd)*x
            enddo
            do k=1,j-1
              nnn1 = 2*k-1
              nnn2 = 2*k
              x1   = zero
              x2   = zero
              do i=1,nn
                x1 = x1 + he(nn1,i,ndd)*he(nnn1,i,ndd)
                x2 = x2 + he(nn2,i,ndd)*he(nnn2,i,ndd)
              enddo
              x1 = x1*dx*two
              x2 = x2*dx*two
              do i=1,nn
                he(nn1,i,ndd) = he(nn1,i,ndd) - x1*he(nnn1,i,ndd)
                he(nn2,i,ndd) = he(nn2,i,ndd) - x2*he(nnn2,i,ndd)
              enddo
            enddo
          endif
          x1 = zero
          x2 = zero
          do i=1,nn
            x1 = x1 + he(nn1,i,ndd)**2
            x2 = x2 + he(nn2,i,ndd)**2
          enddo
          x1 = sqrt(tp5/(dx*x1))
          x2 = sqrt(tp5/(dx*x2))
          do i=1,nn
            he(nn1,i,ndd) = x1*he(nn1,i,ndd)
            he(nn2,i,ndd) = x2*he(nn2,i,ndd)
          enddo
        enddo
      enddo

c........................ building and diagonalization of the hamiltonian
c                         storage of the eigenvectors and the eigenvalues
      ho0 = (hox*hoy*hoz)**(one/tt3)
      ax  = hox/ho0
      ay  = hoy/ho0
      az  = hoz/ho0
      ho0 = ho0*hbm(it)
      ho0 = ho0*cf(it)
      nw  = nwn
      if (it.eq.2) nw = nwp
      np  = npn
      if (it.eq.2) np = npp

c................................................................ nucleus
      x = xk(it)
      y = xmu(it)
      if (np.gt.50) x = xk(it+2)
      if (np.gt.50) y = xmu(it+2)
c..................................................... loop on the blocks
      ia = 0
      do 16 ni=1,nblc
      n  = ns(ni+1)-ns(ni)
      nn = ns(ni)
c............................................... loop on the first vector
      do 17 i=1,n
      nx1 = nx(nn+i)
      ny1 = ny(nn+i)
      nz1 = nz(nn+i)
      nb  = nx1 + ny1 + nz1
      x1  = 1 - 2*mod(nb-nz1,2)
      x2  = 1 - 2*mod(ny1,2)
c.............................................. loop on the second vector
      do 18 j=i,n
      nx2 = nx(nn+j)
      ny2 = ny(nn+j)
      nz2 = nz(nn+j)
c..................................... computation of the matrix elements
      if (j.ne.i) go to 19
      h(i,j) = ax*(nx1+tp5) + ay*(ny1+tp5) + az*(nz1+tp5)
     1        -x*y*((nb*(nb+1))/two -nx1**2 -ny1**2 -nz1**2)
      go to 18
   19 if (nx1.ne.nx2) go to 20
      h(i,j) = zero
      an = ny2
      am = ny1
      if (ny2.eq.ny1+1) h(i,j) = x*sqrt(an*nz1)
      if (ny2.eq.ny1-1) h(i,j) = x*sqrt(am*nz2)
      if (ny2.eq.ny1+2) h(i,j) =-x*y*sqrt(an*(ny2-1)*nz1*(nz1-1))
      if (ny2.eq.ny1-2) h(i,j) =-x*y*sqrt(am*(ny1-1)*nz2*(nz2-1))
      go to 18
   20 if (nz1.ne.nz2) go to 21
      h(i,j) = zero
      am = nx1
      an = nx2
      if (nx2.eq.nx1+1) h(i,j) =-x*x1*sqrt(an*ny1)
      if (nx2.eq.nx1-1) h(i,j) =-x*x1*sqrt(am*ny2)
      if (nx2.eq.nx1+2) h(i,j) =-x*y*sqrt(an*(nx2-1)*ny1*(ny1-1))
      if (nx2.eq.nx1-2) h(i,j) =-x*y*sqrt(am*(nx1-1)*ny2*(ny2-1))
      go to 18
   21 h(i,j) = zero
      if (ny1.ne.ny2) go to 18
      x3 = one-two*mod(nx1,2)
      an = nx2
      am = nx1
      x4 = one-two*mod(nx2,2)
      if (nz2.eq.nz1+1) h(i,j) =-x*x2*x3*sqrt(am*nz2)
      if (nz2.eq.nz1-1) h(i,j) =-x*x2*x4*sqrt(an*nz1)
      if (nz2.eq.nz1+2) h(i,j) = x*y*sqrt(nz2*(nz2-1)*am*(nx1-1))
      if (nz2.eq.nz1-2) h(i,j) = x*y*sqrt(nz1*(nz1-1)*an*(nx2-1))
   18 h(j,i) = h(i,j)
   17 continue
      call diagon (h,ndim,n,s,d,wd)

c.......................storage and shift of the single particle energies
      irep(ni) = ia
      do i=1,n
        do j=1,n
          ia    = ia + 1
          a(ia) = s(i,j)
        enddo
        e(nn+i) = h(i,i)*ho0 - t50
      enddo
   16 continue
      if (it.eq.1) print 101,x,y,cb
      if (it.eq.2) print 102,x,y,ca,cb

c.............. ordering the eigenvalues according to increasing energies
c              nor(i) gives the original position of the i th s.p. energy
      do i=1,nvec
        nor(i) = i
      enddo
      i1 = nvec - 1
      do i=1,i1
        ii = i + 1
        do j=ii,nvec
          if (e(j).lt.e(i)) then
            x      = e(i)
            e(i)   = e(j)
            e(j)   = x
            k      = nor(i)
            nor(i) = nor(j)
            nor(j) = k
          endif
        enddo
      enddo
      do i=1,nvec
        j = nor(i)
        k = nx(j) + ny(j) + nz(j)
        npa(i) = 1 - 2*mod(k,2)
      enddo
      print 103,(i,nor(i),e(i),npa(i),i=1,nvec)

c..................................... selection of the nw wave-functions
c       different filling is obtained by previous change of the array nor
      j = 0
      do i=1,nw
        if (npa(i).ne.-1) then
          j = j + 1
          ntrs(j) = i
        endif
      enddo
      npar(1,it) = j
      do i=1,nw
        if (npa(i).ne.1) then
          j = j + 1
          ntrs(j) = i
        endif
      enddo

      npar(2,it) = j - npar(1,it)
      do 43 iwave=1,nw
        do ix=1,mq
          psi(ix,1,1) = zero
        enddo
        nwave = nwave + 1
      if (iwave.le.npar(1,it)) go to 49
      kparz(nwave) =-1
      i = iwave - npar(1,it)
      go to 50
   49 kparz(nwave) =+1
   50 i = ntrs(iwave)
      esp1(nwave) = e(i)
      j = nor(i)
      do nn=1,nblc
        n  = ns(nn+1) - ns(nn)
        ia = irep(nn)
        do i=1,n
        do ja=1,n
          ia = ia + 1
          s(i,ja) = a(ia)
        enddo
        enddo
        nvv = nn
        if (j.le.ns(nn+1)) go to 45
      enddo
      stop 'random?'
   45 nn  = ns(nvv)
      ny2 = 0
      kk  = 0
      do 46 k=1,4
      if (nsi(nvv,k).eq.0) go to 46
      nx2 = ny2 + 1
      ny2 = ny2 + nsi(nvv,k)
      nz2 = 0
      if (k.eq.1.or.k.eq.3) nz2=3
      do i=nx2,ny2
        nx1 = nx(nn+i) + 1
        ny1 = ny(nn+i) + 1
        nz1 = nz(nn+i) + 1
        xph = s(i,j-nn)
        if (mod(ny1,4).eq.nz2) xph =-xph
        do ix=1,mx
          hex = he(nx1,ix,1)
          do iy=1,my
            hey = he(ny1,iy,2)
            do iz=1,mz
              hez = he(nz1,iz,3)
              psi(ix,iy,kk+iz) = psi(ix,iy,kk+iz) + xph*hex*hey*hez
            enddo
          enddo
        enddo
      enddo
   46 kk = kk + mz
      call scopy (mq,psi,ark(1,nwave))
   43 continue
   15 continue

      return
      end

c______________________________________________________________________________
      subroutine evolve

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0)

      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)

  100 format (//,' *****  final  *****')


c................................ initialisation of several coefficientts
c                                        initial occupation probabilities
c                       orthonormalisation and calculation of the density
c                               calculation of the potential and printout

      call inisp
      call class
      call gaphf
      call ortho
      call densit
      do i=1,mv
        wcd(i) = zero
      enddo
      print 100

      call newpot
      call figaro

      return
      end

c______________________________________________________________________________
      subroutine inisp

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4 afor

      parameter (mmx=mx+mc,mmy=my+mc,mmz=mz+mc,mmv=mmx*mmy*mmz)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (tt8=8.0d0,t12=12.0d0,t16=16.0d0,tp5=0.5d0)
      parameter (epsm3=1.0d-03)
      parameter (e2c=1.43998d0)

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q,b14,b15
     1              ,byt3a,byt3b
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz
      common /mud  / xi (mx,my,mz),yi (mx,my,mz),zi (mx,my,mz)
     1              ,xii(mx,my,mz),yii(mx,my,mz),zii(mx,my,mz)
      common /mudd / xu(mmx+mmy+mmz)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)


c............................................................. constants
c     e2c   = square of the charge of the electron (MeV*Fm)
      e2    = e2c                               ! 1.43998d0
      pi    = tt4*atan2(one,one)
      dv    = tt8*dx*dx*dx
      nwave = nwaven + nwavep

c............................................................. xu /mudd /
      x =-dx/two
      m = max0(mmx,mmy,mmz)
      do i=1,m+1
        x     = x+dx
        xu(i) = x
      enddo

c............................................................... /mud  /
      do k=1,mz
      do j=1,my
      do i=1,mx
        xi (i,j,k) = xu(i)
        xii(i,j,k) = xu(i)**2
        yi (i,j,k) = xu(j)
        yii(i,j,k) = xu(j)**2
        zi (i,j,k) = xu(k)
        zii(i,j,k) = xu(k)**2
      enddo
      enddo
      enddo

c............................................................... /spwf /
      do iw = 1,nwaven
        kiso(iw) = 1
      enddo
      do iw=1,nwavep
        kiso(iw+nwaven) = 2
      enddo

c.... coefficients of the skyrme functional......................./ener /
      b1   = t0*(one+x0/two)/two
      b2   =-t0*( x0+tp5   )/two
      b3   = (t1*(one+x1/two)+t2*(one+x2/two))/tt4
      b4   =-(t1*( x1+tp5   )-t2* (x2+tp5   ))/tt4
      b5   = (tt3*t1*(one+x1/two)-t2*(one+x2/two))/t16
      b6   =-(tt3*t1*( x1+tp5  ) +t2*( x2+tp5   ))/t16
      b7a  = t3a*(one +x3a/two)/t12
      b8a  =-t3a*( x3a+tp5    )/t12
      b7b  = t3b*(one +x3b/two)/t12
      b8b  =-t3b*( x3b+tp5    )/t12
      b9   =-wso /two
      b9q  =-wsoq/two
      b14  =-(t1*x1+t2*x2)/tt8
      b15  = (t1   -t2   )/tt8
      byt3a=yt3a
      byt3b=yt3b

c................................................................ /kfcl /
      coexv = -(tt3/pi)**(one/tt3)*e2
      epscl = epsm3/(dx**16)/(nnx*nny*nnz)/tt4
      e2eff = tt4*pi*e2

c................................................................ /cnstd/
      do i=1,mv
        d = (sqrt(xii(i,1,1)+yii(i,1,1)+zii(i,1,1))-rcut)/acut
        if (d.gt.zero) then
           dd = exp(-d)/(one+exp(-d))
        else
           dd = one/(one+exp(d))
        endif
        cutof2(i) = dd
      enddo

      return
      end

c______________________________________________________________________________
      subroutine class

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)


c......................... ordering according to single particle energies
      nof = 0
      do it=1,2
      do ib=1,2
        nvb = npar(ib,it)
        do iv=1,nvb-1
          iwa = nof + iv
          do jv=iv+1,nvb
            jwa = nof + jv
            if (esp1(jwa).lt.esp1(iwa)) then
              call sswap (mq,a(1,iwa),a(1,jwa))
              x         = esp1(iwa)
              esp1(iwa) = esp1(jwa)
              esp1(jwa) = x
              x         = v2(iwa)
              v2(iwa)   = v2(jwa)
              v2(jwa)   = x
              if (iwa.eq.ntqp) then
                ntqp = jwa
              else
                if (jwa.eq.ntqp) ntqp = iwa
              endif 
           endif
          enddo
        enddo
        nof = nof + nvb
      enddo
      enddo

      return
      end

c______________________________________________________________________________
      subroutine gaphf

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,tbig=1.0d30)

      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)


c........................................................................
      nw=nwaven+nwavep
      do i=1,nw
        eqp(i) = zero
        v2(i)  = zero
      enddo

      do 1 it=1,2
      n1=1     +(it-1)*nwaven
      n2=nwaven+(it-1)*nwavep

c............................................................ hf filling
      do i=n1,n2
        v2(i)   = zero
        delta(i)= zero
      enddo
      npnp = npn*(2-it)+npp*(it-1)
      n = npnp / 2
      do i=1,n
        x = tbig
        k = 0
        do j=n1,n2
          if ( (esp1(j).le.x).and.(v2(j).eq.zero) ) then
            x = esp1(j)
            k = j
          endif
        enddo
        v2(k) = one
      enddo
      delmax(it) = zero
      ambda(it)  = esp1(k)
      epair(it)  = zero
      eproj(it)  = zero
      do i=n1,n2
        eqp(i)=abs(esp1(i)-ambda(it))
      enddo

    1 continue

      xlamb(1)   = zero
      xlamb(2)   = zero
      ntqp       = 0
      eproj(3)   = eproj(1)  + eproj(2)
      disper(3)  = disper(1) + disper(2)

      return
      end

c______________________________________________________________________________
      subroutine ortho

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (one=1.0d0,two=2.0d0)

      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /stor / a(mv,4,mw)

c..................   1) reordering according to single particle energies
c                     2) schimdttt
c                     3) normalisation of the wave function


c........................................................................
      dvs2 = dv/two
      iwa  = 0
      do it=1,2
      do ib=1,2
        nvb = npar(ib,it)
        do iv=1,nvb
          iwa = iwa+1
          if (iv.gt.1) then
            iv1 = iv - 1
            do iv2=1,iv1
              iwb = iwa - iv + iv2
              xr  = sdot(mq,a(1,1,iwb),a(1,1,iwa))
              xr  =-xr*dvs2
              call saxpy (mq,xr,a(1,1,iwb),a(1,1,iwa))
            enddo
          endif
          psint = sdot(mq,a(1,1,iwa),a(1,1,iwa))
          fac   = one/sqrt(psint*dvs2)
          call sscal (mq,fac,a(1,1,iwa))
        enddo
      enddo
      enddo

      return
      end

c______________________________________________________________________________
      subroutine densit

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0)

      common /den  / rho(mv,2)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv),brq(mq)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)

      dimension ta(mv)


c............................................. calculation of the density
      do i=1,mv
        rho(i,1)  = zero
        rho(i,2)  = zero
        vtau(i,1) = zero
        vtau(i,2) = zero
        vdiv(i,1) = zero
        vdiv(i,2) = zero
      enddo

      do iw=1,nwave
        it = kiso(iw)
        iz = kparz(iw)
        oc = v2(iw)
        ocd= oc +oc 
        call scopy (mq,a(1,iw),w1)
        do i=1,mv
          ta(i)=w1(i)**2+w2(i)**2+w3(i)**2+w4(i)**2
        enddo
        call saxpy (mv,oc,ta,rho(1,it))
        call deriv (iz)
        call scopy (3*mq,wx1,da(1,iw))
        do i=1,mv
          vtau(i,it) = vtau(i,it) + oc*(
     1                 wx1(i)**2+wx2(i)**2+wx3(i)**2+wx4(i)**2
     2                +wy1(i)**2+wy2(i)**2+wy3(i)**2+wy4(i)**2
     3                +wz1(i)**2+wz2(i)**2+wz3(i)**2+wz4(i)**2 )
          vdiv(i,it) = vdiv(i,it) + ocd*(
     1                 wx1(i)*wy2(i)+wz4(i)*(wy1(i)+wx2(i))
     2                -wy1(i)*wx2(i)-wz3(i)*(wy2(i)-wx1(i))
     3                -wy4(i)*wx3(i)-wz1(i)*(wy4(i)+wx3(i))
     4                -wx4(i)*wz2(i)+wy3(i)*(wz2(i)+wx4(i)) )
        enddo
      enddo
 
      return
      end

c______________________________________________________________________________
      subroutine newpot

      implicit real*8 (a-h,o-z)


c........................................................................
      call mome
      call vbound
      call vcoul
      call vcal

      return
      end

c______________________________________________________________________________
      subroutine mome

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (epsm4=1.0d-4,t180=180.0d0)

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd
      common /cstn / qxnc,qxcstn,qxfinn,excstn,pentexn,
     1               qync,qycstn,qyfinn,eycstn,penteyn,
     2               qznc,qzcstn,qzfinn,ezcstn,pentezn,
     3               qrnc,qrcstn,qrfinn,ercstn,pentern,
     4               qnc0,gnc0
      common /cstp / qxpc,qxcstp,qxfinp,excstp,pentexp,
     1               qypc,qycstp,qyfinp,eycstp,penteyp,
     2               qzpc,qzcstp,qzfinp,ezcstp,pentezp,
     3               qrpc,qrcstp,qrfinp,ercstp,penterp,
     4               qpc0,gpc0
      common /cstt / qxtc,qxcstt,qxfint,excstt,pentext,
     1               qytc,qycstt,qyfint,eycstt,penteyt,
     2               qztc,qzcstt,qzfint,ezcstt,pentezt,
     3               qrtc,qrcstt,qrfint,ercstt,pentert,
     4               qtc0,gtc0
      common /ddcut/ rhoc(mv),drhoc(mv)
      common /den  / rhon(mv),rhop(mv)
      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /mud  / xi (mx,my,mz),yi (mx,my,mz),zi (mx,my,mz)
     1              ,xii(mx,my,mz),yii(mx,my,mz),zii(mx,my,mz)
      common /nxyz / dx,dv

      dimension rhonc(mv),rhopc(mv)

c     * this subroutine computes the multipoles moments of the densities
c     *  up to the order 4 (neutron, proton and total)
c     * with the proton moments, the coulomb boundary conditions
c     *  are calculated


c........................................................................
      traf = tt4*atan2(one,one)/t180
      s3   = sqrt(tt3)

c................. calculation of the density dependent cut-off functions
      if (rcut.gt.zero) then
        do i=1,mv
          rhoc(i)  = one
          drhoc(i) = one
        end do
      else
        do i=1,mv
          rhot     =(rhon(i)+rhop(i))
          frhoc    = rhot/(-rcut)
          rhoc(i)  = tanh(frhoc)
          drhoc(i) = frhoc*(one-rhoc(i)**2)+rhoc(i)
        enddo
      endif

c.................. calculation of the even moments of the proton density
      ap  = ssum(mv,rhop)
      x2p = sdot(mv,rhop,xii)
      y2p = sdot(mv,rhop,yii)
      z2p = sdot(mv,rhop,zii)
      do i=1,mv
        rhopc(i) = rhop(i)*cutof2(i)*rhoc(i)
      enddo
      x2pc = sdot(mv,rhopc,xii)
      y2pc = sdot(mv,rhopc,yii)
      z2pc = sdot(mv,rhopc,zii)
      ap   = ap *dv
      x2p  = x2p*dv
      y2p  = y2p*dv
      z2p  = z2p*dv
      qxpc = (two*x2pc-y2pc-z2pc)*dv
      qypc = (two*y2pc-x2pc-z2pc)*dv
      qzpc = (two*z2pc-x2pc-y2pc)*dv
      qrpc =    (z2pc+x2pc+y2pc)*dv
      qpc0 = sqrt((two/tt3)*(qxpc**2+qypc**2+qzpc**2))
      if ((abs(qxpc-qypc)+abs(qzpc)).le.epsm4) then
        gpc0 = zero
      else
        gpc0 = atan2(qxpc-qypc,qzpc*s3)/traf
      endif

c......................................... moments of the neutron density
      an  = ssum(mv,rhon)
      x2n = sdot(mv,rhon,xii)
      y2n = sdot(mv,rhon,yii)
      z2n = sdot(mv,rhon,zii)
      do i=1,mv
        rhonc(i) = rhon(i)*cutof2(i)*rhoc(i)
      enddo
      x2nc = sdot(mv,rhonc,xii)
      y2nc = sdot(mv,rhonc,yii)
      z2nc = sdot(mv,rhonc,zii)
      an   = an *dv
      x2n  = x2n*dv
      y2n  = y2n*dv
      z2n  = z2n*dv
      qxnc = (two*x2nc-y2nc-z2nc)*dv
      qync = (two*y2nc-x2nc-z2nc)*dv
      qznc = (two*z2nc-x2nc-y2nc)*dv
      qrnc = (z2nc+x2nc+y2nc)*dv
      qnc0 = sqrt((two/tt3)*(qxnc**2+qync**2+qznc**2))
      if ((abs(qxnc-qync)+abs(qznc)).le.epsm4) then
        gnc0 = zero
      else
        gnc0 = atan2(qxnc-qync,qznc*s3)/traf
      endif

c............................ moments of the total density (with cut-off)
      qxtc = qxpc+qxnc
      qytc = qypc+qync
      qztc = qzpc+qznc
      qrtc = qrpc+qrnc
      qtc0 = sqrt((two/tt3)*(qxtc**2+qytc**2+qztc**2))
      if ((abs(qxtc-qytc)+abs(qztc)).le.epsm4) then
        gtc0 = zero
      else
        gtc0 = atan2(qxtc-qytc,qztc*s3)/traf
      endif

      return
      end

c______________________________________________________________________________
      subroutine vbound

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (mmx=mx+mc,mmy=my+mc,mmz=mz+mc,mmv=mmx*mmy*mmz)
      parameter (one=1.0d0,two=2.0d0,tt8=8.0d0,tp375=0.375d0)

      common /bound/ xbound(mmy,mmz),ybound(mmx,mmz),zbound(mmx,mmy)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz
      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /mudd / xu(mmx+mmy+mmz)
      common /nxyz / dx,dv


c............................ calculation of the coulomb boundary values
      c1 = (two*z2p-x2p-y2p)/tt8
      c2 = tp375*(x2p-y2p)

c................................ rearrangement and multiplication by e2
      c3 =-two * (c2+c1)*e2
      c2 = two * (c2-c1)*e2
      c1 = two * two*c1 *e2
      c0 = two *(ap/two)*e2

c................................................................ xbound
      x2   = xu(nnx+1)**2
      c2x2 = c2*x2
      do j=1,nny
        y2   = xu(j)**2
        c3y2 = c3*y2
        do k=1,nnz
          z2   = xu(k)**2
          c1z2 = c1*z2
          r    = one/sqrt(x2+y2+z2)
          xbound(j,k)=c0*r+(c1z2+c2x2+c3y2)*r**5
        enddo
      enddo
c................................................................ ybound
      y2   = xu(nny+1)**2
      c3y2 = c3*y2
      do i=1,nnx
        x2   = xu(i)**2
        c2x2 = c2*x2
        do k=1,nnz
          z2   = xu(k)**2
          c1z2 = c1*z2
          r    = one/sqrt(x2+y2+z2)
          ybound(i,k)=c0*r+(c1z2+c2x2+c3y2)*r**5
        enddo
      enddo
c................................................................ zbound
      z2   = xu(nnz+1)**2
      c1z2 = c1*z2
      do i=1,nnx
        x2   = xu(i)**2
        c2x2 = c2*x2
        do j=1,nny
          y2   = xu(j)**2
          c3y2 = c3*y2
          r    = one/sqrt(x2+y2+z2)
          zbound(i,j)=c0*r+(c1z2+c2x2+c3y2)*r**5
        enddo
      enddo

      return
      end

c______________________________________________________________________________
      subroutine vcoul

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (mmx=mx+mc,mmy=my+mc,mmz=mz+mc,mmv=mmx*mmy*mmz)
      parameter (zero=0.0d0,tt6=6.0d0)
      parameter (nitmax=100+8*mc)

      common /bound/ xbound(mmy,mmz),ybound(mmx,mmz),zbound(mmx,mmy)
      common /den  / rho(mx,my,mz,2)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mx,my,my),wce(mv)
     1              ,wt3a(mv),wt3b(mv)
      common /wcl  / t(mmx,mmy,mmz),p(mmx,mmy,mmz),z(mmx,mmy,mmz)
     1              ,w(mmx,mmy,mmz),r(mmx,mmy,mmz)

  100 format (' error coulomb nitmax,epscl,zz1,zz2,zt,a,c',/,i5,6e15.8)
  102 format (4x,' coulomb iter=',i5,'  (z,z)=',e15.8)


c......................................................... initalisation
      do i=1,mmv
        z(i,1,1) = zero
        w(i,1,1) = zero
        r(i,1,1) = zero
      enddo
      do k=1,mz
      do j=1,my
      do i=1,mx
        w(i,j,k) = wcd(i,j,k)
        r(i,j,k) = rho(i,j,k,2)
      enddo
      enddo
      enddo

c............................................. calculation of delta(wcd)
      do i=1,mmv
        z(i,1,1) =-tt6*w(i,1,1)
      enddo
      do k=1,nnz
      do j=1,nny
        do i=1,nnx-1
          z(i,j,k)   = z(i,j,k)   + w(i+1,j,k)
        enddo
          z(nnx,j,k) = z(nnx,j,k) + xbound(j,k)
          z(1,j,k)   = z(1,j,k)   + w(1,j,k)
        do i=2,nnx
          z(i,j,k)   = z(i,j,k)   + w(i-1,j,k)
        enddo
      enddo
      enddo
      do k=1,nnz
      do i=1,nnx
        do j=1,nny-1
          z(i,j,k)   = z(i,j,k)   + w(i,j+1,k)
        enddo
          z(i,nny,k) = z(i,nny,k) + ybound(i,k)
          z(i,1,k)   = z(i,1,k)   + w(i,1,k)
        do j=2,nny
          z(i,j,k)   = z(i,j,k)   + w(i,j-1,k)
        enddo
      enddo
      enddo
      do j=1,nny
      do i=1,nnx
        do k=1,nnz-1
          z(i,j,k)   = z(i,j,k)   + w(i,j,k+1)
        enddo
          z(i,j,nnz) = z(i,j,nnz) + zbound(i,j)
          z(i,j,1)   = z(i,j,1)   + w(i,j,1)
        do k=2,nnz
          z(i,j,k)   = z(i,j,k)   + w(i,j,k-1)
        enddo
      enddo
      enddo

      do i=1,mmv
        z(i,1,1) = z(i,1,1) / (dx*dx)
      enddo

c................................ calculation of z=p=(source-delta(wcd))
      do i=1,mmv
        z(i,1,1)=(-e2eff)*r(i,1,1)-z(i,1,1)
        p(i,1,1)=z(i,1,1)
      enddo
      zz1=sdot(mmv,z,z)

c........................................... beginning of the iterations
c--   calculation of tk+1=delta(pk)
c--               of (z.z),ak+1,wcdk+1,zk+1 and (zk+1.zk+1)
c--           and of ck+1 and pk+1
      do niter=1,nitmax
        do i=1,mmv
          t(i,1,1) =-tt6*p(i,1,1)
        enddo
        do k=1,nnz
        do j=1,nny
          do i=1,nnx-1
            t(i,j,k) = t(i,j,k) + p(i+1,j,k)
          enddo
            t(1,j,k) = t(1,j,k) + p(  1,j,k)
          do i=2,nnx
            t(i,j,k) = t(i,j,k) + p(i-1,j,k)
          enddo
        enddo
        enddo
        do k=1,nnz
        do i=1,nnx
          do j=1,nny-1
            t(i,j,k) = t(i,j,k) + p(i,j+1,k)
          enddo
            t(i,1,k) = t(i,1,k) + p(i,  1,k)
          do j=2,nny
            t(i,j,k) = t(i,j,k) + p(i,j-1,k)
          enddo
        enddo
        enddo
        do j=1,nny
        do i=1,nnx
          do k=1,nnz-1
            t(i,j,k) = t(i,j,k) + p(i,j,k+1)
          enddo
            t(i,j,1) = t(i,j,1) + p(i,j,  1)
          do k=2,nnz
            t(i,j,k) = t(i,j,k) + p(i,j,k-1)
          enddo
        enddo
        enddo
        do i=1,mmv
          t(i,1,1) = t(i,1,1) / (dx*dx)
        enddo
        zt = sdot(mmv,z,t)
        a  = zz1/zt
        call saxpy (mmv,a,p,w)
        call saxpy (mmv,-a,t,z)
        zz2 = sdot(mmv,z,z)
        if (zz2.lt.epscl) then
          print 102,niter,zz2
          go to 111
        endif
        c   = zz2/zz1
        zz1 = zz2
        do i=1,mmv
          p(i,1,1)=z(i,1,1)+c*p(i,1,1)
        enddo
      enddo

      print 100,nitmax,epscl,zz1,zz2,zt,a,c

 111  do k=1,mz
      do j=1,my
      do i=1,mx
        wcd(i,j,k) = w(i,j,k)
      enddo
      enddo
      enddo

      return
      end

c______________________________________________________________________________
      subroutine vcal

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (one=1.0d0,tt3=3.0d0)

      common /den  / rho(mv,2)
      common /drho / drhon(mv),drhop(mv)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q,b14,b15
     1              ,byt3a,byt3b
      common /evcen/ r(mq,2)
      common /evpro/ rx(mv,2),ry(mv,2),rz(mv,2)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wtmp / ta(mq)


c..............................................................................
      s13= one/tt3
      call lapla (rho(1,1),drhon,one,one,one)
      call lapla (rho(1,2),drhop,one,one,one)
      do i=1,mv
        wce(i)  = rho(i,2)**s13*coexv
        wt3a(i) = (rho(i,1)+rho(i,2))**byt3a
        wt3b(i) = (rho(i,1)+rho(i,2))**byt3b
      enddo

      return
      end

c______________________________________________________________________________
      subroutine figaro

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,two=2.0d0,tbig=1.0d+10)

      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q,b14,b15
     1              ,byt3a,byt3b
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu,ecm(3)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      dimension etsp(3),idh(mw)

  100 format (/,' ',78('*'))
  101 format (/,' neutron levels')
  102 format (/,'  proton levels')
  103 format (/,'   n',2x,'par',2x,'norm',4x,'eqp',5x,'esp',7x,'jz')
  104 format (2i4,4f8.3)
  105 format (/,' ',78('.'))
  109 format (/,' total energy (from functional)  ',f12.3)


c..............................................................................
      call fprte
      call fprtj
      call fprtz
      
      print 100

      do it=1,2
        if (it.eq.1) then
          na = 1
          nb = nwaven
          print 101
        else
          na = nwaven +1
          nb = nwaven + nwavep
          print 102
        endif
        print 103
        etsp(it) = zero
        do iw=na,nb
          idh(iw) = 0
        enddo
        do iw=na,nb
          eref = tbig
          do jw=na,nb
            if (idh(jw).eq.0) then
              if (esp1(jw).le.eref) then
                eref = esp1(jw)
                kw   = jw
              endif
            endif
          enddo
          idh(kw)  = 1
          iwa    = kw
          etsp(it) = etsp(it)+two*esp1(iwa)*v2(iwa)
          print 104,iwa,kparz(iwa),v2(iwa),eqp(iwa),esp1(iwa),ajzd(iwa)
        enddo
      enddo

      print 105
      call fprtm

      ek(3) = ek(2)+ek(1)
      es    = eso(1)+eso(2)
      etot  = esk+ecoul+ecoex+ek(3)+es

      print 105
      print 109,etot

      return
      end

c______________________________________________________________________________
      subroutine fprte

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,tp5=0.5d0,tp75=0.75d0)

      common /den  / rhon(mv),rhop(mv)
      common /drho / drhon(mv),drhop(mv)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q,b14,b15
     1              ,byt3a,byt3b
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu,ecm(3)
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)


c...................................... calculation of some h.f. energies
      ytt   = zero
      ynp   = zero
      ydtt  = zero
      ydnp  = zero
      y3att = zero
      y3anp = zero
      y3btt = zero
      y3bnp = zero
      ycd   = zero
      yce   = zero
      do i=1,mv
        xn    = rhon(i)
        xp    = rhop(i)
        xt    = xn+xp
        xtt   = xt**2
        xnp   = xn**2+xp**2
        dxn   = drhon(i)
        dxp   = drhop(i)
        ytt   = ytt+xtt
        ynp   = ynp  + xnp
        ydtt  = ydtt + xt*(dxn+dxp)
        ydnp  = ydnp + (xn*dxn+xp*dxp)
        y3att = y3att + xtt*wt3a(i)
        y3anp = y3anp + xnp*wt3a(i)
        y3btt = y3btt + xtt*wt3b(i)
        y3bnp = y3bnp + xnp*wt3b(i)
        ycd   = ycd  + xp*wcd(i)
        yce   = yce  + xp*wce(i)
      enddo
      et0   = dv * (b1*ytt  + b2*ynp)
      etd   =-dv * (b5*ydtt + b6*ydnp)
      et3a  = dv * (b7a*y3att + b8a*y3anp)
      et3b  = dv * (b7b*y3btt + b8b*y3bnp)
      ecoul = dv * ycd * tp5
      ecoex = dv * yce * tp75

      return
      end

c______________________________________________________________________________
      subroutine fprtj

      implicit real*8 (a-h,o-z)
      character*4 afor
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)

      common /den  / rho(mv,2)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q,b14,b15
     1              ,byt3a,byt3b
      common /force/ tsk(14),hbar,hbm(2),xm(3),afor
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu,ecm(3)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /stor / a(mq,mw)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv)
     2              ,p1(mv),p2(mv),p3(mv),p4(mv)


c..............................................................................
      ek(1)  = zero
      ek(2)  = zero
      do iwa=1,nwave
        it = kiso(iwa)
        zp = kparz(iwa)
        oc = v2(iwa)
        call scopy (mq,a(1,iwa),w1)
        call lapla (w1,p1, one, one, zp)
        call lapla (w2,p2,-one,-one, zp)
        call lapla (w3,p3,-one, one,-zp)
        call lapla (w4,p4, one,-one,-zp)
        ek(it) = ek(it) + oc * sdot(mq,w1,p1)
      enddo
      ek(1) = ek(1) * (-hbm(1)/two) * (one-xm(1)/xm(3)) * dv
      ek(2) = ek(2) * (-hbm(2)/two) * (one-xm(2)/xm(3)) * dv

      ett    = zero
      eso(1) = zero
      eso(2) = zero
      do i=1,mv
        xn  = rho(i,1)
        xp  = rho(i,2)
        xt  = xn + xp
        vtn = vtau(i,1)
        vtp = vtau(i,2)
        vtt = vtn+vtp
        ett = ett + b3*xt*vtt + b4*(xn*vtn+xp*vtp)
        eso(1) = eso(1) + (b9*xt+b9q*xn) * vdiv(i,1)
        eso(2) = eso(2) + (b9*xt+b9q*xp) * vdiv(i,2)
      enddo
      ett    = dv * ett
      eso(1) = dv * eso(1)
      eso(2) = dv * eso(2)
      esk    = et0+ett+etd+et3a+et3b

      return
      end

c______________________________________________________________________________
      subroutine fprtz

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,two=2.0d0,tt4=4.0d0)

      common /mud  / xi(mv),yi(mv),zi(mv),xii(mv),yii(mv),zii(mv)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv),brik(mq)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     2              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     3              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /wtmp / wjz1(mv),wjz2(mv),wjz3(mv),wjz4(mv)


c........................ this subroutine computes <Jz> for each s.p.w.f.
      do iwa=1,nwave
        ajzd(iwa)=zero
        call scopy (  mq, a(1,iwa),w1 )
        call scopy (3*mq,da(1,iwa),wx1)
        do i=1,mv
          wjz1(i) = w1(i) + two*(+xi(i)*wy2(i)-yi(i)*wx2(i))
          wjz2(i) = w2(i) + two*(-xi(i)*wy1(i)+yi(i)*wx1(i))
          wjz3(i) =-w3(i) + two*(+xi(i)*wy4(i)-yi(i)*wx4(i))
          wjz4(i) =-w4(i) + two*(-xi(i)*wy3(i)+yi(i)*wx3(i))
        enddo
        ajzd(iwa) = dv*sdot(mq,w1,wjz1)/tt4
      enddo

      return
      end

c______________________________________________________________________________
      subroutine fprtm

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (epsm4=1.0d-4,t180=180.0d0)

      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n

  101 format (' mass ',f12.6,'  charge ',f12.6,'  n ',f12.6)
  102 format (' rt=',f8.4,'Fm   rp=',f8.4,'Fm   rn=',f8.4,'Fm')
  110 format (4x,'q',9x,'x(Fm2)',6x,'y(Fm2)',6x,'z(Fm2)',5x,'q0(Fm2)',
     1        5x,'gamma(dg)')
  111 format (' neutron ',5f12.3)
  112 format (' proton  ',5f12.3)
  113 format (' total   ',5f12.3)


c........................................................................
      traf = tt4*atan2(one,one)/t180
      s3   = sqrt(tt3)

      at = ap + an
      print 101,at,ap,an

      rp = x2p+y2p+z2p
      rn = x2n+y2n+z2n
      rt = sqrt((rn+rp)/at)
      rn = sqrt(rn/an)
      rp = sqrt(rp/ap)
      print 102,rt,rp,rn

      qxp = two*x2p-y2p-z2p
      qyp = two*y2p-z2p-x2p
      qzp = two*z2p-x2p-y2p
      qp0 = sqrt((two/tt3)*(qxp**2+qyp**2+qzp**2))
      if ((abs(qxp-qyp)+abs(qzp)).le.epsm4) then
        gp0 = zero
      else
        gp0 = atan2(qxp-qyp,qzp*s3)/traf
      endif
      qxn = two*x2n-y2n-z2n
      qyn = two*y2n-z2n-x2n
      qzn = two*z2n-x2n-y2n
      qn0 = sqrt((two/tt3)*(qxn**2+qyn**2+qzn**2))
      if ((abs(qxn-qyn)+abs(qzn)).le.epsm4) then
        gn0 = zero
      else
        gn0 = atan2(qxn-qyn,qzn*s3)/traf
      endif
      qxt = qxp+qxn
      qyt = qyp+qyn
      qzt = qzp+qzn
      qt0 = sqrt((two/tt3)*(qxt**2+qyt**2+qzt**2))
      if ((abs(qxt-qyt)+abs(qzt)).le.epsm4) then
        gt0 = zero
      else
        gt0 = atan2(qxt-qyt,qzt*s3)/traf
      endif
      print 110
      print 111,qxn,qyn,qzn,qn0,gn0
      print 112,qxp,qyp,qzp,qp0,gp0
      print 113,qxt,qyt,qzt,qt0,gt0

      return
      end

c______________________________________________________________________________
      subroutine wrini

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,tt5=5.0d0)
      parameter (tp5=0.5d0,tp16=0.16d0,t1250=1250.0d0)

      common /champ/ qcha(12)
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd
      common /cstn / wcstn(22)
      common /cstp / wcstp(22)
      common /cstt / wcstt(22)
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)


c........................................................................
      do i=1,12
        qcha(i) = zero
      enddo

      ral     = zero
      epscst  = zero
      cqr     = zero
      cq2     = zero
      imtd    = 0

      do i = 1,22
        wcstn(i) = zero
        wcstp(i) = zero
        wcstt(i) = zero
      enddo

      delq    = zero
      q1n     = zero
      q1p     = zero
      q1t     = zero
      q2n     = zero
      q2p     = zero
      q2t     = zero

      itert   = 0

      npair   = 0
      gn      = t1250
      gp      = t1250
      encut   = tt5
      epcut   = tt5
      dcut    = tp5
      xcut    = one
      alpha   = one/tp16
      alphap  = one/tp16
      ntqp    = 0

      do iw = 1,mw
        v22(iw) = v2(iw)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine writ8

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4 afor,head

      common /champ/ qxxn,qyyn,qzzn,qrrn
     1              ,qxxp,qyyp,qzzp,qrrp
     2              ,qxxt,qyyt,qzzt,qrrt
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd
      common /cstn / qxnc,qxcstn,qxfinn,excstn,pentexn
     1              ,qync,qycstn,qyfinn,eycstn,penteyn
     2              ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     3              ,qrnc,qrcstn,qrfinn,ercstn,pentern
     4              ,qnc0,gnc0
      common /cstp / qxpc,qxcstp,qxfinp,excstp,pentexp
     1              ,qypc,qycstp,qyfinp,eycstp,penteyp
     2              ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     3              ,qrpc,qrcstp,qrfinp,ercstp,penterp
     4              ,qpc0,gpc0
      common /cstt / qxtc,qxcstt,qxfint,excstt,pentext
     1              ,qytc,qycstt,qyfint,eycstt,penteyt
     2              ,qztc,qzcstt,qzfint,ezcstt,pentezt
     3              ,qrtc,qrcstt,qrfint,ercstt,pentert
     4              ,qtc0,gtc0
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /den  / rhon(mx,my,mz),rhop(mx,my,mz)
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,njmunu,ncm2,nmass
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ifor
      common /pairn/ vn,vp,rangen,rangep
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / p1(mv),p2(mv),p3(mv),p4(mv)
     1              ,w1(mv),w2(mv),w3(mv),w4(mv)


c........................................................................
      nwave = nwaven+nwavep

      if (npair.le.5) iver = 5
      if (npair.ge.6) iver = 6
      write (13) iver
      write (13) head
      write (13) nwaven,nwavep,npn,npp
      write (13) nx,ny,nz

      write (13) itert,dx
      write (13) imtd,ral,epscst,cqr,cq2,rcut,acut
     1          ,delq,q1n,q1p,q1t,q2n,q2p,q2t
      write (13) qxnc,qxcstn,qxfinn,excstn,pentexn
     1          ,qync,qycstn,qyfinn,eycstn,penteyn
     1          ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     1          ,qrnc,qrcstn,qrfinn,ercstn,pentern
     1          ,qnc0,gnc0
      write (13) qxpc,qxcstp,qxfinp,excstp,pentexp
     1          ,qypc,qycstp,qyfinp,eycstp,penteyp
     1          ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     1          ,qrpc,qrcstp,qrfinp,ercstp,penterp
     1          ,qpc0,gpc0
      write (13) qxtc,qxcstt,qxfint,excstt,pentext
     1          ,qytc,qycstt,qyfint,eycstt,penteyt
     1          ,qztc,qzcstt,qzfint,ezcstt,pentezt
     1          ,qrtc,qrcstt,qrfint,ercstt,pentert
     1          ,qtc0,gtc0
      write (13) qxxn,qyyn,qzzn,qrrn
     1          ,qxxp,qyyp,qzzp,qrrp
     1          ,qxxt,qyyt,qzzt,qrrt

      write (13) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                            ,t3b,x3b,yt3b,wso,wsoq

      if (iver.eq.5) then
        write (13) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,ntqp
      else
        write (13) npair,vn,vp,rangen,rangep,ambda,xlamb,ntqp
      endif

      write (13) (eqp(i),i=1,nwave)
      write (13) (delta(i),i=1,nwave)
      write (13) (kparz(i),i=1,nwave)
      write (13) (esp1(i),i=1,nwave)
      write (13) (v2(i),i=1,nwave)
      write (13) (v22(i),i=1,nwave)
      write (13) ((npar(i,it),i=1,2),it=1,2)

      write (13) rhon
      write (13) rhop
      write (13) vtau
      write (13) vdiv

      do iwave=1,nwave
        call scopy (mq,a(1,iwave),w1)
        write (13) w1
        write (13) w2
        write (13) w3
        write (13) w4
      enddo

      return
      end

c______________________________________________________________________________
      subroutine diagon (a,ndim,n,v,d,wd)

c..............................................................................
c  diagonalization of a real symmetric matrix (a(i,j)                         .
c     input : a  n*n matrix           (with declared dimensions ndim*ndim)    .
c             a(i,k) * v(k,j) = v(i,k) * d(k)                                 .
c     output: v block of eigenvectors (with declared dimensions ndim*ndim)    .
c             d eigenvalues in ascending order                                .
c             wd working array                                                .
c             the content of a is lost (actualy, a(i,i) = d(i)                .
c             a(i,j) = v(i,k) * d(k) * v(j,k)                                 .
c                      v is an orthogonal matrix : v(i,k)*v(j,k) = delta_ij   .
c..............................................................................

      implicit real*8 (a-h,o-z)

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (eps=9.0d-12,epsd=1.0d-16,tol=1.0d-36,jstop=30)
      dimension a(ndim,ndim),v(ndim,ndim),d(ndim),wd(ndim)


c........................................................................
      v(1,1) = one
      d(1)   = a(1,1)
      if (n.le.1) return

      do i=1,ndim*ndim
        v(i,1) = a(i,1)
      enddo

      do 9 ii=2,n
        i     = n+2-ii
        l     = i-1
        h     = zero
        scale = zero
        if (l.eq.1) go to 100
        do k=1,l
          scale = scale + abs(v(i,k))
        enddo
        if (scale.gt.tol) go to 3
  100   continue
        wd(i) = v(i,l)
        d(i) = h
        go to 9
    3   do k=1,l
          v(i,k) = v(i,k)/scale
          h      = h + v(i,k)**2
        enddo
        f      = v(i,l)
        g      =-sign(sqrt(h),f)
        wd(i)  = g*scale
        h      = h-f*g
        v(i,l) = f-g
        f      = zero
        do j=1,l
          v(j,i) = v(i,j)/(h*scale)
          g = zero
          do k=1,j
            g = g + v(j,k)*v(i,k)
          enddo
          j1 = j + 1
          if (j1.le.l) then
            do k=j1,l
              g = g + v(k,j)*v(i,k)
            enddo
          endif
          wd(j) = g/h
          f     = f + wd(j)*v(i,j)
        enddo
        hh = f/(h+h)
        do j=1,l
          f     = v(i,j)
          g     = wd(j) - hh*f
          wd(j) = g
          do k=1,j
            v(j,k) = v(j,k) - f*wd(k) - g*v(i,k)
          enddo
        enddo
        do k=1,l
          v(i,k) = scale * v(i,k)
        enddo
        d(i)=h
    9 continue

      d(1)  = zero
      wd(1) = zero
      do i=1,n
        l=i-1
        if ((abs(d(i)).ge.epsd).and.(l.ne.0)) then
          do j=1,l
            g = zero
            do k=1,l
              g = g + v(i,k)*v(k,j)
            enddo
            do k=1,l
              v(k,j) = v(k,j) - g*v(k,i)
            enddo
          enddo
        endif
        d(i)   = v(i,i)
        v(i,i) = one
        if (l.ne.0) then
          do j=1,l
            v(i,j) = zero
            v(j,i) = zero
          enddo
        endif
      enddo

      do i=2,n
        wd(i-1) = wd(i)
      enddo

      wd(n) = zero
      b     = zero
      f     = zero
      do 212 l=1,n
        j = 0
        h = eps * ( abs(d(l)) + abs(wd(l)) )
        if (b.lt.h) b = h
        m = l - 1
  202   m = m + 1
        if (m.gt.n) go to 203
        if (abs(wd(m))-b) 203,203,202
  203   continue
        if (m.eq.l) go to 211
  204   continue
        if (j.eq.jstop) stop 'diagon jstop'
        j = j + 1
        p = (d(l+1)-d(l))/(two*wd(l))
        r = sqrt(p*p+one)
        h =  d(l)-wd(l)/(p+sign(r,p))
        do i=l,n
          d(i) = d(i) - h
        enddo
        f  = f+h
        p  = d(m)
        c  = one
        s  = zero
        m1 = m  - 1
        ml = m1 + l
        do ii=l,m1
          i = ml - ii
          g = c*wd(i)
          h = c*p
          if (abs(p).ge.abs(wd(i))) then
            c       = wd(i)/p
            r       = sqrt(c*c+one)
            wd(i+1) = s*p*r
            s       = c/r
            c       = one/r
          else
            c       = p/wd(i)
            r       = sqrt(c*c+one)
            wd(i+1) = s*wd(i)*r
            s       = one/r
            c       = c/r
          endif
          p      = c*d(i) - s*g
          d(i+1) = h + s*(c*g+s*d(i))
          do k=1,n
            h        = v(k,i+1)
            v(k,i+1) = s*v(k,i) + c*h
            v(k,i)   = c*v(k,i) - s*h
          enddo
        enddo
        wd(l) = s*p
        d(l)  = c*p
        if (abs(wd(l))-b) 211,211,204
  211   continue
        d(l) = d(l)+f
  212 continue

      n1 = n-1
      do i=1,n1
        k  = i
        p  = d(i)
        ii = i + 1
        do j=ii,n
          if (d(j).lt.p) then
            k = j
            p = d(j)
          endif
        enddo
        if (k.ne.i) then
          d(k) = d(i)
          d(i) = p
          do j=1,n
            p      = v(j,i)
            v(j,i) = v(j,k)
            v(j,k) = p
          enddo
        endif
      enddo

      do i=1,n
        a(i,i) = d(i)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine deriv (iz)

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /wave / ps1(mv),ps2(mv),ps3(mv),ps4(mv),brik(mq)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)


c..............................................................................
      call der (ps1,wx1,wy1,wz1, 1, 1, iz)
      call der (ps2,wx2,wy2,wz2,-1,-1, iz)
      call der (ps3,wx3,wy3,wz3,-1, 1,-iz)
      call der (ps4,wx4,wy4,wz4, 1,-1,-iz)

      return
      end

c______________________________________________________________________________
      subroutine der (w,wx,wy,wz,ix,iy,iz)

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (ca=45.0d0,cb=-9.0d0,cc=60.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),wx(mx,my,mz),wy(mx,my,mz),wz(mx,my,mz)


c.......................................................... x derivative
      if (ix.gt.0) then
        do j=1,my*mz
          wx(1,j,1)    = ca*(w(2,j,1)   -w(1,j,1)   )
     1                  +cb*(w(3,j,1)   -w(2,j,1)   )
     2                     +(w(4,j,1)   -w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)   -w(1,j,1)   )
     1                  +cb*(w(4,j,1)   -w(1,j,1)   )
     2                     +(w(5,j,1)   -w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)   -w(2,j,1)   )
     1                  +cb*(w(5,j,1)   -w(1,j,1)   )
     2                     +(w(6,j,1)   -w(1,j,1)   )
          wx(mx-2,j,1) = ca*(w(mx-1,j,1)-w(mx-3,j,1))
     1                  +cb*(w(mx,j,1)  -w(mx-4,j,1))
     2                     +(           -w(mx-5,j,1))
          wx(mx-1,j,1) = ca*(w(mx,j,1)  -w(mx-2,j,1))
     1                  +cb*(           -w(mx-3,j,1))
     2                     +(           -w(mx-4,j,1))
          wx(mx,j,1)   = ca*(           -w(mx-1,j,1))
     1                  +cb*(           -w(mx-2,j,1))
     2                     +(           -w(mx-3,j,1))
        enddo
      else
        do j=1,my*mz
          wx(1,j,1)    = ca*(w(2,j,1)   +w(1,j,1)   )
     1                  +cb*(w(3,j,1)   +w(2,j,1)   )
     2                     +(w(4,j,1)   +w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)   -w(1,j,1)   )
     1                  +cb*(w(4,j,1)   +w(1,j,1)   )
     2                     +(w(5,j,1)   +w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)   -w(2,j,1)   )
     1                  +cb*(w(5,j,1)   -w(1,j,1)   )
     2                     +(w(6,j,1)   +w(1,j,1)   )
          wx(mx-2,j,1) = ca*(w(mx-1,j,1)-w(mx-3,j,1))
     1                  +cb*(w(mx,j,1)  -w(mx-4,j,1))
     2                     +(           -w(mx-5,j,1))
          wx(mx-1,j,1) = ca*(w(mx,j,1)  -w(mx-2,j,1))
     1                  +cb*(           -w(mx-3,j,1))
     2                     +(           -w(mx-4,j,1))
          wx(mx,j,1)   = ca*(           -w(mx-1,j,1))
     1                  +cb*(           -w(mx-2,j,1))
     2                     +(           -w(mx-3,j,1))
        enddo
      endif

      do i=4,mx-3
      do j=1,my*mz
        wx(i,j,1)    = ca*(w(i+1,j,1) -w(i-1,j,1) )
     1                +cb*(w(i+2,j,1) -w(i-2,j,1) )
     2                   +(w(i+3,j,1) -w(i-3,j,1) )
      enddo
      enddo

c.......................................................... y derivative
      if (iy.gt.0) then
        do k=1,mz
        do i=1,mx
          wy(i,1,k)    = ca*(w(i,2,k)   -w(i,1,k)   )
     1                  +cb*(w(i,3,k)   -w(i,2,k)   )
     2                     +(w(i,4,k)   -w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)   -w(i,1,k)   )
     1                  +cb*(w(i,4,k)   -w(i,1,k)   )
     2                     +(w(i,5,k)   -w(i,2,k)   ) 
          wy(i,3,k)    = ca*(w(i,4,k)   -w(i,2,k)   )
     1                  +cb*(w(i,5,k)   -w(i,1,k)   )
     2                     +(w(i,6,k)   -w(i,1,k)   )
          wy(i,my-2,k) = ca*(w(i,my-1,k)-w(i,my-3,k))
     1                  +cb*(w(i,my,k)  -w(i,my-4,k))
     2                     +(           -w(i,my-5,k))
          wy(i,my-1,k) = ca*(w(i,my,k)  -w(i,my-2,k))
     1                  +cb*(           -w(i,my-3,k))
     2                     +(           -w(i,my-4,k))
          wy(i,my,k)   = ca*(           -w(i,my-1,k))
     1                  +cb*(           -w(i,my-2,k))
     2                     +(           -w(i,my-3,k))
        enddo
        enddo
      else
        do k=1,mz
        do i=1,mx
          wy(i,1,k)    = ca*(w(i,2,k)   +w(i,1,k)   )
     1                  +cb*(w(i,3,k)   +w(i,2,k)   )
     2                     +(w(i,4,k)   +w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)   -w(i,1,k)   )
     1                  +cb*(w(i,4,k)   +w(i,1,k)   )
     2                     +(w(i,5,k)   +w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)   -w(i,2,k)   )
     1                  +cb*(w(i,5,k)   -w(i,1,k)   )
     2                     +(w(i,6,k)   +w(i,1,k)   )
          wy(i,my-2,k) = ca*(w(i,my-1,k)-w(i,my-3,k))
     1                  +cb*(w(i,my,k)  -w(i,my-4,k))
     2                     +(           -w(i,my-5,k))
          wy(i,my-1,k) = ca*(w(i,my,k)  -w(i,my-2,k))
     1                  +cb*(           -w(i,my-3,k))
     2                     +(           -w(i,my-4,k))
          wy(i,my,k)   = ca*(           -w(i,my-1,k))
     1                  +cb*(           -w(i,my-2,k))
     2                     +(           -w(i,my-3,k))
        enddo
        enddo
      endif

      do j=4,my-3
      do i=1,mx
      do k=1,mz
        wy(i,j,k)    = ca*(w(i,j+1,k) -w(i,j-1,k) )
     1                +cb*(w(i,j+2,k) -w(i,j-2,k) )
     2                   +(w(i,j+3,k) -w(i,j-3,k) )
      enddo
      enddo
      enddo

c.......................................................... z derivative
      if (iz.gt.0) then
      do i=1,mx*my
        wz(i,1,1)    = ca*(w(i,1,2)   -w(i,1,1)   )
     1                +cb*(w(i,1,3)   -w(i,1,2)   )
     2                   +(w(i,1,4)   -w(i,1,3)   )
        wz(i,1,2)    = ca*(w(i,1,3)   -w(i,1,1)   )
     1                +cb*(w(i,1,4)   -w(i,1,1)   )
     2                   +(w(i,1,5)   -w(i,1,2)   )
        wz(i,1,3)    = ca*(w(i,1,4)   -w(i,1,2)   )
     1                +cb*(w(i,1,5)   -w(i,1,1)   )
     2                   +(w(i,1,6)   -w(i,1,1)   )
        wz(i,1,mz-2) = ca*(w(i,1,mz-1)-w(i,1,mz-3))
     1                +cb*(w(i,1,mz)  -w(i,1,mz-4))
     2                   +(           -w(i,1,mz-5))
        wz(i,1,mz-1) = ca*(w(i,1,mz)  -w(i,1,mz-2))
     1                +cb*(           -w(i,1,mz-3))
     1                   +(           -w(i,1,mz-4))
        wz(i,1,mz)   = ca*(           -w(i,1,mz-1))
     1                +cb*(           -w(i,1,mz-2))
     2                   +(           -w(i,1,mz-3))
      enddo
      else
      do i=1,mx*my
        wz(i,1,1)    = ca*(w(i,1,2)   +w(i,1,1)   )
     1                +cb*(w(i,1,3)   +w(i,1,2)   )
     2                   +(w(i,1,4)   +w(i,1,3)   )
        wz(i,1,2)    = ca*(w(i,1,3)   -w(i,1,1)   )
     1                +cb*(w(i,1,4)   +w(i,1,1)   )
     2                   +(w(i,1,5)   +w(i,1,2)   )
        wz(i,1,3)    = ca*(w(i,1,4)   -w(i,1,2)   )
     1                +cb*(w(i,1,5)   -w(i,1,1)   )
     2                   +(w(i,1,6)   +w(i,1,1)   )
        wz(i,1,mz-2) = ca*(w(i,1,mz-1)-w(i,1,mz-3))
     1                +cb*(w(i,1,mz)  -w(i,1,mz-4))
     2                   +(           -w(i,1,mz-5))
        wz(i,1,mz-1) = ca*(w(i,1,mz)  -w(i,1,mz-2))
     1                +cb*(           -w(i,1,mz-3))
     1                   +(           -w(i,1,mz-4))
        wz(i,1,mz)   = ca*(           -w(i,1,mz-1))
     1                +cb*(           -w(i,1,mz-2))
     2                   +(           -w(i,1,mz-3))
      enddo
      endif

      do k=4,mz-3
      do i=1,mx*my
        wz(i,1,k)    = ca*(w(i,1,k+1) -w(i,1,k-1) )
     1                +cb*(w(i,1,k+2) -w(i,1,k-2) )
     2                   +(w(i,1,k+3) -w(i,1,k-3) )
      enddo
      enddo

      do i=1,mv
        wx(i,1,1) = wx(i,1,1) / (cc*dx)
        wy(i,1,1) = wy(i,1,1) / (cc*dx)
        wz(i,1,1) = wz(i,1,1) / (cc*dx)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine lapla (w,dw,xp,yp,zp)

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (ca=-8064.0d0,cb=1008.0d0,cc=-128.0d0,cd=9.0d0)
      parameter (ce=5040.0d0,cf=43050.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),dw(mx,my,mz)


c................................................ calculation of delta(w)
c                      for a component of the wave function of parity zp

      do i=1,mv
        dw(i,1,1) = cf * w(i,1,1)
      enddo
c                                                                 x sweep
      do k=1,mz
      do j=1,my
      dw(1,j,k)    = dw(1,j,k)   +ca*w(1,j,k)*xp+cb*w(2,j,k)*xp
     1                           +cc*w(3,j,k)*xp+cd*w(4,j,k)*xp
      dw(2,j,k)    = dw(2,j,k)   +ca*w(1,j,k)   +cb*w(1,j,k)*xp
     1                           +cc*w(2,j,k)*xp+cd*w(3,j,k)*xp
      dw(3,j,k)    = dw(3,j,k)   +ca*w(2,j,k)   +cb*w(1,j,k)
     1                           +cc*w(1,j,k)*xp+cd*w(2,j,k)*xp
      dw(4,j,k)    = dw(4,j,k)   +ca*w(3,j,k)   +cb*w(2,j,k)
     1                           +cc*w(1,j,k)   +cd*w(1,j,k)*xp
      dw(mx-1,j,k) = dw(mx-1,j,k)+ca*w(mx,j,k)
      dw(mx-2,j,k) = dw(mx-2,j,k)+ca*w(mx-1,j,k)+cb*w(mx,j,k)
      dw(mx-3,j,k) = dw(mx-3,j,k)+ca*w(mx-2,j,k)+cb*w(mx-1,j,k)
     1                           +cc*w(mx,j,k)     
        do i=1,mx-4
         dw(i,j,k) = dw(i,j,k)   +ca*w(i+1,j,k) +cb*w(i+2,j,k)
     1                           +cc*w(i+3,j,k) +cd*w(i+4,j,k)
        enddo
        do i=5,mx
         dw(i,j,k) = dw(i,j,k)   +ca*w(i-1,j,k) +cb*w(i-2,j,k)
     1                           +cc*w(i-3,j,k) +cd*w(i-4,j,k)
          enddo
      enddo
      enddo
c                                                                 y sweep
      do k=1,mz
      do i=1,mx
      dw(i,1,k)    = dw(i,1,k)   +ca*w(i,1,k)*yp+cb*w(i,2,k)*yp
     1                           +cc*w(i,3,k)*yp+cd*w(i,4,k)*yp
      dw(i,2,k)    = dw(i,2,k)   +ca*w(i,1,k)   +cb*w(i,1,k)*yp
     1                           +cc*w(i,2,k)*yp+cd*w(i,3,k)*yp
      dw(i,3,k)    = dw(i,3,k)   +ca*w(i,2,k)   +cb*w(i,1,k)
     1                           +cc*w(i,1,k)*yp+cd*w(i,2,k)*yp
      dw(i,4,k)    = dw(i,4,k)   +ca*w(i,3,k)   +cb*w(i,2,k)
     1                           +cc*w(i,1,k)   +cd*w(i,1,k)*yp
      dw(i,my-1,k) = dw(i,my-1,k)+ca*w(i,my,k)
      dw(i,my-2,k) = dw(i,my-2,k)+ca*w(i,my-1,k)+cb*w(i,my,k)
      dw(i,my-3,k) = dw(i,my-3,k)+ca*w(i,my-2,k)+cb*w(i,my-1,k)
     1                           +cc*w(i,my,k)     
        do j=1,my-4
         dw(i,j,k) = dw(i,j,k)   +ca*w(i,j+1,k) +cb*w(i,j+2,k)
     1                           +cc*w(i,j+3,k) +cd*w(i,j+4,k)
        enddo
        do j=5,my
         dw(i,j,k) = dw(i,j,k)   +ca*w(i,j-1,k) +cb*w(i,j-2,k)
     1                           +cc*w(i,j-3,k) +cd*w(i,j-4,k)
        enddo
      enddo
      enddo
c                                                                 z sweep
      do i=1,mx*my
      dw(i,1,1)    = dw(i,1,1)   +ca*w(i,1,1)*zp+cb*w(i,1,2)*zp
     1                           +cc*w(i,1,3)*zp+cd*w(i,1,4)*zp
      dw(i,1,2)    = dw(i,1,2)   +ca*w(i,1,1)   +cb*w(i,1,1)*zp
     1                           +cc*w(i,1,2)*zp+cd*w(i,1,3)*zp
      dw(i,1,3)    = dw(i,1,3)   +ca*w(i,1,2)   +cb*w(i,1,1)
     1                           +cc*w(i,1,1)*zp+cd*w(i,1,2)*zp
      dw(i,1,4)    = dw(i,1,4)   +ca*w(i,1,3)   +cb*w(i,1,2)
     1                           +cc*w(i,1,1)   +cd*w(i,1,1)*zp
      dw(i,1,mz-1) = dw(i,1,mz-1)+ca*w(i,1,mz)
      dw(i,1,mz-2) = dw(i,1,mz-2)+ca*w(i,1,mz-1)+cb*w(i,1,mz)
      dw(i,1,mz-3) = dw(i,1,mz-3)+ca*w(i,1,mz-2)+cb*w(i,1,mz-1)
     1                           +cc*w(i,1,mz)     
        do k=1,mz-4
         dw(i,1,k) = dw(i,1,k)   +ca*w(i,1,k+1) +cb*w(i,1,k+2)
     1                           +cc*w(i,1,k+3) +cd*w(i,1,k+4)
        enddo
        do k=5,mz
         dw(i,1,k) = dw(i,1,k)   +ca*w(i,1,k-1) +cb*w(i,1,k-2)
     1                           +cc*w(i,1,k-3) +cd*w(i,1,k-4)
        enddo
       enddo

      do i=1,mv
        dw(i,1,1) = - dw(i,1,1) / (ce*dx*dx)
      enddo

      return
      end

c______________________________________________________________________________
      function ssum (n,a)

      double precision zero,a(n),ssum
      parameter (zero=0.0d0)

      ssum = zero
      do i=1,n
        ssum = ssum + a(i)
      enddo

      return
      end

c______________________________________________________________________________
      function sdot (n,a,b)

      double precision zero,a(n),b(n),sdot
      parameter (zero=0.0d0)

      sdot = zero
      do i=1,n
        sdot = sdot + a(i) * b(i)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine sswap (n,a,b)

      double precision a(n),b(n),tmp

      do i=1,n
        tmp  = b(i)
        b(i) = a(i)
        a(i) = tmp
      enddo

      return
      end

c______________________________________________________________________________
      subroutine scopy (n,a,b)

      double precision a(n),b(n)

      do i=1, n
        b(i) = a(i)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine sscal (n,fac,a)

      double precision a(n),fac

      do i=1, n
        a(i) = fac * a(i)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine saxpy (n,fac,a,b)

      double precision a(n),b(n),fac

      do i=1, n
        b(i) = b(i) + fac * a(i)
      enddo

      return
      end

c______________________________________________________________________________
