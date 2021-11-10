c______________________________________________________________________________
      program int8

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

c    ***********************************************************************
c    *                                                                     *
c    * Interpolation program with the Fourier method                       *
c    *           for a left-right symmetric nucleus                        *
c    *           the fonctions to be interpolated are read from fort.12    *
c    *           the interpolated functions are stored on fort.13          *
c    * The interpolated vtau and vdiv are set to 0.0                       *
c    *                                                                     *
c    ***********************************************************************


c    data  *********************
c          *  head        a4   *
c          *  nx,ny,nz   3i5   *
c          *  dx         f8.3  *
c          *********************


c..............................................................................
      implicit real*8 (a-h,o-z)
      character*4 head
      character*2 num
      character*20 endian
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt4=4.0d0,tt8=8.0d0)
      parameter (tp5=0.5d0,epsm3=1.0d-3)

      common /den  / rho(mx,my,mz,2)
      common /info / bidon(500),head(20),irb,nx1,ny1,nz1
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / p1(mx,my,mz),p2(mx,my,mz),p3(mx,my,mz),p4(mx,my,mz)
     1              ,w1(mx,my,mz),w2(mx,my,mz),w3(mx,my,mz),w4(mx,my,mz)

      dimension cpx(mx,mx),cpy(my,my),cpz(mz,mz),
     1          cmx(mx,mx),cmy(my,my),cmz(mz,mz)


  101 format (20a4)
  102 format (3i5)
  103 format (f8.3)

  200 format (/,'  __________________________________________________ ',
     1        /,' |                                                  |',
     2        /,' |  Program int8     Release 1.0.4                  |',
     3        /,' |  10 March 2014                                   |',
     4        /,' |  Copyright  P. Bonche, H. Flocard, P.H. Heenen   |',
     5        /,' |    in case of trouble, take it easy...           |',
     5        /,' |    and (please) let us know                      |',
     7        /,' |__________________________________________________|')
  201 format (//,' parameters: mx,my,mz ',3i5,
     2         /,'                   mw ',i5)
  203 format (//,' ',78('_'),/,' run information',//,20a4)
  210 format  (/,' ',78('*'),/,'  new parameters')
  211 format                  ('  old values')
  212 format (15x,'nx =',i3,'  ny =',i3,'  nz =',i3,'  dx =',f8.3)
  213 format (/,13x,'number of wave-functions =',i3)
  216 format (//,30X,' THE END',/,' ',78('*'))


c..............................................................................
      pi = tt4*atan2(one,one)
      print 200
      print 201,mx,my,mz,mw

      do i=1,mq
        p1(i,1,1) = zero
        w1(i,1,1) = zero
      enddo

c.................................................... read input wave-functions
      irb    = 12
      iprint =  1
      write (num,'(i2.2)') irb
      
      !When using pgf90 or compilers that don't accept this:
      ! comment the following:
      
      !Checking for the existence of the file and what kind of endianness it is.
      inquire(UNIT=irb, convert=endian, iostat=io) !Comment this!
      
      if(io.ne.0) call stp('I/O error for fort.12 .')  !Comment this!
      
      
      select case(endian)  !Comment this!
      
      case('BIG_ENDIAN')   !Comment this!
            open  (irb,form='unformatted',file='fort.'//num,!Comment this!
     1             convert='big_endian')  
      case default  !Comment this!
            open  (irb,form='unformatted',file='fort.'//num)!Don't comment this!
      end select  !Comment this!
      
      
      call lec8 (iprint)
      close (irb)

      do i=1,mv
        rho(i,1,1,1) = zero
        rho(i,1,1,2) = zero
        vtau(i,1)    = zero
        vtau(i,2)    = zero
        vdiv(i,1)    = zero
        vdiv(i,2)    = zero
      enddo

      dx1   = dx
      nwave = nwaven + nwavep
      read  101,head
      read  102,nx,ny,nz
      read  103,dx
      dv = tt8 *dx**3
      print 203,head
      print 210
      print 212,nx,ny,nz,dx
      print 211
      print 212,nx1,ny1,nz1,dx1
      print 213,nwave
      if (nx.gt.mx) call stp (' === nx > mx ! === ')
      if (ny.gt.my) call stp (' === ny > my ! === ')
      if (nz.gt.mz) call stp (' === nz > mz ! === ')

      dh  = dx/dx1
      fac = tp5/nx1
      ph  = tp5*pi/nx1
      x1  =-tp5*dh
      do i=1,nx
        x1 = x1+dh
        x2 =-tp5
        do j=1,nx1
          x2 = x2 + one
          if (abs(x1-x2).le.epsm3) then
            c= pi/ph
          else
            c  = sin(pi*(x1-x2))/sin(ph*(x1-x2))
          endif
          if (abs((x1+x2)/(2*nx1)-one).le.epsm3) then
            d= pi/ph
          else
            d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
          endif
          cmx(i,j) = fac*(c - d)
          cpx(i,j) = fac*(c + d)
        enddo
      enddo

      fac = tp5/ny1
      ph  = tp5*pi/ny1
      x1  =-tp5*dh
      do i=1,ny
        x1 = x1+dh
        x2 =-tp5
        do j=1,ny1
          x2 = x2 + one
          if (abs(x1-x2).le.epsm3) then
            c= pi/ph
          else
            c  = sin(pi*(x1-x2))/sin(ph*(x1-x2))
          endif
          if (abs((x1+x2)/(2*ny1)-one).le.epsm3) then
            d= pi/ph
          else
            d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
          endif
          cmy(i,j) = fac*(c - d)
          cpy(i,j) = fac*(c + d)
        enddo
      enddo

      fac = tp5/nz1
      ph  = tp5*pi/nz1
      x1  =-tp5*dh
      do i=1,nz
        x1 = x1+dh
        x2 =-tp5
        do j=1,nz1
          x2 = x2 + one
          if (abs(x1-x2).le.epsm3) then
            c= pi/ph
          else
            c  = sin(pi*(x1-x2))/sin(ph*(x1-x2))
          endif
          if (abs((x1+x2)/(2*nz1)-one).le.epsm3) then
            d= pi/ph
          else
            d  = sin(pi*(x1+x2))/sin(ph*(x1+x2))
          endif
          cmz(i,j) = fac*(c - d)
          cpz(i,j) = fac*(c + d)
        enddo
      enddo

      do 11 iwa=1,nwave
      it = 1
      if (iwa.gt.nwaven) it = 2
      call scopy (mq,a(1,iwa),w1)

      do n=1,nz1
      do m=1,ny1
      do i=1,nx
        a1 = zero
        a2 = zero
        a3 = zero
        a4 = zero
        do l=1,nx1
          a1 = a1 + cpx(i,l)*w1(l,m,n)
          a2 = a2 + cmx(i,l)*w2(l,m,n)
          a3 = a3 + cmx(i,l)*w3(l,m,n)
          a4 = a4 + cpx(i,l)*w4(l,m,n)
        enddo
        p1(i,m,n) = a1
        p2(i,m,n) = a2
        p3(i,m,n) = a3
        p4(i,m,n) = a4
      enddo
      enddo
      enddo

      do n=1,nz1
      do j=1,ny
      do i=1,nx
        a1 = zero
        a2 = zero
        a3 = zero
        a4 = zero
        do m=1,ny1
          a1 = a1 + cpy(j,m)*p1(i,m,n)
          a2 = a2 + cmy(j,m)*p2(i,m,n)
          a3 = a3 + cpy(j,m)*p3(i,m,n)
          a4 = a4 + cmy(j,m)*p4(i,m,n)
        enddo
        w1(i,j,n) = a1
        w2(i,j,n) = a2
        w3(i,j,n) = a3
        w4(i,j,n) = a4
      enddo
      enddo
      enddo

      if (kparz(iwa).ge.0 ) then
        do k=1,nz
        do j=1,ny
        do i=1,nx
          a1 = zero
          a2 = zero
          a3 = zero
          a4 = zero
          do n=1,nz1
            a1 = a1 + cpz(k,n)*w1(i,j,n)
            a2 = a2 + cpz(k,n)*w2(i,j,n)
            a3 = a3 + cmz(k,n)*w3(i,j,n)
            a4 = a4 + cmz(k,n)*w4(i,j,n)
          enddo
          p1(i,j,k) = a1
          p2(i,j,k) = a2
          p3(i,j,k) = a3
          p4(i,j,k) = a4
        enddo
        enddo
        enddo
      else
        do k=1,nz
        do j=1,ny
        do i=1,nx
          a1 = zero
          a2 = zero
          a3 = zero
          a4 = zero
          do n=1,nz1
            a1 = a1 + cmz(k,n)*w1(i,j,n)
            a2 = a2 + cmz(k,n)*w2(i,j,n)
            a3 = a3 + cpz(k,n)*w3(i,j,n)
            a4 = a4 + cpz(k,n)*w4(i,j,n)
          enddo
          p1(i,j,k) = a1
          p2(i,j,k) = a2
          p3(i,j,k) = a3
          p4(i,j,k) = a4
        enddo
        enddo
        enddo
      endif

c................................................................ normalisation
c                                          and calculation of the new densities
      x = zero
      do k=1,nz
      do j=1,ny
      do i=1,nx
        x = x + p1(i,j,k)**2+p2(i,j,k)**2+p3(i,j,k)**2+p4(i,j,k)**2
      enddo
      enddo
      enddo

      x = sqrt(two/(x*dv))
      do k=1,nz
      do j=1,ny
      do i=1,nx
        p1(i,j,k) = x*p1(i,j,k)
        p2(i,j,k) = x*p2(i,j,k)
        p3(i,j,k) = x*p3(i,j,k)
        p4(i,j,k) = x*p4(i,j,k)
        ta = p1(i,j,k)**2+p2(i,j,k)**2+p3(i,j,k)**2+p4(i,j,k)**2
        rho(i,j,k,it) = rho(i,j,k,it) + v2(iwa)*ta
      enddo
      enddo
      enddo

      do i=1,mq
        w1(i,1,1) = zero
      enddo

      do k=1,nz
      do j=1,ny
      do i=1,nx
        w1(i,j,k) = p1(i,j,k)
        w2(i,j,k) = p2(i,j,k)
        w3(i,j,k) = p3(i,j,k)
        w4(i,j,k) = p4(i,j,k)
      enddo
      enddo
      enddo

      call scopy (mq,w1,a(1,iwa))

   11 continue

      print 216

      nx1 = nx
      ny1 = ny
      nz1 = nz
      dx1 = dx

      open  (13,form='unformatted',file='fort.13')
      call writ8
      close (13)

      call stp (' **** Bye ! *** ')
      end

c______________________________________________________________________________
      subroutine lec8 (iprint)

c..............................................................................
c     single-particle wave functions and other information read from fort.irb .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character(4):: afor
      character(4):: head
      character(4):: headd

      parameter (zero=0.0d0)

      common /champ/ qxxn,qyyn,qzzn,qrrn
     1              ,qxxp,qyyp,qzzp,qrrp
     2              ,qxxt,qyyt,qzzt,qrrt
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd,imtg
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
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
     4              ,q2fin,g2fin,q2cst,qtc0,gtc0
      common /den  / rhon(mx,my,mz),rhop(mx,my,mz)
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /enert/ c14,c15,t14,t15
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pairn/ vn,vp,rangen,rangep
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)




      common /stor / a(mq,mw)
      common /taudj/ vtau(mx,my,mz,2),vdiv(mx,my,mz,2)
      common /vo8  / delqst,q1tst,q2tst,qrfintst
     1                     ,q1nst,q2nst,qrfinnst
     2                     ,q1pst,q2pst,qrfinpst
      common /wave / p1(mv),p2(mv),p3(mv),p4(mv),w1(mx,my,mz)
     1              ,w2(mx,my,mz),w3(mx,my,mz),w4(mx,my,mz)
      dimension headd(20)

c..............................................................................
  201 format (/,' ',78('_'),/,' tape information (version ',I1,')',/,
     1        /,' tape parameters:',
     1         /,' mx=', i5,' my=', i5
     1          ,' mz=',i5,
     3         /,' Total number of wavefunctions: ', i5  )
  202 format   (20a4)
  203 format   (' Neutron wavefunctions', i5, ' Proton wavefunctions',i5,
     1        /,' mass =',i5,'  charge =',i5,'  n =',i5)
  211 format (/,' nx=',i5,10x,'  ny=',i5,10x,'  nz=',i5,
     1        /,' dx=',f8.3,'Fm  dt=',f8.3,' 10-22s  itert =',i5)
  218 format (/,' constraints data',
     1        /,' cr2           = ',f12.5,
     2        /,' delq iq1, iq2 = ',f12.5,2i5,
     3        /,' cq2           = ',f12.5,
     4        /,' rcut          = ',f12.5)
  259 format   (/,' Pairing on fort.12: ')
  260 format   (' npair  = ',i5,' -- hartree-fock --')
  261 format   (' npair  = ',i5,' -- BCS seniority force --')
  262 format   (' npair  = ',i5,' -- BCS constant delta --')
  263 format   (' npair  = ',i5,' -- BCS+LN  seniority force --')
  264 format   (' npair  = ',i5,' -- BCS delta force --')
  265 format   (' npair  = ',i5,' -- BCS+LN delta force --')
  266 format   (' npair  = ',i5,' -- Finite range pairing --')
  267 format   (' npair  = ',i5,' -- Finite range pairing + LN--')
  271 format   ('  gn =',f8.3,'/(11+n),  encut=',f6.3,
     2        /,'  gp =',f8.3,'/(11+z),  epcut=',f6.3,'  dcut =',f6.3)
  272 format   ('  deltan =',f8.3,'  encut =',f6.3,
     1        /,'  deltap =',f8.3,'  epcut =',f6.3,'  dcut =',f6.3)
  273 format   ('  vn = ',f8.3,'  encut =',f6.3,
     2        /,'  vp = ',f8.3,'  epcut =',f6.3,'  dcut =',f6.3,
     4          '  rhoc =',f6.3)
  275 format   ('    vn = ',f12.6,'  range =',f6.3,
     2        /,'    vp = ',f12.6,'  range =',f6.3)
  281 format   (' cutoff only above the fermi level ')
  282 format   (' cutoff above and below the fermi level ')
  283 format   ('  ntqp = ',i3,' -- BCS ground state --')
  284 format   ('  ntqp = ',i3,' -- 2-qp state --')
  285 format   ('  ntqp = ',i3,' -- 1-qp filling approximation --')
  290 format   (/, 'Constraint data from fort.12: ', / , ' imtd  = ',i5,
     1          ' (0,2)=fixed,(1,3)=readjusted,(0,1)=T=0,(2,3)=T=1')
  293 format   (' delq  = ',e12.4,
     3        /,' q1t   = ',e12.4,'  q2t   = ',e12.4)
  294 format   (' actual isoscalar constraints',
     1        /,' cqr = ',e12.4,'  qrfint = ',e12.4,'  qrcstt = ',e12.4,
     2        /,' cq2 = ',e12.4,'  qxfint = ',e12.4,'  qxcstt = ',e12.4,
     3        /,'       ',12x,  '  qyfint = ',e12.4,'  qycstt = ',e12.4,
     4        /,'       ',12x,  '  qzfint = ',e12.4,'  qzcstt = ',e12.4)
  295 format   (' delq  = ',e15.8,
     3        /,' q1n  = ',e12.4,'  q2n   = ',e12.4,
     4        /,' q1p  = ',e12.4,'  q2p   = ',e12.4)
  296 format   (' actual isovector constraints',
     1        /,' cqr = ',e12.4,'  qrfinn = ',e12.4,'  qrcstn = ',e12.4,
     2        /,'                  qrfinp = ',e12.4,'  qrcstp = ',e12.4,
     3        /,' cq2 = ',e12.4,'  qxfinn = ',e12.4,'  qxcstn = ',e12.4,
     4        /,'       ',12x,  '  qyfinn = ',e12.4,'  qycstn = ',e12.4,
     5        /,'       ',12x,  '  qzfinn = ',e12.4,'  qzcstn = ',e12.4,
     6        /,'       ',12x,  '  qxfinp = ',e12.4,'  qxcstp = ',e12.4,
     7        /,'       ',12x,  '  qyfinp = ',e12.4,'  qycstp = ',e12.4,
     8        /,'       ',12x,  '  qzfinp = ',e12.4,'  qzcstp = ',e12.4)
  240 format   (/, 'Force used on fort.12: ',/,
     1          '  t0 =',f13.6,' x0 =',f13.6,/,
     1          '  t1 =',f13.6,' x1 =',f13.6,/,
     2          '  t2 =',f13.6,' x2 =',f13.6,/,
     3          '  t3a=',f13.6,' x3a=',f13.6,'  sigma =',f10.6,/,
     4          '  t3b=',f13.6,' x3b=',f13.6,'  sigmb =',f10.6,/,
     5          '  w  =',f13.6,' wq =',f13.6,/,
     6          '  te =',f13.6,' to =',f13.6)
  302 format ( /,' Functional parameters used on fort.12')
  311 format (  '  b1  =',f13.6,' b2  =',f12.6)
  312 format (  '  b3  =',f13.6,' b4  =',f12.6)
  313 format (  '  b5  =',f13.6,' b6  =',f12.6)
  314 format (  '  b7  =',f13.6,' b8  =',f12.6,'  sigma =',f10.6)
  315 format (  '  b7a =',f13.6,' b8a =',f12.6,'  sigma =',f10.6,/,
     1          '  b7b =',f13.6,' b8b =',f12.6,'  sigmb =',f10.6)
  316 format (  '  b9  =',f13.6,' b9q =',f12.6)
  317 format (  '  c14 =',f13.6,' c15 =',f12.6)
  318 format (  '  t14 =',f13.6,' t15 =',f12.6)
  319 format (  '  b14 =',f13.6,' b15 =',f12.6)
  320 format (  '  b16 =',f13.6,' b17 =',f12.6)

c..............................................................................
      do i=1,mq
        w1(i,1,1) = zero
        p1(i)     = zero
      enddo

      read  (irb) iver
      read  (irb) headd
      read  (irb) nwaven,nwavep,npn,npp
      read  (irb) nnx,nny,nnz
      nwave = nwaven+nwavep
      if (iprint.eq.1) print 201,iver,nnx,nny,nnz,nwave
                       print 202,headd
      if (iprint.eq.1) print 203,nwaven,nwavep,npn+npp,npp,npn

      if (nwave.gt.mw) then
	call stp (
     1     ' Number of wavefunctions on fort.12 exceeds parameter mw!')
	endif
c     ........................ adapt mesh dimensions if read-in mesh is smaller
c     if    (nx.ne.mx) call stp (' mx !')
c     if    (ny.ne.my) call stp (' my !')
c     if    (nz.ne.mz) call stp (' mz !')
      if   (nnx.gt.mx) call stp (' nnx > mx !')
      if   (nny.gt.my) call stp (' nny > my !')
      if   (nnz.gt.mz) call stp (' nnz > mz !')
      nx = nnx
      ny = nny
      nz = nnz
      if (nnx.ne.mx) nx = mx
      if (nny.ne.my) ny = my
      if (nnz.ne.mz) nz = mz

      if ((iver.lt.2).or.(iver.gt.6)) call stp (' iver !')

      if (iver.eq.2) then
        read (irb) itert,dx
        read (irb) imtd,ral,epscst,cqr,cq2,rcut,acut,
     1             delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst,
     1             qxnc,qxcstn,qxfinn,excstn,pentexn,
     1             qync,qycstn,qyfinn,eycstn,penteyn,
     1             qznc,qzcstn,qzfinn,ezcstn,pentezn,
     1             qrnc,qrcstn,qrfinnst,ercstn,pentern,
     1             qnc0,gnc0,
     1             qxpc,qxcstp,qxfinp,excstp,pentexp,
     1             qypc,qycstp,qyfinp,eycstp,penteyp,
     1             qzpc,qzcstp,qzfinp,ezcstp,pentezp,
     1             qrpc,qrcstp,qrfinpst,ercstp,penterp,
     1             qpc0,gpc0,
     1             qxtc,qxcstt,qxfint,excstt,pentext,
     1             qytc,qycstt,qyfint,eycstt,penteyt,
     1             qztc,qzcstt,qzfint,ezcstt,pentezt,
     1             qrtc,qrcstt,qrfintst,ercstt,pentert,
     1             q2cst,gtc0,
     1             qxxn,qyyn,qzzn,qrrn,
     1             qxxp,qyyp,qzzp,qrrp,
     1             qxxt,qyyt,qzzt,qrrt
        read (irb) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,wso,
     1             npair,gn,gp,encut,epcut,dcut,
     1             ambda,delmax,xlamb,ntqp,(eqp(i),i=1,nwave),
     1             (delta(i),i=1,nwave)
        read (irb) bidon
      endif

      if (iver.eq.3) then
        read (irb) itert,dx
        read (irb) imtd,ral,epscst,cqr,cq2,rcut,acut
     1            ,delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst
        read (irb) qxnc,qxcstn,qxfinn,excstn,pentexn
     1            ,qync,qycstn,qyfinn,eycstn,penteyn
     1            ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     1            ,qrnc,qrcstn,qrfinnst,ercstn,pentern
     1            ,qnc0,gnc0
        read (irb) qxpc,qxcstp,qxfinp,excstp,pentexp
     1            ,qypc,qycstp,qyfinp,eycstp,penteyp
     1            ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     1            ,qrpc,qrcstp,qrfinpst,ercstp,penterp
     1            ,qpc0,gpc0
        read (irb) qxtc,qxcstt,qxfint,excstt,pentext
     1            ,qytc,qycstt,qyfint,eycstt,penteyt
     1            ,qztc,qzcstt,qzfint,ezcstt,pentezt
     1            ,qrtc,qrcstt,qrfintst,ercstt,pentert
     1            ,q2cst,gtc0
        read (irb) qxxn,qyyn,qzzn,qrrn
     1            ,qxxp,qyyp,qzzp,qrrp
     1            ,qxxt,qyyt,qzzt,qrrt
        read (irb) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,wso,wsoq
        read (irb) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)
      endif

      if (iver.eq.4) then
        read (irb) itert,dx
        read (irb) imtd,ral,epscst,cqr,cq2,rcut,acut
     1            ,delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst
        read (irb) qxnc,qxcstn,qxfinn,excstn,pentexn
     1            ,qync,qycstn,qyfinn,eycstn,penteyn
     1            ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     1            ,qrnc,qrcstn,qrfinnst,ercstn,pentern
     1            ,qnc0,gnc0
        read (irb) qxpc,qxcstp,qxfinp,excstp,pentexp
     1            ,qypc,qycstp,qyfinp,eycstp,penteyp
     1            ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     1            ,qrpc,qrcstp,qrfinpst,ercstp,penterp
     1            ,qpc0,gpc0
        read (irb) qxtc,qxcstt,qxfint,excstt,pentext
     1            ,qytc,qycstt,qyfint,eycstt,penteyt
     1            ,qztc,qzcstt,qzfint,ezcstt,pentezt
     1            ,qrtc,qrcstt,qrfintst,ercstt,pentert
     1            ,q2cst,gtc0
        read (irb) qxxn,qyyn,qzzn,qrrn
     1            ,qxxp,qyyp,qzzp,qrrp
     1            ,qxxt,qyyt,qzzt,qrrt
        read (irb) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a,wso,wsoq
        read (irb) npair,vn,vp,rangen,rangep,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)
      endif

      if (iver.eq.5) then
        read (irb) itert,dx
        read (irb) imtd,ral,epscst,cqr,cq2,rcut,acut
     1            ,delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst
        read (irb) qxnc,qxcstn,qxfinn,excstn,pentexn
     1            ,qync,qycstn,qyfinn,eycstn,penteyn
     1            ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     1            ,qrnc,qrcstn,qrfinnst,ercstn,pentern
     1            ,qnc0,gnc0
        read (irb) qxpc,qxcstp,qxfinp,excstp,pentexp
     1            ,qypc,qycstp,qyfinp,eycstp,penteyp
     1            ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     1            ,qrpc,qrcstp,qrfinpst,ercstp,penterp
     1            ,qpc0,gpc0
        read (irb) qxtc,qxcstt,qxfint,excstt,pentext
     1            ,qytc,qycstt,qyfint,eycstt,penteyt
     1            ,qztc,qzcstt,qzfint,ezcstt,pentezt
     1            ,qrtc,qrcstt,qrfintst,ercstt,pentert
     1            ,q2cst,gtc0
        read (irb) qxxn,qyyn,qzzn,qrrn
     1            ,qxxp,qyyp,qzzp,qrrp
     1            ,qxxt,qyyt,qzzt,qrrt
        read (irb, iostat=io) t0,x0,t1,x1,t2,x2,t3a,
     1  x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
        read (irb) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)
      endif

      if(iver.eq.6) then
      	read (irb) itert,dx
        read (irb) imtd,ral,epscst,cqr,cq2,rcut,acut
     1            ,delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst
        read (irb) qxnc,qxcstn,qxfinn,excstn,pentexn
     1            ,qync,qycstn,qyfinn,eycstn,penteyn
     1            ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     1            ,qrnc,qrcstn,qrfinnst,ercstn,pentern
     1            ,qnc0,gnc0
        read (irb) qxpc,qxcstp,qxfinp,excstp,pentexp
     1            ,qypc,qycstp,qyfinp,eycstp,penteyp
     1            ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     1            ,qrpc,qrcstp,qrfinpst,ercstp,penterp
     1            ,qpc0,gpc0
        read (irb) qxtc,qxcstt,qxfint,excstt,pentext
     1            ,qytc,qycstt,qyfint,eycstt,penteyt
     1            ,qztc,qzcstt,qzfint,ezcstt,pentezt
     1            ,qrtc,qrcstt,qrfintst,ercstt,pentert
     1            ,q2cst,gtc0
        read (irb) qxxn,qyyn,qzzn,qrrn
     1            ,qxxp,qyyp,qzzp,qrrp
     1            ,qxxt,qyyt,qzzt,qrrt
c    Note that this line reads two numbers more than iver=5 version:
c    te and to!
        read (irb) t0,x0,t1,x1,t2,x2,t3a,
     1  x3a,yt3a,t3b,x3b,yt3b,wso,wsoq,te,to

        print *, 'read to=', to

c    ..................................reading functional coefficients
c                                      this is a line that was added
c                                      with respect to iver=5    
        read (irb, iostat=io) b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1          ,b14,b15,b16,b17,byt3a,byt3b,c14,c15,t14,t15
c    ..................................reading force options
c                                      this is a line that was added
c                                      with respect to iver=5
	    read (irb, iostat=io) nfunc,njmunu,ncm2,nmass,ndd,ncoex
        read (irb) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)    
  	  endif

      if (iver.eq.20) then
        read (irb) itert,dx
        read (irb) imtd,ral,epscst,cqr,cq2,rcut,acut
     1            ,delqst,q1nst,q1pst,q1tst,q2nst,q2pst,q2tst
        read (irb) qxnc,qxcstn,qxfinn,excstn,pentexn
     1            ,qync,qycstn,qyfinn,eycstn,penteyn
     1            ,qznc,qzcstn,qzfinn,ezcstn,pentezn
     1            ,qrnc,qrcstn,qrfinnst,ercstn,pentern
     1            ,qnc0,gnc0
        read (irb) qxpc,qxcstp,qxfinp,excstp,pentexp
     1            ,qypc,qycstp,qyfinp,eycstp,penteyp
     1            ,qzpc,qzcstp,qzfinp,ezcstp,pentezp
     1            ,qrpc,qrcstp,qrfinpst,ercstp,penterp
     1            ,qpc0,gpc0
        read (irb) qxtc,qxcstt,qxfint,excstt,pentext
     1            ,qytc,qycstt,qyfint,eycstt,penteyt
     1            ,qztc,qzcstt,qzfint,ezcstt,pentezt
     1            ,qrtc,qrcstt,qrfintst,ercstt,pentert
     1            ,q2cst,gtc0
        read (irb) qxxn,qyyn,qzzn,qrrn
     1            ,qxxp,qyyp,qzzp,qrrp
     1            ,qxxt,qyyt,qzzt,qrrt
        read (irb) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                              ,t3b,x3b,yt3b
     2                              ,te,to,wso,wsoq

        read (irb) npair,gn,gp,rangen,rangep
     1            ,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)
      endif

      if (iprint.ne.0) then
        print 211,nx,ny,nz,dx,dt,itert

        if(iver.eq.6) then
c...............................................Printing force coefficients

	        print 240,t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                           ,t3b,x3b,yt3b,
     2                            wso,wsoq,te,to

c............................................Printing functional parameters        	
            print 302
            print 311,b1,b2
            print 312,b3,b4
            print 313,b5,b6
            if (ndd.eq.0) then
              print 314,b7a,b8a,byt3a
            endif
            if (ndd.eq.1) then
              print 315,b7a,b8a,byt3a,b7b,b8b,byt3b
            endif
            print 316,b9,b9q
            if (nfunc.eq.0) then
              print 317,c14,c15
              print 318,t14,t15
            endif
            print 319,b14,b15
            print 320,b16,b17
        endif

        if(iver.ge.5) then
c         ................ overriding numbers that are not properly initialised
          if (cqr .eq. zero ) then
            qrcstt = zero
            qrcstp = zero
            qrcstn = zero
          endif
          if (cq2.eq.0.d0) then
            qxcstn = 0.d0
            qycstn = 0.d0
            qzcstn = 0.d0
            qxcstp = 0.d0
            qycstp = 0.d0
            qzcstp = 0.d0
            qxcstt = 0.d0
            qycstt = 0.d0
            qzcstt = 0.d0
          endif

          print 290,imtd
          if (imtd.le.1) then
            print 293,delqst,q1tst,q2tst
            print 294,cqr,qrfint,qrcstt,
     1                cq2,qxfint,qxcstt,qyfint,qycstt,qzfint,qzcstt
          else
            print 295,delqst,q1nst,q1pst,q2nst,q2pst
            print 296,cqr,qrfinn,qrcstn,qrfinp,qrcstp,
     1                cq2,qxfinn,qxcstn,qyfinn,qycstn,qzfinn,qzcstn,
     2                    qxfinp,qxcstp,qyfinp,qycstp,qzfinp,qzcstp
          endif

        endif

        !Printing pairing information
        print 259
        if (npair.eq.0) print 260,npair
        if (npair.eq.1) print 261,npair
        if (npair.eq.2) print 262,npair
        if (npair.eq.3) print 263,npair
        if (npair.eq.4) print 264,npair
        if (npair.eq.5) print 265,npair
        if (npair.eq.6) print 266,npair
        if (npair.eq.7) print 267,npair

        if (npair.eq.1) print 271,gn,encut,gp,encut,dcut
        if (npair.eq.2) print 272,delmax(1),encut,delmax(2),epcut,dcut
        if (npair.eq.3) print 271,gn,encut,gp,encut,dcut
        if ((npair.eq.4).or.(npair.eq.5)) then
          rhoc = alpha
          if (alpha.ne.zero) rhoc = 1.d0/alpha
          print 273,gn,encut,gp,encut,dcut,rhoc
        endif

        if ((npair.ge.1).and.(npair.le.5)) then
          if (xcut.eq.zero) print 281
          if (xcut.ne.zero) print 282
        endif

        if ((npair.eq.6).or.(npair.eq.7)) then
          print 275,vn,rangen,vp,rangep
        endif

        if (npair.ge.1) then
          if (ntqp.eq.0) print 283,ntqp
          if (ntqp.gt.0) print 284,ntqp
          if (ntqp.lt.0) print 285,ntqp
        endif
      endif

      read (irb) (kparz(i),i=1,nwave)
      read (irb) (esp1 (i),i=1,nwave)
      read (irb) (v2   (i),i=1,nwave)

      if (iver.ge.3) then
        read (irb) (v22(i),i=1,nwave)
      endif

      read (irb) ((npar(i,it),i=1,2),it=1,2)

      if (nnx.eq.mx. and. nny.ne.my .and. nnz.ne.mz) then
c       ........................ read wave function file using the same boxsize
        read (irb) rhon
        read (irb) rhop
        if (iver.ge.3) then
          read (irb) vtau
          read (irb) vdiv
        endif
        do iwave=1,nwave
          read (irb) w1
          read (irb) w2
          read (irb) w3
          read (irb) w4
          call scopy (mq,w1,a(1,iwave))
        enddo
      else
c       ......................... read wave function file using smaller boxsize
        print '(/," read:  nnx = ",i5," nny = ",i5," nnz = ",i5)',
     1         nnx,nny,nnz
        print '(  " here:  mx  = ",i5," my  = ",i5," mz  = ",i5,/)',
     1         mx,my,mz
        rhon(:,:,:)   = 0.0d0
        rhop(:,:,:)   = 0.0d0
        vtau(:,:,:,:) = 0.0d0
        vdiv(:,:,:,:) = 0.0d0
        read (irb) (((rhon(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
        read (irb) (((rhop(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
        read (irb) ((((vtau(i,j,k,it),i=1,nnx),j=1,nny),k=1,nnz),it=1,2)
        read (irb) ((((vdiv(i,j,k,it),i=1,nnx),j=1,nny),k=1,nnz),it=1,2)

c       ddv = 8.0d0*dx*dx*dx
c       testn = ddv * ssum(mv,rhon)
c       testp = ddv * ssum(mv,rhop)
c       print '(" test N Z = ",2d16.6)',testn,testp

        w1(:,:,:) = 0.d0
        w2(:,:,:) = 0.d0
        w3(:,:,:) = 0.d0
        w4(:,:,:) = 0.d0
        do iwave=1,nwave
          read (irb) (((w1(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
          read (irb) (((w2(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
          read (irb) (((w3(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
          read (irb) (((w4(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
c         testnorm = ddv*(  sdot(mv,w1,w1) + sdot(mv,w2,w2)
c    1                    + sdot(mv,w3,w3) + sdot(mv,w4,w4))*0.5d0
c         print '(" iwave = ",i5," testnorm = ",1es16.6)',iwave,testnorm

          call scopy (mq,w1,a(1,iwave))
        enddo

      endif

      return
      end subroutine lec8
c______________________________________________________________________________
      subroutine writ8

c..............................................................................
c     write single-particle wave functions and other information on the       .
c     calculation to tape fort.13                                             .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4 afor,head

      common /champ/ qxxn,qyyn,qzzn,qrrn
     1              ,qxxp,qyyp,qzzp,qrrp
     2              ,qxxt,qyyt,qzzt,qrrt

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),imtd,imtg
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
     4              ,q2fin,g2fin,q2cst,qtc0,gtc0
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /den  / rhon(mv),rhop(mv)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /enert/ c14,c15,t14,t15
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /evpro/ rx(mv,2),ry(mv,2),rz(mv,2)
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pairn/ vn,vp,rangen,rangep
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / p1(mv),p2(mv),p3(mv),p4(mv)
     1              ,w1(mv),w2(mv),w3(mv),w4(mv)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)

c..............................................................................
      nwave = nwaven+nwavep
      one = 1.0d0

      iver=6      
      if(npair.ge.6) iver=20
      
      write (13) iver
      write (13) head
      write (13) nwaven,nwavep,npn,npp
      write (13) nx,ny,nz

      write (13) itert,dx
      write (13) imtd,ral,epscst,cqr,cq2,rcut,acut
     1          ,delq,q1n,q1p,q1t,q2n,q2p,q2t

c     ................. Making sure no erroneous information is written to file
      if (cq2.eq.0.d0) then
        qxcstn = 0.d0
        qycstn = 0.d0
        qzcstn = 0.d0
        qxcstp = 0.d0
        qycstp = 0.d0
        qzcstp = 0.d0
        qxcstt = 0.d0
        qycstt = 0.d0
        qzcstt = 0.d0
      endif
      if (cqr.eq.0.d0) then
        qrcstn = 0.d0
        qrcstp = 0.d0
        qrcstt = 0.d0
      endif

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
     1          ,q2cst,gtc0

      write (13) qxxn,qyyn,qzzn,qrrn
     1          ,qxxp,qyyp,qzzp,qrrp
     1          ,qxxt,qyyt,qzzt,qrrt

      write (13) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                            ,t3b,x3b,yt3b,wso,wsoq
     2                            ,te,to

      write (13) b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1          ,b14,b15,b16,b17,byt3a,byt3b,c14,c15,t14,t15

      !Writing force options
      write (13) nfunc,njmunu,ncm2,nmass,ndd,ncoex

      if (iver.ne.20) then
        write (13) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,ntqp
      else
        write (13) npair,vn,vp,rangen,rangep,ambda,xlamb,ntqp
      endif

      write (13) (eqp  (i),i=1,nwave)
      write (13) (delta(i),i=1,nwave)
      write (13) (kparz(i),i=1,nwave)
      write (13) (esp1 (i),i=1,nwave)
      write (13) (v2   (i),i=1,nwave)
      write (13) (v22  (i),i=1,nwave)
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
      end subroutine writ8

c______________________________________________________________________________
      subroutine scopy (n,a,b)

      double precision a(n),b(n)

      do i=1, n
        b(i) = a(i)
      enddo

      return
      end

c______________________________________________________________________________
      subroutine stp (message)

      character (len = *) :: message
      integer :: ierr = 0

      print *,' '
      print *,' =====  S T O P  ===== '
      print *,' '
      print *,' ',message
      print *,' '

      stop
      end

c__________________________________________________________________________
