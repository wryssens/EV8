c______________________________________________________________________________
      program ev8dens

c    ***********************************************************************
c    *                                                                     *
c    * Copyright  M. Bender and P.H. Heenen                                *
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
c
c            *******************************************************
c            *                                                     *
c            * summary of the data list               -  format -  *
c            *  - nprint,ndir                         -  2i5    -  *
c            *  - irtape,iwtape,idtape                -  3i5    -  *
c            *                                                     *
c            *******************************************************
c
c                            description of data
c
c    ***********************************************************************
c    * nprint                                                              *
c    *                                                                     *
c    * ndir     = 1 : plot cut at x = xpp                                  *
c    *          = 2 : plot cut at y = ypp                                  *
c    *          = 3 : plot cut at z = zpp                                  *
c    *                                                                     *
c    * irtape   fort.irtape the ev8 output is read from                    *
c    * idtape   fort.idtape the cut through the density will be written to *
c    * iwtape   fort.idtape the density in the 1/1 box will be written to  *
c    *                      for later processing (subroutine wridens,      *
c    *                      currently deactivated)                         *
c    *                                                                     *
c    ***********************************************************************
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      character*4 afor,head
      character*2 num

      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv

c..............................................................................
c.... read formats
  101 format (20a4)
  102 format (5e15.8)
  105 format (10i5)

c.... write formats
  200 format (/,'  __________________________________________________ ',
     1        /,' |                                                  |',
     2        /,' |  program ev8dens                 release 1.1.0   |',
     3        /,' |  12 March 2014                                   |',
     4        /,' |  Copyright  M. Bender and P.H. Heenen            |',
     6        /,' |    in case of trouble, take it easy...           |',
     7        /,' |    and (please) let us know                      |',
     8        /,' |__________________________________________________|')
  201 format (//,' parameters: mx,my,mz ',3i5,
     2         /,'                   mw ',1i5)
  203 format ( /,' dimensions: nx,ny,nz ',3i5,
     1         /,'                   dx ',1f8.3,' fm')

c..............................................................................
      print 200

c     ............................................ parameters controling output
      read  105,nprint,ndir
      read  105,irtape,iwtape,idtape

c     ......................... print information on allocated array dimensions
      print 201,mx,my,mz,mw

c     ............................................... read input wave-functions
      irb    = irtape
      write (num,'(i2.2)') irb
      open  (irb,form='unformatted',file='fort.'//num)

      print = 1
      call lec8 (iprint)
      close (irb)

c     ....................... print information on dimensions of read-in arrays
      print 203,nx,ny,nz,dx


c     ................................ print cuts through density to fort.itape
      zero = 0.0d0
      call meshcut2d (ndir,idtape  ,zero)

c     ................................. write entire density in 1/1 box to tape
c     call wridens (iwtape)


      call stp (' ***** THE END ***** ')

      end

c______________________________________________________________________________
      subroutine lec8 (iprint)

c..............................................................................
c
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
      common /den  / rhon(mv),rhop(mv)
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint
     1              ,nfunc,njmunu,ncm2,nmass,ndd
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
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw)
     1              ,eqp(mw),delta(mw),ajzd(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /vo8  / delqst,q1tst,q2tst,qrfintst
     1                     ,q1nst,q2nst,qrfinnst
     2                     ,q1pst,q2pst,qrfinpst
      common /wave / p1(mv),p2(mv),p3(mv),p4(mv)
     1              ,w1(mv),w2(mv),w3(mv),w4(mv)
      dimension headd(20)

c..............................................................................
  201 format (/,' ',78('_'),/,' tape information (version ',I1,')',
     1        /,' parameters: mx,my,mz, mw ',4i5)
  202 format   (20a4)
  203 format   (' nwaven,nwavep =',2i5,
     1        /,' mass =',i5,'  charge =',i5,'  n =',i5)
  211 format (/,' nx=',i5,10x,'  ny=',i5,10x,'  nz=',i5,
     1        /,' dx=',e15.8,'Fm  dt=',e15.8,' 10-22s  itert =',i5)
  218 format (/,' constraints data',
     1        /,' cr2           = ',f12.5,
     2        /,' delq iq1, iq2 = ',f12.5,2i5,
     3        /,' cq2           = ',f12.5,
     4        /,' rcut          = ',f12.5)
  260 format   (' npair  = ',i5,' -- hartree-fock --')
  261 format   (' npair  = ',i5,' -- BCS seniority force --')
  265 format   (' npair  = ',i5,' -- BCS constant delta --')
  267 format   (' npair  = ',i5,' -- BCS+LN  seniority force --')
  268 format   (' npair  = ',i5,' -- BCS delta force --')
  270 format   (' npair  = ',i5,' -- BCS+LN delta force --')
  271 format   (' npair  = ',i5,' -- Gaussian force --')
  272 format   (' npair  = ',i5,' -- Gaussian+LN force --')
  290 format   (' imtd  = ',i5,
     1          ' (0,2)=fixed,(1,3)=readjusted,(0,1)=T=0,(2,3)=T=1')
  293 format   (' delq  = ',e15.8,
     3        /,' q1t   = ',e15.8,'  q2t   = ',e15.8)
  294 format   (' actual isoscalar constraints',
     1        /,' cqr = ',e15.8,'  qrfint = ',e15.8,'  qrcstt = ',e15.8,
     2        /,' cq2 = ',e15.8,'  qxfint = ',e15.8,'  qxcstt = ',e15.8,
     3        /,'       ',15x,  '  qyfint = ',e15.8,'  qycstt = ',e15.8,
     4        /,'       ',15x,  '  qzfint = ',e15.8,'  qzcstt = ',e15.8)
  295 format   (' delq  = ',e15.8,
     3        /,' q1n  = ',e15.8,'  q2n   = ',e15.8,
     4        /,' q1p  = ',e15.8,'  q2p   = ',e15.8)
  296 format   (' actual isovector constraints',
     1        /,' cqr = ',e15.8,'  qrfinn = ',e15.8,'  qrcstn = ',e15.8,
     2        /,'                  qrfinp = ',e15.8,'  qrcstp = ',e15.8,
     3        /,' cq2 = ',e15.8,'  qxfinn = ',e15.8,'  qxcstn = ',e15.8,
     4        /,'       ',15x,  '  qyfinn = ',e15.8,'  qycstn = ',e15.8,
     5        /,'       ',15x,  '  qzfinn = ',e15.8,'  qzcstn = ',e15.8,
     6        /,'       ',15x,  '  qxfinp = ',e15.8,'  qxcstp = ',e15.8,
     7        /,'       ',15x,  '  qyfinp = ',e15.8,'  qycstp = ',e15.8,
     8        /,'       ',15x,  '  qzfinp = ',e15.8,'  qzcstp = ',e15.8)

c.......................... wave-functions of the total system read on tape irb
      do i=1,mq
        w1(i) = zero
        p1(i) = zero
      enddo

      read  (irb) iver
      read  (irb) headd
      read  (irb) nwaven,nwavep,npn,npp
      read  (irb) nx,ny,nz
      nwave = nwaven+nwavep
      if (iprint.eq.1) print 201,iver,nx,ny,nz,nwave
                       print 202,headd
      if (iprint.eq.1) print 203,nwaven,nwavep,npn+npp,npp,npn

      if (nwave.gt.mw) call stp (' mw !')
      if    (nx.ne.mx) call stp (' mx !')
      if    (ny.ne.my) call stp (' my !')
      if    (nz.ne.mz) call stp (' mz !')

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
        read (irb) t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                              ,t3b,x3b,yt3b,wso,wsoq
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
        read (irb) t0,x0,t1,x1,t2,x2,t3a
     1            ,x3a,yt3a,t3b,x3b,yt3b,wso,wsoq,te,to
        read (irb) b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1            ,b14,b15,b16,b17,byt3a,byt3b,c14,c15,t14,t15
        read (irb) nfunc,njmunu,ncm2,nmass,ndd,ncoex
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
     1                              ,t3b,x3b,yt3b,wso,wsoq
        read (irb) npair,gn,gp,rangen,rangep
     1            ,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)
      endif

      if (iprint.ne.0) then
        print 211,nx,ny,nz,dx,dt,itert
        print 290,imtd
        if (imtd.le.1) then
          print 293,delqst,q1tst,q2tst
          print 294,cqr,qrfint,qrcstt,
     1              cq2,qxfint,qxcstt,qyfint,qycstt,qzfint,qzcstt
        else
          print 295,delqst,q1nst,q1pst,q2nst,q2pst
          print 296,cqr,qrfinn,qrcstn,qrfinp,qrcstp,
     1              cq2,qxfinn,qxcstn,qyfinn,qycstn,qzfinn,qzcstn,
     2                  qxfinp,qxcstp,qyfinp,qycstp,qzfinp,qzcstp
        endif
        if (npair.eq.0) print 260,npair
        if (npair.eq.1) print 261,npair
        if (npair.eq.2) print 265,npair
        if (npair.eq.3) print 267,npair
        if (npair.eq.4) print 268,npair
        if (npair.eq.5) print 270,npair
        if (npair.eq.6) print 271,npair
        if (npair.eq.7) print 272,npair
      endif

      read (irb) (kparz(i),i=1,nwave)
      read (irb) (esp1 (i),i=1,nwave)
      read (irb) (v2   (i),i=1,nwave)

      if (iver.ge.3) then
        read (irb) (v22(i),i=1,nwave)
      endif

      read (irb) ((npar(i,it),i=1,2),it=1,2)

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

      return
      end subroutine lec8

c______________________________________________________________________________
      subroutine wridens (iwtape)

c..............................................................................
c     iwtape    : fort.iwtape the density will be written to                  .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      character*4 afor,head
      character*2 num

      common /den     / rhon(mx,my,mz),rhop(mx,my,mz)
      common /info    / bidon(500),head(20),irb,nx,ny,nz
      common /noyau   / nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz    / dx,dv

      dimension rhotb(2*mx,2*my,2*mz)
      dimension rhonb(2*mx,2*my,2*mz),rhopb(2*mx,2*my,2*mz)

c..............................................................................
 101  format (3i5)
 106  format (1f8.4)
 111  format (3i5,3e24.16)

c..............................................................................
      do iz=1,nz
        iz1 = iz + nz
        iz2 = nz + 1 - iz
        do iy=1,ny
          iy1 = iy + ny
          iy2 = ny + 1 - iy
          do ix=1,nx
            ix1 = ix + nx
            ix2 = nx + 1 - ix

            rn = rhon(ix,iy,iz) 
            rp = rhop(ix,iy,iz)
            rt = rn + rp
             
            rhonb(ix1,iy1,iz1) = rn
            rhonb(ix1,iy2,iz1) = rn
            rhonb(ix2,iy1,iz1) = rn
            rhonb(ix2,iy2,iz1) = rn
            rhonb(ix1,iy1,iz2) = rn
            rhonb(ix1,iy2,iz2) = rn
            rhonb(ix2,iy1,iz2) = rn
            rhonb(ix2,iy2,iz2) = rn

            rhopb(ix1,iy1,iz1) = rp
            rhopb(ix1,iy2,iz1) = rp
            rhopb(ix2,iy1,iz1) = rp
            rhopb(ix2,iy2,iz1) = rp
            rhopb(ix1,iy1,iz2) = rp
            rhopb(ix1,iy2,iz2) = rp
            rhopb(ix2,iy1,iz2) = rp
            rhopb(ix2,iy2,iz2) = rp

            rhotb(ix1,iy1,iz1) = rt
            rhotb(ix1,iy2,iz1) = rt
            rhotb(ix2,iy1,iz1) = rt
            rhotb(ix2,iy2,iz1) = rt
            rhotb(ix1,iy1,iz2) = rt
            rhotb(ix1,iy2,iz2) = rt
            rhotb(ix2,iy1,iz2) = rt
            rhotb(ix2,iy2,iz2) = rt
          enddo
        enddo
      enddo

c     .........................................................................
      open  (iwtape)
      write (iwtape,101) 2*nx,2*ny,2*nz
      write (iwtape,106) dx

      do iz=1,2*nz
      do iy=1,2*ny
      do ix=1,2*nx
        write(iwtape,111) ix,iy,iz,rhonb(ix,iy,iz),rhopb(ix,iy,iz)
      enddo
      enddo
      enddo

      close (iwtape)

      return
      end subroutine wridens

c______________________________________________________________________________
      subroutine meshcut2d (icase,iwtape,ppp)

c..............................................................................
c     icase = 1 : plot cut at x = xpp                                         .
c           = 2 : plot cut at y = ypp                                         .
c           = 3 : plot cut at z = zpp                                         .
c     iwtape    : fort.iwtape the density will be written to                  .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      character*4 afor,head
      character*2 num

      common /den     / rhon(mx,my,mz),rhop(mx,my,mz)
      common /info    / bidon(500),head(20),irb,nx,ny,nz
      common /interpol/ cx(2*mx),cy(2*my),cz(2*mz)
      common /noyau   / nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz    / dx,dv

      dimension rhot (2*mx,2*my,2*mz)
      dimension rhotx(2*my,     2*mz)
      dimension rhoty(2*mx,     2*mz)
      dimension rhotz(2*my,2*my     )

c..............................................................................
 101  format('! nx ',1i4,' ny ',1i4,' xmin ',1f6.2,
     1       ' xmax ',1f6.2,' ymin ',1f6.2,' ymax ',1f6.2)
 111  format(3f16.8)

c.......................................................... diagnostic printing
c                              of the interpolation matrices in steps of dx / 4
c                                            from one box boundary to the other
c                                    to be activated for testing purposes only!
c      jcase = 1
c      call primesh (jcase)
c      jcase = 2
c      call primesh (jcase)
c      jcase = 3
c      call primesh (jcase)

c..............................................................................
      open (iwtape)

c     ................................................. calculate total density
      do iz=1,nz
        iz1 = iz + nz
        iz2 = nz + 1 - iz
        do iy=1,ny
          iy1 = iy + ny
          iy2 = ny + 1 - iy
          do ix=1,nx
            ix1 = ix + nx
            ix2 = nx + 1 - ix
            rt = rhon(ix,iy,iz) + rhop(ix,iy,iz)
            rhot(ix1,iy1,iz1) = rt
            rhot(ix1,iy2,iz1) = rt
            rhot(ix2,iy1,iz1) = rt
            rhot(ix2,iy2,iz1) = rt
            rhot(ix1,iy1,iz2) = rt
            rhot(ix1,iy2,iz2) = rt
            rhot(ix2,iy1,iz2) = rt
            rhot(ix2,iy2,iz2) = rt
          enddo
        enddo
      enddo

c     .........................................................................
      select case (icase)
      case (1)
        xpp  = ppp

        zmin = -(nz - 0.5d0) * dx
        zmax =  (nz - 0.5d0) * dx
        ymin = -(ny - 0.5d0) * dx
        ymax =  (ny - 0.5d0) * dx
        write (iwtape,101) 2*nz,2*ny,zmin,zmax,ymin,ymax

        yp = 0.4d0
        zy = 0.4d0
        call setinterpol (1,0,0,xpp,yp,zp)   
        rhotx(:,:) = 0.0d0
        do iz=1,2*nz
        do iy=1,2*ny
        do ix=1,2*nx
          rhotx(iy,iz) = rhotx(iy,iz) + cx(ix)*rhot(ix,iy,iz)
        enddo
        enddo
        enddo

        yp = -0.5d0 * (2*ny+1) * dx
        do iy=1,2*ny
          yp = yp + dx
          zp = -0.5d0 * (2*nz+1) * dx
          do iz=1,2*nz
            zp = zp + dx
c           write (iwtape,111) zp,yp,rhotx(iy,iz)
            write (iwtape,111) rhotx(iy,iz)
          enddo
        enddo

      case (2)
        ypp  = ppp

        zmin = -(nz - 0.5d0) * dx
        zmax =  (nz - 0.5d0) * dx
        xmin = -(nx - 0.5d0) * dx
        xmax =  (nx - 0.5d0) * dx
        write (iwtape,101) 2*nz,2*nx,zmin,zmax,xmin,xmax

        zy = 0.4d0
        xp = 0.4d0
        call setinterpol (0,1,0,xp,ypp,zp)
        rhoty(:,:) = 0.0d0
        do iz=1,2*nz
        do iy=1,2*ny
        do ix=1,2*nx
          rhoty(ix,iz) = rhoty(ix,iz) + cy(iy)*rhot(ix,iy,iz)
        enddo
        enddo
        enddo

        xp = -0.5d0 * (2*nx+1) * dx
        do ix=1,2*nx
          xp = xp + dx
          zp = -0.5d0 * (2*nz+1) * dx
          do iz=1,2*nz
            zp = zp + dx
c           write (iwtape,111) zp,yp,rhoty(ix,iz)
            write (iwtape,111) rhoty(ix,iz)
          enddo
        enddo

      case (3)
        zpp  = ppp

        xmin = -(nx   - 0.5d0) * dx
        xmax =  (nx   - 0.5d0) * dx
        ymin = -(ny   - 0.5d0) * dx
        ymax =  (ny   - 0.5d0) * dx
        write (iwtape,101) 2*ny,2*nx,ymin,ymax,xmin,xmax

        xy = 0.4d0
        yp = 0.4d0
        call setinterpol (0,0,1,xp,yp,zpp)
        rhotz(:,:) = 0.0d0
        do iz=1,2*nz
        do iy=1,2*ny
        do ix=1,2*nx
          rhotz(ix,iy) = rhotz(ix,iy) + cz(iz)*rhot(ix,iy,iz)
        enddo
        enddo
        enddo

        xp = -0.5d0 * (2*nx+1) * dx
        do ix=1,2*nx
          xp = xp + dx
          yp = -0.5d0 * (2*ny+1) * dx
          do iy=1,2*ny
            yp = yp + dx
c           write (iwtape,111) yp,xp,rhotz(ix,iy)
            write (iwtape,111) rhotz(ix,iy)
          enddo
        enddo

      case default
        call stp(' meshcut2d: icase not known!')

      end select

      close (iwtape)

      return
      end subroutine meshcut2d

c______________________________________________________________________________
      subroutine setinterpol (ix,iy,iz,x,y,z)

c..............................................................................
c     set up interpolation matrices to calculate functions at (x,y,z)         .
c     ix # 0 : set up interpolation to x                                      .
c     iy # 0 : set up interpolation to y                                      .
c     iz # 0 : set up interpolation to z                                      .
c                                                                             .
c     ATTENTION: check if xm/ym/zm have to be  = (y1-y) or (y-y1) !!!!!       .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4 head

      parameter (epsm3=1.0d-3)

      common /info    / bidon(500),head(20),irb,nx,ny,nz
      common /interpol/ cx(2*mx),cy(2*my),cz(2*mz)
      common /nxyz    / dx,dv

c     .......................................................... initialization
      pi    = 4.0d0 * atan2(1.0d0,1.0d0)
      cx(:) = 0.0d0
      cy(:) = 0.0d0
      cz(:) = 0.0d0

c............................... set up the interpolation matrix in x direction
      if (ix.ne.0) then
        pnx  =  pi  / (2*nx)
        fac  =  1.0d0 / (2*nx)
        x1   = -0.5d0 * (2*nx+1) * dx
        do i=1,2*nx
          x1 = x1+dx
          xm = (x1-x)/dx
          if (abs(xm).le.epsm3) then
            c = 1.0d0 / fac
          else
            c = sin(pi * xm)/sin(pnx * xm)
          endif
          cx(i) = fac * c
        enddo
      endif

c............................... set up the interpolation matrix in y direction
      if (iy.ne.0) then
        pny  =  pi  / (2*ny)
        fac  =  1.0d0 / (2*ny)
        y1   = -0.5d0 * (2*ny+1) * dx
        do i=1,2*ny
          y1 = y1+dx 
          ym = (y1-y)/dx
          if (abs(ym).le.epsm3) then
            c = 1.0d0 / fac
          else
            c = sin(pi * ym)/sin(pny * ym)
          endif
          cy(i) = fac * c
        enddo
      endif

c............................... set up the interpolation matrix in z direction
      if (iz.ne.0) then
        pnz  =  pi  / (2*nz)
        fac  =  1.0d0 / (2*nz)
        z1   = -0.5d0 * (2*nz+1) * dx
        do i=1,2*nz
          z1 = z1+dx
          zm = (z1-z)/dx
          if (abs(zm).le.epsm3) then
            c = 1.d0 / fac
          else
            c = sin(pi * zm)/sin(pnz * zm)
          endif
          cz(i) = fac * c
        enddo
      endif

      return
      end subroutine setinterpol

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
      subroutine primesh (icase)

c..............................................................................
c     diagnostic routine to print interpolation functions                     .
c     usually de-activated                                                    .
c     activate for debugging purposes only                                    .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4 head

      parameter (epsm3=1.0d-3)

      common /info    / bidon(500),head(20),irb,nx,ny,nz
      common /interpol/ cx(2*mx),cy(2*my),cz(2*mz)
      common /nxyz    / dx,dv

c..............................................................................
 81   format (/,1x,180('='),/)
 82   format (' interpolation matrix for x direction',/)
 84   format (' interpolation matrix for y direction',/)
 86   format (' interpolation matrix for z direction',/)
 91   format (1f6.2,25es9.1)

c..............................................................................
      xp = 0.5d0 * dx
      yp = 0.5d0 * dx
      zp = 0.5d0 * dx

      print 81
      if (icase.eq.1) print 82
      if (icase.eq.2) print 84
      if (icase.eq.3) print 86

      select case (icase)
      case (1)
        if (2*nx.gt.25) call stp('primesh: nx > 25!')
        xp = - ( 0.5d0 * (2*nx+1 ) - 0.25d0) * dx
        do ix = 1,2*nx
          xp = xp + 0.25d0 * dx
          call setinterpol (1,0,0,xp,yp,zp)
          print 91,xp,(cx(jx),jx=1,2*nx)
          xp = xp + 0.25d0 * dx
          call setinterpol (1,0,0,xp,yp,zp)
          print 91,xp,(cx(jx),jx=1,2*nx)
          xp = xp + 0.25d0 * dx
          call setinterpol (1,0,0,xp,yp,zp)
          print 91,xp,(cx(jx),jx=1,2*nx)
          xp = xp + 0.25d0 * dx
          call setinterpol (1,0,0,xp,yp,zp)
          print 91,xp,(cx(jx),jx=1,2*nx)
        enddo
        xp = xp + 0.25d0 * dx
        call setinterpol (1,0,0,xp,yp,zp)
        print 91,xp,(cx(jx),jx=1,2*nx)

      case (2)
        if (2*ny.gt.25) call stp('primesh: ny > 25!')
        yp = - ( 0.5d0 * (2*ny+1) - 0.25d0 ) * dx
        do iy = 1,2*ny
          yp = yp + 0.25d0 * dx
          call setinterpol (0,1,0,xp,yp,zp)
          print 91,yp,(cy(jy),jy=1,2*ny)
          yp = yp + 0.25d0 * dx
          call setinterpol (0,1,0,xp,yp,zp)
          print 91,yp,(cy(jy),jy=1,2*ny)
          yp = yp + 0.25d0 * dx
          call setinterpol (0,1,0,xp,yp,zp)
          print 91,yp,(cy(jy),jy=1,2*ny)
          yp = yp + 0.25d0 * dx
          call setinterpol (0,1,0,xp,yp,zp)
          print 91,yp,(cy(jy),jy=1,2*ny)
        enddo
        yp = yp + 0.25d0 * dx
        call setinterpol (0,1,0,xp,yp,zp)
        print 91,yp,(cy(jy),jy=1,2*ny)

      case (3)
        if (2*nz.gt.25) call stp('primesh: nz > 25!')
        zp = - ( 0.5d0 * (2*nz+1)  - 0.25d0 ) * dx
        do iz = 1,2*nz
          zp = zp + 0.25d0 * dx
          call setinterpol (0,0,1,xp,yp,zp)
          print 91,zp,(cz(jz),jz=1,2*nz)
          zp = zp + 0.25d0 * dx
          call setinterpol (0,0,1,xp,yp,zp)
          print 91,zp,(cz(jz),jz=1,2*nz)
          zp = zp + 0.25d0 * dx
          call setinterpol (0,0,1,xp,yp,zp)
          print 91,zp,(cz(jz),jz=1,2*nz)
          zp = zp + 0.25d0 * dx
          call setinterpol (0,0,1,xp,yp,zp)
          print 91,zp,(cz(jz),jz=1,2*nz)
        enddo
        zp = zp + 0.25d0 * dx
        call setinterpol (0,0,1,xp,yp,zp)
        print 91,zp,(cz(jz),jz=1,2*nz)

      case default
        call stp(' primesh: icase not known!')

      end select

      print 81

      return
      end subroutine primesh

c______________________________________________________________________________
      subroutine scopy (n,a,b)

      double precision a(n),b(n)

      do i=1, n
        b(i) = a(i)
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
      subroutine sscal (n,fac,a)

      double precision a(n),fac

      do i=1, n
        a(i) = fac * a(i)
      enddo

      return
      end subroutine sscal

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
      end subroutine stp

c__________________________________________________________________________
      double precision function smaxf (n,a)

      double precision a(n),tmp

      tmp = -1.d20
      do i=1,n
        tmp = max(a(i),tmp)
      enddo

      smaxf = tmp

      return
      end

c______________________________________________________________________________
