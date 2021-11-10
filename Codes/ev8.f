c______________________________________________________________________________
      program ev8     

c    ***********************************************************************
c    *                                                                     *
c    * static hartree-fock + b.c.s + lipkin-nogami (ln).                   *
c    * time-reversal symmetry                                              *
c    * triaxial nuclei                                                     *
c    * symmetry with respect to x=0, y=0 and z=0 planes                    *
c    * monopole and quadrupole constraints                                 *
c    * isoscalar and isovector constraints                                 *
c    *                                                                     *
c    * each s.p. wave-function which describes a nucleon                   *
c    * and its t-reversed (t-r) has four 3d components:                    *
c    * 1= plus (t-r plus)  the real part of the spin up (t-r down)         *
c    * 2= plus (t-r minus) the imaginary part of the spin up (t-r down)    *
c    * 3= plus (t-r minus) the real part of the spin down (t-r up)         *
c    * 4= plus (t-r plus)  the imaginary part of the spin down (t-r up)    *
c    *                                                                     *
c    * the parities of the four components with respect to the             *
c    *  x=0, y=0 and z=0 planes are:                                       *
c    *  1; spin  up  real         (+,+, kparz)                             *
c    *  2; spin  up  imaginary    (-,-, kparz)                             *
c    *  3; spin down real         (-,+,-kparz)                             *
c    *  4; spin down imaginary    (+,-,-kparz)                             *
c    *                                                                     *
c    *  central + spin-orbit + tensor skyrme interaction                   *
c    *                                                                     *
c    *  either from the anti-symmetrized density-dependent vertex          *
c    *     t0*(1+x0*ps)*del                                                *
c    *   + (t1/2)*(1+x1*ps)*(k'.k'*del+del*k.k)                            *
c    *   + t2*(1+x2*ps)*(k'.del*k)                                         *
c    *   + (t3a/6)*(1+x3a*ps)*(rho**yt3a)*del                              *
c    *   + (t3b/6)*(1+x3b*ps)*(rho**yt3b)*del                              *
c    *   + i*wso*(sig1+sig2)*(k'^del*k)                                    *
c    *   + (te/2)*{[3*(sig1.k')*(sig2.k') - (sig1.sig2)*(k'.k')]*del       *
c    *            +del*[3*(sig1.k)*(sig2.k) - (sig1.sig2)*(k.k)]}          *
c    *   + to*[3*(sig1.k')*del*(sig2.k) - (sig1.sig2) k'.del*k]            *
c    *                                                                     *
c    *  or from a functional containing the same terms, but with the       *
c    *  freedom to decouple the coupling constants of spin-orbit           *
c    *  (wso != wsoq) and tensor terms (b14, b15, b16 and b17 decoupled    *
c    *  from the interrelations between t1, x1, t2, x2, te and to)         *
c    *                                                                     *
c    *  For the Coulomb interaction, the exact central term and the        *
c    *  exchange term in local Slater approximation are used, in both      *
c    *  cases using the point proton density as charge density             *
c    *                                                                     *
c    *  pairing options:                                                   *
c    *   0) hartree-fock                                                   *
c    *   1) bcs    seniority force                                         *
c    *   2) bcs    fixed pairing gaps                                      *
c    *   3) bcs+ln pairing force                                           *
c    *   4) bcs    delta force: v*(1-ps)*del/2                             *
c    *   5) bcs+ln delta force: v*(1-ps)*del/2                             *
c    *        active-pairing-space defined with respect to fermi level     *
c    *          the cut off may done either only above the fermi level     *
c    *          or above and below                                         *
c    *      bcs ground state or two quasiparticle state                    *
c    *   6) finite range pairing                                           *
c    *   7) finite range pairing + ln                                      *
c    *                                                                     *
c    *  when the lipkin-nogami option is used, constraints are on          *
c    *  the ln-corrected moments                                           *
c    *                                                                     *
c    *  decorelated two-quasiparticle states (for even-even nuclei)        *
c    *  filling approximation for odd-A nuclei                             *
c    *                                                                     *
c    *  two-body centre of mass correction and/or the square of the        *
c    *  spin-current tensor density Jmunu can be included selfconsistently *
c    *  (see njmunu and ncm2 below)                                        *
c    *                                                                     *
c    * units: length, fm; time, 10-22 sec; energy; mev                     *
c    *                                                                     *
c    ***********************************************************************

c    ****************************************************************
c    *                                                              *
c    * to compile properly the program requires a param8.h file     *
c    *   containing the folowing information:                       *
c    * parameter (mx=14,my=14,mz=14,mc=8,mv=mx*my*mz,mq=4*mv,mw=40) *
c    *   where                                                      *
c    *   mx, my,mz are : numbers of points in each directions       *
c    *                                                              *
c    ****************************************************************

c    **************************************************************
c    * Program logic                                              *
c    *   - ev8                                                    *
c    *     - lec8             reads the initial spwf on fort.12   *
c    *     - lecfor           skyrme force                        *
c    *       - blockdata case                                     *
c    *     - lecmas           proton and neutron masses ...       *
c    *     - evolve           evolves the spwf                    *
c    *       - inisp                                              *
c    *       - class                                              *
c    *       - ortho          orthogonalization of the spwf       *
c    *       - gap            pairing                             *
c    *         - gaphf                HF                          *
c    *         - gapbcs               BCS                         *
c    *           - gform                                          *
c    *         - gapln                Lipkin Nogami               *
c    *           - gform      density dependent form factor       *
c    *         - gfrbcs               BCS finite range            *
c    *           - gapwf                                          *
c    *         - gfrln                BCS + LN finite range       *
c    *           - gapwf                                          *
c    *       - densit         density                             *
c    *       - newpot         H.F. potential                      *
c    *         - mome         multipole moments                   *
c    *           - distance   damping of constraint at large radii*
c    *         - vbound       coulomb boundary conditions         *
c    *         - vcoul        coulomb potential                   *
c    *         - vcal         skyrme contributions                *
c    *         - vcalg        skyrme contributions                *
c    *         - vcalte       skyrme contributions  (tensor)      *
c    *         - vcalq        constraint contributions            *
c    *       - conver         convergence check                   *
c    *       - fsum           prints a summary                    *
c    *         - cpling       pring force parameters & options    *
c    *         - ReanalyseLag Reanalyses with lagrange derivatives*
c    *       - figaro         printout                            *
c    *         - fprte                                            *
c    *         - fprtj                                            *
c    *         - fprtz                                            *
c    *         - fprtm                                            *
c    *           - momehigh                                       *
c    *       - reord          diagonalization if ndiag is not 0   *
c    *          - hpsi        calculates h|psi>                   *
c    *             - vcen                                         *
c    *       - loop on the H.F. iterative scheme...............   *
c    *           - pro        iterates the spwf               .   *
c    *              - hpsi, vcen                              .   *
c    *           - class                                      .   *
c    *           - ortho      orthogonalization of the spwf   .   *
c    *           - gap                                        .   *
c    *           - densit     density                         .   *
c    *           - newpot     H.F. potential                  .   *
c    *           - fsum       prints a summary                .   *
c    *           - figaro     printout                        .   *
c    .       - end loop........................................   *
c    *     - jdeux            average angular momentum J**2       *
c    *     - jmunu            so tensor energy in perturbation    *
c    *     - pipj             2body centre of mass in perturbation*
c    *     - writ8            storage of the spwf on fort.13      *
c    * and also :                                                 *
c    *     - diagon           real sym. matrix diaginalization    *
c    *     - matin1           real sym. matrix inversion          *
c    *     - stp              print error message from stop to    *
c    *                        standard out and standard error     *
c    *     - deriv, der, derx, dery, derz, lapla                  *
c    *     - inilag, derlag, laplalag                             *
c    *     - scopy, sdot, saxpy, sscal, sswap                     *
c    **************************************************************
c
c    ************************************************************
c    *                                                          *
c    * data(per nucleus)                           -  format -  *
c    *  - head                                     -  20a4   -  *
c    *  - npx,npy,npz,iCoul                        -  4i5    -  *
c    *  - dt                                       -   e15.8 -  *
c    *  - nitert,nxmu,ndiag                        -  3i5    -  *
c    *  - nprint, iverb                            -  2i5    -  *
c    *  - npn,npp                                  -  2i5    -  *
c    *  - afor                                     -   a20   -  *
c    *      if (afor.eq.'XXXX')                              -  *
c    *       - ncm2,nfunc,nmass,ncoex              -  4i5    -  *
c    *       - t0,x0,                              -  2e15.8 -  *
c    *       - t1,x1,t2,x2                         -  4e15.8 -  *
c    *       - t3a,x3a,yt3a                        -  3e15.8 -  *
c    *       - t3b,x3b,yt3b                        -  3e15.8 -  *
c    *         if (nfunc.eq.0)                                  *
c    *         - wso                               -   e15.8 -  *
c    *         - te,to                             -  2e15.8 -  *
c    *         if (nfunc.eq.1)                                  *
c    *         - wso,wsoq                          -  2e15.8 -  *
c    *         - b14,b15,b16,b17                   -  4e15.8 -  *
c    *         if (nmass.eq.2)                                  *
c    *         - hbm(1),hbm(2)                     -  2e15.8 -  *
c    *  - npair,icut                               -  2i5    -  *
c    *  - gn,encut,delmax(1),alpha                 -  4e15.8 -  *
c    *  - gp,epcut,delmax(2)                       -  3e15.8 -  *
c    *  - epse,epsdh,epsq,epsf                     -  4e15.8 -  *
c    *  - imtd,ifrt,imtg, icutq                    -  4i5    -  *
c    *  - ral,epscst,rcut                          -  3e15.8 -  *
c    *  - cqr,rtcst or cqr,rncst,rpcst             -  3e15.8 -  *
c    *  - cq2,delq                                 -  2e15.8 -  *
c    *  - q1t,q2t or q1n,q2n,q1p,q2p               -  4e15.8 -  *
c    *                                                          *
c    ************************************************************

c                            description of data

c    ***********************************************************************
c    * head                                                                *
c    *        title of calculation                                         *
c    ***********************************************************************

c    ***********************************************************************
c    * npx,npy,npz                                                         *
c    *        number of extra points to calculate the coulomb field        *
c    * iCoul                                                               *
c    *        Order of the discretisation of the Coulomb problem (1 or 2)  *
c    ***********************************************************************

c    ***********************************************************************
c    * dt                                                                  *
c    *        imaginary time step in 10**(-22)s                            *
c    ***********************************************************************

c    ***********************************************************************
c    * nitert                                                              *
c    *        number of iterations = abs(nitert)                           *
c    *        nitert<0: iteration index is reset to zero                   *
c    * nxmu                                                                *
c    *        damping factor for the mean-field potential (evolve)         *
c    *        xmu=nxmu/100.; w(n+1)=xmu*w(n+1)+(1.-xmu)*w(n)               *
c    *        if nxmu is read 0, xmu is set to 0.25                        *
c    * ndiag                                                               *
c    *       diagonalization of the HF hamiltonian in the spwf subspace    *
c    *       ndiag = 1 diagonalization at time 0.0                         *
c    ***********************************************************************

c    ***********************************************************************
c    *                                                                     *
c    * nprint                                                              *
c    *        number of time steps between two consecutive printouts       *
c    ***********************************************************************
c
c    ***********************************************************************
c    * npn, npp                                                            *
c    *       number of neutrons and protons npn,npp                        *
c    *       npn should be less than 2*nwaven, npp less than 2*nwavep      *
c    ***********************************************************************

c    ***********************************************************************
c    * afor                                                                *
c    *   available forces:                                                 *
c    *        SkM* SIII Ska SGII Sly4 Sly5 Sly6 Sly7 SkP                   *
c    *        SkM  SKI4 T22 T24 T26 T42 T44 T46 T62 T64 T65 T66            *
c    *        SkX  Sly4T Sly5+T SKI3 SVMIN UNEDF0 UNEDF1 UNEDF1SO          *
c    *        SlyIIIK07 SLYIIIK08 SLYIIIK09 SLYIIIK10                      *
c    *                                                                     *
c    * to read a new force, set afor to a unrecongnised string             *
c    * and read t0,... accordingly:                                        *
c    *     ncm2,nfunc,nmass,ncoex                                          *
c    *     t0,x0                                                           *
c    *     t1,x1,t2,x2                                                     *
c    *     t3a,x3a,yt3a                                                    *
c    *     t3b,x3b,yt3b                                                    *
c    *     if (nfunc.eq.0) then                                            *
c    *       wso                                                           *
c    *       te,to                                                         *
c    *     else                                                            *
c    *       wso,wsoq                                                      *
c    *       b14,b15,b16,b17                                               *
c    *     endif                                                           *
c    *     if (nmass.eq.2) then                                            *
c    *       hbm(1),hbm(2)                                                 *
c    *     endif                                                           *
c    ***********************************************************************

c    ***********************************************************************
c    * ncm2,nfunc,nmass,ncoex                                              *
c    *        options of the interaction. For pre-defined parametrizations *
c    *        they are no longer read from data. This gave rise to too many*
c    *        mistakes. If you want to use a standard interaction with an  *
c    *        non-standard option, you have to read the entire set of      *
c    *        coupling constants                                           *
c    *                                                                     *
c    * ncm2                                                                *
c    *        centre of mass correction                                    *
c    *        = -2 : 2-body  and 1-body Centre-of-mass correction included *
c    *               perturbatively: meaning included in the total energy  *
c    *               but not in the single-particle hamiltonian            *
c    *        = -1 : no c.m. correction at all                             *
c    *        =  0 : only 1-body c.m. correction is considered in the      *
c    *               variational equations and contained in the energies   *
c    *               printed at the end.                                   *
c    *               2-body c.m. correction is calculated perturbatively   *
c    *               after the last iteration; it is not included in the   *
c    *               printed total energy                                  *
c    *        =  1 : 2-body c.m. correction included self consistently     *
c    *               and occurs then in the detail of the energy terms     *
c    * nfunc                                                               *
c    *        antisymmetric Skyrme force or Skyrme energy functional       *
c    *        = 0 :  coupling constants of the energy functional (with the *
c    *               exception of b9q so far for compatibility reasons)    *
c    *               calculated from the t,x and W of the central, tensor, *
c    *               and spin-orbit Skyrme force (as in previous versions  *
c    *               of this code)                                         *
c    *        = 1 :  coupling constants of the J^2 terms not related to    *
c    *               central + tensor Skyrme vertices. Additional coupling *
c    *               constants of the functional b14, b15, b16 and b17     *
c    *               read from extra line of data. c14 and c15 set to zero.*
c    * nmass                                                               *
c    *        proton and neutron masses                                    *
c    *        = 0 : both masses are equal to the average (m_n+m_p)/2       *
c    *        = 1 : masses are different and taken from the CERN booklet   *
c    *        = 2 : read hbar^2/2m from data                               *
c    *              and set hbar to experimental value. m_n and m_p are    *
c    *              then defined accordingly                               *
c    *              (which plays some role for cm correction)
c    *        The main use of the nuclon masses is to determine hbar^2/2m, *
c    *        which enters the calculation of the kinetic energy and       *
c    *        centre-of-mass-correction. The value of hbar is also used    *
c    *        in the imaginary time step. The mass itself seems not to be  *
c    *        used in ev8 for another purpose.                             *
c    *        For consistency with legacy calculations the default options *
c    *        are 0 and 1 as in the earlier versions of the codes.         *
c    *        To use default parameterizations with the hbar^2/2m values   *
c    *        used in the respective fit, for the moment their parameters  *
c    *        have to be read from data (afor=XXXX) with nmass = 2. In     *
c    *        this case, an additional line with h^2/2m is read.           *
c    *                                                                     *
c    * ncoex                                                               *
c    *        treatment of Coulomb exchange                                *
c    *        = 0 : Slater approximation                                   *
c    *        = 1 : no Coulomb exchange                                    *
c    ***********************************************************************

c    ***********************************************************************
c    * npair                                                               *
c    *       method used to define the occupation of the orbitals          *
c    *           0   hf                                                    *
c    *           1   bcs       with seniority pairing force                *
c    *           2   bcs       with fixed pairing gaps                     *
c    *           3   bcs + ln  with seniority pairing force                *
c    *           4   bcs       with delta pairing force                    *
c    *           5   bcs + ln  with delta pairing force                    *
c    * icut                                                                *
c    *       0 cutoff only above the fermi level                           *
c    *       1 cutoff above and below                                      *
c    ***********************************************************************
c
c    ***********************************************************************
c    * gn,gp                                                               *
c    *       intensity of the neutron and proton T=1 pairing forces        *
c    *        seniority force                                              *
c    *           strengthes: vn=gn/(11+n), vp=gp/(11+p)                    *
c    *        delta force (L=0, S=0, T=1): (v/4)*(1-sigma1*sigma2)*delta   *
c    *           strengthes:  vn=gn,  vp=gp                                *
c    * alpha                                                               *
c    *        density dependence of the strength                           *
c    *           vn=gn*(1-alpha*rho/rhosat)                                *
c    *           vp=gp*(1-alpha*rho/rhosat)                                *
c    *           strengthes:  vn=gn,  vp=gp                                *
c    * encut,epcut                                                         *
c    *       cutoff parameters of pairing phase space                      *
c    *                 cut-off is defined relative to the fermi level      *
c    * delmax(1), delmax(2)                                                *
c    *       initial n and p pairing gaps                                  *
c    *       enforced values when npair = 2                                *
c    * alpha                                                               *
c    *       coeff of the density dependent term, if any                   *
c    ***********************************************************************
c
c    ***********************************************************************
c    *  epse,epsdh,epsq,epsf                                               *
c    *    convergence criteria for termination of iteration                *
c    *    The iteration is ended when the following criteria for the       *
c    *    relative change of the total energy                              *
c    *            (E(i) - E(i-n))/E(i) < epse                              *
c    *    weighted sum of dispersions of the s.p. Hamiltonian              *
c    *            sum_i v_i^2 [ (i|h^2|i) - (i|h|i)^2 ] / A < epsdh        *
c    *    change in proton and neutron Fermi energies                      *
c    *            lambda(i) - lambda(i-n) < epsl                           *
c    *    and deviation from the value of the quadrupole constraint        *
c    *            (Q(i-n) - Qcon ) / Qcon < epsq                           *
c    *    are simultaneously met for the last 7 iterations.                *
c    *    For unconstrained calculations, epsq has to be set to a very     *
c    *    large value                                                      *
c    ***********************************************************************

c    ***********************************************************************
c    * imtd                                                                *
c    *       type of constraint                                            *
c    *         0 isoscalar constraint, parameters fixed                    *
c    *         1 isoscalar constraint, parameters readjusted               *
c    *         2 isovector constraint, parameters fixed                    *
c    *         3 isovector constraint, parameters readjusted               *
c    * ifrt                                                                *
c    *        0 constraint strengths read in data                          *
c    *        1 constraint strengths from file fort.12                     *
c    * imtg                                                                *
c    *        0 constraint on q, gamma                                     *
c    *        1 constraint on q only                                       *
c    * icutq                                                               *
c    *        cut-off type.                                                *
c    *                                                                     *
c    *        0 density-dependent cutoff a la Klemens Rutz                 *
c    *          cutoff radius is rcut fm above the nuclear surface         *
c    *          with a width of acut = 0.4 fm                              *
c    *        1 density-dependent cut-off function = tan(rho/rhoc)         *
c    *          with rhoc=-rcut (fm-3)                                     *
c    *        2 spherical fermi cut-off function of radius rcut (fm)       *
c    *          if rcut is read 0 it is set to 100                         *
c    ***********************************************************************
c
c    ***********************************************************************
c    * ral                                                                 *
c    *        damping factor for the constraint potential (mome)           *
c    *        0. <= ral <= 1.                                              *
c    *        if  ral is read 0., it is set to 0.1                         *
c    * epscst                                                              *
c    *        slow-down factor for readjustment of constraint (imtd=1 or 3)*
c    *        0. <= epscst <= 1.                                           *
c    *        if epscst is read 0., it is set to 0.02                      *
c    * rcut                                                                *
c    *        cut-off radius                                               *
c    ***********************************************************************
c
c    ***********************************************************************
c    * cqr                                                                 *
c    *        intensity of monopole constraint                             *
c    * rtcst or rncst,rpcst                                                *
c    *        imtd=0,1                                                     *
c    *          rtcst = radius of mass    monopole constraint              *
c    *        imtd=2,3                                                     *
c    *          rncst = radius of neutron monopole constraint              *
c    *          rpcst = radius of proton  monopole constraint              *
c    ***********************************************************************

c    ***********************************************************************
c    * delq                                                                *
c    *        unit for quadrupole deformation mesh; strictly positive      *
c    * cq2                                                                 *
c    *        intensity of the quadrupole constraint                       *
c    *        cq2 is positive or zero                                      *
c    ***********************************************************************

c    ***********************************************************************
c    * q1t,q2t or q1n,q2n,q1p,q2p                                          *
c    *        quadrupole constraint (in units of delq) and triaxility      *
c    *        q1 and q2 can be positive, zero or negative                  *
c    *        imtd=0,1                                                     *
c    *          q1t,q2t = mass    quadrupole constraint                    *
c    *        imtd=2,3                                                     *
c    *          q1n,q2n = neutron quadrupole constraint                    *
c    *          q1p,q2p = proton  quadrupole constraint                    *
c    *        relations between data and moments                           *
c    *        qx  = -delq*(q1-q2)/2    = -q0*cos(gam+60)                   *
c    *        qy  = -delq*(q1+2*q2)/2  = -q0*cos(gam-60)                   *
c    *        qz  =  delq*(2*q1+q2)/2  =  q0*cos(gam)                      *
c    *        q0  =  delq*sqrt(q1**2+q2**2+q1*iq2)                         *
c    *            =  sqrt(2*(qx**2+qy**2+qz**2)/3)                         *
c    *        gam = atan2(q2*sqrt(3.),2*q1+q2) = atan2(qx-qy,sqrt(3.)*qz)  *
c    *            = 0 at the spherical point (q1=q2=0)                     *
c    ***********************************************************************

c    ***********************************************************************
c    *                                                                     *
c    *                     y2=z2<x2           y2<x2=z2                     *
c    *                     prolate            oblate                       *
c    *                           \             /                           *
c    *                            \           /                            *
c    *                             \         / IQ2 axis, gamma = 60        *
c    *                              \       /                              *
c    *                               \     /                               *
c    *                                \   /                                *
c    *                                 \ /                                 *
c    *     z2<x2=y2     _______________ v_______________     x2=y2<z2      *
c    *     oblate                       /\      IQ1 axis     prolate       *
c    *                                 /  \                  gamma = 0     *
c    *                                /    \                               *
c    *                               /      \                              *
c    *                              /        \                             *
c    *                             /          \                            *
c    *                            /            \                           *
c    *                     prolate            oblate                       *
c    *                     x2=z2<y2           x2<y2=z2                     *
c    *                                                                     *
c    ***********************************************************************
c    *                                                                     *
c    *  Symmetry axis,  z axis iq2 = 0        gamma =   0 degrees          *
c    *                  y axis iq1 = 0        gamma =  60 degrees          *
c    *                  x axis iq1+1q2 = 0    gamma = 120 degrees          *
c    *                                                                     *
c    ***********************************************************************

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4   head
      character*20  afor
      character(2)  num
      logical ok

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (epsm8=1.0d-8,epsm4=1.0d-4,epsm3=1.0d-3)
      parameter (tp01=0.01d0,tp02=0.02d0,tp1=0.1d0,tp5=0.5d0)
      parameter (t100=100.0d0,t180=180.0d0)
      parameter (t1000=1000.d0,tp4=0.4d0)

      common /Lag  / ilag, iAnaLag
      common /conv / econv1,econv2,iconv
      common /condat/ epsq,epse,epsdh,epsf,epsm,iprintconver
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),
     1               imtd,imtg,icutq
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
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /enert/ c14,c15,t14,t15
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     3              ,hbar,hbm(2),xm(3),afor
      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pairn/ vn,vp,rangen,rangep
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /vo8  / delqst,q1tst,q2tst,qrfintst
     1                     ,q1nst,q2nst,qrfinnst
     2                     ,q1pst,q2pst,qrfinpst

c................................................................. read formats
  101 format (20a4)
 1010 format (20a20) 
  102 format (5e15.8)
  103 format (6i5)

c................................................................ write formats
  200 format (/,'  __________________________________________________ ',
     1        /,' |                                                  |',
     2        /,' |  Program ev8      CPC version 2                  |',
     3        /,' |  August 2014                                     |',
     7        /,' |__________________________________________________|')
  201 format (//,' compilation parameters: ', / , ' mx=', i5,' my=', i5
     1          ,' mz=',i5,
     2         /,' Maximum number of extra coulomb points: ',i5,
     3         /,' Total number of wavefunctions: ', i5  )
  203 format (//,' ',78('_'),/,' run information',//,20a4)
  204 format (/,' nx=',i5,6x,'  ny=',i5,6x,'  nz=',i5,
     1        /,' dx=',f8.3,' Fm',
     2        /,' dt=',f8.3,' 10-22s,  number of iterations =',i5)
 2041 format   (' number of extra points for coulomb ',
     1        /,' nx=',i5,6x,'  ny=',i5,6x,'  nz=',i5)
 2042 format (' Order of the laplacian discretisation for Coulomb: '
     &        , i5/)
 2043 format (' Convergence output: ', A20)
 2044 format (' Invalid value for iverb, reverting to Benders output. ')
  205 format (/,' nxmu   = ',i5,
     1        /,' ndiag  = ',i5)
  206 format (  ' nprint = ',i5)
  207 format (  ' EDF parameters  from a Skyrme force, nfunc = ',i5)
  208 format (  ' njmunu = ',i5,' -- spin tensor is self-consistent --',
     1        /,22x,' -- taken from central Skyrme force only --')
 2080 format (/,' nfunc  = ',i5,' -- EDF option --',/)
 2081 format (  ' njmunu = ',i5,' -- self-consistent tensor force --',
     1        /,15x,' -- taken from central+tensor Skyrme force --')
 2082 format (  ' njmunu = ',i5,' -- spin tensor is self-consistent --',
     1        /,15x,'  -- Jij Jij term from functional --')
 2083 format (  ' njmunu = ',i5,' -- self-consistent tensor force --',
     1        /,15x,' -- Jij Jij + Jij Jj terms from functional --')
 2084 format (  ' njmunu = ',i5,' no tensor terms in the functional')
  209 format (  ' ncm2   = ',i5,' -- 1-body c.m. correction only --')
 2090 format (  ' ncm2   = ',i5,' -- no c.m. correction at all --')
 2091 format (  ' ncm2   = ',i5,' -- perturbative 1-& 2-body c.m. --') 
  210 format (  ' ncm2   = ',i5,' -- 2-body c.m. is self-consistent --')
  211 format (  ' nmass  = ',i5,' m_n = m_p')
  212 format (  ' nmass  = ',i5,' experimental nucleon masses')
 2012 format (  ' nmass  = ',i5,' h^2/2m read from data',/,
     1          ' h^2/2m_n = ',f10.6,' h^2/2m_p = ',f10.6,/)
 2013 format (  ' ncoex  = ',i5,' Coulomb exchange',
     1          ' in Slater approximation')
 2014 format (  ' ncoex  = ',i5,' no Coulomb exchange')
  220 format (/,' number of wave-functions n=',i5,'  p=',i5)
  221 format   (' nucleus n =',i5,'  z =',i5,'  a =',i5)
  224 format (/,' blocking or filling of level ',i5,' attempted',/)
  225 format (/,' filling of neutron level ',i5,' attempted',
     1          ' but neutron number is even N = ',i5,/)
  226 format (/,' filling of proton level ',i5,' attempted',
     1          ' but proton number is even Z = ',i5,/)
  227 format (/,' blocking of neutron level ',i5,' attempted',
     1          ' but neutron number is odd N = ',i5,/)
  228 format (/,' blocking of proton level ',i5,' attempted',
     1          ' but proton number is odd Z = ',i5,/)
  230 format (/,' Skyrme force: ',a20,'  (new set of parameters)')
  231 format (/,' Skyrme force: ',a20)
  240 format   ('  t0 =',f13.6,' x0 =',f13.6,/,
     1          '  t1 =',f13.6,' x1 =',f13.6,/,
     2          '  t2 =',f13.6,' x2 =',f13.6,/,
     3          '  t3a=',f13.6,' x3a=',f13.6,'  et3a =',f10.6,/,
     4          '  t3b=',f13.6,' x3b=',f13.6,'  et3b =',f10.6,/,
     5          '  te =',f13.6,' to =',f13.6,/,
     6          '  w  =',f13.6,' wq =',f13.6,/)
  250 format (/,' npair  = ',i5,' -- hartree-fock --')
  251 format (/,' npair  = ',i5,' -- BCS seniority force --')
  252 format   (' gn  =',f8.3,'/(11+N) MeV,  encut=',f8.3,'MeV',
     2        /,' gp  =',f8.3,'/(11+Z) MeV,  epcut=',f8.3,'MeV',
     3        /,' dcut=',f8.3,'MeV')
  253 format   (' -- BCS ground state --')
  256 format   (' index of blocked state taken from data')
  257 format   (' index of blocked state taken from tape,',
     1           ' value on data is ',i5)
  258 format (/,' no blocked state found on tape,',
     1           ' will use index from data instead')
  265 format (/,' npair  = ',i5,' -- BCS constant delta --')
  266 format (  ' deltan=',f8.3,'  encut=',f8.3,
     1        /,' deltap=',f8.3,'  epcut=',f8.3,
     2        /,' dcut  =',f8.3)
  267 format (/,' npair  = ',i5,' -- BCS+LN  seniority force --')
  268 format (/,' npair  = ',i5,' -- BCS delta force --')
  269 format (  ' vn   =',f8.3,'MeV Fm3,  encut=',f8.3,'MeV',
     2        /,' vp   =',f8.3,'MeV Fm3,  epcut=',f8.3,'MeV',
     3        /,' dcut =',f8.3,'MeV',
     4        /,' rhoc =',f8.3)
  270 format (/,' npair  = ',i5,' -- BCS+LN delta force --')
  273 format (  ' vn     = ',f12.6,' MeV    range =',f8.3,' Fm',
     2        /,' vp     = ',f12.6,' MeV    range =',f8.3,' Fm')
  275 format (/,' cutoff only above the fermi level ')
  276 format (/,' cutoff above and below the fermi level ')
  280 format (/,' imtd   = ',i5,
     1          ' (0,2)=fixed,(1,3)=readjusted,(0,1)=T=0,(2,3)=T=1',
     2        /,' ifrt   = ',i5,' strength in data: 0,  on tape: 1',
     3        /,' imtg   = ',i5,' constraint on (q,gamma): 0',
     4          ' on q only: 1')
  281 format   (' ral    =',f5.2,' slow-down of constraint potential',
     1        /,' epscst =',f5.2,' constraint readjustment')
  282 format   (' cut-off:(1+exp((rcut-r)/acut)',
     1        /,' rcut =',f12.3,'fm, acut =',f12.3,'fm')
  283 format   (' cut-off:tanh(rho/rhoc) rho c=',f12.4,'fm-3')
  284 format   (' cut-off:(1+exp((Delta R-rcut)/acut)^-1'
     1        /,' rcut =',f12.3,'fm, acut =',f12.3,'fm')
  290 format (/,' constraints data',
     1        /,'  cr2  = ',e12.4,' rms   = ',e12.4,
     2        /,'  cq2  = ',e12.4,' delq  = ',e12.4,
     3        /,'  q1t  = ',e12.4,' q2t   = ',e12.4)
  291 format (/,' actual constraints',
     1  /,'  cqr =',e12.4,
     2          '  Requested r   = ',e12.4,'  Readjusted r   = ',e12.4,
     3  /,'  cq2 =',e12.4,
     4          '  Requested Qxx = ',e12.4,'  Readjusted Qxx = ',e12.4,  
     5  /,19x  ,'  Requested Qyy = ',e12.4,'  Readjusted Qyy = ',e12.4,
     6  /,19x  ,'  Requested Qzz = ',e12.4,'  Readjusted Qzz = ',e12.4)
  292 format (/,'  constraints data',
     1        /,'  cr2  = ',e12.4,' rmsn  = ',e12.4,' rmsp  = ',e12.4,
     2        /,'  cq2  = ',e12.4,' delq  = ',e12.4,
     3        /,'  q1n  = ',e12.4,' q2n   = ',e12.4,
     4        /,'  q1p  = ',e12.4,' q2p   = ',e12.4)
  293 format (/,' actual constraints',
     1  /,      '  cqr =',e12.4,
     2          '  Requested N r   = ',e12.4,
     2          '  Readj. N r   = ',e12.4,
     3  /,19x  ,
     4          '  Requested P r   = ',e12.4,
     4          '  Readj. P r   = ',e12.4,
     5  /,      '  cq2 =',e12.4,
     6          '  Requested N Qxx = ',e12.4,
     7          '  Readj. N Qxx = ',e12.4,
     8  /,19x  ,'  Requested N Qyy = ',e12.4,
     9          '  Readj. N Qyy = ',e12.4,
     1  /,19x  ,'  Requested N Qzz = ',e12.4,
     2          '  Readj. N Qzz = ',e12.4,
     3  /,19x  ,'  Requested P Qxx = ',e12.4,
     4          '  Readj. P Qxx = ',e12.4,
     5  /,19x  ,'  Requested P Qyy = ',e12.4,
     6          '  Readj. P Qyy = ',e12.4,
     7  /,19x  ,'  Requested P Qzz = ',e12.4,
     8          '  Readj. P Qzz = ',e12.4)
  294 format (/,' Convergence tests : '
     1        /,'     (delta e)/e                < ',es10.2,
     2        /,'     sum v^2 (<h^2>-<h>^2)  / a < ',es10.2,
     3        /,'     (q2-q2con)/q2con           < ',es10.2,
     4        /,'     delta e_fermi              < ',es10.2)
  295 format (
     1 /,' actual constraints on intrinsic moment only',
     2 /,'  cqr =',e12.4,'  Requested R  = ',e12.4,
     3                   '  Readj. R   = ',e12.4,
     4 /,'  cq2 =',e12.4,'  Requested Q2 = ',e12.4,
     5                   '  Readj. Q2  = ',e12.4)
  300 format (  ' EDF coefficients read: '       , /,
     1          '  b14=',f13.6,' b15=', f13.6 , /,
     1          '  b16=',f13.6,' b18=', f13.6  )

c..............................................................................
      pi   = tt4*atan2(one,one)
      traf = pi / t180

c     ........................................ data -- title of the calculation
      read  101,head
      print 200
      print 201,mx,my,mz,mc,mw

c     ............................................... read input wave-functions
c                                        checking for the existence of the file
      irb    = 12
      iprint =  1
      write (num,'(i2.2)') irb
      inquire(exist= ok , file= 'fort.'//num)
      if (.not. ok) then
        call stp("EV8 can't find fort.12!")
      endif
      open (irb,form='unformatted',file='fort.'//num)
      call lec8 (iprint)
      close (irb)

c     ........... characteristics of the mesh, the time step and the total time
      print 203,head
      read  103,npx,npy,npz,iCoul
      read  102,dt

      iLag=0

c     ..................................... characteristics of iteration scheme
      nxmu = 0
      read  103,nitert,nxmu,ndiag

c     ............................................ parameters controling output
      read  103,nprint, iverb
      print 204,nx,ny,nz,dx,dt,nitert

      select case(iverb)
        case (1) !Benders output
          print 2043, trim('Bender')
          iprintCoul=1
          iprintconver=1

        case (0) !Public output
          print 2043, trim('Normal')
          iprintCoul=0
          iprintconver=0

        case DEFAULT !
          print 2044
          iprintCoul=0
          iprintconver=0

      end select

c     .................................................. print and control data
      if (npx.gt.mc) then
        call stp (
     1    'Extra points in the x direction for Coulomb
     2     more than maximum!')
      endif
      if (npy.gt.mc) then
        call stp (
     1  'Extra points in the y direction for Coulomb
     2  more than maximum!')
      endif
      if (npz.gt.mc)  then
        call stp (
     1  'Extra points in the z direction for Coulomb
     2  more than maximum!')
      endif

      if (iCoul.lt.0 .or. iCoul.gt.2) then
        call stp ('Invalid value for icoul. It can only be 0,1 or 2!')
      endif

      if (iCoul.eq.0) iCoul=2 !Defaulting to second order

      nnx = mx + npx
      nny = my + npy
      nnz = mz + npz
      print 2041,npx,npy,npz
      print 2042, iCoul

      if (dt.lt.epsm8) then
        call stp ('Too small value for dt !')
      endif
      if (nitert.lt.0) itert = 0
      nitert = abs(nitert)
      if (nxmu.eq.0) nxmu = 25
      print 205,nxmu,ndiag
      if (nxmu.lt.0.or.nxmu.gt.100) then
        call stp (
     1  'Invalid value of nxmu! (Values lie between 0 and 100')
      endif
      print 206,nprint
      if (nprint.le.0) then
        call stp (
     1  'Invalid value for nprint! (It should be bigger than 0)')
      endif

c     ......................................... number of particle and orbitals
      read  103,npn,npp
      xaf = npn + npp
      print 220,nwaven,nwavep
      print 221,npn,npp,npn+npp
      if (npp.gt.2*nwavep) then
        call stp (
     1  'Insufficient wavefunctions to accomodate all protons!')
      endif
      if (npn.gt.2*nwaven) then
        call stp (
     1  'Insufficient wavefunctions to accomodate all neutrons!')
      endif
c     ................................................. skyrme force parameters
      read  1010,afor
      print 231, trim(afor)
      call lecfor (kfor)
      if (kfor.eq.0) then
        read 103,ncm2,nfunc,nmass,ncoex
        read 102,t0,x0
        read 102,t1,x1,t2,x2
        read 102,t3a,x3a,yt3a
        read 102,t3b,x3b,yt3b
        ndd = 1
        if (t3b.eq.0.0) ndd = 0
        njmunu = 2
        if (nfunc.eq.0) then
          read 102,wso
          read 102,te,to
          wsoq = wso
          if (te.eq.0.0.and.to.eq.0.0) njmunu = 1 
        else
          read 102,wso,wsoq
          read 102,b14,b15,b16,b17
          if (b16.eq.0.0.and.b17.eq.0.0) then
            if (b14.eq.0.0.and.b15.eq.0.0) then
              njmunu = 0
            else
              njmunu = 1
            endif
          endif  
        endif
        if (nmass.eq.2) then
          read 102,hbm(1),hbm(2)
        endif
        !print 230,afor
      else
        !print 231,afor
      endif
      call lecmas
      print 240,t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                           ,t3b,x3b,yt3b,
     2          te,to,wso,wsoq
      if(nfunc .eq. 1) print 300, b14,b15,b16,b17

      if (nfunc .eq.0) then         
        print 207,nfunc
        if (njmunu.eq.1) print 208 ,njmunu
        if (njmunu.eq.2) print 2081,njmunu
      else
        print 2080,nfunc
        if (njmunu.eq.0) print 2084,njmunu
        if (njmunu.eq.1) print 2082,njmunu
        if (njmunu.eq.2) print 2083,njmunu
      endif
      if (ncm2  .eq.-2)               print 2091,ncm2
      if (ncm2  .eq.-1)               print 2090,ncm2
      if (ncm2  .eq. 0)               print 209 ,ncm2
      if (ncm2  .eq. 1)               print 210 ,ncm2
      if (nmass .eq.0)                print 211 ,nmass
      if (nmass .eq.1)                print 212 ,nmass
      if (nmass .eq.2)                print 2012,nmass,
     1                                0.5d0*hbm(1),0.5d0*hbm(2)
      if (ncoex .eq.0)                print 2013,ncoex
      if (ncoex .eq.1)                print 2014,ncoex

c     .............................................................. test input
      if (nfunc .lt. 0.or.nfunc .gt.1) then
        call stp(' Invalid value of nfunc! It should be 0 or 1.')
      endif
      if (njmunu.lt. 0.or.njmunu.gt.2) then
        call stp(' Invalid value of njmunu! It should be 0, 1 or 2.')
      endif
      if (ncm2  .lt.-2.or.ncm2  .gt.1)then
        call stp(' Invalid value of ncm2! It should be -2,-1, 0 or 1.')
      endif
      if (ndd   .lt. 0.or.ndd   .gt.1) then
        call stp(' Invalid value of ndd! It should be 0 or 1.')
      endif
      if (nmass .lt. 0.or.nmass .gt.2) then
        call stp(' Invalid value of nmass! It should be 0, 1 or 2.')
      endif
      if (ncoex .lt. 0.or.ncoex .gt.1) then
        call stp(' Invalid value of ncoex! It should be 0 or 1.')
      endif

c     ..................... occupation of the orbitals pairing force parameters
      read 103,npair,icut
      
      icut = 0
      ifor = 0 
      mtqp = 0 
      ilqp = 0
      
      if (npair.le.5) then
        read 102,gn,encut,delmax(1),alpha
        read 102,gp,epcut,delmax(2)
        dcut = tp5
      endif
      if (npair.eq.0) then
      	!Reset the Lipkin-Nogami occupations
      	v22 = v2
        print 250,npair
      endif
      xcut = zero
      if (icut.ne.0) xcut = one
      if (npair.eq.1) then
        !Reset the Lipkin-Nogami occupations
      	v22 = v2
        print 251,npair
        print 252,gn,encut,gp,epcut,dcut
      endif
      if (npair.eq.2) then
      	!Reset the Lipkin-Nogami occupations
      	v22 = v2
        print 265,npair
        print 266,delmax(1),encut,delmax(2),epcut,dcut
      endif
      if (npair.eq.3) then
        print 267,npair
        print 252,gn,encut,gp,epcut,dcut
      endif
      if (npair.eq.4) then
        !Reset the Lipkin-Nogami occupations
      	v22 = v2        
        print 268,npair
        print 269,gn,encut,gp,epcut,dcut,alpha
      endif
      if (npair.eq.5) then
        print 270,npair
        print 269,gn,encut,gp,epcut,dcut,alpha
      endif

      if (npair.le.5) then
         alpha = alpha/0.16d0
        if (icut.eq.0) print 275
        if (icut.ne.0) print 276
      endif

      if (npair.gt.5) call stp ('npair.gt.5 !')

c............................................................. convergency test
      read  102,epse,epsdh,epsq,epsf
      print 294,epse,epsdh,epsq,epsf

c     ............................................................. constraints
      imtdvo = imtd
      read  103,imtd,ifrt,imtg,icutq
      if ((imtd.ge.2).and.(imtdvo.le.1)) ifrt = 0
      print 280,imtd,ifrt,imtg
      if ((imtd.lt.0).or.(imtd.gt.3)) then
        call stp(' Invalid value of imtd! It should be 0,1,2 or 3.')
      endif
      if ((ifrt.lt.0).or.(ifrt.gt.1)) then
        call stp(' Invalid value of ifrt! It should be 0 or 1.')
      endif
      read 102,ral,epscst,rcut
      if    (ral.eq.zero)    ral = tp1
      if (epscst.eq.zero) epscst = tp02
      print 281,ral,epscst
      if (( ral.lt.zero).or.(   ral.gt.one)) then
        call stp(
     1     ' Invalid value of ral! It should be between 0 and 1.')
      endif
      if ((epscst.lt.zero).or.(epscst.gt.one))then
        call stp(
     1     ' Invalid value of epcst! It should be between 0 and 1.')
      endif
      if (icutq.lt.0 .or. icutq.gt.2) then
        call stp(
     1     ' Invalid value for icutq! It should be either 0,1 or 2')
      endif
      if(rcut.lt.0.0d0) then
        call stp(
     1     ' Invalid value for rcut! It should be strictly positive!')
      endif

      if(rcut.eq.0.0d0) rcut=100.0d0

      select case(icutq) 
      case(0)
        acut=0.4d0
        print 284, rcut, acut
      case(1)
        print 283, rcut
      case(2)
        acut = 2.0d0 * dx
        print 282, rcut,acut 
      end select

c     ............................................ data -- isoscalar constraint
      cqrvo = cqr
      cq2vo = cq2
      if (imtd.le.1) then
        read 102,cqr,rms
        read 102,cq2,delq
        read 102,q1t,q2t
        qrfint = float(npn+npp)*rms**2
        if (abs(qrfint-qrfintst).gt.epsm8) itert=0
        if       (abs(q1t-q1tst).gt.epsm4) itert=0
        if       (abs(q2t-q2tst).gt.epsm4) itert=0
        if     (abs(delq-delqst).gt.epsm4) itert=0
        print 290,cqr,rms,cq2,delq,q1t,q2t
        if (cqr.lt.zero)  call stp ('cqr !')
        if (rms.lt.zero)  call stp ('rms !')
        if (cq2.lt.zero)  call stp ('cq2 !')
        if (delq.lt.zero) call stp ('delq !')
        qxfint =-delq*(q1t-  q2t)/two
        qyfint =-delq*(q1t+2*q2t)/two
        qzfint = delq*(2*q1t+q2t)/two
        q2fin  = delq*sqrt(q1t**2+q2t**2+q1t*q2t)
        g2fin  = zero
        if ((abs(q1t).gt.epsm3).or.(abs(q2t).gt.epsm3))
     1    gtc0 = atan2(q2t*sqrt(tt3),two*q1t+q2t)/traf
        if (cqr.ne.zero) then
          if (abs(qrcstt).lt.tp01) qrcstt = qrfint
        endif
        if (ifrt.ne.1) then
          qrcstt = qrfint
          qxcstt = qxfint
          qycstt = qyfint
          qzcstt = qzfint
          q2cst  = q2fin
          g2cst  = g2fin
        endif
        if (cqrvo.eq.zero) then
          qrcstt = qrfint
        endif
        if (cq2vo.eq.zero) then
          qxcstt = qxfint
          qycstt = qyfint
          qzcstt = qzfint
        endif
        if (imtg.eq.0) then
          print 291,cqr,sqrt(qrfint/(npn+npp)),sqrt(qrcstt/(npn+npp)),
     1            cq2,qxfint,qxcstt,qyfint,qycstt,qzfint,qzcstt
        else
          print 295,cqr,sqrt(qrfint/(npn+npp)),sqrt(qrcstt/(npn+npp)),
     1            cq2,q2fin,q2cst
        endif

c     ............................................ data -- isovector constraint
      else
        read 102,cqr,rmsn,rmsp
        read 102,cq2,delq
        read 102,q1n,q2n,q1p,q2p
        qrfinn = float(npn)*rmsn**2
        qrfinp = float(npp)*rmsp**2
        if (abs(qrfinn-qrfinnst).gt.epsm8) itert=0
        if       (abs(q1n-q1nst).gt.epsm4) itert=0
        if       (abs(q2n-q2nst).gt.epsm4) itert=0
        if (abs(qrfinp-qrfinpst).gt.epsm8) itert=0
        if       (abs(q1p-q1pst).gt.epsm4) itert=0
        if       (abs(q2p-q2pst).gt.epsm4) itert=0
        if     (abs(delq-delqst).gt.epsm4) itert=0
        print 292,cqr,rmsn,rmsp,cq2,delq,q1n,q2n,q1p,q2p
        if  (cqr.lt.zero) call stp ('cqr !')
        if (rmsn.lt.zero) call stp ('rmsn !')
        if (rmsp.lt.zero) call stp ('rmsp !')
        if  (cq2.lt.zero) call stp ('cq2 !')
        if (delq.lt.zero) call stp ('delq !')
        qxfinn =-delq*(q1n-  q2n)/two
        qyfinn =-delq*(q1n+2*q2n)/two
        qzfinn = delq*(2*q1n+q2n)/two
        qnc0   = delq*sqrt(q1n**2+q2n**2+q1n*q2n)
        gnc0   = zero
        if ((abs(q1n).gt.epsm3).or.(abs(q2n).gt.epsm3))
     1    gnc0 = atan2(q2n*sqrt(tt3),two*q1n+q2n)/traf
        qxfinp =-delq*(q1p-  q2p)/two
        qyfinp =-delq*(q1p+2*q2p)/two
        qzfinp = delq*(2*q1p+q2p)/two
        qpc0   = delq*sqrt(q1p**2+q2p**2+q1p*q2p)
        gpc0   = zero
        if ((abs(q1p).gt.epsm3).or.(abs(q2p).gt.epsm3))
     1    gpc0 = atan2(q2p*sqrt(tt3),two*q1p+q2p)/traf
        if (cqr.ne.zero) then
          if (abs(qrcstn).lt.tp01) qrcstn = qrfinp
          if (abs(qrcstp).lt.tp01) qrcstn = qrfinp
        endif
        if (ifrt.ne.1) then
          qrcstn = qrfinn
          qxcstn = qxfinn
          qycstn = qyfinn
          qzcstn = qzfinn
          qrcstp = qrfinp
          qxcstp = qxfinp
          qycstp = qyfinp
          qzcstp = qzfinp
        endif
        if (cqrvo.eq.zero) then
          qrcstn = qrfinn
          qrcstp = qrfinp
        endif
        if (cq2vo.eq.zero) then
          qxcstn = qxfinn
          qycstn = qyfinn
          qzcstn = qzfinn
          qxcstp = qxfinp
          qycstp = qyfinp
          qzcstp = qzfinp
        endif

        print 293,cqr,sqrt(qrfinn/(npn)),sqrt(qrcstn/(npn))
     1           ,sqrt(qrfinp/(npp)),sqrt(qrcstp/(npp)),
     1            cq2,qxfinn,qxcstn,qyfinn,qycstn,qzfinn,qzcstn,
     2                qxfinp,qxcstp,qyfinp,qycstp,qzfinp,qzcstp
      endif

c     .......................................................... time evolution
      call evolve

c     ..................... storage of the wave functions at the last time step
      if (iconv.ge.0) then
        open  (13,form='unformatted',file='fort.13')
        call writ8
        close (13)
      endif

      call stp (' ***** THE END ***** ')

      end

c______________________________________________________________________________
      subroutine lec8 (iprint)

c..............................................................................
c     single-particle wave functions and other information read from fort.irb .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character(20):: afor
      character(4):: head
      character(4):: headd

      parameter (zero=0.0d0)

      common /champ/ qxxn,qyyn,qzzn,qrrn
     1              ,qxxp,qyyp,qzzp,qrrp
     2              ,qxxt,qyyt,qzzt,qrrt
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),
     1               imtd,imtg,icutq
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
  271 format   ('  gn =',f8.3,'/(11+n),  encut=',f6.3,
     2        /,'  gp =',f8.3,'/(11+z),  epcut=',f6.3,'  dcut =',f6.3)
  272 format   ('  deltan =',f8.3,'  encut =',f6.3,
     1        /,'  deltap =',f8.3,'  epcut =',f6.3,'  dcut =',f6.3)
  273 format   ('  vn = ',f8.3,'  encut =',f6.3,
     2        /,'  vp = ',f8.3,'  epcut =',f6.3,'  dcut =',f6.3,
     4          '  alpha =',f6.3)
  275 format   ('    vn = ',f12.6,'  range =',f6.3,
     2        /,'    vp = ',f12.6,'  range =',f6.3)
  281 format   (' cutoff only above the fermi level ')
  282 format   (' cutoff above and below the fermi level ')
  290 format   (/, 'Constraint data from fort.12: ', / , ' imtd  = ',i5,
     1          ' (0,2)=fixed,(1,3)=readjusted,(0,1)=T=0,(2,3)=T=1')
  293 format   (' delq  = ',e12.4,
     3        /,' q1t   = ',e12.4,'  q2t   = ',e12.4)

  294 format (/,' actual constraints',
     1  /,'  cqr =',e12.4,
     2          '  Requested r   = ',e12.4,'  Readjusted r   = ',e12.4,
     3  /,'  cq2 =',e12.4,
     4          '  Requested Qxx = ',e12.4,'  Readjusted Qxx = ',e12.4,  
     5  /,19x  ,'  Requested Qyy = ',e12.4,'  Readjusted Qyy = ',e12.4,
     6  /,19x  ,'  Requested Qzz = ',e12.4,'  Readjusted Qzz = ',e12.4)
  295 format   (' delq  = ',e15.8,
     3        /,' q1n  = ',e12.4,'  q2n   = ',e12.4,
     4        /,' q1p  = ',e12.4,'  q2p   = ',e12.4)
  296 format (/,' actual constraints',
     1  /,      '  cqr =',e12.4,
     2          '  Requested N r   = ',e12.4,
     2          '  Readj.  N r   = ',e12.4,
     3  /,19x  ,
     4          '  Requested P  r  = ',e12.4,
     4          '  Readj.  P  r  = ',e12.4,
     5  /,      '  cq2 =',e12.4,
     6          '  Requested N Qxx = ',e12.4,
     7          '  Readj.  N Qxx = ',e12.4,
     8  /,19x  ,'  Requested N Qyy = ',e12.4,
     9          '  Readj.  N Qyy = ',e12.4,
     1  /,19x  ,'  Requested N Qzz = ',e12.4,
     2          '  Readj.  N Qzz = ',e12.4,
     3  /,19x  ,'  Requested P Qxx = ',e12.4,
     4          '  Readj.  P Qxx = ',e12.4,
     5  /,19x  ,'  Requested P Qyy = ',e12.4,
     6          '  Readj.  P Qyy = ',e12.4,
     7  /,19x  ,'  Requested P Qzz = ',e12.4,
     8          '  Readj.  P Qzz = ',e12.4)
  240 format   (/, 'Force used on fort.12: ',/,
     1          '  t0 =',f13.6,' x0 =',f13.6,/,
     1          '  t1 =',f13.6,' x1 =',f13.6,/,
     2          '  t2 =',f13.6,' x2 =',f13.6,/,
     3          '  t3a=',f13.6,' x3a=',f13.6,'  et3a  =',f10.6,/,
     4          '  t3b=',f13.6,' x3b=',f13.6,'  et3b  =',f10.6,/,
     5          '  w  =',f13.6,' wq =',f13.6,/,
     6          '  te =',f13.6,' to =',f13.6)
  302 format ( /,' Functional parameters used on fort.12')
  311 format (  '  b1  =',f13.6,' b2  =',f12.6)
  312 format (  '  b3  =',f13.6,' b4  =',f12.6)
  313 format (  '  b5  =',f13.6,' b6  =',f12.6)
  314 format (  '  b7  =',f13.6,' b8  =',f12.6,'  et3a  =',f10.6)
  315 format (  '  b7a =',f13.6,' b8a =',f12.6,'  et3a  =',f10.6,/,
     1          '  b7b =',f13.6,' b8b =',f12.6,'  et3b  =',f10.6)
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
      if   (nnx.gt.mx) call stp (' nnx > mx !')
      if   (nny.gt.my) call stp (' nny > my !')
      if   (nnz.gt.mz) call stp (' nnz > mz !')
      nx = nnx
      ny = nny
      nz = nnz
      if (nnx.ne.mx) nx = mx
      if (nny.ne.my) ny = my
      if (nnz.ne.mz) nz = mz

      if (((iver.lt.2).or.(iver.gt.6)).and.(iver.ne.20)) then
       call stp (' iver !')
      endif
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
     1             x3a,yt3a,t3b,x3b,yt3b,wso,wsoq
        read (irb) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,ntqp
        read (irb) (eqp(i),i=1,nwave)
        read (irb) (delta(i),i=1,nwave)
      endif

      if (iver.eq.6) then
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
c              Note that this line reads two numbers more than iver=5 version:
c              te and to!
        read (irb) t0,x0,t1,x1,t2,x2,t3a,
     1              x3a,yt3a,t3b,x3b,yt3b,wso,wsoq,te,to

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
c             .............................................ATTENTION
c             B5 & B6 as used in the code are not equal to the ones
c             discussed in all of the articles; they carry a 
c             relative minus sign in this code!
c             ......................................................           
            print 313,-b5,-b6
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
          if (cqr .eq. 0.d0 ) then
            qrcstt = 0.d0
            qrcstp = 0.d0
            qrcstn = 0.d0
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
            print 294,cqr,sqrt(qrfint/(npn+npp)),sqrt(qrcstt/(npn+npp)),
     1                cq2,qxfint,qxcstt,qyfint,qycstt,qzfint,qzcstt
          else
            print 295,delqst,q1nst,q1pst,q2nst,q2pst
            print 296,cqr,sqrt(qrfinn/(npn)),sqrt(qrcstn/(npn))
     1                   ,sqrt(qrfinp/(npp)),sqrt(qrcstp/(npp)),
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

        if (npair.eq.1) print 271,gn,encut,gp,encut,dcut
        if (npair.eq.2) print 272,delmax(1),encut,delmax(2),epcut,dcut
        if (npair.eq.3) print 271,gn,encut,gp,encut,dcut
        if ((npair.eq.4).or.(npair.eq.5)) then
c         Note that at this point in the program alpha has not yet been divided
c         by 0.16.        
          print 273,gn,encut,gp,encut,dcut,1.0d0/alpha
        endif

        if ((npair.ge.1).and.(npair.le.5)) then
          if (xcut.eq.zero) print 281
          if (xcut.ne.zero) print 282
        endif

      endif

      read (irb) (kparz(i),i=1,nwave)
      read (irb) (esp1 (i),i=1,nwave)
      read (irb) (v2   (i),i=1,nwave)

      if (iver.ge.3) then
        read (irb) (v22(i),i=1,nwave)
      endif

      read (irb) ((npar(i,it),i=1,2),it=1,2)

      if (nnx.eq.mx. and. nny.eq.my .and. nnz.eq.mz) then
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


        w1(:,:,:) = 0.d0
        w2(:,:,:) = 0.d0
        w3(:,:,:) = 0.d0
        w4(:,:,:) = 0.d0
        do iwave=1,nwave
          read (irb) (((w1(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
          read (irb) (((w2(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
          read (irb) (((w3(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)
          read (irb) (((w4(i,j,k),i=1,nnx),j=1,nny),k=1,nnz)

          call scopy (mq,w1,a(1,iwave))
        enddo

      endif

      return
      end subroutine lec8
c______________________________________________________________________________
      subroutine lecfor (kfor)

c..............................................................................
c      This subroutine tries to recognize the input force name from among the 
c      forces Ev8 knows.
c      The blockdata for a force has the following structure:
c      -------------------------------------------------------
c       t0 , t1  ,   t2  ,   t3a, t3b
c       x0 , x1  ,   x2  ,  x3a , x3b
c       na , da  ,   n   ,  db  , 0.0
c       wso, wsoq, ncm2  , nfunc, nmass
c       te , to  , hbm(1), hbm(2),ncoex
c       b14, b15 , b16   , b17
c      -------------------------------------------------------
c       Where:
c       * hbm(i) = 0 when nmass != 2
c       * na/da = yt3a and nb/db = yt3b
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      character*20 afor
      character*20 UPafor

      parameter (zero=0.0d0)

      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /skm  / prm(30,34)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b

      kfor = 0
      call to_upper(afor, UPafor)
      if(trim(UPafor).eq. 'SKM*') kfor=1
      if(trim(UPafor).eq. 'SKM ') kfor=2
      if(trim(UPafor).eq. 'SIII') kfor=3
      if(trim(UPafor).eq. 'SKA ') kfor=4
      if(trim(UPafor).eq. 'SGII') kfor=5
      if(trim(UPafor).eq. 'SLYIIIK10') kfor=6
      if(trim(UPafor).eq. 'SLY4') kfor=7
      if(trim(UPafor).eq. 'SLY5') kfor=8
      if(trim(UPafor).eq. 'SLY6') kfor=9
      if(trim(UPafor).eq. 'SLY7') kfor=10
      if(trim(UPafor).eq. 'SKP ') kfor=11
      if(trim(UPafor).eq. 'SKI4') kfor=12
      if(trim(UPafor).eq. 'T22 ') kfor=13
      if(trim(UPafor).eq. 'T24 ') kfor=14
      if(trim(UPafor).eq. 'T26 ') kfor=15
      if(trim(UPafor).eq. 'T42 ') kfor=16
      if(trim(UPafor).eq. 'T44 ') kfor=17
      if(trim(UPafor).eq. 'T46 ') kfor=18
      if(trim(UPafor).eq. 'T62 ') kfor=19
      if(trim(UPafor).eq. 'T64 ') kfor=20
      if(trim(UPafor).eq. 'T65 ') kfor=21
      if(trim(UPafor).eq. 'T66 ') kfor=21
      if(trim(UPafor).eq. 'SKX ') kfor=23
      if(trim(UPafor).eq. 'SLY4T') kfor=24
      if(trim(UPafor).eq. 'SLY5+T') kfor=25
      if(trim(UPafor).eq. 'SKI3') kfor=26
      if(trim(UPafor).eq. 'SVMIN ') kfor=27
      if(trim(UPafor).eq. 'UNEDF0 ') kfor=28
      if(trim(UPafor).eq. 'UNEDF1 ') kfor=29
      if(trim(UPafor).eq. 'UNEDF1SO') kfor=30
      if(trim(UPafor).eq. 'SLYIIIK07') kfor=32
      if(trim(UPafor).eq. 'SLYIIIK08') kfor=33
      if(trim(UPafor).eq. 'SLYIIIK09') kfor=34
  
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

      xncm2  = prm (18, kfor)
      
      xnfunc = prm (19, kfor)
      xnmass = prm (20, kfor)

      te     = prm (21, kfor)
      to     = prm (22, kfor)

      if(xnmass.eq.2.0d0) then
        nmass=2
        hbm(1) = prm( 23, kfor)
        hbm(2) = prm( 24, kfor)
      elseif(xnmass.eq.1.0d0) then
        nmass=1
      else
        nmass=0
      endif
      
      scoex = prm( 25,kfor)
      if(scoex.eq.1.0d0) then
        ncoex=1
      else
        ncoex=0
      endif

c.................................... Determining other force parameters
      if(xncm2 .eq. -1.0d0) then
        ncm2=-1
      elseif(xncm2 .eq. 0.0d0) then
        ncm2=0
      elseif(xncm2 .eq. -2.0d0) then
        ncm2 =-2
      else
        ncm2=1
      endif

      if(xnfunc .eq. 1.0d0) then
        nfunc=1
        b14 = prm( 26, kfor)
        b15 = prm( 27, kfor)
        b16 = prm( 28, kfor)
        b17 = prm( 29, kfor)
      else
        nfunc=0
      endif

      ndd=1
      if(t3b.eq.0.0d0) ndd=0

      njmunu=2

      if(nfunc.eq. 0) then
        if(te.eq.0.0d0 .and. to .eq. 0.0d0) then
          njmunu=1
        endif
      else
        if(b16.eq.0.0d0 .and. b17.eq.0.0d0) then
          if(b14 .eq. 0.0d0 .and. b15 .eq. 0.0d0) then
            njmunu=0
          else
            njmunu=1
          endif
        endif
      endif
      return
      end subroutine lecfor

c______________________________________________________________________________
      subroutine lecmas
c..............................................................................
c   Constants used in the program are declared in this routine.               .
c                                                                             .
c   It concerns:                                                              .
c       - hbar in units of MeV * 10^{-22} s                                   .
c         builtin value:                                                      .
c             hbar   = 6.58211928                                             .
c       - m_n and m_p in units of MeV * c^{-2}                                .
c         builtin values:                                                     .
c             m_n    = 939.565379                                             .
c             m_p    = 938.272046                                             .
c       - speed of light in units of fm/(10^{-22} s)                          .
c         builtin value:                                                      .
c             c      = 29.9792458                                             .
c       - hbar multiplied by C in units of MeV * fm                           .
c         builtin value:                                                      .
c             hbar*c = 197.326 9718                                           .
c                                                                             .
c   Three treatments are used, based on the value of nmass:                   .
c                                                                             .
c   nmass = 0 : Use the builtin constants, taken from the 2010 NIST           .
c               recommendation. Nucleon masses are taken to be equal to the   .
c               average of proton and neutron mass.                           .
c         = 1 : Use the builtin constants, taken from the 2010 NIST           .
c               recommendation. Neutron and proton masses are different.      .                      .
c         = 2 : read (hbar c)^2/2mc of protons and neutrons from data, set    .
c               hbar to the builtin value and calculate the corresponding     .
c               nucleon masses.                                               .                                                                            .
c                                                                             .
c  For the rest of program, the variables xm(1-3) contain the total mass of   .
c  neutrons (1), protons (2) and protons plus neutrons(3).                    .
c                                                                             .
c..............................................................................
c Reference for the newer NIST values:                                        .
c..............................................................................
c P.J. Mohr et al., Rev. Mod. Phys. 84, 2012, CODATA recommended values of    . 
c the fundamental physical constants constants:2010                           .
c                                                                             .
c Or more practically on: http://physics.nist.gov/cuu/Constants/index.html    .
c..............................................................................


      implicit real*8 (a-h,o-z)
      character*20 afor

      parameter (hc    = 197.3269718d0)
      parameter (clum  =  29.9792458d0) 
      parameter (hhbar =   6.58211928d0)
      parameter (xmasp = 938.272046d0)
      parameter (xmasn = 939.565379d0)

c      ........................................................ Old parameters
c      parameter (hhbar=6.58218d0,xxmn =1.044673d0)
c      parameter (hc= 197.327053d0,xmasn=939.565360d0,xmasp=938.272029)
c      parameter (clum=29.9792458d0)

      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     3              ,hbar,hbm(2),xm(3),afor
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)

c..............................................................................
      select case (nmass)

      case (0)
        hbar   = hhbar                    
        hbm(1) = 2.0d0*hc*hc/(xmasn + xmasp)
        hbm(2) = 2.0d0*hc*hc/(xmasn + xmasp)

        xmasrd = 0.5d0*(xmasn+xmasp)
        xm(1)  = xmasrd
        xm(2)  = xmasrd

      case (1)
        hbarc  = hc                 
        hbar   = hc/clum            
        xm(1)  = xmasn              
        xm(2)  = xmasp              
        hbm(1) = hc*hc / xmasn
        hbm(2) = hc*hc / xmasp

      case (2)
        hbm(1) = hbm(1) * 2.0d0       
        hbm(2) = hbm(2) * 2.0d0
        hbar   = hhbar
        xm(1)  = hc*hc / hbm(1)
        xm(2)  = hc*hc / hbm(2)

      case default
        print '(/," nmass = ",i8,/)',nmass
        call stp(' lecmas : nmass!')

      end select

c     .................................. sum of nucleon masses for this nucleus
      xm(3) = npn*xm(1) + npp*xm(2)

      return
      end subroutine lecmas
c______________________________________________________________________________
      blockdata case

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)

      common /skm  / prm (30,34)

      data prm
c...................................................................... SkM*
     1 / -2645.0d0,    410.0d0,   -135.0d0,   15595.0d0,       0.0d0,
     1       0.090d0,    0.000d0,    0.000d0,   0.000d0,       0.0d0,
     1       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     1     130.0d0,    130.0d0,      0.0d0,       1.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c..................................................................... SkM 
     1 -2645.000000d0, 385.000000d0, -120.000000d0, 15595.0d0, 0.0d0, 
     1 0.090000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 130.000000d0, 130.000000d0, 0.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.730000d0, 20.730000d0, 0.000000d0, 
     1 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0d0,  
c...................................................................... SIII  
     2   -1128.75d0,   395.0d0,    -95.0d0,   14000.0d0,       0.0d0,
     2       0.450d0,    0.000d0,    0.000d0,   1.000d0,       0.0d0,
     2       1.0d0,      1.0d0,      0.0d0,       1.0d0,       0.0d0,
     2     120.0d0,    120.0d0,      0.0d0,       1.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,     
c...................................................................... Ska   
     3   -1602.78d0,   570.88d0,   -67.7d0,    8000.0d0,       0.0d0,
     3      -0.020d0,    0.000d0,    0.000d0,    -0.286d0,     0.0d0,
     3       1.0d0,      3.0d0,      0.0d0,       1.0d0,       0.0d0,
     3     125.0d0,    125.0d0,      0.0d0,       1.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c...................................................................... SGII  
     4   -2645.00d0,   340.0d0,    -41.9d0,   15595.0d0,       0.0d0,
     4       0.090d0,   -0.0588d0,   1.425d0,     0.0604d0,    0.0d0,
     4       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     4     105.0d0,    105.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c...................................................................... Sly3         
     1 -1066.975870d0, 245.430973d0, -245.314386d0, 16026.086d0, 0.0d0, 
     1 0.525497d0, 0.603399d0, -0.500115d0, 0.366056d0, 0.000000d0, 
     1 1.000000d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 97.977319d0, 97.977319d0, 0.000000d0, 0.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 0.000000d0, 
     1 -33.847261d0, 61.343170d0, 0.000000d0, 0.000000d0, 0d0, 
c...................................................................... Sly4  
     7   -2488.913d0,  486.818d0, -546.395d0, 13777.0d0,       0.0d0,
     7       0.834d0,   -0.344d0,   -1.0d0,       1.354d0,       0.0d0,
     7       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     7     123.0d0,    123.0d0,      0.0d0,       1.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c...................................................................... Sly5  
     8   -2484.880d0,  483.130d0, -549.400d0, 13763.0d0,       0.0d0,
     8       0.778d0,   -0.328d0,   -1.0d0,       1.267d0,     0.0d0,
     8       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     8     126.0d0,    126.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c...................................................................... Sly6  
     9   -2479.500d0,  462.180d0, -448.610d0, 13673.0d0,       0.0d0,
     9       0.825d0,   -0.465d0,   -1.0d0,       1.355d0,     0.0d0,
     9       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     9     122.0d0,    122.0d0,      1.0d0,       1.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c...................................................................... Sly7 
     A   -2482.410d0,  457.970d0, -419.850d0, 13677.0d0,       0.0d0,
     A       0.846d0,   -0.511d0,   -1.0d0,       1.391d0,     0.0d0,
     A       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     A     126.0d0,    126.0d0,      1.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c...................................................................... SkP  
     B   -2931.70d0,   320.62d0,  -337.41d0,  18708.97d0,      0.0d0,
     B       0.29215d0,  0.65318d0, -0.53732d0,   0.18103d0,   0.0d0,
     B       1.0d0,      6.0d0,      0.0d0,       1.0d0,       0.0d0,
     B     100.0d0,    100.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
     1       0.0d0,      0.0d0,      0.0d0,       0.0d0,       0.0d0,
c..................................................................... ski4 
     1 -1855.827000d0, 473.829000d0, 1006.855000d0, 9703.607d0, 0.0d0, 
     1 0.405082d0, -2.889148d0, -1.325150d0, 1.145203d0, 0.000000d0, 
     1 0.250000d0, 1.000000d0, 1.000000d0, 1.000000d0, 0.000000d0, 
     1 366.194000d0, -360.702000d0,-1.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.752500d0, 20.752500d0, 0.000000d0, 
     1 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0d0,    
c..................................................................... T22 
     1 -2484.397160d0, 484.494607d0, -471.453764d0, 13786.9739d0, 0.0d0, 
     1 0.730120d0, -0.442635d0, -0.944655d0, 1.188194d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 123.224775d0, 123.224775d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 118.684490d0, -72.504183d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 -17.318240d0, 71.696378d0, -17.317615d0, 71.695752d0, 0d0,  
c..................................................................... T24     
     1 -2482.931250d0, 484.346259d0, -433.184957d0, 13768.5605d0, 0.0d0, 
     1 0.729639d0, -0.503889d0, -0.921044d0, 1.190192d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 139.272483d0, 139.272483d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 11.245703d0, 19.739460d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 -11.619436d0, 116.814841d0, -11.619436d0, -3.185159d0, 0d0,        
c..................................................................... T26     
     1 -2476.673460d0, 484.489962d0, -482.591486d0, 13699.0388d0, 0.0d0, 
     1 0.767612d0, -0.434554d0, -0.962725d0, 1.254753d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 156.145858d0, 156.145858d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 -69.885228d0, 120.698480d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 -19.054971d0, 168.531108d0, -19.054970d0, -71.468891d0, 0d0,   
c..................................................................... T42     
     1 -2492.153070d0, 494.635402d0, -251.272272d0, 13869.0593d0, 0.0d0, 
     1 0.690625d0, -0.785802d0, -0.630399d0, 1.121129d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 145.089434d0, 145.089434d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 243.562380d0, -97.619156d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 65.271292d0, 7.943075d0, -54.728709d0, 127.943076d0, 0d0,      
c..................................................................... T44     
     1 -2485.670090d0, 494.477066d0, -337.960664d0, 13794.7464d0, 0.0d0, 
     1 0.721557d0, -0.661848d0, -0.803184d0, 1.175908d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 161.367047d0, 161.367047d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 173.661370d0, 7.173828d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 52.186799d0, 62.432831d0, -67.813199d0, 62.432828d0, 0d0,   
c..................................................................... T46     
     1 -2484.405150d0, 495.224612d0, -356.434894d0, 13769.0702d0, 0.0d0, 
     1 0.735176d0, -0.639443d0, -0.833399d0, 1.201318d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 176.278588d0, 176.278588d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 83.204427d0, 104.872530d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 49.471144d0, 111.874464d0, -70.528859d0, -8.125539d0, 0d0, 
c..................................................................... T62     
     1 -2495.047960d0, 499.980676d0, -197.374135d0, 13901.2355d0, 0.0d0, 
     1 0.690739d0, -0.868510d0, -0.431559d0, 1.117413d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 162.687937d0, 162.687937d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 418.829530d0, -104.641430d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 122.179458d0, -43.698389d0, -117.820537d0, 196.301610d0, 0d0,     
c..................................................................... T64     
     1 -2487.323230d0, 501.095502d0, -284.539013d0, 13818.0333d0, 0.0d0, 
     1 0.705320d0, -0.746420d0, -0.694782d0, 1.148322d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 180.134825d0, 180.134825d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 348.929960d0, -0.196973d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 109.225149d0, 10.922581d0, -130.774870d0, 130.922600d0, 0d0,
c..................................................................... T65     
     1 -2489.412570d0, 497.528219d0, -194.991823d0, 13841.0429d0, 0.0d0, 
     1 0.699857d0, -0.875605d0, -0.446927d0, 1.137559d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 1.000000d0, 1.000000d0, 0.000000d0, 
     1 183.698044d0, 183.698044d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 274.402910d0, 39.898898d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 122.136824d0, 27.939002d0, -117.863178d0, 87.939005d0, 0d0,       
c..................................................................... T66     
     1 -2485.362700d0, 500.799266d0, -228.479017d0, 13794.5593d0, 0.0d0, 
     1 0.715164d0, -0.832653d0, -0.566420d0, 1.165944d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 195.348576d0, 195.348576d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 236.170150d0, 90.314490d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 117.568262d0, 54.695870d0, -122.431740d0, 54.695873d0, 0d0,       
c..................................................................... SkX               
     1 -1445.300000d0, 246.900000d0, -131.800000d0, 12103.90d0, 0.0d0, 
     1 0.340000d0, 0.580000d0, 0.127000d0, 0.030000d0, 0.000000d0, 
     1 0.500000d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 148.600000d0, 0.000000d0, -1.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 1.000000d0, 
     1 -15.807925d0, 47.337500d0, 0.000000d0, 0.000000d0, 0d0,  
c..................................................................... ZAL1    
     1 -2488.913000d0, 486.818000d0, -546.395000d0, 13777.0d0, 0.0d0, 
     1 0.834000d0, -0.344000d0, -1.000000d0, 1.354000d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 80.000000d0, 80.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 -47.366201d0, 129.151625d0, 0.000000d0, 0.000000d0, 0d0,
c................................................................... SLy5colo   
     1 -2484.880000d0, 483.130000d0, -549.400000d0, 13763.0d0, 0.0d0, 
     1 0.778000d0, -0.328000d0, -1.000000d0, 1.267000d0, 0.000000d0, 
     1 0.166667d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 126.000000d0, 126.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 296.000000d0, -136.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 
     1 -8.866670d0, 21.066250d0, -60.000000d0, 162.000000d0, 0d0,  
c.................................................................... ski3     
     1 -1762.880000d0, 561.608000d0, -227.090000d0, 8106.200d0, 0.0d0, 
     1 0.308300d0, -1.172200d0, -1.090700d0, 1.292600d0, 0.000000d0, 
     1 0.250000d0, 1.000000d0, 1.000000d0, 1.000000d0, 0.000000d0, 
     1 188.508000d0, 0.00d0,-1.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.752500d0, 20.752500d0, 0.000000d0, 
     1 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0d0,    
c..................................................................... SVMIN   
     1 -2112.248000d0, 295.781000d0, 142.268000d0, 13988.567d0, 0.0d0, 
     1 0.243886d0, -1.434926d0, -2.625899d0, 0.258070d0, 0.000000d0, 
     1 0.255368d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 -45.936150d0, 157.227150d0, -2.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.721260d0, 20.749821d0, 0.000000d0, 
     1 0.000000d0, 0.000000d0, 0.000000d0, 0.000000d0, 0d0, 
c..................................................................... UNEDF0 
     1 -1883.687810d0, 277.500212d0, 608.430905d0, 13901.9483d0, 0.0d0, 
     1 0.009744d0, -1.777844d0, -1.676990d0, -0.380790d0, 0.000000d0, 
     1 0.321956d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 250.322000d0, -182.520800d0, 0.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 0.000000d0, 
     1 189.210604d0, -41.366337d0, -0.000000d0, 0.000000d0, 0d0, 
c..................................................................... UNEDF1 
     1 -2078.328020d0, 239.400812d0, 1575.11954d0, 14263.6462d0, 0.0d0, 
     1 0.053757d0, -5.077232d0, -1.366506d0, -0.162491d0, 0.000000d0, 
     1 0.270018d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 76.736144d0, 142.633040d0, -1.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 0.000000d0, 
     1 0.000000d0, 0.000000d0, -0.000000d0, 0.000000d0, 0d0,  
c................................................................... UNEDFE1so
     1 -2078.32802d0, 239.400812d0, 1575.11954d0, 14263.6462d0, 0.0d0, 
     1 0.053757d0, -5.077232d0, -1.366506d0, -0.162491d0, 0.000000d0, 
     1 0.270018d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 193.016000d0, -33.832000d0, -1.000000d0, 1.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 0.000000d0,
     1 420.987905d0, -166.964841d0, -0.000000d0, 0.000000d0, 0d0,       
c..................................................................... UNEDF2
     1 -1735.456852d0, 262.814456d0, 1183.31263d0, 12293.2433d0, 0.0d0, 
     1 0.172247d0, -3.766927d0, -1.383498d0, 0.055286d0, 0.000000d0, 
     1 0.351455d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 51.317335d0, 154.600779d0, -1.000000d0, 1.000000d0, 2.000000d0, 
     1 -240.140598d0, -266.930662d0, 20.735530d0, 20.735530d0, 0.00d0,  
     1 201.621390d0, -121.759788d0, 190.151723d0, 10.046274d0, 0d0,
c.................................................................... slyiii_0.7          
     1 -1122.408110d0, 440.572081d0, -197.528490d0, 11906.2992d0, 0.0d0, 
     1 0.394119d0, 0.068384d0, -0.752728d0, 0.946945d0, 0.000000d0, 
     1 1.000000d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 119.125079d0, 119.125079d0, 0.000000d0, 0.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 1.000000d0, 
     1 -22.351679d0, 79.762571d0, 0.000000d0, 0.000000d0, 0d0, 
c.................................................................... slyiii_0.8          
     1 -1100.271870d0, 359.568376d0, -210.839819d0, 13653.8453d0, 0.0d0, 
     1 0.445280d0, 0.224693d0, -0.615015d0, 0.639947d0, 0.000000d0, 
     1 1.000000d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 110.827690d0, 110.827690d0, 0.000000d0, 0.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 1.000000d0, 
     1 -26.307798d0, 71.301024d0, 0.000000d0, 0.000000d0, 0d0,         
c..................................................................... slyiii_0.9          
     1 -1082.609230d0, 295.998896d0, -240.652615d0, 15003.161d0, 0.0d0, 
     1 0.491775d0, 0.389884d0, -0.579284d0, 0.512106d0, 0.000000d0, 
     1 1.000000d0, 1.000000d0, 0.000000d0, 1.000000d0, 0.000000d0, 
     1 103.515557d0, 103.515557d0, 0.000000d0, 0.000000d0, 2.000000d0, 
     1 0.000000d0, 0.000000d0, 20.735530d0, 20.735530d0, 1.000000d0, 
     1 -31.851421d0, 67.081439d0, 0.000000d0, 0.000000d0, 0d0/     
 
      end
c______________________________________________________________________________
      subroutine evolve

c..............................................................................
c     this subroutine commands the time evolution of the system
c itpri = 0 no printout
c       = 1    printout  itprii = 0 summary of the single particle spectrum
c                               = 1           full single particle spectrum
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,t100=100.0d0)

      common /Lag  / ilag, iAnaLag
      common /clou / rvst(2*mv,3),rvjj(2*mv*9)
      common /conv / econv1,econv2,iconv
      common /den  / rho(2*mv)
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /taudj/ vtau(2*mv),vdiv(2*mv)
      common /wj2  / vjj (2*mv*9)

c..............................................................................
  102 format (/)
  103 format ( ' ', 78('_'))

c............................................................... initialisation
      print 102

      itpri  = 1
      itprii = 1

      xmu = float(nxmu)/t100
      ymu = one-xmu

      call inisp

c     ..... ordering of the level indices according to single-particle energies
      call class

c     .... Gram-Schmidt orthonomalization of the single-particle wave functions
      call ortho

c     ...... solve pairing equations according to pairing model and interaction
      call gap (itpri)

c     .................................. calculate local densities and currents
      call densit

c     .................................... calculation of the initial potential
c                                          initialisation of the coulomb solver
c                                          printing the initial energies
      do i=1,mv
        wcd(i) = zero
      enddo
      do iwave=1,nwave
        esp3(iwave) = zero
        esp2(iwave) = zero
      enddo

c     .............................. calculate mean fields from local densities
      call newpot

      if (nitert.eq.0) then

            call fsum(itprii)
            print 103
c     ..................................................... diagnostic printing
      else
            call figaro (itprii,0,0)
            !print 103
      endif
      print 102

c     ....................................
      if (nitert.eq.0) return

c     ........................... save local densities and currents for damping
      do i=1,2*mv
        rvst(i,1) = rho (i)
        rvst(i,2) = vtau(i)
        rvst(i,3) = vdiv(i)
      enddo

c     ................................. diagonalize single-particle Hamiltonian
      if (ndiag.ne.0) call reord

c     .......................................................... time evolution
      npack = nitert/nprint
      do ipack=1,npack
        itprii = ipack/npack
        do ielem=1,nprint
          itpri = ielem/nprint
          itert = itert + 1

c         ............................... gradient step of s. p. wave functions
          call pro

c         ................... reorder level indices according to s. p. energies
          call class

c         .......... Gram-Schmidt orthonomalization of the s. p. wave functions
          call ortho

c         .. solve pairing equations according to pairing model and interaction
          call gap (itpri)

c         .............................. calculate local densities and currents
          call densit

c         ... of damping the local densities entering the mean-field potentials
c                            MB: this might be moved to a subroutine of its own
          do i=1,2*mv
            rho (i) = xmu* rho(i) + ymu*rvst(i,1)
            vtau(i) = xmu*vtau(i) + ymu*rvst(i,2)
            vdiv(i) = xmu*vdiv(i) + ymu*rvst(i,3)
            rvst(i,1) = rho (i)
            rvst(i,2) = vtau(i)
            rvst(i,3) = vdiv(i)
          enddo

c         ................................................ calculate potentials
          call newpot

c         ............................ convergence test and diagnostic printing
          iprint = 1
          iterm  = 0
          call conver (iprint,iterm)
          if (iterm.gt.0) then
            itprii = 1
            call densit
            call fsum(itprii)
            go to 12
          endif
          if (itpri.eq.1) then
            if (itprii.eq.1) then
                  !Recalculate the final densities without smoothing!
                  call densit
                  ! print *, 'Recalculated densities 2'
                  call fsum(itprii)
            else
                  call figaro (itprii,0,0)
            endif
            print 102
            if (iconv.eq.-1) go to 12
            if (itprii.eq.1) go to 12
          endif
        enddo
      enddo

   12 continue

      return
      end subroutine evolve
      
c______________________________________________________________________________
      subroutine ReanalyseLag(itprii)

c..............................................................................
c     Reanalyse the wavefunction with lagrangian derivatives                  .
c..............................................................................
      implicit real*8 (a-h,o-z)

      common /Lag  / ilag, iAnaLag

c..............................................................................
    1 format (' Energies with Lagrange derivatives')

c..............................................................................
      !Signalling to all routines that lagrange derivatives should be used.
      iLag = 1

      !Calculating the coefficients for the lagrangian derivatives
      call iniderlag

      !Orthogonalising the wavefunctions for safety
      call ortho

      !Recalculate the densities (some are calculated using derivatives)
      call densit
      
      !Recalculate potentials. This is not needed, but the calculation of 
      !the laplacian of rho is hidden in here somewhere. In order to not have
      !more unpleasant surprises involving derivatives, just calculate all 
      !potentials again.
      call newpot
      
      print 1
      call figaro (itprii,iLag,0)

      return
      end subroutine ReanalyseLag
c______________________________________________________________________________
      subroutine inisp

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (mpx=mx+mc,mpy=my+mc,mpz=mz+mc,mpv=mpx*mpy*mpz)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (tt8=8.0d0,t10=10.0d0,t12=12.0d0,t16=16.0d0)
      parameter (tp2=0.2d0,tp5=0.5d0,t1p75=1.75d0)
      parameter (ee2= 1.43996446d0) 
      parameter (y44rnm=0.1057855d0,epsm3=1.0d-3)
      parameter (t1000=1000.d0)

      common /conv / econv1,econv2,iconv
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),
     1               imtd,imtg,icutq
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /enert/ c14,c15,t14,t15
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     3              ,hbar,hbm(2),xm(3),afor
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /kfmom/ tfac,tfac1,tfac2,s3
      common /kfpro/ tn
      common /mud  / xi (mx,my,mz),yi (mx,my,mz),zi (mx,my,mz)
     1              ,xii(mx,my,mz),yii(mx,my,mz),zii(mx,my,mz)
      common /mudd / xu(mpx+mpy+mpz)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)

c.................................................... initialize some constants
c                  e2 = square of the charge of the electron = 1.43998d0 MeV*Fm
      e2   =  ee2
      pi    = tt4*atan2(one,one)
      dv    = tt8*dx*dx*dx
      nwave = nwaven + nwavep

c     ...................................................... initialize /mud  /
      x =-tp5*dx
      m = max0(mpx,mpy,mpz)
      do i=1,m+1
        x     = x + dx
        xu(i) = x
      enddo

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

c     ..................................................... initialize  /spwf /
      do iw = 1,nwaven
        kiso(iw) = 1
      enddo
      do iw=1,nwavep
        kiso(iw+nwaven) = 2
      enddo

c     ................................... reset non-convergence indices /conv /
c       iconv is set to -1 if
c         a) the total energy is positive
c         b) the total energy oscillates in a divergent way
      iconv  = 0
      econv1 = zero
      econv2 = zero

c     ................ initialize coefficients of the skyrme functional /ener /
      b1    = t0*(one+x0/two)/two
      b2    =-t0*( x0+tp5   )/two
      b3    = (t1*(one+x1/two)+t2*(one+x2/two))/tt4
      b4    =-(t1*( x1+tp5   )-t2* (x2+tp5   ))/tt4
c             .............................................ATTENTION
c             B5 & B6 as used in the code are not equal to the ones
c             discussed in all of the articles; they carry a 
c             relative minus sign in this code!
c             ...................................................... 
      b5    = (tt3*t1*(one+x1/two)-t2*(one+x2/two))/t16
      b6    =-(tt3*t1*( x1+tp5  ) +t2*( x2+tp5   ))/t16
      b7a   = t3a*(one +x3a/two)/t12
      b8a   =-t3a*( x3a+tp5    )/t12
      b7b   = t3b*(one +x3b/two)/t12
      b8b   =-t3b*( x3b+tp5    )/t12
      byt3a = yt3a
      byt3b = yt3b
c                               njmunu overrides set values for te and to.
c                               having finite values for te and to override
c                               njmunu might be more intelligent, but is
c                               difficult to envision for the contribution
c                               from central terms.
      if (njmunu.lt.2) then
        te  = zero
        to  = zero
      endif

      select case (nfunc)
      case (0)

        select case (njmunu) 

        case (1)
          c14 =-(t1*x1+t2*x2)/tt8
          c15 = (t1   -t2   )/tt8
          t14 = zero
          t15 = zero
          b14 = c14
          b15 = c15
          b16 = zero
          b17 = zero

        case (2)
          c14 =-(t1*x1+t2*x2)/tt8
          c15 = (t1   -t2   )/tt8
          t14 =  (te+to)/tt4
          t15 = -(te-to)/tt4
          b14 = c14+t14
          b15 = c15+t15
          b16 = -tt3*(te+to)/tt8
          b17 =  tt3*(te-to)/tt8

        case default
          print '(/," njmunu = ",i8,/)',njmunu
          call stp(' setfor : njmunu!')
        end select

      case (1)
        if (njmunu.eq.0) then
        te  = zero
        to  = zero
        c14 = zero
        c15 = zero
        t14 = zero
        t15 = zero
        b14 = zero
        b15 = zero
        b16 = zero
        b17 = zero
        endif
        if (njmunu.eq.1) then
          b16 = zero
          b17 = zero
        endif

      case default
        print '(/," nfunc = ",i8)',nfunc
        call stp(' inisp : nfunc!')
      end select

c     ..........................b9/b9q fall in principle into the same category
c                     but remain treated in the tradional manner for the moment
      b9    =-wso /two
      b9q   =-wsoq/two

c     ................................................................. /kfpro/
      tn = - dt/hbar

c     ................................................................. /kfcl /
      coexv = 0.0d0
      if (ncoex.eq.0) coexv = -(tt3/pi)**(one/tt3)*e2

      epscl = epsm3*epsm3/(dx**3)/ ( nnx*nny*nnz)
      e2eff = tt4*pi*e2

c     ................................................................. /kfmom/
      tfac  = sqrt(1/(pi)) * 3.d0/16.d0 * dv
      tfac1 = tfac *sqrt(t10)
      tfac2 = tfac1*sqrt(t1p75)
      s3    = sqrt(tt3)

c     ................................................................. /cst  /
      do i=1,mv
        dd = one
        if(icutq.eq. 2) then
          d = (sqrt(xii(i,1,1)+yii(i,1,1)+zii(i,1,1))-rcut)/acut
          if (d.gt.zero) dd = exp(-d)/(one+exp(-d))
          if (d.le.zero) dd =    one /(one+exp(+d))
        endif
        cutof2(i) = dd
      enddo

c     ..................... initialisation of the fermi levels and of the gaps
c                                                                       /pair /
      do iw=1,nwave
        if (abs(delta(iw)).le.1.d-3) delta(iw)=0.01d0
      enddo

      if ((npair.eq.3).or.(npair.eq.5)) then
        if (xlamb(1).eq.zero) xlamb(1) = tp2
        if (xlamb(2).eq.zero) xlamb(2) = tp2
      else
        xlamb(1)=zero
        xlamb(2)=zero
      endif

      if (npair.eq.2) then
        do it=1,2
          n1=1     +(it-1)*nwaven
          n2=nwaven+(it-1)*nwavep
          do i=n1,n2
            delta(i)=delmax(it)
          enddo
        enddo
      endif

      return
      end subroutine inisp
c______________________________________________________________________________
      subroutine cpling

c..............................................................................
c     diagnostic printing of coupling constants                               .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (tt8=8.0d0,t16=16.0d0)
      parameter (tp2=0.2d0,tp5=0.5d0)
      parameter (tt5=5.0d0,tt9=9.0d0)
      parameter (t24=24.0d0,t32=32.0d0,t48=48.0d0,t64=64.0d0)

      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /enert/ c14,c15,t14,t15
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     3              ,hbar,hbm(2),xm(3),afor
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)

c..............................................................................
  101 format (/,' ',78('_'))
  103 format (/,' Skyrme interaction: ',a20)
  105 format (/,' ncm2    =',i5,/,
     1          ' ndd     =',i5,/,
     1          ' njmunu  =',i5,/,
     2          ' nfunc   =',i5,/,
     9          ' ncoex   =',i5,/,
     3          ' nmass   =',i5,/,/,
     4          ' m_n         =',f14.7,/,
     5          ' m_p         =',f14.7,/,
     6          ' <m>         =',f14.7,/,
     7          ' hbar^2/2m_n =',f14.7,/,
     8          ' hbar^2/2m_p =',f14.7)
  110 format (/,' Coupling constants of the Skyrme interaction ',/)
  111 format (  '  t0 =',f13.6,' x0 =',f12.6)
  112 format (  '  t1 =',f13.6,' x1 =',f12.6)
  113 format (  '  t2 =',f13.6,' x2 =',f12.6)
  114 format (  '  t3 =',f13.6,' x3 =',f12.6,'  sigma =',f10.6)
  115 format (  '  t3a=',f13.6,' x3a=',f12.6,'  sigma =',f10.6,/,
     1          '  t3b=',f13.6,' x3b=',f12.6,'  sigmb =',f10.6)
  116 format (  '   w =',f13.6,' wq =',f12.6)
  117 format (  '  te =',f13.6,' to =',f12.6,
     1          '   T =',f12.6,'   U =',f13.6)
  210 format (/,' Coupling constants of the Skyrme functional ',/)
  211 format (  '  b1  =',f13.6,' b2  =',f12.6)
  212 format (  '  b3  =',f13.6,' b4  =',f12.6)
  213 format (  '  b5  =',f13.6,' b6  =',f12.6)
  214 format (  '  b7  =',f13.6,' b8  =',f12.6,'  sigma =',f10.6)
  215 format (  '  b7a =',f13.6,' b8a =',f12.6,'  sigma =',f10.6,/,
     1          '  b7b =',f13.6,' b8b =',f12.6,'  sigmb =',f10.6)
  216 format (  '  b9  =',f13.6,' b9q =',f12.6)
  217 format (  '  c14 =',f13.6,' c15 =',f12.6)
  218 format (  '  t14 =',f13.6,' t15 =',f12.6)
  219 format (  '  b14 =',f13.6,' b15 =',f12.6)
  220 format (  '  b16 =',f13.6,' b17 =',f12.6)
  310 format (/,' Coupling constants of the Skyrme functional ',/,/,
     1          ' C^rho_0[_0_] =',f13.6,' C^rho_1[_0_] =',f13.6,/,
     2          ' C^rho_0[sat] =',f13.6,' C^rho_1[sat] =',f13.6,/,
     3          ' C^tau_0      =',f13.6,' C^tau_1      =',f13.6,/,
     4          ' C^Drho_0     =',f13.6,' C^Drho_1     =',f13.6,/,
     5          ' C^divJ_0     =',f13.6,' C^divJ_1     =',f13.6,/,
     6          ' C^J(1)_0     =',f13.6,' C^J(1)_1     =',f13.6,/,
     7          ' C^J(2)_0     =',f13.6,' C^J(2)_1     =',f13.6)
 320  format (/,' Decomposition of tensor coupling constants ',/,/,
     1          ' alpha(centr) =',f13.6,' beta(centr)  =',f13.6,/,
     1          ' alpha(tens)  =',f13.6,' beta(tens)   =',f13.6,/,
     1          ' alpha(total) =',f13.6,' beta(total)  =',f13.6,/,
     6          ' A^J_0        =',f13.6,' A^J_1        =',f13.6,/,
     6          ' B^J_0        =',f13.6,' B^J_1        =',f13.6,/,
     6          ' C^J_0        =',f13.6,' C^J_1        =',f13.6,/,
     6          ' A^J(1)_0     =',f13.6,' A^J(1)_1     =',f13.6,/,
     7          ' A^J(2)_0     =',f13.6,' A^J(2)_1     =',f13.6,/,
     6          ' B^J(1)_0     =',f13.6,' B^J(1)_1     =',f13.6,/,
     7          ' B^J(2)_0     =',f13.6,' B^J(2)_1     =',f13.6,/,
     1          ' A^T_0        =',f13.6,' A^T_1        =',f13.6,/,
     9          ' B^T_0        =',f13.6,' B^T_1        =',f13.6,/,
     8          ' C^T_0        =',f13.6,' C^T_1        =',f13.6,/,
     7          ' C^F_0        =',f13.6,' C^F_1        =',f13.6,/)
 321  format (/,' Decomposition of tensor coupling constants ',/,/,
     1          ' alpha(total) =',f13.6,' beta(total)  =',f13.6,/,
     2          ' C^J_0        =',f13.6,' C^J_1        =',f13.6,/,
     3          ' C^T_0        =',f13.6,' C^T_1        =',f13.6,/,
     4          ' C^F_0        =',f13.6,' C^F_1        =',f13.6,/)

c............................... coupling constants in isovector representation
      rhosat  = 0.16d0           ! should be self-consistent saturation density
      rhosata = rhosat**byt3a
      if (ndd.eq.1) then
        rhosatb = rhosat**byt3b
      else
        rhosatb = zero
      endif

c     ................................................... rho^2 term at rho = 0
      ccrho00 =  tt3/tt8 * t0
      ccrho01 = -one/tt4 * t0 * (tp5 + x0)

c     ................................................ rho^2 term at rho = 0.16
      ccrhos0 =  tt3/tt8 * t0
     1          +tt3/t48 * t3a               * rhosata
      ccrhos1 = -one/tt4 * t0  * (tp5 + x0 )
     1          -one/t24 * t3a * (tp5 + x3a) * rhosata

      if (ndd.eq.1) then
        ccrhos0 = ccrhos0 +tt3/t48 * t3b               * rhosatb
        ccrhos1 = ccrhos1 -one/t24 * t3b * (tp5 + x3b) * rhosatb
      endif

c     ............................................................ rho tau term
      cctau0  =  tt3/t16 * t1 + one/tt4 * t2*(tt5/tt4 + x2)
      cctau1  = -one/tt8 * ( t1*(tp5 + x1) - t2*(tp5 + x2) )

c     ...................................................... rho Delta rho term
      cclrho0 = -tt9/t64 * t1 + one/t16 * t2 * (tt5/tt4 + x2)
      cclrho1 =  one/t32 * ( tt3*t1*(tp5+x1) + t2*(tp5+x2) )

c     .......................................................... rho div J term
      ccdivj0 = b9 + tp5*b9q
      ccdivj1 =      tp5*b9q

c     .......................................................... J_ij J_ij term
      if (nfunc.eq.0) then
          act0 = -one/tt8 * (t1*(tp5-x1) - t2*(tp5+x2))
          act1 = -one/t16 * (t1 - t2)
          bct0 = -one/tt8 * (te + tt3*to)
          bct1 =  one/tt8 * (te -     to)
          cct0 =  act0 + bct0
          cct1 =  act1 + bct1
      else
        cct0  = -(b14 + tp5*b15)               ! additional sign as the C refer
        cct1  = -       tp5*b15                ! to s*T, the b to J_ij J_ij
      endif

c     .......................................................... J_ij J_ji term
      if (nfunc.eq.0) then
        if (njmunu.gt.1) then
          ccf0  =  tt3/tt8*(te + tt3*to)
          ccf1  = -tt3/tt8*(te -     to)
        else
          ccf0  = zero
          ccf1  = zero
        endif
      else
        ccf0  = -two * (b16 + tp5*b17)     ! additional factor -2 sign as the C
        ccf1  = -two *        tp5*b17      ! refer to s*F, the b to J_ij J_ij
      endif

c     ................................................... C^(1) and C^(2) terms
c            note that ccJ1x is 1/2 the coupling constant of the [J^(1)]^2 term
c            note also that the separation into "central" and "tensor pieces
c            only makes sense when there are anti-symmetrized vertices

      if (nfunc.eq.0) then
        acJ10  = -tp5*act0
        acJ11  = -tp5*act1
        acJ20  = -act0
        acJ21  = -act1
        acJs0  =  two*acJ10
        acJs1  =  two*acJ11
        aalpha =  acJs0 + acJs1
        abeta  =  acJs0 - acJs1
        bcJ10  = -tp5*(bct0 - tp5*ccf0)
        bcJ11  = -tp5*(bct1 - tp5*ccf1)
        bcJ20  = -    (bct0 + tp5*ccf0)
        bcJ21  = -    (bct1 + tp5*ccf1)
        bcJs0  =  two*bcJ10
        bcJs1  =  two*bcJ11
        balpha =  bcJs0 + bcJs1
        bbeta  =  bcJs0 - bcJs1
        ccJ10  = acJ10 + bcJ10
        ccJ11  = acJ11 + bcJ11
        ccJ20  = acJ20 + bcJ20
        ccJ21  = acJ21 + bcJ21
        ccJs0  = acJs0 + bcJs0
        ccJs1  = acJs1 + bcJs1
        calpha = ccJs0 + ccJs1
        cbeta  = ccJs0 - ccJs1
      else
        acJ10  = zero
        acJ11  = zero
        acJ20  = zero
        acJ21  = zero
        acJs0  = zero
        acJs1  = zero
        aalpha = zero
        abeta  = zero
        bcJ10  = zero
        bcJ11  = zero
        bcJ20  = zero
        bcJ21  = zero
        bcJs0  = zero
        bcJs1  = zero
        balpha = zero
        bbeta  = zero
        ccJ10  = -tp5*(cct0 - tp5*ccf0)
        ccJ11  = -tp5*(cct1 - tp5*ccf1)
        ccJ20  = -    (cct0 + tp5*ccf0)
        ccJ21  = -    (cct1 + tp5*ccf1)
        ccJs0  = two*ccJ10
        ccJs1  = two*ccJ11
        calpha = ccJs0 + ccJs1
        cbeta  = ccJs0 - ccJs1
      endif

c..............................................................................
      print 101
      print 103,afor
      print 105,ncm2,ndd,njmunu,nfunc,ncoex,nmass,
     1          xm(1),xm(2),
     2          xm(3)/(1.0d0*(npn+npp)),
     3          hbm(1)/2.0d0,hbm(2)/2.0d0

c     ................................................... Skyrme representation
      print 110
      print 111,t0,x0
      print 112,t1,x1
      print 113,t2,x2
      if (ndd.eq.0) then
        print 114,t3a,x3a,yt3a
      endif
      if (ndd.eq.1) then
        print 115,t3a,x3a,yt3a,t3b,x3b,yt3b
      endif
      print 116,wso,wsoq
      if (nfunc.eq.0) print 117,te,to,tt3*te,tt3*to

c     .................................... Bonche-Flocard-Heenen representation
      print 210
      print 211,b1,b2
      print 212,b3,b4
c             .............................................ATTENTION
c             B5 & B6 as used in the code are not equal to the ones
c             discussed in all of the articles; they carry a 
c             relative minus sign in this code!
c             ...................................................... 
      print 213,-b5,-b6
      if (ndd.eq.0) then
        print 214,b7a,b8a,byt3a
      endif
      if (ndd.eq.1) then
        print 215,b7a,b8a,byt3a,b7b,b8b,byt3b
      endif
      print 216,b9,b9q
      if (nfunc.eq.0) then
        print 217,c14,c15
        print 218,t14,t15
      endif
      print 219,b14,b15
      print 220,b16,b17

c     .................................................. isospin representation
      print 310,ccrho00,ccrho01,ccrhos0,ccrhos1,
     1          cctau0 ,cctau1 ,cclrho0,cclrho1,
     2          ccdivj0,ccdivj1,
     3          ccJ10,ccJ11,ccJ20,ccJ21

c     ............................................................. tensor part
      if (njmunu.ne.0) then
        if (nfunc.eq.0) then
          print 320,aalpha,abeta,balpha,bbeta,calpha,cbeta,
     1              acJs0,acJs1,bcJs0,bcJs1,ccJs0,ccJs1,
     2              acJ10,acJ11,acJ20,acJ21,
     3              bcJ10,bcJ11,bcJ20,bcJ21,
     4              act0,act1,bct0,bct1,cct0,cct1,
     5              ccf0,ccf1
        else
          print 321,calpha,cbeta,ccJs0,ccJs1,cct0,cct1,
     1              ccf0,ccf1
        endif
      else
        print *,' '
      endif

      return
      end subroutine cpling

c______________________________________________________________________________

c______________________________________________________________________________
      subroutine class

c..............................................................................
c     ordering of the level indices according to single-particle energies     .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)

c..............................................................................

c..............................................................................
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
            endif
          enddo
        enddo
        nof = nof + nvb
      enddo
      enddo

      return
      end subroutine class

c______________________________________________________________________________
      subroutine ortho

c..............................................................................
c     1) Gram-Schmidt orthonomalization                                       .
c     2) normalisation of the wave functions to 2                             .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (one=1.0d0,two=2.0d0)

      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /stor / a(mq,mw)
      common /nxyz / dx,dv

c..............................................................................
      dvs2 = dv/two
      nof = 0
      do it=1,2
      do ib=1,2
        nvb = npar(ib,it)
        do iv=1,nvb
          iwa = nof + iv
          if (iv.gt.1) then
            do jv=1,iv-1
              jwa = nof + jv
              xr  = sdot(mq,a(1,jwa),a(1,iwa))
              xr  =-xr * dvs2
              call saxpy (mq,xr,a(1,jwa),a(1,iwa))
            enddo
          endif
          psint = sdot(mq,a(1,iwa),a(1,iwa))
          fac   = one/sqrt(psint*dvs2)
          call sscal (mq,fac,a(1,iwa))
        enddo
        nof = nof + nvb
      enddo
      enddo

      return
      end subroutine ortho
c______________________________________________________________________________
      subroutine gap (itpri)

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor

c..............................................................................
      if (npair.eq.0) call gaphf
      if (npair.eq.1) call gapbcs
      if (npair.eq.2) call gapbcs
      if (npair.eq.3) call gapln  (itpri)
      if (npair.eq.4) call gapbcs
      if (npair.eq.5) call gapln  (itpri)
      return
      end subroutine gap

c______________________________________________________________________________
      subroutine gaphf

c..............................................................................
c     set various pairing-related quantities to trivial HF values             .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0)
      parameter (tbig=1.0d+30)

      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)

c.......................................................... some initialization
      do i=1,nwave
        eqp(i) = zero
        v2 (i) = zero
      enddo

      do 1 it=1,2
      n1=1     +(it-1)*nwaven
      n2=nwaven+(it-1)*nwavep

c     .............................................................. hf filling
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
      ambda (it) = esp1(k)
      epair (it) = zero
      eproj (it) = zero
      do i=n1,n2
        eqp(i) = abs(esp1(i)-ambda(it))
      enddo

    1 continue

      epair (3) = epair(1)  + epair(2)
      eproj (3) = eproj(1)  + eproj(2)
      disper(3) = disper(1) + disper(2)

      return
      end subroutine gaphf

c______________________________________________________________________________
      subroutine gapbcs

c..............................................................................
c     solve BCS pairing equations for schematic and zero-range pairing forces .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt8=8.0d0,t11=0.0d0)
      parameter (epsm5=1.0d-5)
      parameter (prbcs=1.d-7,prdel=1.d-3)
      parameter (niter=200)

      common /gmatr/ phi2(mv,mw),v(mw,mw,2)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gnp(2),delmax(2),dcut,enpcut(2),xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pfrm / pform(mv,2)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /wave / wf1(mv),wf2(mv),wf3(mv),wf4(mv),psi(mq)
      dimension fi(mw),cc(mw),fs(mw)

c..............................................................................
  100 format(/,' WARNING ',/
     1 ,' BCS iterations did not converge for isopin ',i2,/
     2 ,'  Previous lambda    : ',e20.10, ' Last Lambda: ', e20.10,/
     3 ,'  Largest pairing gap: ',e20.10, ' at orbital : ', i4,/)

c............................ initialize occupations and quasiparticle energies
      do i = 1,nwave
        eqp(i) = zero
        v2 (i) = zero
      enddo


c     ...................................... calculate density-dependent factor
      call gform

c     ................................ loop over isospin to solve BCS equations
c                                            for given single-particle spectrum
      
      do it=1,2
        n1 = 1     +(it-1)*nwaven
        n2 = nwaven+(it-1)*nwavep
        n3 = n1-1

        npnp = npn*(2-it) + npp*(it-1)
        dn = npnp

c       ...................... matrix elements of the schematic seniority force
c                    delta force involves half of matrix element --> factor 1/2
        if (npair.le.3) then
          do i=n1,n2
            i3 = i - n3
            do j=i,n2
              j3 = j - n3
              x           = gnp(it) / (t11+dn) / two
              v(i3,j3,it) = x
              v(j3,i3,it) = x
            enddo
          enddo
        endif

c       .......... calculation of pairing matrix elements of a zero-range force
c                v12(it) = 0.5 * gnp(it) * (1-alpha*rho) * (1-ps12)* delta(r12)
c                                               see Nucl. Phys. A517 (1990) 275
c
c                  integration over one octant only      --> factor 8*dx**3=dv
c                  4 wf's normalized to 2                --> factor 1/4
c                  delta involves half of matrix element --> factor 1/2
c
c                  This could be done differently from a local pair potential
c                  instead of matrix elements

        if (npair.gt.3) then
          do i=n1,n2
            call scopy (mq,a(1,i),wf1)
            do  j=1,mv
              phi2(j,i) = wf1(j)**2+wf2(j)**2+wf3(j)**2+wf4(j)**2
            enddo
          enddo
          do i=n1,n2
            i3 = i - n3
            do j=i,n2
              j3 = j - n3
              x  = zero
              do k=1,mv
                x  = x + phi2(k,i)*pform(k,it)*phi2(k,j)
              enddo
              x           = gnp(it) * x * dv / tt8
              v(i3,j3,it) = x
              v(j3,i3,it) = x
            enddo
          enddo
        endif

c      ........... reinitialization to avoid missing a non-trivial bcs solution
c                  delmax(it) is the largest gap for this isospin from the last
c                  HF iteration. If it was smaller than 10^-5, reinitialize the
c                  pairing gaps. 1 is the typical size of pairing gaps.

        if (delmax(it).lt.epsm5) then
          do i=n1,n2
            delta(i) = one
          enddo
        endif
        
        cutoff = enpcut(it)

c       ......................................................... bcs iteration
        nit = 0
    7   continue
          nit = nit+1
          do i = n1,n2
            if (i.eq.-ntqp) delta(i) = zero
            cc(i) = delta(i)
          enddo

c       ....................................... calculation of the new lambda's
c         fi = squared cut-off factor, see Nucl. Phys. A517 (1990) 275, Eq.(5)
c         mind the misprint in Nucl. Phys. A443 (1985) 39 eq. (14)

          ambdai = ambda(it)
          soer   = zero
          sor    = zero
          xfill  = zero
          do i = n1,n2
            x = esp1(i)
c           .................................................... cutoff factors
            xup   = x - (ambdai + cutoff)
            xdo   = x - (ambdai - cutoff)
            f     = one/(one+exp(xup/dcut))
            f     = f*(one-xcut+xcut/(one+exp(-xdo/dcut)))
            
            fi(i) = f
            x     = esp1(i)
            xml   = ambdai - x

c           ............................................ quasiparticle energies
c                           usually positive, but negative for blocked articles
c                                       (which swaps the roles of u^2 and v^2)
            y = sqrt(xml**2 + f*cc(i)**2)
            if (i.eq.ntqp) y = -y
            eqp(i) = y

c           ............ sum up quantities needed to calculate the Fermi energy
c                        filled orbitals have to be taken out of the sum
            if(i.ne.-ntqp) then
             soer = soer + (one-x/y)
             sor  = sor  +  one/y
            else
              xfill = one
            end if
          enddo


c         ................ calculate fermi energy for given distribution of v^2
c                      N = 2 sum v2
c                        = sum 1 + {(e-lambda)/sqrt[(e-lambda)^2 + f delta^2]}
c                        = soer + lambda * sor
c              => lambda = (N-soer)/sor
c              In case of filling, the filled levels are taken out of the sums
c              soer and sor, but their contribution has to be subtracted from N
c
           ambda(it) = (dn-xfill-soer)/sor

c         ......................................... calculation of the new gaps
          k = 0
          x = zero
          if (npair.ne.2) then
            do  i=n1,n2
              i3 = i - n3
              y  = zero
              do j = n1,n2
                if (j.ne.-ntqp) then
                  y = y + v(i3,j-n3,it)*fi(j)*cc(j)/eqp(j)
                endif
              enddo
              if (i.eq.-ntqp) y = zero    ! set gap to zero for filled orbitals
              delta(i) = y

c             ......... after the loop, the local variable x is the largest gap
c                   for this isospin and k the index of the coresponding level
c                   The common block variable delmax(it) doubles this, which
c                   seems unnecessary
c
              if (y.gt.delmax(it)) delmax(it) = y
              y = abs(cc(i)-y)
              if (y.gt.x) then
                x = y
                k = i
              endif
            enddo
          else
c           ............................. with (npair.eq.2) all gaps are equal
            delmax(it) = delta(1)
          endif


c       ..................................................... check convergence
        if (nit.gt.niter) go to 9
        if ( abs(ambda(it)-ambdai).ge.prbcs ) go to 7
        if ((x.ge.prbcs).and.(delmax(it).ge.prdel)) go to 7
        if ((x.ge.prdel).and.(delmax(it).lt.prdel)) go to 7
    9   if (nit.gt.niter) print 100,it*2-3,ambdai,ambda(it),x,k

c....................... - recalculate everything with the final Fermi energy
c                        - sum up the dispersion of particle number <N^2>-<N>^2
c                          + a decorrelated 2-qp state contributes to the
c                            dispersion as in the BCS vacuum
c                          + in case of filling, the filled particle should not
c                            contribute to the dispersion
        ss2 = zero
        do i = n1,n2
          if(i.ne.-ntqp) then

c           ............................................ single-particle energy
            e = esp1(i)

c           ..................................................... cutoff factor
            xup   = e - (ambda(it) + cutoff)
            xdo   = e - (ambda(it) - cutoff)
            f     = one/(one+exp(xup/dcut))
            f     = f*(one-xcut+xcut/(one+exp(-xdo/dcut)))
            fi(i) = f

c           .............................................. quasiparticle energy
c                                    in case of blocked 2-qp, the quasiparticle
c                                    energy is negative which swaps u^2 and v^2
            xml   = ambda(it)-e
            y     = sqrt(xml**2+f*delta(i)**2)
            if (i.eq.ntqp) y = -y
            eqp(i) = y

c           ............................
            c      = xml/eqp(i)
            cc (i) = c
            v2 (i) = (one+c)/two
            v2 (i) = max(v2(i),zero)
            v22(i) = v2(i)

c           ................. ............ 2 |u v| calculated in an unusual way
            f      = sqrt(fi(i))
            s      = f*delta(i)/eqp(i)
            fs (i) = f*s

c           ....................................... contribution to <N^2>-<N>^2
            s2     = s*s
            ss2    = ss2 + s2
          else
c          ............................... filling approximation for this level
            eqp(i) = zero
            v2 (i) = one/two
            v22(i) = one/two
            fs (i) = zero
            cc (i) = zero
          endif
        enddo

c       .............................................. calculate pairing energy
c                         Eq. (7) in Krieger et al, Nucl. Phys. A517 (1990) 275
        sgs2  = zero
        sgs2c = zero
        do i=n1,n2
          i3 = i-n3
          x  = zero
          if (npair.ne.2) then
            do j=n1,n2
              z = v(i3,j-n3,it)*fs(j)
              x = x+z
            enddo
          else
c           ............................. with (npair.eq.2) all gaps are equal,
c                                                 so we can use the largest one
            x = delmax(it)
          endif
          sgs2  = sgs2  + fs(i)*x
          sgs2c = sgs2c + fs(i)*cc(i)*x                                 ! ?????
        enddo

c       ......................... pairing energy and particle-number dispersion
        epair (it) = -sgs2/two
        disper(it) =  ss2

c       ...................... this is BCS, so the LN correction energy is zero
        eproj (it) =  zero

      enddo
c     ...... total pairing energy, LN correction and particle-number dispersion
      epair (3) = epair (1) + epair (2)
      eproj (3) = eproj (1) + eproj (2)
      disper(3) = disper(1) + disper(2)

      return
      end subroutine gapbcs

c______________________________________________________________________________
      subroutine gapln (itpri)

c..............................................................................
c     solve BCS + LN pairing equations for schematic and zero-range pairing   .
c     forces                                                                  .
c                                                                             .
c     determination of                                                        .
c        1) hfbcs occupation of the orbitals:     v2                          .
c        2) fermi levels:                         ambda                       .
c           LN lambda2                            xlamb                       .
c        3) pairing gaps:                         delta                       .
c        4) bcs pairing energies:                 epair                       .
c        5) bcs quasiparticle energies            eqp                         .
c        6) bcs nucleon dispersions:              disper                      .
c        7) lipkin-nogami level occupations       v22                         .
c        8) lipkin-nogami energy corrections      eproj                       .
c        9) lipkin-nogami radii and q2            rln,qxln,qyln,qzln          .
c..............................................................................

c..............................................................................
c     We start with the old set of occupations, lambda and lambda2 and a new  .
c     set of single-particle wave functions and energies. This routine        .
c     determines occupations, lambda, and lambda2 that correspond to the new  .
c     single-particle wave functions and energies.                            .
c                                                                             .
c     The iteration scheme is as follows                                      .
c                                                                             .
c     update gaps from old u and v and new wave functions
c     11 loop over values of lambda2 --------------------------------\
c     22   loop over sets of gaps gaps for given lambda2 -------\    |
c     33     loop over lambda for given lambda2 and gaps ---\   |    |
c              update lambda until new v2 give N = sum v2   |   |    |
c            go to 33 --------------------------------------/   |    |
c     44   update gaps with new v2 until they are stable        |    |
c          go to 22 --------------------------------------------/    |
c        update lambda2 with new gaps and v2 until it is stable      |
c        go to 11 ---------------------------------------------------/
c     66 calculate observables
c                                                                             .
c     The density entering the density-dependent form factor of the pairing   .
c     energy functional is left untouched during the iterative scheme         .
c..............................................................................

c..............................................................................
c     The gaps are calculated using the wave functions stored in /gmatr/      .
c     A less storage intensive alternative would be to reformulate the method .
c     using a local pair density, which, however, requires to recalculate     .
c     what is stored in /gmatr/ over and over when summing up the pair density.
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt4=4.0d0,tt8=8.0d0)
      parameter (tp1=0.1d0,tp15=0.15d0,tp85=0.85d0,tp2=0.2d0)
      parameter (tp5=0.5d0,tp25=0.25d0,tp75=0.75d0)
      parameter (epsm5=1.0d-5)
      parameter (prbcs=1.d-6,priini=1.d-3,gapmin=1.d-3)
      parameter (niter1=300,niter2=200,niter3=100)

      logical lpair

      common /gmatr/ phi2(mv,mw),v(mw,mw,2)
      common /lnog / rln(3),qxln(3),qyln(3),qzln(3)
      common /mud  / xi(mv),yi(mv),zi(mv),xii(mv),yii(mv),zii(mv)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gnp(2),delmax(2),dcut,enpcut(2),xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pfrm / pform(mv,2)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /wave / wf1(mv),wf2(mv),wf3(mv),wf4(mv),psi(mq)

      dimension fi(mw),ss(mw),cc(mw),fs(mw),prip(mx,my,mz,mw)
      equivalence (prip,phi2)

c..............................................................................
  101 format(' bcs does not converge, it=',i2,'  niter ',i4,/,
     1       '  lambda ',2f14.7)
  102 format(' bcs does not converge, it=',i2,2e16.9,/,
     1       '  gap orbital, err, old, new ',i4,3e13.6)
  104 format(' ln does not catch',i2,5e10.3)
  110 format(' xlamb divergence in ln ',i2,' alpha ',2e15.8,/,
     1       '                        ',2x,' ambda ',2e15.8,/,
     2       '                        ',2x,' delta ',2e15.8)

c..............................................................................

c     ...................................... calculate density-dependent factor
      call gform

c     ................................. loop over isospin to solve LN equations
c                                            for given single-particle spectrum
      do it=1,2
        lpair = .true.
        if (abs(gnp(it)).lt.tp1) lpair = .false.
        n1   = 1     +(it-1)*nwaven
        n2   = nwaven+(it-1)*nwavep
        n3   = n1-1
        npnp = npn*(2-it)+npp*(it-1)
        dn   = npnp

c       ...................... matrix elements of the schematic seniority force
c                          delta involves half of matrix element --> factor 1/2
        if (npair.le.3) then
          do i=n1,n2
            i3 = i - n3
            do j=i,n2
              j3 = j - n3
              x           = gnp(it) / (11+dn) / two
              v(i3,j3,it) = x
              v(j3,i3,it) = x
            enddo
          enddo
        endif

c       ............................. calculate the two-body pair wave function
c               and calculate the pairing matrix elements of a zero-range force
c                v12(it) = 0.5 * gnp(it) * (1-alpha*rho) * (1-ps12)* delta(r12)
c                                               see Nucl. Phys. A517 (1990) 275
c
c            integration over one octant only      --> factor 8*dx**3=dv
c            4 wf's normalized to 2                --> factor 1/4
c            delta involves half of matrix element --> factor 1/2
c
c            phi2 is the two-body wave function entering the pairing matrix
c            element. It is stored for later recalculation of the gaps with
c            new values for the u and v, and also can be used to calculate
c            the density with LN occupation numbers after convergence.

        if (npair.gt.3) then
          do i=n1,n2
            call scopy (mq,a(1,i),wf1)
            do j=1,mv
              phi2(j,i) =   wf1(j)*wf1(j) + wf2(j)*wf2(j)
     1                    + wf3(j)*wf3(j) + wf4(j)*wf4(j)
            enddo
          enddo
          do i=n1,n2
            i3=i-n3
            do j=i,n2
              j3=j-n3
              x  = zero
              do k=1,mv
                x  = x + phi2(k,i)*pform(k,it)*phi2(k,j)
              enddo
              x = gnp(it)*x*dv / tt8
              v(i3,j3,it) = x
              v(j3,i3,it) = x
            enddo
          enddo
        endif

c       .......... reinitialization of the gaps in case of trivial BCS solution
c                                        found in the preceeding HF+BCS+LN step
c                                delmax(it) is the largest gap found previously
        if (delmax(it).lt.epsm5) then
          do i=n1,n2
            delta(i) = one
          enddo
        endif
        cutoff = enpcut(it)
        nfail  = 0
        pritmp = priini

c       ....................................... LN iteration: loop over lambda2
c                          ambda0 : Fermi energy from preceeding HF+BCS+LN step
c                          alph2i : lambda2 from preceding HF+BCS+LN step
c                          dalph2 : change of lambda2 during the LN iteration
c                          delta0 : largest gap from preceeding HF+BCS+LN step
c                       in case of lpair = .false. there is no point to iterate
c
        nit1 = 0
   11   continue
          nit1   = nit1 + 1
          ambda0 = ambda (it)
          alph2i = xlamb (it)
          dalph2 = zero
          delta0 = delmax(it)

c         ............................... BCS iteration: loop over sets of gaps
          nit2 = 0
   22     continue
            nit2 = nit2 + 1
            do i=n1,n2
              if (i.eq.-ntqp) delta(i) = zero
              cc(i) = delta(i)
            enddo

c           .................................................. loop over lambda
c                                             for given lambda2 and set of gaps
c             fi = squared cut-off factor, Nucl. Phys. A517 (1990) 275, Eq. (5)
c             mind the misprint in Nucl. Phys. A443 (1985) 39 Eq. (14)
c           ...................................................................
            nit3 = 0
   33       continue

              nit3 = nit3 + 1
              ambdai = ambda(it)
              soer   = zero
              sor    = zero
              xfill  = zero
              do i=n1,n2
                e     = esp1(i)

c               ................................................ cutoff factors
                xup   = e - (ambdai + cutoff)
                xdo   = e - (ambdai - cutoff)
                f     = one/(one+exp(xup/dcut))
                f     = f*(one-xcut+xcut/(one+exp(-xdo/dcut)))
                fi(i) = f
                x     = e + tt4*alph2i * (v2(i)-tp5)
                xml   = ambdai - x

c               ........................................ quasiparticle energies
c                          usually positive, but negative for blocked particles
c                                       (which swaps the roles of u^2 and v^2)
                y = sqrt(xml*xml + f*cc(i)*cc(i))
                if (i.eq.ntqp) y =-y
                eqp(i) = y

c               ........ sum up quantities needed to calculate the Fermi energy
c                               filled orbitals have to be taken out of the sum
                if(i.ne.-ntqp) then
                  soer = soer + (one-x/y)
                  sor  = sor  +  one/y
                else
                  xfill = one
                end if

                v2f   = (one+xml/y)/two
                v2(i) = tp25*v2f + tp75*v2(i)
                v2(i) = max(v2(i),zero)

                if (i.eq.-ntqp) v2(i) = one/two
              enddo

c             ............ calculate Fermi energy for given distribution of v^2
c                      N = 2 sum v2
c                        = sum 1 + {(e-lambda)/sqrt[(e-lambda)^2 + f delta^2]}
c                        = soer + lambda * sor
c              => lambda = (N-soer)/sor
c              In case of filling, the filled levels are taken out of the sums
c              soer and sor, but their contribution has to be subtracted from N
c
              ambda(it) = (dn-xfill-soer)/sor

c           ................................................. convergence check
            if (abs(ambda(it)-ambdai).lt.pritmp) go to 44
            if (nit3.lt.niter3) go to 33                ! next lambda iteration

c           ................. calculation new gaps and store them for later use
   44       x = zero
            k = 0
            delmax(it) = zero
            do i=n1,n2
              i3 = i-n3
              y  = zero
              do j=n1,n2
                if (j.ne.-ntqp) then
                  y = y + v(i3,j-n3,it)*fi(j)*cc(j)/eqp(j)
                endif
              enddo

              if(i.eq.-ntqp) y = zero
              y        = abs(y)
              delta(i) = y

c             ......... after the loop, the local variable x is the largest
c                       difference between any old and new gaps
c                       and k the index of the corresponding level

              y = fi(i)*y
              if (y.gt.delmax(it)) delmax(it) = y
              y = fi(i) * abs( cc(i) - delta(i) )
              if (y.gt.x) then
                x = y
                k = i
              endif
            enddo
          if (x.   lt.pritmp) go to 55                 ! trivial solution found
          if (nit2.lt.niter2) go to 22                 ! next BCS iteration

c         .............................while in ln mode bcs degenerated into hf
c                           which should not happen for finite pairing strength
   55     if (delmax(it).lt.gapmin.and.lpair) then
            nfail = nfail + 1
            if (nfail.gt.4) then
              print 104,it,ambda(it),delmax(it),alph2i
              call stp (' nfail !')
            endif

c           ................................................. increasing xlambi
            alph2i     = tp2 + nfail * tp1
            do i=n1,n2
              delta(i) = one + nfail * tp1
            enddo
            go to 22                                      ! next BCS iteration
          endif

c         ............ moments of Delta N that enter the calculation of lambda2
c                          for given single-particle spectrum, gaps, and lambda
          ss2   = zero
          ss2c  = zero
          ss4   = zero
          ses2  = zero
          ses2c = zero
          do i=n1,n2
            if (i.ne.-ntqp) then
              e      = esp1(i)

c             ......................... NOTE: the u and v factors are expressed
c                in terms of epsilon', lambda and Eqp throughout the subroutine
c
c                    1  [     eps' - lambda ]
c             v^2 = --- | 1 - ------------- |
c                    2  [           Eqp     ]
c
c                    1  [     eps' - lambda ]
c             u^2 = --- | 1 + ------------- |
c                    2  [           Eqp     ]
c
c             where Eqp = sqrt[ (eps' - lambda)^2 + f^2 Delta^2 ]
c
c             From this follows 
c
c             (lambda - eps') / Eqp = 2v^2 - 1 = v^2 - u^2
c
c             and, using the shorthand x = (eps' - lambda)/Eqp
c
c             (2uv)^2 = 4 u^2 v^2 = (1+x)*(1-x) = 1 - x^2 
c
c                                   (eps' - lambda)^2
c                       = 1 - -------------------------------
c                             (eps' - lambda)^2 + f^2 Delta^2
c
c                                   f^2 Delta^2 
c                       =  -------------------------------
c                          (eps' - lambda)^2 + f^2 Delta^2
c
c             or
c
c             2uv = f Delta / Eqp   
c
c             which is stored as variable "s" in the context of LN corrections.
c
c             .......................................................... cutoff
              xup    = e - (ambda(it) + cutoff)
              xdo    = e - (ambda(it) - cutoff)
              f      = one/(one+exp(xup/dcut))                  
              f      = f*(one-xcut+xcut/(one+exp(-xdo/dcut)))
              fi (i) = f                    ! ATTENTION "f" is here still f^2 !
              xml    = ambda(it) - e - tt4*alph2i * (v2(i)-tp5) ! lambda - eps'
              y      = sqrt(xml**2 + f*delta(i)**2)             ! Eqp
              if (i.eq.ntqp) y = -y         ! 2qp state: exchanges u and v
              eqp(i) = y                    ! Eqp saved for later use
              c      = xml/eqp(i)           ! 2 v^2 - 1 = v^2 - u^2
              cc (i) = c                    ! save c for later use
              v2 (i) = (one+c)/two          ! v^2 = (1 + 2 v^2 - 1)/2 
              v2 (i) = max(v2(i),zero)      ! avoid negative v2 from num. noise
              v22(i) = v2(i)                ! initialize LN occupations
              f      = sqrt(fi(i))          ! ATTENTION: now f is really f
              s      = f*delta(i)/eqp(i)    ! 2 |u v|
              ss (i) = s*s                  ! save s^2 for later use
              fs (i) = f*s                  ! 2 |u v| * f
            else                            ! filling: this orbit should not
              fi (i) = zero                 ! be occupied
              e      = zero
              eqp(i) = zero
              c      = zero
              cc (i) = zero
              v2 (i) = one/two
              fi (i) = zero
              s      = zero
              ss (i) = zero
              fs (i) = zero
            endif
            s2    = s*s                     ! 4 u^2 v^2
            ss2   = ss2   + s2              ! sum |2uv|^2
            ss2c  = ss2c  + s2*c            ! sum |2uv|^2 (v^2 - u^2)
            ss4   = ss4   + s2*s2           ! sum |2uv|^4
            ses2  = ses2  + e*s2            ! sum |2uv|^2 epsilon 
            ses2c = ses2c + e*s2*c          ! sum |2uv|^2 epsilon (v^2 - u^2)
          enddo

c         ................ recalculate lambda2 using the new values for u and v
c                           only the pairing contribution if taken into account
          sgs2   = zero
          sgs2c  = zero
          sgs2c2 = zero
          do i=n1,n2
            i3 = i-n3
            x = zero
            y = zero
            do j=n1,n2
              z = v(i3,j-n3,it)*fs(j)
              x = x+z
              y = y+z*cc(j)
            enddo
            sgs2   = sgs2   + fs(i)*x
            sgs2c  = sgs2c  + fs(i)*cc(i)*x
            sgs2c2 = sgs2c2 + fs(i)*( (2*cc(i)*cc(i)-1)*x + cc(i)*y )
          enddo
c                                       if lipkin-nogami: redefinition of xlamb
c                        remember that there is a factor 1/2 included in v(i,j)
c                                            Z. Phys. A348 (1994) 183, Eq. (18)
          if (lpair) then
            hd2n      = -(two*ses2c+sgs2c2)
            hdn       = ses2+sgs2c
            x         = ss2*hd2n+2*ss2c*hdn
            deno      = ss2*2*(ss2*ss2-3*ss4+2*ss2)-4*ss2c*ss2c
            alph2f    = x/deno
          else
            hd2n      = zero
            hdn       = zero
            x         = zero
            deno      = one
            alph2f    = x/deno
          endif

c         ...................... mixing of old (alph2i) and new (alph2f) lambda
c                                                        to stabilize iteration
          xlamb(it) = tp15*alph2f + tp85*alph2i

c         ..................................
          dalph2    = abs(alph2f-alph2i)
          dambda    = abs(ambda (it)-ambda0)
          ddelta    = abs(delmax(it)-delta0)
          errmax    = max(dalph2,dambda,ddelta)

c         .................... no converged within allowed number of iterations
c                          (stopping the program is a very drastic measure ...)
          if (nit1.gt.niter1) then
            print 110,it,alph2i,alph2f,ambda0,ambda(it),
     1                   delta0,delmax(it)
            call stp (' nit1 !')
          endif

c         .................................................... convergence test
c                                           maximum change of - lambda
c                                                             - lambda2
c                                                             - the largest gap
c                                          has to be smaller than prbcs = 10^-6
          if (errmax.lt.prbcs) go to 66
          if (.not.lpair) go to 66
          pritmp = min(errmax*tp1,priini)
        go to 11                        ! if not, repeat with the new lambda2

c       .......................................................................
c       .                  end of LN loop for given it                        .
c       .......................................................................

c       ...................................................... some observables
   66   epair (it) =-sgs2/two                 ! pairing energy
        disper(it) = ss2                      ! <N^2> - <N>^2
        eproj (it) = xlamb(it)*ss2            ! LN correction to pairing energy

c       ....... calculate effective LN occupations v22 for one-body observables
c       This is Eq. (3.3) of Bennour et al, PRC 40 (1989) 2834, although in a
c       massaged form in order to save a few microseconds of CPU time and to 
c       waste the time of everybody trying to understand this subroutine.
c
c  Using the definitions
c 
c       s = 2 |u v|
c       c = 2 v^2 - 1 = v^2 - u^2
c
c the Bennour expression is             
c
c v22 = v2 
c
c   1             s^2 * { c * [sum s2] - [sum s2c]} * [sum s2] 
c + - * -----------------------------------------------------------------------
c   4 [sum s2 c2][sum s2] - [sum s2c]^2 + 1/4[sum s2]^3 - 1/2[sum s4][sum s2]
c
c where the sums run over ALL single-particle states, as stressed in the paper.
c
c Reducing the summations to half the single-particle states gives
c a global factor 4 in the nominator, but different factors in the denominator
c
c                       4   s^2 * { c * [sum s^2] - [sum s^2c]} * [sum s^2] 
c           v22 = v^2 + - * -----------------------------------------------
c                       4                      deno
c       where 
c
c       deno =  4 [sum s2][sum s2c2]
c             - 4 [sum s2c]^2 
c             + 2 [sum s2]^3 
c             - 2 [sum s4][sum s2]
c
c       Using that
c
c           [sum s2c2]  = sum 4 u^2 v^2 (2 v^2 - 1)^2
c                       = sum 4 u^2 v^2 (4 v^4 - 4 v^2 + 1)
c                       = sum 4 u^2 v^2 [4 v^2 (1 - u^2) - 4 v^2 + 1]
c                       = sum 4 u^2 v^2 [1 - 4 u^2 v^2)]
c                       = [sum s2] - [sum s4]
c
c       'deno' can be rewritten as
c
c       deno  =   4 [sum s2][sum s2] - 4 [sum s2][sum s4]
c               - 4 [sum s2c]^2 
c               + 2 [sum s2]^3 
c               - 2 [sum s4][sum s2]
c
c       Taking out [sum s2] where possible and combining the [sum s2][sum s4] 
c       terms
c 
c       deno  =   2 * [sum s2] * { 2 [sum s2] - 3 [sum s4] + [sum s2]^2 }
c               - 4 * [sum s2c]^2 
c             = [sum s2] * 2 * {\A0[sum s2]^2 - 3 [sum s4] + 2 [sum s2] }
c               - 4 * [sum s2c]^2 
c
c       which is deno as defined above: ss2*2*(ss2*ss2-3*ss4+2*ss2)-4*ss2c*ss2c
c       .......................................................................
c       sv2  = zero
c       sv22 = zero
        do i=n1,n2
c         s      = fs (i)                                   ! why 2 |u v| f ???
c         v22(i) = v22(i) + s*s*ss2*(cc(i)*ss2-ss2c)/deno
          v22(i) = v22(i) + ss(i)*ss2*(cc(i)*ss2-ss2c)/deno
          if (i.eq.-ntqp) v22(i) = one/two
c         sv2  = sv2  + 2.0d0*v2 (i)
c         sv22 = sv22 + 2.0d0*v22(i)
        enddo
c       print '(/," it = ",i2," sum v2 = ",e15.8,
c    1            " sum v22 = ",e15.8,/)',it,sv2,sv22

c       ......... recalculate radius and quadrupole moments with LN occupations
c                   these values could be used for the constraints, but are not
        if (itpri.eq.1) then
          rln (it) = zero
          qxln(it) = zero
          qyln(it) = zero
          qzln(it) = zero
          do  i=n1,n2
            xx2 = v22(i) * sdot(mv,phi2(1,i),xii)
            yy2 = v22(i) * sdot(mv,phi2(1,i),yii)
            zz2 = v22(i) * sdot(mv,phi2(1,i),zii)
            rln (it) = rln (it) + (xx2+yy2+zz2)
            qxln(it) = qxln(it) + (two*xx2-yy2-zz2)
            qyln(it) = qyln(it) + (two*yy2-zz2-xx2)
            qzln(it) = qzln(it) + (two*zz2-xx2-yy2)
          enddo
          rln (it) = rln (it) * dv
          qxln(it) = qxln(it) * dv
          qyln(it) = qyln(it) * dv
          qzln(it) = qzln(it) * dv
        endif
      enddo                                               ! end of isospin loop

c     .............................................. calculate some observables
      epair(3)  = epair(1)  + epair(2)
      eproj(3)  = eproj(1)  + eproj(2)
      disper(3) = disper(1) + disper(2)

      if (itpri.eq.1) then
        rln(3)  = sqrt((rln(1)+rln(2))/(npn+npp))
        rln(1)  = sqrt(rln(1)/npn)
        rln(2)  = sqrt(rln(2)/npp)
        qxln(3) = qxln(1)+qxln(2)
        qyln(3) = qyln(1)+qyln(2)
        qzln(3) = qzln(1)+qzln(2)
      endif

      return
      end subroutine gapln

c______________________________________________________________________________
      subroutine gform

c..............................................................................
c     density dependent "form factor" of the pairing force                    .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (one=1.0d0)

      common /den  / rhon(mv),rhop(mv)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pairf/ gnp(2),delmax(2),dcut,enpcut(2),xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pfrm / pform(mv,2)

c..............................................................................
      do  j=1,mv
        pform(j,1) = one - alpha*(rhon(j)+rhop(j))
        pform(j,2) = one - alpha*(rhon(j)+rhop(j))
      enddo

      return
      end subroutine gform
c______________________________________________________________________________
      subroutine densit

c..............................................................................
c     calculate                                                               .
c       - density rho                                                         .
c       - kinetic density tau                                                 .
c       - divergence of spin-orbit current nabla.J                            .
c       - spin-tensor density J_ij                 (njmunu > 0 only)          .
c       - matrix elements of the nabla operator    (ncm    > 0 only)          .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (epsm5=1.0d-5)

      common /Lag  / ilag, iAnaLag iAnaLag
      common /den  / rho(mv,2)
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv)
     1              ,p1(mv),p2(mv),p3(mv),p4(mv)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /wj2  / vjxx(mv,2),vjyx(mv,2),vjzx(mv,2)
     1              ,vjxy(mv,2),vjyy(mv,2),vjzy(mv,2)
     2              ,vjxz(mv,2),vjyz(mv,2),vjzz(mv,2)
      dimension ta(mv)

c..............................................................................
      do i=1,mv
        rho (i,1) = zero
        rho (i,2) = zero
        vtau(i,1) = zero
        vtau(i,2) = zero
        vdiv(i,1) = zero
        vdiv(i,2) = zero
      enddo

      do iw=1,nwave
        it  = kiso (iw)
        iz  = kparz(iw)
        oc  = v2   (iw)
        ocd = oc + oc
        call scopy (mq,a(1,iw),w1)
        do i=1,mv
          ta(i) = w1(i)**2 + w2(i)**2 + w3(i)**2 + w4(i)**2
        enddo
        call saxpy (mv,oc,ta(1),rho(1,it))
        if(ilag.eq.0) then
            call deriv (iz)
        else
            call derivlag(iz)
        endif
        call scopy (3*mq,wx1(1),da(1,iw))
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

      do it=1,2
        do i=1,mv
          rho(i,it) = max(rho(i,it),zero)
        enddo
      enddo

c     .................................................. full spin-orbit tensor
      if (njmunu.ge.1) then
        do i=1,9*mv*2
          vjxx(i,1) = zero
        enddo
        do iwave=1,nwave
          it = kiso(iwave)
          call scopy (mq,a(1,iwave),w1)
          call scopy (3*mq,da(1,iwave),wx1)
          v2i=v2(iwave)
          do i=1,mv
            vjxx(i,it) = vjxx(i,it)+v2i*( w1(i)*wx4(i)-w2(i)*wx3(i)
     1                                   +w3(i)*wx2(i)-w4(i)*wx1(i))
            vjyx(i,it) = vjyx(i,it)+v2i*( w1(i)*wy4(i)-w2(i)*wy3(i)
     1                                   +w3(i)*wy2(i)-w4(i)*wy1(i))
            vjzx(i,it) = vjzx(i,it)+v2i*( w1(i)*wz4(i)-w2(i)*wz3(i)
     1                                   +w3(i)*wz2(i)-w4(i)*wz1(i))
            vjxy(i,it) = vjxy(i,it)+v2i*(-w1(i)*wx3(i)-w2(i)*wx4(i)
     1                                   +w3(i)*wx1(i)+w4(i)*wx2(i))
            vjyy(i,it) = vjyy(i,it)+v2i*(-w1(i)*wy3(i)-w2(i)*wy4(i)
     1                                   +w3(i)*wy1(i)+w4(i)*wy2(i))
            vjzy(i,it) = vjzy(i,it)+v2i*(-w1(i)*wz3(i)-w2(i)*wz4(i)
     1                                   +w3(i)*wz1(i)+w4(i)*wz2(i))
            vjxz(i,it) = vjxz(i,it)+v2i*(+w1(i)*wx2(i)-w2(i)*wx1(i)
     1                                   -w3(i)*wx4(i)+w4(i)*wx3(i))
            vjyz(i,it) = vjyz(i,it)+v2i*(+w1(i)*wy2(i)-w2(i)*wy1(i)
     1                                   -w3(i)*wy4(i)+w4(i)*wy3(i))
            vjzz(i,it) = vjzz(i,it)+v2i*(+w1(i)*wz2(i)-w2(i)*wz1(i)
     1                                   -w3(i)*wz4(i)+w4(i)*wz3(i))
          enddo

        enddo
      endif

c     ................... matrix elements of the nabla for 2-body cm correction
      if (ncm2.eq.1) call nabxyz

      return
      end subroutine densit
c______________________________________________________________________________
      subroutine nabxyz

c..............................................................................
c     calculate matrix elements of the nabla operator                         .
c..............................................................................
c     Note the following properties of the matrix elements of the nabla       .
c     [cf. MB's notes on promesse where their calculation and symmetries are  .
c      discussed in detail]                                                   .
c                                                                             .
c       nablaX_kl = nablaX_-l-k = 0     nablaX_l-k = -nablaX_-lk real         .
c       nablaY_kl = nablaY_-l-k = 0     nablaY_l-k =  nablaY_-lk imaginary    .
c       nablaZ_kl = nablaY_-l-k real    nablaZ_l-k =  nablaZ_-lk = 0          .
c                                                                             .
c     meaning that tx(ibar,iket) and ty(ibar,iket) are in fact the matrix     .
c     elements between the states ibar and the conjugate state of iket that   .
c     has to be constructed from time-reversal of iket.                       .
c                                                                             .
c     what is stored as ty(ibar,iket) is i*nablaY_ibar,-iket, i.e. MINUS the  .
c     imaginary part of ty(ibar,iket)                                         .
c                                                                             .
c     Also, we have the symmetry        nabla_ij = -nabla_ji^*                .
c     which, in principle, could be used to reduce the number of matrix       .
c     elements to be calculated.                                              .
c                                                                             .
c     Also, the nabla only has matrix elements between single-particle states .
c     of opposite parity.                                                     .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (epsm5=1.0d-5)

      common /Lag  / ilag, iAnaLag iAnaLag
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv)
     1              ,p1(mv),p2(mv),p3(mv),p4(mv)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /wcm2 / tx(mw,mw),ty(mw,mw),tz(mw,mw)

c..............................................................................
      do it=1,2
        i1 = 1 + nwaven*(it-1)
        i2 =     nwaven        + nwavep*(it-1)
        do ibar=i1,i2
          do iket=i1,i2
            tx (ibar,iket) = zero
            ty (ibar,iket) = zero
            tz (ibar,iket) = zero
          enddo
        enddo
      enddo

      do it=1,2
        i1 = 1 + nwaven*(it-1)
        i2 =     nwaven        + nwavep*(it-1)
        do iket=i1,i2
          if (v2(iket).gt.epsm5) then
            call scopy (3*mq,da(1,iket),wx1)
            do ibar=i1,i2
              if (v2(ibar).gt.epsm5) then
                if (kparz(iket).ne.kparz(ibar)) then
                  call scopy (mq,a(1,ibar),p1)
                  tx(ibar,iket) = sdot(mv,p1,wx3)-sdot(mv,p2,wx4)
     1                           -sdot(mv,p3,wx1)+sdot(mv,p4,wx2)
                  ty(ibar,iket) = sdot(mv,p1,wy4)+sdot(mv,p2,wy3)
     1                           -sdot(mv,p3,wy2)-sdot(mv,p4,wy1)
                  tz(ibar,iket) = sdot(mv,p1,wz1)+sdot(mv,p2,wz2)
     1                           +sdot(mv,p3,wz3)+sdot(mv,p4,wz4)
                  tx(ibar,iket) = tx(ibar,iket)*(dv/two)
                  ty(ibar,iket) = ty(ibar,iket)*(dv/two)
                  tz(ibar,iket) = tz(ibar,iket)*(dv/two)
                endif
              endif
            enddo
          endif
        enddo
      enddo

      return
      end subroutine nabxyz

c______________________________________________________________________________
      subroutine newpot

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)

      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
c..............................................................................
      call mome
      if(ICoul.eq.1) then
        !Constructing Coulomb with first order discretisation
        call vbound
        call vcoul
      elseif(ICoul.eq.2) then
        !Using second order discretisation
        call vboundOrder2
        call vcoulOrder2
      endif

      call vcal
      call vcalg
      call vcalte
      call vcalq

      return
      end subroutine newpot

c______________________________________________________________________________
      subroutine mome

c..............................................................................
c     compute even multipole moments of the density up to the order 4         .
c     calculate multipole moments with cutoff used for the constraint.        .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (epsm4=1.0d-4,epsm8=1.0d-8,t180=180.0d0)
      parameter (t1000=1000.d0,t10=10.d0,t1pm12=1.0d-12)

      common /champ/ qxxn,qyyn,qzzn,qrrn,
     1               qxxp,qyyp,qzzp,qrrp,
     2               qxxt,qyyt,qzzt,qrrt
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),
     1               imtd,imtg,icutq
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
     4               q2fin,g2fin,q2cst,qtc0,gtc0
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /ddcut/ rhoc(mv),drhoc(mv)
      common /den  / rhon(mv),rhop(mv)
      common /nxyz / dx,dv
      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /mud  / xi (mx,my,mz),yi (mx,my,mz),zi (mx,my,mz)
     1              ,xii(mx,my,mz),yii(mx,my,mz),zii(mx,my,mz)
      common /kfmom/ tfac,tfac1,tfac2,s3

      dimension rhonc(mv),rhopc(mv)
      dimension densc(mv)

c..............................................................................
      traf = tt4*atan2(one,one)/t180                    ! 3.141592653d0/180.0d0

c     .................. calculation of the density dependent cut-off functions
      select case(icutq)
      case(0)
        if (imtd.eq.0.or.imtd.eq.1) then
          do i=1,mv
            densc(i) = rhon(i)+rhop(i)
          enddo
        else
          do i=1,mv
            densc(i) = rhop(i)
          enddo
        endif
          do i=1,mv
           rhoc (i) = one
           drhoc(i) = one
          enddo
          racut = smaxf(mv,densc)/t10
          call distance(densc,racut)
       case(1)
          do i=1,mv
            rhot     =(rhon(i)+rhop(i))
            frhoc    = rhot/(-rcut)
            rhoc(i)  = tanh(frhoc)
            drhoc(i) = frhoc*(one-rhoc(i)**2)+rhoc(i)
          enddo
        case(2)
          do i=1,mv
            rhoc (i) = one
            drhoc(i) = one
          enddo
      end select

c     ........................................... moments of the proton density
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
      qrpc =     (z2pc+x2pc+y2pc)*dv
      qpc0 = sqrt((two/tt3)*(qxpc**2+qypc**2+qzpc**2))
      if ((abs(qxpc-qypc)+abs(qzpc)).le.epsm4) then
        gpc0 = zero
      else
        gpc0 = atan2(qxpc-qypc,qzpc*s3)/traf
      endif

      qxp = (two*x2p-y2p-z2p)
      qyp = (two*y2p-x2p-z2p)
      qzp = (two*z2p-x2p-y2p)
      qrp = (z2p+x2p+y2p)

c     .......................................... moments of the neutron density
      an  = ssum(mv,rhon)
      x2n = sdot(mv,rhon,xii)
      y2n = sdot(mv,rhon,yii)
      z2n = sdot(mv,rhon,zii)

      if (imtd.eq.2.or.imtd.eq.3) then
        do i=1,mv
          densc(i) = rhon(i)
        enddo
        racut = smaxf(mv,densc)/t10
        call distance(densc,racut)
      endif
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
      if ((abs(qxnc-qync)+abs(qznc)).le.epsm4) then
        gnc0 = zero
      else
        gnc0 = atan2(qxnc-qync,qznc*s3)/traf
      endif

      qnc0 = sqrt((two/tt3)*(qxnc**2+qync**2+qznc**2))
      qxn = (two*x2n-y2n-z2n)
      qyn = (two*y2n-x2n-z2n)
      qzn = (two*z2n-x2n-y2n)
      qrn = (z2n+x2n+y2n)

c     ............................. moments of the total density (with cut-off)
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

      qxt = qxp+qxn
      qyt = qyp+qyn
      qzt = qzp+qzn
      qrt = qrp+qrn
      qt0 = sqrt((two/tt3)*(qxt**2+qyt**2+qzt**2))

c     .......................................... calculation of the constraints
      if (imtd.eq.1) then
        if (imtg.eq.0) then
          if(cq2.ne.0.0d0) then
            qxcstt = qxcstt-epscst*(qxt-qxfint)
            qycstt = qycstt-epscst*(qyt-qyfint)
            qzcstt = qzcstt-epscst*(qzt-qzfint)
          endif
          if(cqr.ne.0.0d0) then
            qrcstt = qrcstt-epscst*(qrt-qrfint)
          endif
        else
          q2cst  = q2cst -epscst*(qt0-q2fin )
        endif
      endif
      if (imtd.le.1) then
        if (imtg.eq.0) then
          pentext = two*cq2*(qxtc-qxcstt)
          penteyt = two*cq2*(qytc-qycstt)
          pentezt = two*cq2*(qztc-qzcstt)
          pentert = two*cqr*(qrtc-qrcstt)
          qqxx    = (two*pentext-penteyt-pentezt)*ral
          qqyy    = (two*penteyt-pentext-pentezt)*ral
          qqzz    = (two*pentezt-pentext-penteyt)*ral
          qqrr    =      pentert                 *ral
          excstt  = cq2*(qxtc-qxcstt)**2
          eycstt  = cq2*(qytc-qycstt)**2
          ezcstt  = cq2*(qztc-qzcstt)**2
          ercstt  = cqr*(qrtc-qrcstt)**2
        else
          qx      = tt4*cq2*(one-q2cst/(qtc0+t1pm12))
          pentext = two*cq2*(qtc0-q2cst)
          penteyt = zero
          pentezt = zero
          pentert = zero
          excstt  = cq2*(qtc0-q2cst)**2
          eycstt  = zero
          ezcstt  = zero
          ercstt  = zero
          qqxx    = qx*qxtc*ral
          qqyy    = qx*qytc*ral
          qqzz    = qx*qztc*ral
          qqrr    = zero
        endif
        if (abs(qxxt).lt.epsm8) qqxx = qqxx/ral
        if (abs(qyyt).lt.epsm8) qqyy = qqyy/ral
        if (abs(qzzt).lt.epsm8) qqzz = qqzz/ral
        if (abs(qrrt).lt.epsm8) qqrr = qqrr/ral
        qxxt   = qqxx+(one-ral)*qxxt
        qyyt   = qqyy+(one-ral)*qyyt
        qzzt   = qqzz+(one-ral)*qzzt
        qrrt   = qqrr+(one-ral)*qrrt
        qxxn   = qxxt
        qyyn   = qyyt
        qzzn   = qzzt
        qrrn   = qrrt
        qxxp   = qxxt
        qyyp   = qyyt
        qzzp   = qzzt
        qrrp   = qrrt
      endif

      if (imtd.eq.3) then
        if(cq2.ne.0.0d0) then
          qxcstn = qxcstn-epscst*(qxn-qxfinn)
          qycstn = qycstn-epscst*(qyn-qyfinn)
          qzcstn = qzcstn-epscst*(qzn-qzfinn)
          qxcstp = qxcstp-epscst*(qxp-qxfinp)
          qycstp = qycstp-epscst*(qyp-qyfinp)
          qzcstp = qzcstp-epscst*(qzp-qzfinp)
        endif
        if(cqr.ne.0.0d0) then
          qrcstn = qrcstn-epscst*(qrn-qrfinn)
          qrcstp = qrcstp-epscst*(qrp-qrfinp)
        endif
      endif
      if (imtd.ge.2) then
        cq2n=cq2*(an+ap)/an
        cqrn=cqr*(an+ap)/an
        cq2p=cq2*(an+ap)/ap
        cqrp=cqr*(an+ap)/ap
        pentexn = two*cq2n*(qxnc-qxcstn)
        penteyn = two*cq2n*(qync-qycstn)
        pentezn = two*cq2n*(qznc-qzcstn)
        pentern = two*cqrn*(qrnc-qrcstn)
        pentexp = two*cq2p*(qxpc-qxcstp)
        penteyp = two*cq2p*(qypc-qycstp)
        pentezp = two*cq2p*(qzpc-qzcstp)
        penterp = two*cqrp*(qrpc-qrcstp)
        qqxxn   = (two*pentexn-penteyn-pentezn)*ral
        qqyyn   = (two*penteyn-pentexn-pentezn)*ral
        qqzzn   = (two*pentezn-pentexn-penteyn)*ral
        qqrrn   =      pentern                 *ral
        qqxxp   = (two*pentexp-penteyp-pentezp)*ral
        qqyyp   = (two*penteyp-pentexp-pentezp)*ral
        qqzzp   = (two*pentezp-pentexp-penteyp)*ral
        qqrrp   =      penterp                 *ral
        excstn  = cq2n*(qxnc-qxcstn)**2
        eycstn  = cq2n*(qync-qycstn)**2
        ezcstn  = cq2n*(qznc-qzcstn)**2
        ercstn  = cqrn*(qrnc-qrcstn)**2
        excstp  = cq2p*(qxpc-qxcstp)**2
        eycstp  = cq2p*(qypc-qycstp)**2
        ezcstp  = cq2p*(qzpc-qzcstp)**2
        ercstp  = cqrp*(qrpc-qrcstp)**2
        if (abs(qxxn).lt.epsm8) qqxxn = qqxxn/ral
        if (abs(qyyn).lt.epsm8) qqyyn = qqyyn/ral
        if (abs(qzzn).lt.epsm8) qqzzn = qqzzn/ral
        if (abs(qrrn).lt.epsm8) qqrrn = qqrrn/ral
        if (abs(qxxp).lt.epsm8) qqxxp = qqxxp/ral
        if (abs(qyyp).lt.epsm8) qqyyp = qqyyp/ral
        if (abs(qzzp).lt.epsm8) qqzzp = qqzzp/ral
        if (abs(qrrp).lt.epsm8) qqrrp = qqrrp/ral
        qxxn    = qqxxn+(one-ral)*qxxn
        qyyn    = qqyyn+(one-ral)*qyyn
        qzzn    = qqzzn+(one-ral)*qzzn
        qrrn    = qqrrn+(one-ral)*qrrn
        qxxp    = qqxxp+(one-ral)*qxxp
        qyyp    = qqyyp+(one-ral)*qyyp
        qzzp    = qqzzp+(one-ral)*qzzp
        qrrp    = qqrrp+(one-ral)*qrrp
      endif

      return
      end subroutine mome

c______________________________________________________________________________
      subroutine vbound

c..............................................................................
c     calculation of the coulomb boundary values                              .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (mpx=mx+mc,mpy=my+mc,mpz=mz+mc,mpv=mpx*mpy*mpz)
      parameter (one=1.0d0,two=2.0d0,tt8=8.0d0,tp375=0.375d0)

      common /bound/ xbound(mpy,mpz),ybound(mpx,mpz),zbound(mpx,mpy)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /mudd / xu(mpx+mpy+mpz)
      common /nxyz / dx,dv

c     .........................................................................
      c1 = (two*z2p-x2p-y2p)/tt8
      c2 = tp375*(x2p-y2p)

c     .................................. rearrangement and multiplication by e2
      c3 =-two * (c2+c1)*e2
      c2 = two * (c2-c1)*e2
      c1 = two * two*c1 *e2
      c0 = two *(ap/two)*e2

c     .................................................................. xbound
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

c     .................................................................. ybound
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

c     .................................................................. zbound
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
      end subroutine vbound

c______________________________________________________________________________
      subroutine vboundOrder2

c..............................................................................
c     calculation of the coulomb boundary values, for use in a second order   .
c     laplacian. Concept is completely analogous to the older subroutine      .
c     vbound, only more points calculated.                                    .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (mpx=mx+mc,mpy=my+mc,mpz=mz+mc,mpv=mpx*mpy*mpz)
      parameter (one=1.0d0,two=2.0d0,tt8=8.0d0,tp375=0.375d0)

      common /bound2/xbound(mpy,mpz,2),ybound(mpx,mpz,2),
     1               zbound(mpx,mpy,2)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /mudd / xu(mpx+mpy+mpz)
      common /nxyz / dx,dv

c     Numerical constants, to put in front of the multipole factors
      c1 = (two*z2p-x2p-y2p)/tt8
      c2 = tp375*(x2p-y2p)

c     .................................. rearrangement and multiplication by e2
      c3 =-two * (c2+c1)*e2
      c2 = two * (c2-c1)*e2
      c1 = two * two*c1 *e2
      c0 = two *(ap/two)*e2

c     .................................................................. xbound

      do n=1,2
          ! xu only stores coordinates up to one outside the box. We need one
          ! more and do it this way to not disturb the other code.
        x2   = (xu(nnx+1) + (n-1)*dx)**2
        c2x2 = c2*x2
        do j=1,nny
          y2   = xu(j)**2
          c3y2 = c3*y2
          do k=1,nnz
            z2   = xu(k)**2
            c1z2 = c1*z2
            r    = one/sqrt(x2+y2+z2)
            xbound(j,k,n)=c0*r+(c1z2+c2x2+c3y2)*r**5
          enddo
        enddo
      enddo

c     .................................................................. ybound
      do n=1,2
          ! xu only stores coordinates up to one outside the box. We need one
          ! more and do it this way to not disturb the other code.
        y2   = (xu(nny+1) + (n-1)*dx)**2
        c3y2 = c3*y2
        do i=1,nnx
          x2   = xu(i)**2
          c2x2 = c2*x2
          do k=1,nnz
            z2   = xu(k)**2
            c1z2 = c1*z2
            r    = one/sqrt(x2+y2+z2)
            ybound(i,k,n)=c0*r+(c1z2+c2x2+c3y2)*r**5
          enddo
        enddo
      enddo
c     .................................................................. zbound
      do n=1,2
          ! xu only stores coordinates up to one outside the box. We need one
          ! more and do it this way to not disturb the other code.
        z2   = (xu(nnz+1) + (n-1)*dx)**2
        c1z2 = c1*z2
        do i=1,nnx
          x2   = xu(i)**2
          c2x2 = c2*x2
          do j=1,nny
            y2   = xu(j)**2
            c3y2 = c3*y2
            r    = one/sqrt(x2+y2+z2)
            zbound(i,j,n)=c0*r+(c1z2+c2x2+c3y2)*r**5
          enddo
        enddo
      enddo

      return
      end subroutine vboundOrder2

c______________________________________________________________________________
      subroutine vcoulOrder2

c..............................................................................
c Subroutine that uses the conjugate gradients algorithm to solve the Coulomb
c problem in the nucleus.
c Discretisation of the laplacian of order 2 and boundary conditions are
c supplied by the common block /bound2/ and are calculated in vboundOrder2
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (mpx=mx+mc,mpy=my+mc,mpz=mz+mc,mpv=mpx*mpy*mpz)
      parameter (zero=0.0d0,tt3=3.0d0)
      parameter (nitmax=300+ 8*mc)

      parameter (Coef1=-15.0d0/2.0d0, Coef2=4.0d0/3.0d0)
      parameter (Coef3=-1.0d0/12.0d0)

      common /bound2/ xbound(mpy,mpz,2),ybound(mpx,mpz,2),
     1                zbound(mpx,mpy,2)
      common /den  / rho(mx,my,mz,2)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mx,my,mz),wce(mv)
     1              ,wt3a(mv),wt3b(mv)
      common /wcl  / t(mpx,mpy,mpz),p(mpx,mpy,mpz),z(mpx,mpy,mpz)
     1              ,w(mpx,mpy,mpz),r(mpx,mpy,mpz)

      data ikal /0/

c..............................................................................
 100  format (' error coulomb nitmax,epscl,zz1,zz2,zt,a,c',/,i5,6e15.8)
 102  format (4x,' coulomb iter=',i5,'  (z,z)=',e15.8)

c  ........................................................... initalisation
      if (ikal.eq.0) then
c      If this is the first time this subroutine is run, we guess Zero for the
c      Coulomb potential.
        ikal = 1
        do i=1,mpv
          w(i,1,1) = zero
          r(i,1,1) = zero
        enddo
      endif

      do k=1,mz
      do j=1,my
      do i=1,mx
        r(i,j,k) = rho(i,j,k,2) !Only take the proton value
      enddo
      enddo
      enddo

c     ............................................... calculation of delta(wcd)
c     using second order Laplacian discretisation. (It is here that this
c     routine differs from the old vcoul.)
      do i=1,mpv
        z(i,1,1) =  Coef1*w(i,1,1)
      enddo

      do k=1,nnz
      do j=1,nny

        !Forward part of the finite difference
        !Points at the start of the line
        do i=1,nnx-2
          z(i,j,k)   = z(i,j,k)   + Coef2*w(i+1,j,k)
     1                            + Coef3*w(i+2,j,k)
        enddo
        !Points at (end-of-the-line -1)
        z(nnx-1,j,k) = z(nnx-1,j,k) + Coef2*w(nnx,j,k)
     1                              + Coef3*xbound(j,k,1)
        !Points at the end of the line
        z(nnx,j,k)   = z(nnx,j,k) + Coef2*xbound(j,k,1)
     1                            + Coef3*xbound(j,k,2)

        !Backward part of the finite difference
        !Points at the start of the line
          z(1,j,k)   = z(1,j,k)   + Coef2*w(1,j,k)
     1                            + Coef3*w(2,j,k)
          z(2,j,k)   = z(2,j,k)   + Coef2*w(1,j,k)
     1                            + Coef3*w(1,j,k)
        !Points not at the start of the line
        do i=3,nnx
          z(i,j,k)   = z(i,j,k)   + Coef2*w(i-1,j,k)
     1                            + Coef3*w(i-2,j,k)
        enddo
      enddo
      enddo

      do k=1,nnz
      do i=1,nnx

        !Forward part of the finite difference
        !Points at the start of the line
        do j=1,nny-2
          z(i,j,k)   = z(i,j,k)   + Coef2*w(i,j+1,k)
     1                            + Coef3*w(i,j+2,k)
        enddo
        !Points at (end-of-the-line -1)
        z(i,nny-1,k) = z(i,nny-1,k) + Coef2*w(i,nny,k)
     1                              + Coef3*ybound(i,k,1)
        !Points at the end of the line
        z(i,nny,k)   = z(i,nny,k)   + Coef2*ybound(i,k,1)
     1                              + Coef3*ybound(i,k,2)

        !Backward part of the finite difference
        !Points at the start of the line
          z(i,1,k)   = z(i,1,k)   + Coef2*w(i,1,k)
     1                            + Coef3*w(i,2,k)
          z(i,2,k)   = z(i,2,k)   + Coef2*w(i,1,k)
     1                            + Coef3*w(i,1,k)
        !Points not at the start of the line
        do j=3,nny
          z(i,j,k)   = z(i,j,k)   + Coef2*w(i,j-1,k)
     1                            + Coef3*w(i,j-2,k)
        enddo
      enddo
      enddo

      do j=1,nny
      do i=1,nnx

        !Forward part of the finite difference
        !Points at the start of the line
        do k=1,nnz-2
          z(i,j,k)   = z(i,j,k)   + Coef2*w(i,j,k+1)
     1                            + Coef3*w(i,j,k+2)
        enddo
        !Points at (end-of-the-line -1)
        z(i,j,nnz-1) = z(i,j,nnz-1) + Coef2*w(i,j,nnz)
     1                              + Coef3*zbound(i,j,1)
        !Points at the end of the line
        z(i,j,nnz)   = z(i,j,nnz)   + Coef2*zbound(i,j,1)
     1                              + Coef3*zbound(i,j,2)

        !Backward part of the finite difference
        !Points at the start of the line
          z(i,j,1)   = z(i,j,1)   + Coef2*w(i,j,1)
     1                            + Coef3*w(i,j,2)
          z(i,j,2)   = z(i,j,2)   + Coef2*w(i,j,1)
     1                            + Coef3*w(i,j,1)
        !Points not at the start of the line
        do k=3,nnz
          z(i,j,k)   = z(i,j,k)   + Coef2*w(i,j,k-1)
     1                            + Coef3*w(i,j,k-2)
        enddo
      enddo
      enddo

      fac = 1.d0/(dx*dx)
      do i=1,mpv
        z(i,1,1) = z(i,1,1) *fac
      enddo

c............................................. Calculating the Poisson Equation
      do i=1,mpv
        z(i,1,1)=(-e2eff)*r(i,1,1)-z(i,1,1)
        p(i,1,1)=z(i,1,1)
      enddo

c............ Starting the iterations........................
      zz1=sdot(mpv,z,z)

      do niter=1,nitmax
c............ Step 1: Calculate the laplacian of p...........
        do i=1,mpv
          t(i,1,1) = Coef1*p(i,1,1)
        enddo

          do k=1,nnz
          do j=1,nny
            !Forward part of the finite difference
            !Points at the start of the line
            do i=1,nnx-2
              t(i,j,k)   = t(i,j,k)   + Coef2*p(i+1,j,k)
     1                                + Coef3*p(i+2,j,k)
            enddo
            !Notice that the boundary conditions are absent: this
            !is because they already solve poissons equation
            !Points at (end-of-the-line -1)
            t(nnx-1,j,k) = t(nnx-1,j,k) + Coef2*p(nnx,j,k)


            !Backward part of the finite difference
            !Points at the start of the line
              t(1,j,k)   = t(1,j,k)   + Coef2*p(1,j,k)
     1                                + Coef3*p(2,j,k)
              t(2,j,k)   = t(2,j,k)   + Coef2*p(1,j,k)
     1                                + Coef3*p(1,j,k)
            !Points not at the start of the line
            do i=3,nnx
              t(i,j,k)   = t(i,j,k)   + Coef2*p(i-1,j,k)
     1                                + Coef3*p(i-2,j,k)
            enddo
          enddo
          enddo

          do k=1,nnz
          do i=1,nnx

            !Forward part of the finite difference
            !Points at the start of the line
            do j=1,nny-2
              t(i,j,k)   = t(i,j,k)   + Coef2*p(i,j+1,k)
     1                                + Coef3*p(i,j+2,k)
            enddo
            !Points at (end-of-the-line -1)
            t(i,nny-1,k) = t(i,nny-1,k) + Coef2*p(i,nny,k)
            !Backward part of the finite difference
            !Points at the start of the line
              t(i,1,k)   = t(i,1,k)   + Coef2*p(i,1,k)
     1                                + Coef3*p(i,2,k)
              t(i,2,k)   = t(i,2,k)   + Coef2*p(i,1,k)
     1                                + Coef3*p(i,1,k)
            !Points not at the start of the line
            do j=3,nny
              t(i,j,k)   = t(i,j,k)   + Coef2*p(i,j-1,k)
     1                                + Coef3*p(i,j-2,k)
            enddo
          enddo
          enddo

          do j=1,nny
          do i=1,nnx

            !Forward part of the finite difference
            !Points at the start of the line
            do k=1,nnz-2
              t(i,j,k)   = t(i,j,k)   + Coef2*p(i,j,k+1)
     1                                + Coef3*p(i,j,k+2)
            enddo
            !Points at (end-of-the-line -1)
            t(i,j,nnz-1) = t(i,j,nnz-1) + Coef2*p(i,j,nnz)

            !Backward part of the finite difference
            !Points at the start of the line
            t(i,j,1)   = t(i,j,1)   + Coef2*p(i,j,1)
     1                              + Coef3*p(i,j,2)
            t(i,j,2)   = t(i,j,2)   + Coef2*p(i,j,1)
     1                              + Coef3*p(i,j,1)
            !Points not at the start of the line
            do k=3,nnz
              t(i,j,k) = t(i,j,k)   + Coef2*p(i,j,k-1)
     1                              + Coef3*p(i,j,k-2)
            enddo
          enddo
      enddo

      do i=1,mpv
        t(i,1,1) = (fac)*t(i,1,1)
      enddo

c................. Step 2: Calculate the various inproducts
        zt = sdot(mpv,z,t)
        a  = zz1/zt
        call saxpy (mpv,a,p,w)
        call saxpy (mpv,-a,t,z)
        zz2 = sdot(mpv,z,z)
c................. Step 3: Check for convergence
        if (zz2.lt.epscl) then
          if(iprintcoul.eq.1) print 102,niter,zz2
          go to 111
        endif
        c   = zz2/zz1
        zz1 = zz2
c................. Step 4: adjust the potential
        do i=1,mpv
          p(i,1,1)=z(i,1,1)+c*p(i,1,1)
        enddo
      enddo

      if(iprintcoul.eq.1) print 100,nitmax,epscl,zz1,zz2,zt,a,c
c................. Copy the potential to the common block
 111  do k=1,mz
      do j=1,my
      do i=1,mx
        wcd(i,j,k) = w(i,j,k)
      enddo
      enddo
      enddo

      return
      end subroutine vcoulOrder2

c______________________________________________________________________________
      subroutine vcoul

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (mpx=mx+mc,mpy=my+mc,mpz=mz+mc,mpv=mpx*mpy*mpz)
      parameter (zero=0.0d0,tt6=6.0d0)
      parameter (nitmax=100+8*mc)

      common /bound/ xbound(mpy,mpz),ybound(mpx,mpz),zbound(mpx,mpy)
      common /den  / rho(mx,my,mz,2)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mx,my,mz),wce(mv)
     1              ,wt3a(mv),wt3b(mv)
      common /wcl  / t(mpx,mpy,mpz),p(mpx,mpy,mpz),z(mpx,mpy,mpz)
     1              ,w(mpx,mpy,mpz),r(mpx,mpy,mpz)
      data ikal /0/

c..............................................................................
 100  format (' error coulomb nitmax,epscl,zz1,zz2,zt,a,c',/,i5,6e15.8)
 102  format (4x,' coulomb iter=',i5,'  (z,z)=',e15.8)

c     ........................................................... initalisation
      if (ikal.eq.0) then
        ikal = 1
        do i=1,mpv
          w(i,1,1) = zero
          r(i,1,1) = zero
        enddo
      endif

      do k=1,mz
      do j=1,my
      do i=1,mx
        r(i,j,k) = rho(i,j,k,2)
      enddo
      enddo
      enddo

c     ............................................... calculation of delta(wcd)
      do i=1,mpv
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

      do i=1,mpv
        z(i,1,1) = z(i,1,1) / (dx*dx)
      enddo

c     .................................. calculation of z=p=(source-delta(wcd))
      do i=1,mpv
        z(i,1,1)=(-e2eff)*r(i,1,1)-z(i,1,1)
        p(i,1,1)=z(i,1,1)
      enddo
      zz1=sdot(mpv,z,z)

c     ............................................. beginning of the iterations
c                         calculation of tk+1=delta(pk)
c                                     of (z.z),ak+1,wcdk+1,zk+1 and (zk+1.zk+1)
c                                 and of ck+1 and pk+1
c     .........................................................................
      do niter=1,nitmax

        do i=1,mpv
          t(i,1,1) = -tt6*p(i,1,1)
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

        do i=1,mpv
          t(i,1,1) = t(i,1,1) / (dx*dx)
        enddo

        zt = sdot(mpv,z,t)
        a  = zz1/zt
        call saxpy (mpv,a,p,w)
        call saxpy (mpv,-a,t,z)
        zz2 = sdot(mpv,z,z)

        if (zz2.lt.epscl) then
          if(iprintcoul.eq.1) print 102,niter,zz2
          go to 111
        endif
        c   = zz2/zz1
        zz1 = zz2
        do i=1,mpv
          p(i,1,1)=z(i,1,1)+c*p(i,1,1)
        enddo
      enddo

       if(iprintcoul.eq.1)  print 100,nitmax,epscl,zz1,zz2,zt,a,c

 111  do k=1,mz
      do j=1,my
      do i=1,mx
        wcd(i,j,k) = w(i,j,k)
      enddo
      enddo
      enddo

      return
      end subroutine vcoul

c______________________________________________________________________________
      subroutine vcal

c..............................................................................
c     calculation of the central potential U                                  .
c     including the Coulomb contribution for protons                          .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (one=1.0d0,two=2.0d0,tt3=3.0d0,epsm20=1.0d-20)

      common /Lag   / ilag, iAnaLag
      common /den  / rho(mv,2)
      common /drho / drhon(mv),drhop(mv)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /taudj/ vtau(mv,2),vdiv(mv,2)

c..............................................................................
      if (ilag.eq.0) then
        call lapla (rho(1,1),drhon,one,one,one)
        call lapla (rho(1,2),drhop,one,one,one)
      else
        call laplalag (rho(1,1),drhon,1,1,1)
        call laplalag (rho(1,2),drhop,1,1,1)
      endif

      do i=1,mv
        wt3a(i) = (rho(i,1)+rho(i,2))**byt3a
      enddo

      if (ndd.eq.1) then
        do i=1,mv
          wt3b(i) = (rho(i,1)+rho(i,2))**byt3b
        enddo
      else
        do i=1,mv
          wt3b(i) = 0.0d0
        enddo
      endif

      if (ncoex.eq.0) then
        do i=1,mv
          wce(i) = coexv * rho(i,2)**(one/tt3)
        enddo
      else
        do i=1,mv
          wce(i) = 0.0d0
        enddo
      endif

      eps = epsm20

      do i=1,mv
        xn   = rho(i,1)
        xp   = rho(i,2)
        xt   = xn+xp
        xnp  = xn**2+xp**2
        dxn  = drhon(i)
        dxp  = drhop(i)
        dxt  = dxn+dxp
        vtnp = vtau(i,1)+vtau(i,2)
        vdnp = vdiv(i,1)+vdiv(i,2)
c             .............................................ATTENTION
c             B5 & B6 as used in the code are not equal to the ones
c             discussed in all of the articles; they carry a 
c             relative minus sign in this code!
c             ...................................................... 
        www  = two*(b1*xt-b5*dxt)
     1       + (b7a*(two+byt3a)*(xt+eps)+b8a*byt3a*xnp/(xt+eps))*wt3a(i)
     2       + (b7b*(two+byt3b)*(xt+eps)+b8b*byt3b*xnp/(xt+eps))*wt3b(i)
     3       + b3*vtnp+b9*vdnp
        wnn(i) = www + two*(b2*xn-b6*dxn+(b8a*wt3a(i)+b8b*wt3b(i))*xn)
     1               + b4*vtau(i,1)+b9q*vdiv(i,1)
        wpp(i) = www + two*(b2*xp-b6*dxp+(b8a*wt3a(i)+b8b*wt3b(i))*xp)
     1               + b4*vtau(i,2)+b9q*vdiv(i,2)
     2               + wcd(i)+wce(i)
      enddo

      return
      end subroutine vcal

c______________________________________________________________________________
      subroutine vcalg

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /Lag   / ilag, iAnaLag
      common /den  / rho(mv,2)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /evcen/ r(mq,2)
      common /evpro/ rx(mv,2),ry(mv,2),rz(mv,2)
      common /wtmp / ta(mq)

c..............................................................................
      do it=1,2
        do i=1,mv
          r(i,it)       = b3*(rho(i,1)+rho(i,2)) + b4*rho(i,it)
          r(i+mv,it)    = r(i,it)
          r(i+2*mv,it)  = r(i,it)
          r(i+3*mv,it)  = r(i,it)
          ta(i)         = b9*(rho(i,1)+rho(i,2)) + b9q*rho(i,it)
        enddo

        if(ilag.eq.0) then
            call der (ta,rx(1,it),ry(1,it),rz(1,it),1,1,1)
        else
            call derlagx(1,ta,rx(1,it))
            call derlagy(1,ta,ry(1,it))
            call derlagz(1,ta,rz(1,it))
        endif
      enddo

      return
      end subroutine vcalg

c______________________________________________________________________________
      subroutine vcalte

c..............................................................................
c     tensor contribution to the spin-orbit potential                         .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (two=2.d0)

      common /den  / rho(mv,2)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /wmunu/ wjxx(mv,2),wjyx(mv,2),wjzx(mv,2)
     1              ,wjxy(mv,2),wjyy(mv,2),wjzy(mv,2)
     2              ,wjxz(mv,2),wjyz(mv,2),wjzz(mv,2)
      common /wj2  / vjxx(mv,2),vjyx(mv,2),vjzx(mv,2)
     1              ,vjxy(mv,2),vjyy(mv,2),vjzy(mv,2)
     2              ,vjxz(mv,2),vjyz(mv,2),vjzz(mv,2)
      common /wtmp / wt1(mv),wt2(mv),wt3(mv),wt4(mv)

c........................................ contribution of the spin-orbit tensor
c                                                 from the central Skyrme force
c                                            to the single-particle Hamiltonian
      if (njmunu.eq.1) then
        bb1 = two * (b14+b15)
        bb2 = two *  b14
        do it=1,2
          do i=1,mv
            wjxx(i,it) = bb1 * vjxx(i,it) + bb2 * vjxx(i,3-it)
            wjyx(i,it) = bb1 * vjyx(i,it) + bb2 * vjyx(i,3-it)
            wjzx(i,it) = bb1 * vjzx(i,it) + bb2 * vjzx(i,3-it)
            wjxy(i,it) = bb1 * vjxy(i,it) + bb2 * vjxy(i,3-it)
            wjyy(i,it) = bb1 * vjyy(i,it) + bb2 * vjyy(i,3-it)
            wjzy(i,it) = bb1 * vjzy(i,it) + bb2 * vjzy(i,3-it)
            wjxz(i,it) = bb1 * vjxz(i,it) + bb2 * vjxz(i,3-it)
            wjyz(i,it) = bb1 * vjyz(i,it) + bb2 * vjyz(i,3-it)
            wjzz(i,it) = bb1 * vjzz(i,it) + bb2 * vjzz(i,3-it)
          enddo

        enddo
      endif

c     ................................... contribution of the spin-orbit tensor
c                               from the central + genuine tensor Skyrme forces
c                                            to the single-particle Hamiltonian
      if (njmunu.eq.2) then
        bb1 = two * (b14 + b15)
        bb2 = two *  b14
        bb3 = two * (b16 + b17)
        bb4 = two *  b16
        cc1 = bb1 + bb3
        cc2 = bb2 + bb4

        do it=1,2
          do i=1,mv
            wjxx(i,it) = cc1 * vjxx(i,it) + cc2 * vjxx(i,3-it)
            wjyy(i,it) = cc1 * vjyy(i,it) + cc2 * vjyy(i,3-it)
            wjzz(i,it) = cc1 * vjzz(i,it) + cc2 * vjzz(i,3-it)
          enddo
          do i=1,mv
            wjyx(i,it) = bb1 * vjyx(i,it) + bb2 * vjyx(i,3-it)
     1                  +bb3 * vjxy(i,it) + bb4 * vjxy(i,3-it)
            wjxy(i,it) = bb1 * vjxy(i,it) + bb2 * vjxy(i,3-it)
     1                  +bb3 * vjyx(i,it) + bb4 * vjyx(i,3-it)
            wjzx(i,it) = bb1 * vjzx(i,it) + bb2 * vjzx(i,3-it)
     1                  +bb3 * vjxz(i,it) + bb4 * vjxz(i,3-it)
            wjxz(i,it) = bb1 * vjxz(i,it) + bb2 * vjxz(i,3-it)
     1                  +bb3 * vjzx(i,it) + bb4 * vjzx(i,3-it)
            wjzy(i,it) = bb1 * vjzy(i,it) + bb2 * vjzy(i,3-it)
     1                  +bb3 * vjyz(i,it) + bb4 * vjyz(i,3-it)
            wjyz(i,it) = bb1 * vjyz(i,it) + bb2 * vjyz(i,3-it)
     1                  +bb3 * vjzy(i,it) + bb4 * vjzy(i,3-it)
          enddo

        enddo
      endif

      return
      end subroutine vcalte

c______________________________________________________________________________
      subroutine vcalq

c..............................................................................
c     contribution of the constraints to the potentials                       .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),
     1               imtd,imtg,icutq
      common /champ/ qxxn,qyyn,qzzn,qrrn,
     1               qxxp,qyyp,qzzp,qrrp,
     2               qxxt,qyyt,qzzt,qrrt
      common /ddcut/ rhoc(mv),drhoc(mv)
      common /den  / rhon(mv),rhop(mv)
      common /mud  / xi(mv),yi(mv),zi(mv),xii(mv),yii(mv),zii(mv)
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)

      qxn = qxxn + qrrn
      qyn = qyyn + qrrn
      qzn = qzzn + qrrn
      qxp = qxxp + qrrp
      qyp = qyyp + qrrp
      qzp = qzzp + qrrp
      do i=1,mv
        qct    = cutof2(i)*drhoc(i)
        wnn(i) = wnn(i) + qct * (qxn*xii(i)+qyn*yii(i)+qzn*zii(i))
        wpp(i) = wpp(i) + qct * (qxp*xii(i)+qyp*yii(i)+qzp*zii(i))
      enddo

      return
      end subroutine vcalq

c______________________________________________________________________________
      subroutine fsum(itprii)

c..............................................................................
c     Subroutine that prints out all info at the end of the iteration process.
c
c     1) Print out all parameters (Coupling constants, ....)
c     2) Recalculate energies with lagrangian derivatives, and print
c     3) Compute all multipole moments up to l=10
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor
      character*4  head

      parameter (zero=0.0d0,one=1.0d0)

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv),
     1               imtd,imtg,icutq
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
     4               q2fin,g2fin,q2cst,qtc0,gtc0
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ t0,x0,t1,x1,t2,x2,t3a,x3a,yt3a
     1                                ,t3b,x3b,yt3b
     2              ,te,to,wso,wsoq
     2              ,hbar,hbm(2),xm(3),afor
      common /info / bidon(500),head(20),irb,nx,ny,nz
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /pairn/ vn,vp,rangen,rangep
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)

c..............................................................................

  201 format (/,' ',78('_'),
     1        /,'   **final**  '
     3        /,' parameters: mx,my,mz (',3i3,' ) mw (',i4
     4         ,' ) mc (', i3 , ' )' ,
     5        /,' Neutron & proton wavefunctions ', 2i3,' ' )
  203 format   (' ',78('_'),/,20a4)
  204 format   (' dx=',f8.3,' dt=',f8.3,'  number of iterations =',i5)
  205 format   (' nprint = ',i5,' ndiag  = ',i5)
  260 format   (' npair  = ',i5,' -- hartree-fock --')
  261 format   (' npair  = ',i5,' -- BCS seniority force --')
  262 format   (' npair  = ',i5,' -- BCS constant delta --')
  263 format   (' npair  = ',i5,' -- BCS+LN  seniority force --')
  264 format   (' npair  = ',i5,' -- BCS delta force --')
  265 format   (' npair  = ',i5,' -- BCS+LN delta force --')
  271 format   ('  gn =',f8.3,'/(11+n),  encut=',f6.3,
     2        /,'  gp =',f8.3,'/(11+z),  epcut=',f6.3,'  dcut =',f6.3)
  272 format   ('  deltan =',f8.3,'  encut =',f6.3,
     1        /,'  deltap =',f8.3,'  epcut =',f6.3,'  dcut =',f6.3)
  273 format   ('  vn = ',f8.3,'  encut =',f6.3,
     2        /,'  vp = ',f8.3,'  epcut =',f6.3,'  dcut =',f6.3,
     4          '  alpha =',f6.3)
  275 format   ('    vn = ',f12.6,'  range =',f6.3,
     2        /,'    vp = ',f12.6,'  range =',f6.3)
  281 format   (' cutoff only above the fermi level ')
  282 format   (' cutoff above and below the fermi level ')
  291 format (/,' actual monopole constraints',
     1  /,'  cqr =',e12.4,
     2          '  Requested r   = ',e12.4,'  Readjusted r   = ',e12.4)
  292 format (/,' actual quadrupole constraints',
     3  /,'  cq2 =',e12.4,
     4          '  Requested Qxx = ',e12.4,'  Readjusted Qxx = ',e12.4,  
     5  /,19x  ,'  Requested Qyy = ',e12.4,'  Readjusted Qyy = ',e12.4,
     6  /,19x  ,'  Requested Qzz = ',e12.4,'  Readjusted Qzz = ',e12.4)
  293 format (/,' actual monopole constraints',
     1  /,      '  cqr =',e12.4,
     2          '  Requested N  r  = ',e12.4,
     2          '  Readj. N  r  = ',e12.4,
     3  /,19x  ,
     4          '  Requested P  r  = ',e12.4,
     4          '  Readj. P  r  = ',e12.4)
  294 format (/, ' actual quadrupole constraints',    
     5  /,      '  cq2 =',e12.4,
     6          '  Requested N Qxx = ',e12.4,
     7          '  Readj. N Qxx = ',e12.4,
     8  /,19x  ,'  Requested N Qyy = ',e12.4,
     9          '  Readj. N Qyy = ',e12.4,
     1  /,19x  ,'  Requested N Qzz = ',e12.4,
     2          '  Readj. N Qzz = ',e12.4,
     3  /,19x  ,'  Requested P Qxx = ',e12.4,
     4          '  Readj. P Qxx = ',e12.4,
     5  /,19x  ,'  Requested P Qyy = ',e12.4,
     6          '  Readj. P Qyy = ',e12.4,
     7  /,19x  ,'  Requested P Qzz = ',e12.4,
     8          '  Readj. P Qzz = ',e12.4)

c..............................................................................
      print 201,mx,my,mz,mw,mc, nwaven, nwavep
      print 203,head
      print 204,dx,dt,nitert
      print 205,nprint,ndiag
      call cpling

      if (npair.eq.0) print 260,npair
      if (npair.eq.1) print 261,npair
      if (npair.eq.2) print 262,npair
      if (npair.eq.3) print 263,npair
      if (npair.eq.4) print 264,npair
      if (npair.eq.5) print 265,npair

      if (npair.eq.1) print 271,gn,encut,gp,encut,dcut
      if (npair.eq.2) print 272,delmax(1),encut,delmax(2),epcut,dcut
      if (npair.eq.3) print 271,gn,encut,gp,encut,dcut
      if ((npair.eq.4).or.(npair.eq.5)) then
        rhoc = alpha
        if (alpha.ne.zero) rhoc = one/alpha
        print 273,gn,encut,gp,encut,dcut,rhoc
      endif

      if ((npair.ge.1).and.(npair.le.5)) then
        if (xcut.eq.zero) print 281
        if (xcut.ne.zero) print 282
      endif

      if (imtd.le.1) then
        if (cqr.ne.zero) then
          print 291,cqr,sqrt(qrfint/(npp+npn)),sqrt(qrcstt/(npn+npp))
        endif
        if (cq2.ne.zero) then
          print 292,cq2,qxfint,qxcstt,qyfint,qycstt,qzfint,qzcstt
        endif
      else
        if (cqr.ne.zero) then
          print 293,cqr,sqrt(qrfinn/npn),sqrt(qrcstn/npn)
     1             ,sqrt(qrfinp/npp),sqrt(qrcstp/npp)
        endif
        if (cq2.ne.zero) then
          print 294,cq2,qxfinn,qxcstn,qyfinn,qycstn,qzfinn,qzcstn
     1                 ,qxfinp,qxcstp,qyfinp,qycstp,qzfinp,qzcstp
        endif
      endif

      !All normal output
      call figaro(1, 0,1)

      !Reanalysing the energies with exact derivatives
      call ReanalyseLag(itprii)

      return
      end subroutine fsum

c______________________________________________________________________________
c______________________________________________________________________________
      subroutine figaro (itprii,iprintlag,ifinal)

c..............................................................................
c     calculate observables and print them to standard out                    .
c..............................................................................
c     iprintlag determines whether output is printed after a reanalysis with  .
c     lagrange derivatives (1) or not(0).                                     .
c..............................................................................
c     this subroutine prints the energies computed in                         .
c          fprte: (et0,etdecoul,ejmunu)                                       .
c          fprtj: (esmf,espo,zsp,ek,eso,esk)                                  .
c          pro  : (esp1,esp2)                                                 .
c      and gap  : (v2,v22,eqp,delta,ambda,epair,eproj,disper)                 .
c..............................................................................
c      if itprii=0   only the s.p. levels near the fermi surface are printed  .
c                    and reduced printout in fprtm                            .
c      if itprii=1   all the levels are printed.                              .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0)
      parameter (epsm8=1.0d-8,esplm=5.0d0,t50=50.0d0,tbig=1.0d+10)
      parameter (epsm4=1.0d-4)

      common /conv / econv1,econv2,iconv
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv)
     1              ,imtd,imtg,icutq
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
      common /epente/ epente
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu(3),ecm(2,3),ecmp(3)
      common /lnog / rln(3),qxln(3),qyln(3),qzln(3)
      common /mmt  / ap,an,at
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      dimension etsp(3),idh(mw),etmp(3)
      dimension v2sum(2),uvsum(2),uvDelta(2),v2Delta(2)

c..............................................................................
  100 format (/,' ',78('_'),/,'   iteration : ',i5)
  101 format (/,' neutron levels')
  102 format (/,' proton levels')
  103 format   ('   i   n',2x,'par',2x,'vbcs',5x,'vln',6x,'delta',4x,
     1          'E_qp',4x,'E_sp',4x,'jz',6x,'j',6x,'d2(h)')
  104 format   (3i4,2f9.5,5f8.3,1pe10.2)
  105 format (/,' ',78('_'))
  110 format (  ' pairing : ',
     1        /,43x,'neutron',6x,'proton',
     2        /,' lambda  (MeV)                        ', 2f12.4,
     3        /,' lambda2 (MeV)                        ', 2f12.4,
     4        /,' dispersion                           ', 2f12.4)
  111 format (/,' Sum ( UV )                           ', 2f12.4,
     1        /,' Sum ( UV  Delta)/[Sum ( UV ) ] (MeV) ', 2f12.4,
     2        /,' Sum ( V^2 Delta)/[Sum ( V^2) ] (MeV) ', 2f12.4)
  120 format (  ' total energies  (MeV)')
  121 format (/,32x,'neutron',6x,'proton',6x,'total',
     1        /,' kinetic                   ',3f12.3)
 1221 format (  ' 1-body c.m.               ',3f12.3)
  122 format (  ' 2-body c.m.               ',3f12.3)
  123 format   (' spin orbit                ',3f12.3)
  124 format   (' spin tensor (central     )',24x,f12.3,/,
     1          ' spin tensor (tensor, sym.)',24x,f12.3,/,
     2          ' spin tensor (tensor, asy.)',24x,f12.3)
 1124 format   (' spin tensor sym           ',24x,f12.3,/,
     1          ' spin tensor asy           ',24x,f12.3)
  125 format   (' single particle           ',3f12.3)
  126 format   (' pairing                   ',3f12.3)
  127 format   (' Lipkin-Nogami             ',3f12.3)
  128 format (/,' Coulomb direct            ',12x,f12.3,
     1        /,' Coulomb exchange          ',12x,f12.3)
  143 format (/,' Skyrme energies',/,
     1          ' E_[rho rho]      =',f12.3,
     2          ' E_[ rho tau ]    =',f12.3,/,
     3          ' E_[rho^(2+alpha)]=',f12.3,
     4          ' E_[rho^(2+beta)] =',f12.3,/
     5          ' E_[rho Delta rho]=',f12.3,
     6          ' E_[rho Nabla J]  =',f12.3,/,
     7          ' E_Skyrme         =',f12.3)
  144 format (/,' Energy without LN correction    ',f15.6)
  145 format (/,' total energy (from functional)  ',f15.6,/,
     1          '              (from sp energies) ',f15.6)
  146 format (/,' total energy (Lagrange)         ',f15.6,/
     1          '              (Non-Lagrange)     ',f15.6)
  155 format   ('                at t-  dt*nprint ',f15.6)
  156 format   ('                at t-2*dt*nprint ',f15.6)
  150 format (/,' ',78('*'),/,' ',78('*'),//,
     1          '  program stops because of positive total energy',//,
     2          ' ',78('*'),/,' ',78('*'))
  151 format (/,' ',78('*'),/,' ',78('*'),//,
     1          '  program stops because the total energy diverges',//,
     2          ' ',78('*'),/,' ',78('*'))

c................................................... calculate various energies
      call fprte(itprii)
      call fprtj

c     ..................... single-particle matrix elements of angular momentum
      call fprtz

c     ..................................................... some initialization
      exmax = 0.0d0
      xaf   = npn + npp

      if (iprintlag.ne.1) then
        print 100,itert
      else

      endif

c     ............................................. sum up average pairing gaps
      v2sum  (:) = zero
      uvsum  (:) = zero
      v2Delta(:) = zero
      uvDelta(:) = zero
      if (iprintlag.ne.1) then
        do it=1,2
          if (it.eq.1) then
            na = 1
            nb = nwaven
          else
            na = nwaven + 1
            nb = nwaven + nwavep
          endif
          do iwave=na,nb
            u2v2        = (1.0d0 - v2(iwave))*v2(iwave)
            uv          = sqrt(u2v2)
            v2sum  (it) = v2sum(it) + 2.0d0*v2(iwave)
            uvsum  (it) = uvsum(it) + 2.0d0*uv
            v2Delta(it) = v2Delta(it) + 2.0d0*delta(iwave)*v2(iwave)
            uvDelta(it) = uvDelta(it) + 2.0d0*delta(iwave)*uv
          enddo
          if (uvsum(it).lt.epsm4) then
            uvsum  (it) = 0.0d0
            uvDelta(it) = 0.0d0
            v2Delta(it) = 0.0d0
          else
            uvDelta(it) = uvDelta(it)/uvsum(it)
            v2Delta(it) = v2Delta(it)/v2sum(it)
          endif
        enddo
      endif

c     ........................................... print single-particle spectra
      if (iprintlag.ne.1) then
        do it=1,2
          if (it.eq.1) then
            na = 1
            nb = nwaven
            print 101
          else
            na = nwaven + 1
            nb = nwaven + nwavep
            print 102
          endif
          print 103
          etsp(it) = zero
          do iw=na,nb
            idh(iw) = 0
          enddo
          do iw=na,nb
            iww = 2*(iw-na+1)
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
            iwave    = kw
            etsp(it) = etsp(it) + two*esp1(iwave)*v2(iwave)
            esp4     = esp2(iwave) - esp1(iwave)**2
            esp4     = abs(esp4)
            if (esp4.lt.esp3(iwave)) esp4 = -esp4
            if ((itprii.ge.1) .or.
     1          (abs(esp1(iwave)-ambda(it)).le.esplm)) then
               print 104,iwave,iww,kparz(iwave),v2(iwave),v22(iwave)
     1                  ,delta(iwave),eqp(iwave)
     2                  ,esp1(iwave),ajzd(iwave),aj2d(iwave),esp4
            endif
            esp4        = abs(esp4)
            esp3(iwave) = esp4
            if (esp4.gt.exmax) exmax = esp4
          enddo
        enddo
      endif

c     .............................. calculate and print moments of the density
      if (iprintlag.ne.1) then
        print 105
        call fprtm (itprii,iprintlag,ifinal)
      else
        call fprtm (0     ,iprintlag, 0)
      endif

c     ................................... calculate mean value of j(x,y,z)(**2)
      if (iprintlag.ne.1) then
        print 105
        call pipj
      endif

c     ................................... calculate mean value of j(x,y,z)(**2)
      if (iprintlag.ne.1) then
        print 105
        call jdeux
      endif

c     ...................................................... pairing properties
      if (iprintlag.ne.1 .and. npair.ge.1) then
        print 105
        print 110,ambda(1:2),xlamb  (1:2),disper (1:2)
        print 111,uvsum(1:2),uvDelta(1:2),v2Delta(1:2)
      endif

c     ................................................................ energies
      epair(3) = zero
      if (npair.ne.0) epair(3) = epair(1) + epair(2)
      ek(3)   = ek  (1) + ek  (2)
      es      = eso (1) + eso (2)
      etsp(3) = etsp(1) + etsp(2)
      etot    = esk  + ecoul + ecoex + ek(3)
      etot    = etot + epair(3)
      if (ncm2.eq.-2) etot = etot + ecm(1,3) + ecm(2,3) 
      if (ncm2.eq. 1) etot = etot            + ecm(2,3)
      if (iprintlag.ne.1) then
        print 105
        print 120
      endif

c     ................................ this now separates kinetic and cm energy
c            the kinetic energy in ek(:) might contain the 1_body cm correction
      if (ncm2.ge.0) then
        print 121,ek(1)-ecm(1,1),ek(2)-ecm(1,2),ek(3)-ecm(1,3)
      else
        print 121,ek(1:3)
      endif
      print 1221,ecm(1,1:3)
      print 122 ,ecm(2,1:3)
      if (njmunu.ge.1) then
        if (nfunc.eq.0) print  124,ejmunu(1),ejmunu(2),ejmunu(3)
        if (nfunc.eq.1) print 1124,ejmunu(1),ejmunu(3)
      endif
      if (iprintlag.ne.1) print 125,etsp
      if (npair.ne.0) print 126,epair
      if (npair.eq.3 .or. npair.eq.5 ) print 127, -eproj
      print 128,ecoul,ecoex

c     ............................ energy from summing single-particle energies
      eto2 = etsp(3) + ek(3) - byt3a*et3a + epente
      if (ndd .eq .1) eto2 = eto2 - byt3b*et3b
      eto2 = eto2/two + ecoex/tt3 + epair(3)
      if (ncm2.eq.-2) eto2 = eto2 + ecm(1,3) + ecm(2,3) 
      if (ncm2.eq. 1) eto2 = eto2 + ecmp (3) ! pairing part absent in spes

      et3  = et3a + et3b
      print 143,et0,ett,et3a,et3b,etd,es,esk

      eprojt = etot - (eproj(1)+eproj(2))

      if (iprintlag.ne.1) then
        if (npair.eq.3 .or. npair.eq.5) then
          print 144,etot
          print 145,eprojt,eto2 - (eproj(1) + eproj(2))
        else
          print 145,etot,eto2
        endif
      else
        if (npair.eq.3 .or. npair.eq.5) then
          print 144,etot
          print 146,eprojt,econv2
        else
          print 146,etot,econv2
        endif
      endif

      if (iprintlag.ne.1) then
        iconv = iconv + 1
        if (etot.ge.zero) then
          print 150
          print 105
          iconv = -1
          return
        endif
        if (iconv.le.1) then
          econv1 = etot
          if (npair.eq.3 .or. npair.eq.5)
     1      econv1 = econv1 - eproj(3)
          print 105
          return
        endif
        if (iconv.le.2) then
          econv2 = etot
          if (npair.eq.3 .or. npair.eq.5)
     1      econv2 = econv2 - eproj(3)
          print 155,econv1
          print 105
          return
        endif
        print 155,econv2
        print 156,econv1

        if(npair.eq.3 .or. npair.eq.5) then
          es = abs(etot -eproj(3) - econv2)
        else
          es = abs(etot-econv2)
        endif
        if (es.lt.t50) go to 12
        es = es/(abs(econv2-econv1)+epsm8)
        if (es.lt.one) go to 12
        print 151
        print 105
        iconv = -1
        return
   12   econv1 = econv2
        econv2 = etot
        if(npair.eq.3 .or. npair.eq.5) then
          econv2 = econv2 - eproj(3)
        endif
      endif

      print 105

      return
      end subroutine figaro
c______________________________________________________________________________
      subroutine fprte (ipri)

c..............................................................................
c     calculation of                                                          .
c     - the Skyrme rho rho           energy                                   .
c     - the Skyrme rho rho rho^alpha energy                                   .
c     - the Skyrme rho Delta rho     energy                                   .
c     - the Skyrme J J               energy                                   .
c     - the Coulomb                  energy                                   .
c     - the two-body part of the centre-of-mass correction                    .
c..............................................................................
c     the routine fprtj calculates the remaining energies                     .
c..............................................................................
      include 'param8.h'
      implicit real*8 (a-h,o-z)
      character*20 afor

      parameter (zero=0.0d0,two=2.0d0,tp5=0.5d0,tp75=0.75d0)

      common /den  / rhon(mv),rhop(mv)
      common /drho / drhon(mv),drhop(mv)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /enert/ c14,c15,t14,t15
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu(3),ecm(2,3),ecmp(3)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /wcm2 / tx(mw,mw),ty(mw,mw),tz(mw,mw)
      common /wj2  / vjxx(mv,2),vjyx(mv,2),vjzx(mv,2)
     1              ,vjxy(mv,2),vjyy(mv,2),vjzy(mv,2)
     2              ,vjxz(mv,2),vjyz(mv,2),vjzz(mv,2)

c..............................................................................
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
        ytt   = ytt   + xtt
        ynp   = ynp   + xnp
        ydtt  = ydtt  + xt*(dxn+dxp)
        ydnp  = ydnp  + (xn*dxn+xp*dxp)
        y3att = y3att + xtt*wt3a(i)
        y3anp = y3anp + xnp*wt3a(i)
        y3btt = y3btt + xtt*wt3b(i)
        y3bnp = y3bnp + xnp*wt3b(i)
        ycd   = ycd   + xp*wcd(i)
        yce   = yce   + xp*wce(i)
      enddo
      et0   = dv * (b1*ytt  + b2*ynp)
c             .............................................ATTENTION
c             B5 & B6 as used in the code are not equal to the ones
c             discussed in all of the articles; they carry a 
c             relative minus sign in this code!
c             ...................................................... 
      etd   =-dv * (b5*ydtt + b6*ydnp)
      et3a  = dv * (b7a*y3att + b8a*y3anp)
      et3b  = dv * (b7b*y3btt + b8b*y3bnp)

c     ................................................................. Coulomb
      ecoul = dv * ycd * tp5
      ecoex = 0.0d0
      if (ncoex.eq.0) ecoex = dv * yce * tp75

c     .......................... two-body part of the centre-of-mass correction
c                                the one-body part if contained in the kinetic
c                                energy ek(:) calculated in subroutine fprtj

      if (ncm2.ne.-2) then       ! don't reset ecm to zero for ncm = -2
        do it=1,3                ! to save cpu time, the value from the 
          ecm(2,it) = zero       ! last intermediate printout is used
          ecmp (it) = zero
        enddo
      endif

      if (ncm2.eq.1 .or. (ncm2.eq.-2.and.ipri.gt.0)) then
        if (ipri.gt.0) call nabxyz
        do it=1,2
          ecm(2,it) = zero
          ecmp (it) = zero
          i1 = 1      + nwaven*(it-1)
          i2 = nwaven + nwavep*(it-1)
          do ia = i1,i2
          do ib = i1,i2
            if (kparz(ia).ne.kparz(ib)) then
              u2ia = 1.0d0 - v2(ia)
              u2ib = 1.0d0 - v2(ib)
              v2v2 = v2(ia) * v2(ib)
              uvuv = sqrt(v2(ia) * u2ia * v2(ib) * u2ib)
              txyz = tx(ia,ib)**2 + ty(ia,ib)**2 + tz(ia,ib)**2
              ecm(2,it) = ecm(2,it) + v2v2 * txyz
              ecmp (it) = ecmp (it) + uvuv * txyz
            endif
          enddo
          enddo
          ecm(2,it) = hbm(it)*(xm(it)/xm(3))*ecm(2,it)
          ecmp (it) = hbm(it)*(xm(it)/xm(3))*ecmp (it)
        enddo
        ecm(2,3) = ecm(2,1) + ecm(2,2)
        ecmp (3) = ecmp (1) + ecmp (2)
      endif

      ecm(2,:) = ecm(2,:) + ecmp (:)

c     ................................. the energy from spin-orbit tensor terms
c                              for non-functionals, the energy is splitted into
c                                    ejmunu(1)         central force, J_mn J_mn
c                                    ejmunu(2)         tensor  force, J_mn J_mn
c                                    ejmunu(3)         tensor  force, J_mn J_nm
c                              for functionals, the energy is splitted into
c                                    ejmunu(1)                        J_mn J_mn
c                                    ejmunu(2)                        zero
c                                    ejmunu(3)                        J_mn J_nm
      ejmunu(:) = zero

c     ..... the symmetric combination from the central force and a tensor force
c     the contributions from central and tensor forces to the total binding
c     energy are separated via the coupling constants c14/15 and t14/15. In
c     case of a functional (nfunc.eq.1) this separation makes no sense anymore.
      x = 0.0d0
      y = 0.0d0
      if (njmunu.ge.1) then
        x = sdot(18*mv,vjxx(1,1),vjxx(1,1))
        y =    sdot(mv,vjxx(1,1),vjxx(1,2))
     2       + sdot(mv,vjyx(1,1),vjyx(1,2))
     3       + sdot(mv,vjzx(1,1),vjzx(1,2))
     4       + sdot(mv,vjxy(1,1),vjxy(1,2))
     5       + sdot(mv,vjyy(1,1),vjyy(1,2))
     6       + sdot(mv,vjzy(1,1),vjzy(1,2))
     7       + sdot(mv,vjxz(1,1),vjxz(1,2))
     8       + sdot(mv,vjyz(1,1),vjyz(1,2))
     9       + sdot(mv,vjzz(1,1),vjzz(1,2))
        x = x * dv
        y = y * dv
        if (nfunc.eq.0) ejmunu(1) = (c14+c15)*x + two*c14*y
        if (nfunc.eq.1) ejmunu(1) = (b14+b15)*x + two*b14*y
      endif

c     ................................ asymmetric combination from tensor force
c           the symmetric combination is splitted into a central part ejmunu(1)
c                                                   and a tensor part ejmunu(2)
c                 the asymmetric combination from the tensor force is ejmunu(3)
c                the first line uses "x" and "y" from the previous if statement
      if (njmunu.eq.2) then
        if (nfunc.eq.0) ejmunu(2) = (t14+t15)*x + two*t14*y
        xa =         sdot(mv,vjxx(1,1),vjxx(1,1))
     1       +       sdot(mv,vjyy(1,1),vjyy(1,1))
     2       +       sdot(mv,vjzz(1,1),vjzz(1,1))
     3       + two * sdot(mv,vjxy(1,1),vjyx(1,1))
     4       + two * sdot(mv,vjxz(1,1),vjzx(1,1))
     5       + two * sdot(mv,vjyz(1,1),vjzy(1,1))
     6       +       sdot(mv,vjxx(1,2),vjxx(1,2))
     7       +       sdot(mv,vjyy(1,2),vjyy(1,2))
     8       +       sdot(mv,vjzz(1,2),vjzz(1,2))
     9       + two * sdot(mv,vjxy(1,2),vjyx(1,2))
     8       + two * sdot(mv,vjxz(1,2),vjzx(1,2))
     7       + two * sdot(mv,vjyz(1,2),vjzy(1,2))
        ya =         sdot(mv,vjxx(1,1),vjxx(1,2))
     1       +       sdot(mv,vjxy(1,1),vjyx(1,2))
     2       +       sdot(mv,vjxz(1,1),vjzx(1,2))
     3       +       sdot(mv,vjyx(1,1),vjxy(1,2))
     4       +       sdot(mv,vjyy(1,1),vjyy(1,2))
     5       +       sdot(mv,vjyz(1,1),vjzy(1,2))
     6       +       sdot(mv,vjzx(1,1),vjxz(1,2))
     7       +       sdot(mv,vjzy(1,1),vjyz(1,2))
     8       +       sdot(mv,vjzz(1,1),vjzz(1,2))
        xa = xa * dv
        ya = ya * dv
        ejmunu(3) = (b16+b17)*xa + two*b16*ya
      endif

      return
      end subroutine fprte

c______________________________________________________________________________
      subroutine fprtj

c..............................................................................
c     calculation of                                                          .
c     - kinetic energy (+ cm. correction)                                     .
c     - the Skyrme rho tau energy                                             .
c     - the skyrme spin-orbit energy                                          .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)

      common /lag  / ilag, iAnalag
      common /den  / rho(mv,2)
      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu(3),ecm(2,3),ecmp(3)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pot  / wnn(mv),wpp(mv),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /stor / a(mq,mw)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /taudj/ vtau(mv,2),vdiv(mv,2)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv)
     1              ,p1(mv),p2(mv),p3(mv),p4(mv)

c..............................................................................
      ek(1)  = zero
      ek(2)  = zero
      do iwa=1,nwave
        it = kiso(iwa)
        zp = kparz(iwa)
        oc = v2(iwa)
        call scopy (mq,a(1,iwa),w1)
        if(ilag.eq.0) then
          call lapla (w1,p1, one, one, zp)
          call lapla (w2,p2,-one,-one, zp)
          call lapla (w3,p3,-one, one,-zp)
          call lapla (w4,p4, one,-one,-zp)
        else
          call laplalag (w1,p1, 1, 1, kparz(iwa))
          call laplalag (w2,p2,-1,-1, kparz(iwa))
          call laplalag (w3,p3,-1, 1,-kparz(iwa))
          call laplalag (w4,p4, 1,-1,-kparz(iwa))
        endif

        ek(it) = ek(it) + oc * sdot(mq,w1,p1)
      enddo

c     .................................................. one-body cm correction
c               attention: this is to be added to the energy only for ncm2 = -2
c               for all other values of ncm2 this has only diagnostic purposes
      ecm(1,:) = 0.0d0
      if (ncm2.ne.-1) then
        ecm(1,1) = dv * 0.5d0 * hbm(1) * ek(1) * xm(1)/xm(3)
        ecm(1,2) = dv * 0.5d0 * hbm(2) * ek(2) * xm(2)/xm(3)
        ecm(1,3) = ecm(1,1) + ecm(1,2)
      endif

c     ............................ kinetic energy with or without cm correction
      if (ncm2.ge.0) then
        ek(1) = -dv * 0.5d0 * hbm(1) * ek(1) * (1.0d0 - xm(1)/xm(3))
        ek(2) = -dv * 0.5d0 * hbm(2) * ek(2) * (1.0d0 - xm(2)/xm(3))
      endif
      if (ncm2.lt.0) then
        ek(1) = -dv * 0.5d0 * hbm(1) * ek(1)
        ek(2) = -dv * 0.5d0 * hbm(2) * ek(2)
      endif

c     .................................................... some skyrme energies
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
      esk    = et0 + ett + etd + et3a + et3b + sum(eso) + sum(ejmunu)

      return
      end subroutine fprtj

c______________________________________________________________________________
c______________________________________________________________________________
      subroutine fprtz

c..............................................................................
c     compute expectation value of <jz> for each single-particle wave function.
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,two=2.0d0,tt4=4.0d0)

      common /mud  / xi(mv),yi(mv),zi(mv),xii(mv),yii(mv),zii(mv)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv),
     1               p1(mv),p2(mv),p3(mv),p4(mv)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /wtmp / wjz1(mv),wjz2(mv),wjz3(mv),wjz4(mv)

c..............................................................................
      do iwa=1,nwave
        call scopy (  mq, a(1,iwa),w1 (1))
        call scopy (3*mq,da(1,iwa),wx1(1))

c       .................................................................. <jz>
        ajzd(iwa) = zero
        do i=1,mv
          wjz1(i) = w1(i) + two*(+xi(i)*wy2(i)-yi(i)*wx2(i))
          wjz2(i) = w2(i) + two*(-xi(i)*wy1(i)+yi(i)*wx1(i))
          wjz3(i) =-w3(i) + two*(+xi(i)*wy4(i)-yi(i)*wx4(i))
          wjz4(i) =-w4(i) + two*(-xi(i)*wy3(i)+yi(i)*wx3(i))
        enddo
        ajzd(iwa) = dv*sdot(mq,w1(1),wjz1(1))/tt4

c       .................................................................. <j2>
c       the factor 1/8 comes from 1) |psi) being normalized to 2 (psi|psi) = 2
c                                 2) wjz1 etc being 2 * jz |psi)
        xjz2 = dv*sdot(mq,wjz1(1),wjz1(1)) / 8.0d0                      ! <jz2>

        do i=1,mv
          p1(i) =  w4(i) + two*(+zi(i)*wx2(i) - xi(i)*wz2(i))
          p2(i) = -w3(i) + two*(-zi(i)*wx1(i) + xi(i)*wz1(i))
          p3(i) = -w2(i) + two*(+zi(i)*wx4(i) - xi(i)*wz4(i))
          p4(i) =  w1(i) + two*(-zi(i)*wx3(i) + xi(i)*wz3(i))
        enddo
        xjy2 = dv*sdot(mq,p1(1),p1(1)) / 8.0d0                          ! <jy2>

        do i=1,mv
          p1(i) = w3(i) + two*(+yi(i)*wz2(i) - zi(i)*wy2(i))
          p2(i) = w4(i) + two*(-yi(i)*wz1(i) + zi(i)*wy1(i))
          p3(i) = w1(i) + two*(+yi(i)*wz4(i) - zi(i)*wy4(i))
          p4(i) = w2(i) + two*(-yi(i)*wz3(i) + zi(i)*wy3(i))
        enddo
        xjx2 = dv*sdot(mq,p1(1),p1(1)) / 8.0d0                         ! <jx2>

        xj2  = xjx2 + xjy2 + xjz2                ! <j2> = <jx2> + <jy2> + <jz2>
        xj   = -0.5d0 + sqrt(0.25d0 + xj2)       !         j from <j2> = j(j+1)
        aj2d(iwa) = xj
      enddo

      return
      end subroutine fprtz
      
c______________________________________________________________________________
      subroutine fprtm (itprii,iprintlag,ifinal)

c..............................................................................
c     calculate and print                                                     .
c     - multipole moments of the density                                      .
c     - information on constraints                                            .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0)
      parameter (tt4=4.0d0,tt5=5.0d0,tt6=6.0d0,t24=24.0d0)
      parameter (tt8=8.0d0,t180=180.0d0)
      parameter (epsm3=1.0d-3,epsm4=1.0d-4)

      parameter (tt7=7.0d0)
      parameter (tt9=9.0d0)
      parameter (t60=60.0d0,t36=36.0d0,t77=77.0d0)
      parameter (t300=300.0d0,t729=729.0d0,t1001=1001.0d0)
      parameter (clum  =  29.9792458d0)  !Speed of light added 

      common /champ/ qxxn,qyyn,qzzn,qrrn
     1              ,qxxp,qyyp,qzzp,qrrp
     2              ,qxxt,qyyt,qzzt,qrrt
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv)
     1               ,imtd,imtg,icutq
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
      common /epente/ epente
      common /cstw / delq,q1n,q1p,q1t,q2n,q2p,q2t
      common /den  / rhon(mv),rhop(mv)
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /kfmom/ tfac,tfac1,tfac2,s3
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /mmt  / ap,an,at
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /mud  / xi (mx,my,mz),yi (mx,my,mz),zi (mx,my,mz)
     1              ,xii(mx,my,mz),yii(mx,my,mz),zii(mx,my,mz)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /lnog / rln(3),qxln(3),qyln(3),qzln(3)
      dimension qln0(3),gln0(3)

c..............................................................................
   99 format (/,' ',78('_'),/, ' Multipole constraints')
  100 format (  ' Moments of the density',/,
     1        /,25x,'neutron',7x,'proton',8x,'total')
  101 format (  ' particle number  : ',3(f12.6,1x))
  102 format (/,' Radii',
     1        /,' ms  radius (fm2) : ',3(f12.4,1x),
     2        /,' rms radius (fm)  : ',3(f12.4,1x),
     3        /,' skin       (fm)  : ',26x,1f12.4)
!  104 format (/,' Radius Constraint'             )
!  105 format (/,'                    ', 5x, 'N', 5x, 'P',/,
!     1          ' Mono. Constraint :',2(f12.4,1x))    
!  106 format (/,' Mono. Constraint :', (f12.4,1x))  
  103 format (/,' Lipkin-Nogami corrected radii: ',
     1        /,' ms  radius (fm2) : ',3(f12.4,1x),
     2        /,' rms radius (fm)  : ',3(f12.4,1x),
     3        /,' skin       (fm)  : ',26x,1f12.4)
  110 format (/,' Cartesian multipole moments',
     1        /,' Qx         (fm2) : ',3(f12.3,1x),
     2        /,' Qy         (fm2) : ',3(f12.3,1x),
     3        /,' Qz         (fm2) : ',3(f12.3,1x),
     4        /,' Q0         (fm2) : ',3(f12.3,1x),
     5        /,' gamma      (deg) : ',3(f12.3,1x))
  114 format (/,' Quadrupole moment as used in constraints' 
     1         ,'(with cutoff)',
     1        /,' q1=iq1*delq(fm2) : ',3(f12.3,1x),
     2        /,' q2=iq1*delq(fm2) : ',3(f12.3,1x),
     3        /,' q0         (fm2) : ',3(f12.3,1x),
     4        /,' gamma      (deg) : ',3(f12.3,1x))
  115 format (/,' Lipkin-Nogami corrected cartesian multipole moments',
     1        /,' Qx         (fm2) : ',3(f12.3,1x),
     2        /,' Qy         (fm2) : ',3(f12.3,1x),
     3        /,' Qz         (fm2) : ',3(f12.3,1x),
     4        /,' Q0         (fm2) : ',3(f12.3,1x),
     5        /,' gamma      (deg) : ',3(f12.3,1x))
  121 format (/,' Equivalent Linear Constraint Intensity',
     1          ' (MeV/fm and MeV/fm^2)',/,
     1           8x,'r',8x,'Qxx',7x,'Qyy',7x,'Qzz')
  122 format ( /, ' Readjusted Constraints (fm and fm^2)'
     1        ,/,10x,'r',10x,'Qxx',9x,'Qyy',
     1           9x,'Qzz',/,' T',4f12.3)
 1122 format ( ' Readjusted Constraint:',/,4x,f8.4,'Fm',
     1          2( 6x, f8.4,'Fm'),/)
  123 format (/,' Equivalent Linear Constraint Intensity',
     1          ' (MeV/fm and MeV/(fm^2)) ', 
     1         /,10x,'r',10x,'Qxx',9x,'Qyy',9x,'Qzz'
     1        ,/,' N',4f12.3,
     1         /,' P',4f12.3)
  124 format ( /, ' Readjusted Constraints (fm and fm^2)'
     1        ,/,10x,'r',10x,'Qxx',9x,'Qyy',
     1           9x,'Qzz',/,' N',4f12.3,/,' P',4f12.3)
  150 format (/,' Con.Ener.,x,y,z',4f12.3)
  151 format (  ' slopes  r,x,y,z',4e13.5,' (MeV/fm2)')
  152 format (  ' cstrn eqr,x,y,z',4f12.3,
     1        /,' cstrp eqr,x,y,z',4f12.3)
  153 format (/,' slopesn r,x,y,z',4e12.3,' (MeV/fm2)',
     1        /,' slopesp r,x,y,z',4e12.3,' (MeV/fm2)')
  154 format (/,' Constraint Energies (MeV)',/,
     1          10x,'r',10x,'Qxx',9x,'Qyy',9x,'Qzz')   
  155 format (' N' , 4f12.3)
  156 format (' P' , 4f12.3)
  157 format (' T' , 4f12.3)
  158 format (/,' Eqv. Linear Constraint (MeV/fm and MeV/(fm^2))  ',/,
     1          10x,'r',10x,'Qxx',9x,'Qyy',9x,'Qzz')       

c..............................................................................
      pi   = tt4*atan2(one,one)
      traf = pi/t180

      ccc =  sqrt(  tt5*pi)/tt3
      ddd =  sqrt(      pi)
      c22 =  tt6/(  tt7*pi)
      c44 =  t60/(  t77*pi)
      c24 =  t36/(  tt7*pi*sqrt(tt5))
      d22 =  tt9/(  tt7*pi)
      d44 = t729/( t1001*pi)
      d24 = t300/(   t77*pi*sqrt(tt5))

c     ........................................................ particle numbers
      at = ap + an
      if (iprintlag.ne.1) then
        print 100
        print 101,an,ap,at
      endif

c     ................................................................... radii
      rp = x2p+y2p+z2p
      rn = x2n+y2n+z2n
      rr = (tt5/tt3)*(rn+rp)/at
      rt = sqrt((rn+rp)/at)
      rn = sqrt(rn/an)
      rp = sqrt(rp/ap)
      if (iprintlag.ne.1) then
        print 102,rn*rn,rp*rp,rt*rt,
     1            rn,rp,rt,
     2            rn-rp
        if (npair.eq.3 .or. npair.eq.5) then
          print 103,rln(1)*rln(1),rln(2)*rln(2),rln(3)*rln(3),
     1              rln(1),rln(2),rln(3),
     2              rln(1)-rln(2)
        endif
      endif

c     ............................................ cartesian quadrupole moment
      qxp = two*x2p - y2p - z2p
      qyp = two*y2p - z2p - x2p
      qzp = two*z2p - x2p - y2p
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
      if(iprintlag.ne.1) then
        print 110,qxn,qxp,qxt,
     1            qyn,qyp,qyt,
     2            qzn,qzp,qzt,
     3            qn0,qp0,qt0,
     4            gn0,gp0,gt0
      endif

c     .......................... cartesian quadrupole moment with LN correction
      if ( (npair.eq.3 .or. npair.eq.5) .and. iprintlag.ne.1 ) then
        do itt=1,3
          fln = (two/tt3)
          qln0(itt) = sqrt(fln*(qxln(itt)**2+qyln(itt)**2+qzln(itt)**2))
          if ((abs(qxln(itt)-qyln(itt))+abs(qzln(itt))).le.epsm4) then
            gln0(itt) = zero
          else
            gln0(itt) = atan2(qxln(itt)-qyln(itt),qzln(itt)*s3)/traf
          endif
        enddo
        print 115,qxln(1),qxln(2),qxln(3),
     1            qyln(1),qyln(2),qyln(3),
     2            qzln(1),qzln(2),qzln(3),
     3            qln0(1),qln0(2),qln0(3),
     4            gln0(1),gln0(2),gln0(3)
      endif

c     ................................................... quadrupole constraint
      if (iprintlag.ne.1) then
        qq1n = (two/tt3) * (qznc-qxnc)
        qq2n = (two/tt3) * (qxnc-qync)
        qq1p = (two/tt3) * (qzpc-qxpc)
        qq2p = (two/tt3) * (qxpc-qypc)
        qq1t = (two/tt3) * (qztc-qxtc)
        qq2t = (two/tt3) * (qxtc-qytc)
        print 114,qq1n,qq1p,qq1t,
     1            qq2n,qq2p,qq2t,
     2            qnc0,qpc0,qtc0,
     3            gnc0,gpc0,gtc0
      endif

c     ............................................. spherical multipole moments
      if ( iprintlag.ne.1 .and. itprii.eq.1) then
        call momehigh
      endif

c     .......................................... diagnostics of the constraints
                       ! MB: NEEDS FURTHER WORK:
                       ! the printed slope contains also the radius constraint,
                       ! not just the quadrupole. To be kept like this or not?
      if (iprintlag.ne.1) then
        if (imtd.le.1) then
          if (imtg.eq.0) then
            f  = excstt/two-cq2*qxcstt*(qxtc-qxcstt)
     1          +eycstt/two-cq2*qycstt*(qytc-qycstt)
     2          +ezcstt/two-cq2*qzcstt*(qztc-qzcstt)
     3          +ercstt/two-cqr*qrcstt*(qrtc-qrcstt)
          else
            f  = excstt/two-cq2*q2cst*(qtc0-q2cst)
          endif
          if (imtg.eq.0) then
            epente= -(pentext*qxtc+penteyt*qytc+pentezt*qztc)
     1              -(pentert*qrtc)
          else
            epente= -pentext*qtc0
          endif
          print 154
          print 157,ercstt,excstt,eycstt,ezcstt
          if(ifinal.eq.1) then
            print 158
            print 157,pentert,pentext,penteyt,pentezt
          endif
        else
          cq2n = cq2*(at/an)
          cq2p = cq2*(at/ap)
          cqrn = cqr*(at/an)
          cqrp = cqr*(at/ap)
          f  = excstn/two - cq2n*qxcstn*(qxnc-qxcstn)
     1        +eycstn/two - cq2n*qycstn*(qync-qycstn)
     2        +ezcstn/two - cq2n*qzcstn*(qznc-qzcstn)
     3        +ercstn/two - cqrn*qrcstn*(qrnc-qrcstn)
     4        +excstp/two - cq2p*qxcstp*(qxpc-qxcstp)
     1        +eycstp/two - cq2p*qycstp*(qypc-qycstp)
     2        +ezcstp/two - cq2p*qzcstp*(qzpc-qzcstp)
     3        +ercstp/two - cqrp*qrcstp*(qrpc-qrcstp)
          epente= -(pentexn*qxnc+penteyn*qync+pentezn*qznc)
     1            -(pentexp*qxpc+penteyp*qypc+pentezp*qzpc)
     2            -(pentern*qrnc+penterp*qrpc)
          print 154
          print 155,ercstn,excstn,eycstn,ezcstn
          print 156,ercstp,excstp,eycstp,ezcstp
					if(ifinal.eq.1) then
            print 158
            print 155,pentern,pentexn,penteyn,pentezn
            print 156,penterp,pentexp,penteyp,pentezp
          endif
        endif
      endif

      if(iprintlag.ne.1 .and. ifinal.ne.1) then
        if (imtd.le.1) then
          if (cq2 .ne. zero ) then
            print 122,sqrt(qrcstt/(npn+npp)),qxcstt,qycstt,qzcstt
          endif
        else
          print 124,sqrt(qrcstn/(npn)),qxcstn,qycstn,qzcstn
     1             ,sqrt(qrcstp/(npp)),qxcstp,qycstp,qzcstp
        endif        
      endif

      return
      end subroutine fprtm
c______________________________________________________________________________      
      subroutine momehigh

c..............................................................................
c     calculate multipole moments of the densities                            .
c..............................................................................
c     adapted from standard.1.3.4 by WR                                       .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (epsm3=1.0d-3)

      common /nxyz  / dx,dv
      common /mud   / x1(mv),y1(mv),z1(mv),
     1                x2(mv),y2(mv),z2(mv)
      common /den   / rhot(mv,2)
      dimension anu(2),r2(2)
      dimension q20   (2),q22r  (2)
      dimension q40   (2),q42r  (2),q44r  (2)
      dimension q60   (2),q62r  (2),q64r  (2),q66r  (2)
      dimension q80   (2),q82r  (2),q84r  (2),q86r  (2),q88r(2)
      dimension q1000 (2),q1002r(2),q1004r(2),q1006r(2),
     1          q1008r(2),q1010r(2)
      dimension beta20(3),beta22(3)
      dimension beta40(3),beta42(3),beta44(3)
      dimension beta60(3),beta62(3),beta64(3),beta66(3)
      dimension beta80(3),beta82(3),beta84(3),beta86(3),beta88(3)
      dimension beta1000(3),beta1002(3),beta1004(3),
     1          beta1006(3),beta1008(3),beta1010(3)

c..............................................................................
  100 format (/,' ',78('_'))
  101 format (/,' Spherical multipole moments')
  102 format (/,25x,'neutron',7x,'proton',8x,'total')
  105 format (  ' Q20        (fm2) : ',3(f12.3,1x),
     1        /,' Q22        (fm2) : ',3(f12.3,1x))
  115 format (/,' Q40        (fm4) : ',3(f12.3,1x),
     1        /,' Q42        (fm4) : ',3(f12.3,1x),
     2        /,' Q44        (fm4) : ',3(f12.3,1x))
  125 format (/,' Q60        (fm6) : ',3(f12.3,1x),
     1        /,' Q62        (fm6) : ',3(f12.3,1x),
     2        /,' Q64        (fm6) : ',3(f12.3,1x),
     3        /,' Q66        (fm6) : ',3(f12.3,1x))
  106 format (/,' Deformation parameters',
     1        /,' beta20           : ',3(f12.4,1x),
     2        /,' beta22           : ',3(f12.4,1x))
  116 format (/,' beta40           : ',3(f12.4,1x),
     1        /,' beta42           : ',3(f12.4,1x),
     2        /,' beta44           : ',3(f12.4,1x))
  126 format (/,' beta60           : ',3(f12.4,1x),
     1        /,' beta62           : ',3(f12.4,1x),
     2        /,' beta64           : ',3(f12.4,1x),
     3        /,' beta66           : ',3(f12.4,1x))
  136 format (/,' beta80           : ',3(f12.4,1x),
     1        /,' beta82           : ',3(f12.4,1x),
     2        /,' beta84           : ',3(f12.4,1x),
     3        /,' beta86           : ',3(f12.4,1x),
     4        /,' beta88           : ',3(f12.4,1x))
  146 format (/,' beta1000         : ',3(f12.4,1x),
     1        /,' beta1002         : ',3(f12.4,1x),
     2        /,' beta1004         : ',3(f12.4,1x),
     3        /,' beta1006         : ',3(f12.4,1x),
     4        /,' beta1008         : ',3(f12.4,1x),
     5        /,' beta1010         : ',3(f12.4,1x))

c..............................................................................
c     print 100
      print 101

c     ................................................ some useful coefficients
      pi      = 4.0d0*atan2(1.0d0,1.0d0)
      cof20   = -sqrt(     5.0d0/pi)         /    4.0d0
      cof22   =  sqrt(    30.0d0/pi)         /    8.0d0
      cof40   =  sqrt(     1.0d0/pi) * 3.0d0 /   16.0d0
      cof42   = -sqrt(    10.0d0/pi) * 3.0d0 /   16.0d0
      cof44   =  sqrt(    70.0d0/pi) * 3.0d0 /   32.0d0
      cof60   = -sqrt(    13.0d0/pi)         /   32.0d0
      cof62   =  sqrt(  1365.0d0/pi)         /   64.0d0
      cof64   = -sqrt(   182.0d0/pi) * 3.0d0 /   64.0d0
      cof66   =  sqrt(  3003.0d0/pi)         /   64.0d0
      cof80   =  sqrt(    17.0d0/pi)         /  256.0d0
      cof82   = -sqrt(   595.0d0/pi) * 3.0d0 /  128.0d0
      cof84   =  sqrt(  2618.0d0/pi) * 3.0d0 /  256.0d0
      cof86   = -sqrt(  7293.0d0/pi)         /  128.0d0
      cof88   =  sqrt( 24310.0d0/pi) * 3.0d0 /  512.0d0
      cof1000 = -sqrt(    21.0d0/pi)         /  512.0d0
      cof1002 =  sqrt(   770.0d0/pi) * 3.0d0 / 1024.0d0
      cof1004 = -sqrt( 10010.0d0/pi) * 3.0d0 /  512.0d0
      cof1006 =  sqrt(  5005.0d0/pi) * 3.0d0 / 1024.0d0
      cof1008 = -sqrt(510510.0d0/pi)         / 1024.0d0
      cof1010 =  sqrt(969969.0d0/pi)         / 1024.0d0

c     ............................................. calculate multipole moments
      anu (:) = 0.0d0
      r2  (:) = 0.0d0
      q20 (:) = 0.0d0
      q22r(:) = 0.0d0
      q40 (:) = 0.0d0
      q42r(:) = 0.0d0
      q44r(:) = 0.0d0
      q60 (:) = 0.0d0
      q62r(:) = 0.0d0
      q64r(:) = 0.0d0
      q66r(:) = 0.0d0
      q80 (:) = 0.0d0
      q82r(:) = 0.0d0
      q84r(:) = 0.0d0
      q86r(:) = 0.0d0
      q88r(:) = 0.0d0
      q1000 (:) = 0.0d0
      q1002r(:) = 0.0d0
      q1004r(:) = 0.0d0
      q1006r(:) = 0.0d0
      q1008r(:) = 0.0d0
      q1010r(:) = 0.0d0

      do it=1,2
        do ir=1,mv
          xx2    = x2(ir)
          yy2    = y2(ir)
          zz2    = z2(ir)

          xx4    = xx2 * xx2
          xx6    = xx4 * xx2
          xx8    = xx6 * xx2
          xx10   = xx8 * xx2

          yy4    = yy2 * yy2
          yy6    = yy4 * yy2
          yy8    = yy6 * yy2
          yy10   = yy8 * yy2

          zz4    = zz2 * zz2
          zz6    = zz4 * zz2
          zz8    = zz6 * zz2
          zz10   = zz8 * zz2

          rr2    = xx2 + yy2 + zz2
          rr4    = rr2 * rr2
          rr6    = rr4 * rr2
          rr8    = rr6 * rr2
          rr10   = rr8 * rr2

c         .......................................................... (x + iy)^n
          xpiy2  = xx2 -           yy2
          xpiy4  = xx4 -  6.d0*xx2*yy2 +           yy4
          xpiy6  = xx6 - 15.d0*xx4*yy2 + 15.d0*xx2*yy4 -           yy6
          xpiy8  = xx8 - 28.d0*xx6*yy2 + 70.d0*xx4*yy4 - 28.d0*xx2*yy6
     1                 +           yy8
          xpiy10 = xx10- 45.d0*xx8*yy2 +210.d0*xx6*yy4 -210.d0*xx4*yy6
     1                 + 45.d0*xx2*yy8 -           yy10

          rho = rhot(ir,it)

c         .............................. contribution of this point to the sums
          dnu    = rho
          dr2    = (xx2 + yy2 + zz2) * rho
          d20    = cof20   *(-     3.0d0 * zz2
     1                       +                   rr2 )         * rho
          d22r   = cof22                               * xpiy2 * rho
          d40    = cof40   *(     35.0d0 * zz4
     1                       -    30.0d0 * zz2 * rr2
     2                       +     3.0d0 *       rr4 )         * rho
          d42r   = cof42   *(-     7.0d0 * zz2
     1                       +                   rr2 ) * xpiy2 * rho
          d44r   = cof44                               * xpiy4 * rho
          d60    = cof60   *(-   231.0d0 * zz6
     1                       +   315.0d0 * zz4 * rr2
     2                       -   105.0d0 * zz2 * rr4
     2                       +     5.0d0 *       rr6 )         * rho
          d62r   = cof62   *(     33.0d0 * zz4
     1                       -    18.0d0 * zz2 * rr2
     2                       +                   rr4 ) * xpiy2 * rho
          d64r   = cof64   *(-    11.0d0 * zz2
     1                       +                   rr2 ) * xpiy4 * rho
          d66r   = cof66                               * xpiy6 * rho
          d80    = cof80   *(   6435.0d0 * zz8
     1                       - 12012.0d0 * zz6 * rr2
     2                       +  6930.0d0 * zz4 * rr4
     3                       -  1260.0d0 * zz2 * rr6
     4                       +    35.0d0 *       rr8)          * rho
          d82r   = cof82   *(-   143.0d0 * zz6
     1                       +   143.0d0 * zz4 * rr2
     2                       -    33.0d0 * zz2 * rr4
     3                       +                   rr6 ) * xpiy2 * rho
          d84r   = cof84   *(     65.0d0 * zz4
     1                       -    26.0d0 * zz2 * rr2
     2                       +                   rr4 ) * xpiy4 * rho
          d86r   = cof86   *(-    15.0d0 * zz2
     1                       +                   rr2 ) * xpiy6 * rho
          d88r   = cof88                               * xpiy8 * rho
          d1000  = cof1000 *(- 46189.0d0 * zz10
     1                       +109395.0d0 * zz8 * rr2
     2                       - 90090.0d0 * zz6 * rr4
     3                       + 30030.0d0 * zz4 * rr6
     4                       -  3465.0d0 * zz2 * rr8
     5                       +    63.0d0       * rr10)         * rho
          d1002r = cof1002 *(   4199.0d0 * zz8
     1                       -  6188.0d0 * zz6 * rr2
     2                       +  2730.0d0 * zz4 * rr4
     3                       -   364.0d0 * zz2 * rr6
     4                       +     7.0d0       * rr8 ) * xpiy2 * rho
          d1004r = cof1004 *(    323.0d0 * zz6
     1                       +   255.0d0 * zz4 * rr2
     2                       -    45.0d0 * zz2 * rr4
     3                       +                   rr6 ) * xpiy4 * rho
          d1006r = cof1006 *(    323.0d0 * zz4
     1                       -   102.0d0 * zz2 * rr2
     2                       +     3.0d0       * rr4 ) * xpiy6 * rho
          d1008r = cof1008 *(-    19.0d0 * zz2
     1                       +                   rr2 ) * xpiy8 * rho
          d1010r = cof1010                             * xpiy10* rho

c         ................................................ add increment to sum
          anu   (it) = anu   (it) + dv * dnu
          r2    (it) = r2    (it) + dv * dr2
          q20   (it) = q20   (it) + dv * d20
          q22r  (it) = q22r  (it) + dv * d22r
          q40   (it) = q40   (it) + dv * d40
          q42r  (it) = q42r  (it) + dv * d42r
          q44r  (it) = q44r  (it) + dv * d44r
          q60   (it) = q60   (it) + dv * d60
          q62r  (it) = q62r  (it) + dv * d62r
          q64r  (it) = q64r  (it) + dv * d64r
          q66r  (it) = q66r  (it) + dv * d66r
          q80   (it) = q80   (it) + dv * d80
          q82r  (it) = q82r  (it) + dv * d82r
          q84r  (it) = q84r  (it) + dv * d84r
          q86r  (it) = q86r  (it) + dv * d86r
          q88r  (it) = q88r  (it) + dv * d88r
          q1000 (it) = q1000 (it) + dv * d1000
          q1002r(it) = q1002r(it) + dv * d1002r
          q1004r(it) = q1004r(it) + dv * d1004r
          q1006r(it) = q1006r(it) + dv * d1006r
          q1008r(it) = q1008r(it) + dv * d1008r
          q1010r(it) = q1010r(it) + dv * d1010r
        enddo
      enddo

c     .................................... dimensionless deformation parameters
      s13 = 1.0d0/3.0d0
      r   = 1.2d0 * (anu(1)+anu(2))**s13

      fac2 = 4.0d0 * pi / (3.0d0*r*r)
      beta20(1) = fac2 * q20 (1) / anu(1)
      beta22(1) = fac2 * q22r(1) / anu(1)
      beta20(2) = fac2 * q20 (2) / anu(2)
      beta22(2) = fac2 * q22r(2) / anu(2)
      beta20(3) = fac2 * (q20 (1) + q20 (2)) / (anu(1) + anu(2))
      beta22(3) = fac2 * (q22r(1) + q22r(2)) / (anu(1) + anu(2))

      fac4 = 4.0d0 * pi / (3.0d0*r*r*r*r)
      beta40(1) = fac4 * q40 (1) / anu(1)
      beta42(1) = fac4 * q42r(1) / anu(1)
      beta44(1) = fac4 * q44r(1) / anu(1)
      beta40(2) = fac4 * q40 (2) / anu(2)
      beta42(2) = fac4 * q42r(2) / anu(2)
      beta44(2) = fac4 * q44r(2) / anu(2)
      beta40(3) = fac4 * (q40 (1) + q40 (2)) / (anu(1) + anu(2))
      beta42(3) = fac4 * (q42r(1) + q42r(2)) / (anu(1) + anu(2))
      beta44(3) = fac4 * (q44r(1) + q44r(2)) / (anu(1) + anu(2))

      fac6 = 4.0d0 * pi / (3.0d0*r*r*r*r*r*r)
      beta60(1) = fac6 * q60 (1) / anu(1)
      beta62(1) = fac6 * q62r(1) / anu(1)
      beta64(1) = fac6 * q64r(1) / anu(1)
      beta66(1) = fac6 * q66r(1) / anu(1)
      beta60(2) = fac6 * q60 (2) / anu(2)
      beta62(2) = fac6 * q62r(2) / anu(2)
      beta64(2) = fac6 * q64r(2) / anu(2)
      beta66(2) = fac6 * q66r(2) / anu(2)
      beta60(3) = fac6 * (q60 (1) + q60 (2)) / (anu(1) + anu(2))
      beta62(3) = fac6 * (q62r(1) + q62r(2)) / (anu(1) + anu(2))
      beta64(3) = fac6 * (q64r(1) + q64r(2)) / (anu(1) + anu(2))
      beta66(3) = fac6 * (q66r(1) + q66r(2)) / (anu(1) + anu(2))

      fac8 = 4.0d0 * pi / (3.0d0*r*r*r*r*r*r*r*r)
      beta80(1) = fac8 * q80 (1) / anu(1)
      beta82(1) = fac8 * q82r(1) / anu(1)
      beta84(1) = fac8 * q84r(1) / anu(1)
      beta86(1) = fac8 * q86r(1) / anu(1)
      beta88(1) = fac8 * q88r(1) / anu(1)
      beta80(2) = fac8 * q80 (2) / anu(2)
      beta82(2) = fac8 * q82r(2) / anu(2)
      beta84(2) = fac8 * q84r(2) / anu(2)
      beta86(2) = fac8 * q86r(2) / anu(2)
      beta88(2) = fac8 * q88r(2) / anu(2)
      beta80(3) = fac8 * (q80 (1) + q80 (2)) / (anu(1) + anu(2))
      beta82(3) = fac8 * (q82r(1) + q82r(2)) / (anu(1) + anu(2))
      beta84(3) = fac8 * (q84r(1) + q84r(2)) / (anu(1) + anu(2))
      beta86(3) = fac8 * (q86r(1) + q86r(2)) / (anu(1) + anu(2))
      beta88(3) = fac8 * (q88r(1) + q88r(2)) / (anu(1) + anu(2))

      fac10 = 4.0d0 * pi / (3.0d0*r*r*r*r*r*r*r*r*r*r)
      beta1000(1) = fac10 *  q1000 (1) / anu(1)
      beta1002(1) = fac10 *  q1002r(1) / anu(1)
      beta1004(1) = fac10 *  q1004r(1) / anu(1)
      beta1006(1) = fac10 *  q1006r(1) / anu(1)
      beta1008(1) = fac10 *  q1008r(1) / anu(1)
      beta1010(1) = fac10 *  q1010r(1) / anu(1)
      beta1000(2) = fac10 *  q1000 (2) / anu(2)
      beta1002(2) = fac10 *  q1002r(2) / anu(2)
      beta1004(2) = fac10 *  q1004r(2) / anu(2)
      beta1006(2) = fac10 *  q1006r(2) / anu(2)
      beta1008(2) = fac10 *  q1008r(2) / anu(2)
      beta1010(2) = fac10 *  q1010r(2) / anu(2)
      beta1000(3) = fac10 * (q1000 (1) + q1000 (2)) / (anu(1) + anu(2))
      beta1002(3) = fac10 * (q1002r(1) + q1002r(2)) / (anu(1) + anu(2))
      beta1004(3) = fac10 * (q1004r(1) + q1004r(2)) / (anu(1) + anu(2))
      beta1006(3) = fac10 * (q1006r(1) + q1006r(2)) / (anu(1) + anu(2))
      beta1008(3) = fac10 * (q1008r(1) + q1008r(2)) / (anu(1) + anu(2))
      beta1010(3) = fac10 * (q1010r(1) + q1010r(2)) / (anu(1) + anu(2))

c     ..................................................... diagnostic printing
      xnnn = anu(1)
      xnnz = anu(2)
      xnna = xnnz + xnnn
      rmsn = 0.0d0
      rmsp = 0.0d0
      rmsa = 0.0d0
      if (xnnn.gt.1d-14) rmsn =  r2(1)        / xnnn
      if (xnnz.gt.1d-14) rmsp =  r2(2)        / xnnz
      if (xnna.gt.1d-14) rmsa = (r2(1)+r2(2)) / xnna
      rrmsn = 0.0d0
      rrmsp = 0.0d0
      rrmsa = 0.0d0
      if (rmsn.gt.0.0d0) rrmsn = sqrt(rmsn)
      if (rmsp.gt.0.0d0) rrmsp = sqrt(rmsp)
      if (rmsa.gt.0.0d0) rrmsa = sqrt(rmsa)

      q20t = q20 (1) + q20 (2)
      q22t = q22r(1) + q22r(2)

      q40t = q40 (1) + q40 (2)
      q42t = q42r(1) + q42r(2)
      q44t = q44r(1) + q44r(2)

      q60t = q60 (1) + q60 (2)
      q62t = q62r(1) + q62r(2)
      q64t = q64r(1) + q64r(2)
      q66t = q66r(1) + q66r(2)

      q80t = q80 (1) + q80 (2)
      q82t = q82r(1) + q82r(2)
      q84t = q84r(1) + q84r(2)
      q86t = q86r(1) + q86r(2)
      q88t = q88r(1) + q88r(2)

      q1000t = q1000 (1) + q1000 (2)
      q1002t = q1002r(1) + q1002r(2)
      q1004t = q1004r(1) + q1004r(2)
      q1006t = q1006r(1) + q1006r(2)
      q1008t = q1008r(1) + q1008r(2)
      q1010t = q1010r(1) + q1010r(2)

      print 105,q20   (1),q20   (2),q20t  ,
     1          q22r  (1),q22r  (2),q22t
      print 115,q40   (1),q40   (2),q40t  ,
     1          q42r  (1),q42r  (2),q42t  ,
     2          q44r  (1),q44r  (2),q44t
      print 125,q60   (1),q60   (2),q60t  ,
     1          q62r  (1),q62r  (2),q62t  ,
     2          q64r  (1),q64r  (2),q64t  ,
     3          q66r  (1),q66r  (2),q66t
      print 106,beta20(1),beta20(2),beta20(3),
     1          beta22(1),beta22(2),beta22(3)
      print 116,beta40(1),beta40(2),beta40(3),
     1          beta42(1),beta42(2),beta42(3),
     2          beta44(1),beta44(2),beta44(3)
      print 126,beta60(1),beta60(2),beta60(3),
     1          beta62(1),beta62(2),beta62(3),
     2          beta64(1),beta64(2),beta64(3),
     3          beta66(1),beta66(2),beta66(3)
      print 136,beta80(1),beta80(2),beta80(3),
     1          beta82(1),beta82(2),beta82(3),
     2          beta84(1),beta84(2),beta84(3),
     3          beta86(1),beta86(2),beta86(3),
     4          beta88(1),beta88(2),beta88(3)
      print 146,beta1000(1),beta1000(2),beta1000(3),
     1          beta1002(1),beta1002(2),beta1002(3),
     2          beta1004(1),beta1004(2),beta1004(3),
     3          beta1006(1),beta1006(2),beta1006(3),
     4          beta1008(1),beta1008(2),beta1008(3),
     5          beta1010(1),beta1010(2),beta1010(3)

c     .................................. call routine for LDM shape parameters
c                                            only makes sense for axial shapes
c     ........................................................................
c                for historical reasons, the input is the cartesian quadrupole
c                moment Q0 and the spherical hexadecapole moment Y40
c     ........................................................................
      if ( abs(beta22(3)) .lt. 2.0d-2 ) then
        print 100
        qt0 = -q20t / cof20
        rrr = 5.0d0/3.0d0 * rmsa
        callLDMshape (rrr,xnna,qt0,q40t)
      endif

      print 100

      return
      end subroutine momehigh
      
c______________________________________________________________________________
      subroutine LDMshape (rr,at,qt0,y40)

c..............................................................................
c     calculate the equivalent multipole expansion parameters of a liquid drop.
c     assuming an axial shape, symmetric around the z axis                    .
c..............................................................................
c     The formulas for these parameters b2 & b4 can be found in               .
c     S. Cwiok et al., Nuclear Physics A 611, (1996) pages 211-246            .
c                                                                             .
c     They are a second order development of equations 6.53 & 6.54 in         .
c     R. W. Hasse & W. D. Myers, Geometrical relationships of Macroscopic     .
c     Nuclear Physics, Springer, Berlin, 1988                                 .
c                                                                             .
c     These equations are solved iteratively.                                 .
c..............................................................................
      implicit real*8 (a-h,o-z)
      parameter (epsm3=1.0d-3)

c..............................................................................
  100 format (  ' Equivalent LDM multipole deformation parameters',
     1        /,' alpha2 : ',f6.3,
     2        /,' alpha4 : ',f6.3)

c............................................................... some constants
      pi  = 4.0d0*atan2(1.0d0,1.0d0)
      ccc = sqrt(5.0d0*pi)/3.0d0
      ddd = sqrt(      pi)
      c22 =   6.0d0 / (   7.0d0*pi            )
      c44 =  60.0d0 / (  77.0d0*pi            )
      c24 =  36.0d0 / (   7.0d0*pi*sqrt(5.0d0))
      d22 =   9.0d0 / (   7.0d0*pi            )
      d44 = 729.0d0 / (1001.0d0*pi            )
      d24 = 300.0d0 / (  77.0d0*pi*sqrt(5.0d0))

c..............................................................................
      b2  = 0.0d0
      b4  = 0.0d0
      do n24=1,60
        bb2 = b2
        bb4 = b4
        q2m = qt0/(at*rr)
        q4m = y40/(at*rr*rr)
        cmm = c22*bb2*bb2+c44*bb4*bb4+c24*bb2*bb4
        dmm = d22*bb2*bb2+d44*bb4*bb4+d24*bb2*bb4
        b2  = ( bb2 + ccc*(q2m-cmm) ) * 0.5d0
        b4  = ( bb4 + ddd*(q4m-dmm) ) * 0.5d0
        if ((abs(b2-bb2).le.epsm3).and.
     1      (abs(b4-bb4).le.epsm3)) exit
      enddo

      print 100,b2,b4

      return
      end subroutine LDMshape
c______________________________________________________________________________
      subroutine reord

c..............................................................................
c     diagonalize the matrix of the h.f. hamiltonian                          .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,ndim=100)

      common /Lag  / ilag, iAnaLag
      common /big  / h(ndim,ndim),s(ndim,ndim),d(ndim),wd(ndim)
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / wf(mq),phi(mq)
      common /waved/ wx1(3*mq)

      dimension index(ndim)

c..............................................................................
  100 format (/,
     1    ' ***danger in reord*** too many orbits for dh size it=',i2)

c..............................................................................
      do it=1,2
        if (max0(npar(1,it),npar(2,it)).gt.ndim) then
          print 100,it
          return
        endif
      enddo

      nof = 0
      do it=1,2
      do ip=1,2
        nvb=npar(ip,it)
        do iv=1,nvb
          iw = nof + iv
          iz = kparz(iw)
          call scopy (mq,a(1,iw),wf)
          call scopy (3*mq,da(1,iw),wx1)
          call hpsi (iw,iz,it)
          do jv=1,iv
            jw = nof + jv
            c  = sdot(mq,phi,a(1,jw))
            h(iv,jv) = c*dv
            h(jv,iv) = c*dv
          enddo
        enddo
        n=nvb
        call diagon (h,ndim,nvb,s,d,wd)
        do iv=1,nvb
          iw = nof+iv
          esp1(iw) = d(iv)
        enddo

c       ...................... before phi(i) =   sum(1 to nvb) of s(k,i)*psi(k)
c                              after  phi(i) =   sum(1 to i-1) of s(i,k)*phi(k)
c                                              + sum(i to nvb) of s(i,k)*psi(k)
c              note the permutation of the role of the indices in the matrix s
        do i=1,nvb
          s(1,i)=s(i,1)
        enddo
        do iv=2,nvb
          id=iv-1
          n=nvb-id
          do j=1,n
          do i=1,n
            h(i,j)=s(j+id,i+id)
          enddo
          enddo
          call matin1 (h,ndim,n,0,index,nerror,determ)
          do i=iv,nvb
            s(iv,i)=h(i-id,1)
          enddo
          do i=1,id
            ts = zero
            do j=iv,nvb
              ts=ts-s(iv,j)*s(j,i)
            enddo
            s(iv,i)=ts
          enddo
        enddo

c       ............................................................ new orbits
        do  iv=1,nvb
          do i=1,mq
            wf(i) = zero
          enddo
          do jv=1,nvb
            call scopy (mq,a(1,nof+jv),phi)
            call saxpy (mq,s(iv,jv),phi,wf)
          enddo
          iz = kparz(nof+iv)
          call scopy (mq,wf,a(1,nof+iv))

          if(ilag.eq.0) then
            call deriv (iz)
          else
            call derivlag(iz)
          endif

          call scopy (3*mq,wx1,da(1,nof+iv))
        enddo
        nof = nof + nvb
      enddo
      enddo

      return
      end subroutine reord
      
c______________________________________________________________________________
      subroutine pro

c..............................................................................
c     this routine computes p=(1.-dt*h/hbar)*p                                .
c     the initial wave-function p is destroyed                                .
c     esp1 and esp2 are divided by 2 as the spwf are normalized to 2          .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (two=2.0d0)

      common /kfpro/ tn
      common /nxyz / dx,dv
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / wf(mq),psi(mq)
      common /waved/ wx(mq),wy(mq),wz(mq)

c..............................................................................
      do iwa=1,nwave
        it = kiso(iwa)
        iz = kparz(iwa)
        call scopy (mq,a(1,iwa),wf)
        call scopy (3*mq,da(1,iwa),wx)
        call hpsi  (iwa,iz,it)
        call saxpy (mq,tn,psi,wf)
        call scopy (mq,wf,a(1,iwa))
        x  = dv * sdot(mq,psi,psi) / two
        xx = dv * sdot(mq,psi,wf)  / two
        esp1(iwa) = xx - tn * x
        esp2(iwa) = x
      enddo

      return
      end subroutine pro

c______________________________________________________________________________
      subroutine hpsi (iw,iz,it)

c..............................................................................
c     this routine computes |p> = h |w>                                       .
c     vcen computes the action of the off-diagonal part of the kinetic        .
c     energy and effective mass operator                                      .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (half=0.5d0)

      common /ener / b1,b2,b3,b4,b5,b6,b7a,b8a,b7b,b8b,b9,b9q
     1              ,b14,b15,b16,b17,byt3a,byt3b
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /evpro/ rx(mv,2),ry(mv,2),rz(mv,2)
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /pot  / wnp(mv,2),wcd(mv),wce(mv),wt3a(mv),wt3b(mv)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stord/ wd(mv,4,3,mw)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv)
     1              ,p1(mv),p2(mv),p3(mv),p4(mv)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /wcm2 / tx(mw,mw),ty(mw,mw),tz(mw,mw)
      common /wj2  / vjxx(mv,2),vjyx(mv,2),vjzx(mv,2)
     1              ,vjxy(mv,2),vjyy(mv,2),vjzy(mv,2)
     2              ,vjxz(mv,2),vjyz(mv,2),vjzz(mv,2)
      common /wmunu/ wjxx(mv,2),wjyx(mv,2),wjzx(mv,2)
     1              ,wjxy(mv,2),wjyy(mv,2),wjzy(mv,2)
     2              ,wjxz(mv,2),wjyz(mv,2),wjzz(mv,2)

c........................ contribution of the off-diagonal matrix elements of h
c                                    and from the spin-orbit one-body potential

      do i=1,mv
        p1(i) = wnp(i,it)*w1(i) - rx(i,it) * (wy2(i)+wz3(i))
     1                          - ry(i,it) * (wz4(i)-wx2(i))
     2                          + rz(i,it) * (wx3(i)+wy4(i))
        p2(i) = wnp(i,it)*w2(i) + rx(i,it) * (wy1(i)-wz4(i))
     1                          + ry(i,it) * (wz3(i)-wx1(i))
     2                          + rz(i,it) * (wx4(i)-wy3(i))
        p3(i) = wnp(i,it)*w3(i) + rx(i,it) * (wy4(i)+wz1(i))
     1                          - ry(i,it) * (wz2(i)+wx4(i))
     2                          - rz(i,it) * (wx1(i)-wy2(i))
        p4(i) = wnp(i,it)*w4(i) - rx(i,it) * (wy3(i)-wz2(i))
     1                          + ry(i,it) * (wz1(i)+wx3(i))
     2                          - rz(i,it) * (wx2(i)+wy1(i))
      enddo

c     .........................................................................
      call vcen (iz,it)

c     .................... contribution of the 2-body centre of mass correction
c                                                     to the one-body potential
      if (ncm2.eq.1) then
        i1 = 1 + nwaven*(it-1)
        i2 =     nwaven        + nwavep*(it-1)
        do jw=i1,i2
          fac = hbm(it)*(xm(it)/xm(3))*v2(jw)
          do i=1,mv
            p1(i) = p1(i) + fac * ( tx(iw,jw)*wd(i,3,1,jw)
     1                             +ty(iw,jw)*wd(i,4,2,jw)
     2                             +tz(iw,jw)*wd(i,1,3,jw) )
            p2(i) = p2(i) + fac * (-tx(iw,jw)*wd(i,4,1,jw)
     1                             +ty(iw,jw)*wd(i,3,2,jw)
     2                             +tz(iw,jw)*wd(i,2,3,jw) )
            p3(i) = p3(i) + fac * (-tx(iw,jw)*wd(i,1,1,jw)
     1                             -ty(iw,jw)*wd(i,2,2,jw)
     2                             +tz(iw,jw)*wd(i,3,3,jw) )
            p4(i) = p4(i) + fac * ( tx(iw,jw)*wd(i,2,1,jw)
     1                             -ty(iw,jw)*wd(i,1,2,jw)
     2                             +tz(iw,jw)*wd(i,4,3,jw) )
          enddo
        enddo
      endif

c     ................................... contribution of the spin-orbit tensor
      if (njmunu.ge.1) then
        do i=1,mv
          t1 = wjxx(i,it)
          t2 = wjyx(i,it)
          t3 = wjzx(i,it)
          t4 = wjxy(i,it)
          t5 = wjyy(i,it)
          t6 = wjzy(i,it)
          t7 = wjxz(i,it)
          t8 = wjyz(i,it)
          t9 = wjzz(i,it)
          ux = 0.0d0
          uy = 0.0d0
          uz = 0.0d0
          p1(i) = p1(i) + (t1*wx4(i)+t2*wy4(i)+t3*wz4(i)+ux*w4(i))
     1                  - (t4*wx3(i)+t5*wy3(i)+t6*wz3(i)+uy*w3(i))
     2                  + (t7*wx2(i)+t8*wy2(i)+t9*wz2(i)+uz*w2(i))
          p2(i) = p2(i) - (t1*wx3(i)+t2*wy3(i)+t3*wz3(i)+ux*w3(i))
     1                  - (t4*wx4(i)+t5*wy4(i)+t6*wz4(i)+uy*w4(i))
     2                  - (t7*wx1(i)+t8*wy1(i)+t9*wz1(i)+uz*w1(i))
          p3(i) = p3(i) + (t1*wx2(i)+t2*wy2(i)+t3*wz2(i)+ux*w2(i))
     1                  + (t4*wx1(i)+t5*wy1(i)+t6*wz1(i)+uy*w1(i))
     2                  - (t7*wx4(i)+t8*wy4(i)+t9*wz4(i)+uz*w4(i))
          p4(i) = p4(i) - (t1*wx1(i)+t2*wy1(i)+t3*wz1(i)+ux*w1(i))
     1                  + (t4*wx2(i)+t5*wy2(i)+t6*wz2(i)+uy*w2(i))
     2                  + (t7*wx3(i)+t8*wy3(i)+t9*wz3(i)+uz*w3(i))
        enddo
      endif

      return
      end subroutine hpsi

c______________________________________________________________________________
      subroutine vcen (iz,it)

c..............................................................................
c
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (one=1.0d0,two=2.0d0)

      common /lag  / ilag, iAnaLag
      common /evcen/ r(mq,2)
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /nxyz / dx,dv
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv),ps(mq)
      common /waved/ wx(mq),wy(mq),wz(mq)
      common /work / rx1(mv),rx2(mv),rx3(mv),rx4(mv)
     1              ,ry1(mv),ry2(mv),ry3(mv),ry4(mv)
     2              ,rz1(mv),rz2(mv),rz3(mv),rz4(mv)
     3              ,tx1(mv),tx2(mv),tx3(mv),tx4(mv)
     4              ,ty1(mv),ty2(mv),ty3(mv),ty4(mv)
     5              ,tz1(mv),tz2(mv),tz3(mv),tz4(mv)
     6              ,t1(mv),t2(mv),t3(mv),t4(mv)
      common /wtmp / ta(mq)

c..............................................................................
      zp = iz

      if(ilag.eq.0) then
        call lapla (w1,t1, one, one, zp)
        call lapla (w2,t2,-one,-one, zp)
        call lapla (w3,t3,-one, one,-zp)
        call lapla (w4,t4, one,-one,-zp)
      else
        call laplalag (w1,t1, 1, 1, iz)
        call laplalag (w2,t2,-1,-1, iz)
        call laplalag (w3,t3,-1, 1,-iz)
        call laplalag (w4,t4, 1,-1,-iz)
      endif

      do i=1,mq
        rx1(i) = r(i,it) * wx(i)
        ry1(i) = r(i,it) * wy(i)
        rz1(i) = r(i,it) * wz(i)
      enddo

      call derx (rx1,tx1,-1)
      call derx (rx2,tx2, 1)
      call derx (rx3,tx3, 1)
      call derx (rx4,tx4,-1)

      call dery (ry1,ty1,-1)
      call dery (ry2,ty2, 1)
      call dery (ry3,ty3,-1)
      call dery (ry4,ty4, 1)

      call derz (rz1,tz1,-iz)
      call derz (rz2,tz2,-iz)
      call derz (rz3,tz3, iz)
      call derz (rz4,tz4, iz)

      cmfac = 1.0d0
      if (ncm2.ge.0) cmfac = 1.0d0 - xm(it)/xm(3)
      fac = 0.5d0 * hbm(it) * cmfac
      do i=1,mq
        ta(i) = fac * t1(i) + (tx1(i)+ty1(i)+tz1(i))
        ps(i) = ps(i) - ta(i)
      enddo

      return
      end subroutine vcen
c______________________________________________________________________________      
      subroutine jdeux

c..............................................................................
c     compute <Jx**2>, <Jy**2> and <Jz**2>                                    .
c     the Belyaev moments of inertia, the rigid moments of inertia and        .
c     the rotational correction for the b.c.s. state (not printed!)           .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (half=0.5d0)
      parameter (epst=0.01d0)
      parameter (clum  =  29.9792458d0)  !Speed of light added 

      common /mud  / xi(mv),yi(mv),zi(mv),xii(mv),yii(mv),zii(mv)
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / w1(mv),w2(mv),w3(mv),w4(mv),brik(mq)
      common /waved/ wx1(mv),wx2(mv),wx3(mv),wx4(mv)
     1              ,wy1(mv),wy2(mv),wy3(mv),wy4(mv)
     2              ,wz1(mv),wz2(mv),wz3(mv),wz4(mv)
      common /work / wjx1(mv),wjx2(mv),wjx3(mv),wjx4(mv)
     1              ,wjy1(mv),wjy2(mv),wjy3(mv),wjy4(mv)
     2              ,wjz1(mv),wjz2(mv),wjz3(mv),wjz4(mv)
     3              ,wks(4*mq)

c..............................................................................
  100 format (/,' ',78('x'),/,13x,'jx2',5x,'jy2',5x,'jz2',4x,'total')
  101 format (' is=',i1,'ip=',i1,4f8.2)
  102 format (' is=',i1,4x,4f8.2)
  103 format (9x,4f8.2)
  104 format (' ',78('x'))
  110 format (/,14x,'x',9x,'y',9x,'z')
  111 format (  ' Rotational properties',
     1        /,35x,'x',12x,'y',12x,'z',11x,'total',
     2        /,' <J2>        (hbar^2)    : ',4(f12.3,1x))
  112 format (  ' Belyaev     (hbar^2/MeV): ',3(f12.3,1x))
  120 format (  ' Rigid rotor (hbar^2/MeV): ',3(f12.3,1x))
  
c..............................................................................
      ipridebug = 0

      if (ipridebug.eq.1) print 100

      iwa = 0
      dvs2  = dv/two         ! volume element for integrals over wave functions
                             ! canceling the norm factor 2
      do i=1,mq
        wjx1(i) = zero
        wjy1(i) = zero
        wjz1(i) = zero
      enddo

      xjt    = zero
      yjt    = zero
      zjt    = zero
      xtheta = zero
      ytheta = zero
      ztheta = zero

      do it=1,2
        xjis = zero
        yjis = zero
        zjis = zero
        do ib=1,2
          xj2 = zero
          yj2 = zero
          zj2 = zero
          nvb = npar(ib,it)

c         .......................................... first loop on the vectors
          do iv=1,nvb
            iwa=iwa+1
            call scopy (mq,a(1,iwa),w1)
            v2i = v2(iwa)
            asq = v2i*(one-v2i)
            uvi = zero
            if (asq.gt.zero) uvi = sqrt(asq)
            call scopy (3*mq,da(1,iwa),wx1)

c           ......................... calculation of (i|j_mu|j) and (i|jmu^2|i)
            do i=1,mv
              wjx1(i) = half*w3(i) + (+yi(i)*wz2(i)-zi(i)*wy2(i))
              wjx2(i) = half*w4(i) + (-yi(i)*wz1(i)+zi(i)*wy1(i))
              wjx3(i) = half*w1(i) + (+yi(i)*wz4(i)-zi(i)*wy4(i))
              wjx4(i) = half*w2(i) + (-yi(i)*wz3(i)+zi(i)*wy3(i))
              wjy1(i) = half*w4(i) + (+zi(i)*wx2(i)-xi(i)*wz2(i))
              wjy2(i) =-half*w3(i) + (-zi(i)*wx1(i)+xi(i)*wz1(i))
              wjy3(i) =-half*w2(i) + (+zi(i)*wx4(i)-xi(i)*wz4(i))
              wjy4(i) = half*w1(i) + (-zi(i)*wx3(i)+xi(i)*wz3(i))
              wjz1(i) = half*w1(i) + (+xi(i)*wy2(i)-yi(i)*wx2(i))
              wjz2(i) = half*w2(i) + (-xi(i)*wy1(i)+yi(i)*wx1(i))
              wjz3(i) =-half*w3(i) + (+xi(i)*wy4(i)-yi(i)*wx4(i))
              wjz4(i) =-half*w4(i) + (-xi(i)*wy3(i)+yi(i)*wx3(i))
            enddo
            xxj = dvs2 * sdot(mq,wjx1,wjx1)
            yyj = dvs2 * sdot(mq,wjy1,wjy1)
            zzj = dvs2 * sdot(mq,wjz1,wjz1)
            xj  = dvs2 *(-sdot(mv,w3,wjx1) + sdot(mv,w4,wjx2)
     1                   +sdot(mv,w1,wjx3) - sdot(mv,w2,wjx4))
            yj  = dvs2 *(-sdot(mv,w4,wjy1) - sdot(mv,w3,wjy2)
     1                   +sdot(mv,w2,wjy3) + sdot(mv,w1,wjy4))
            zj  = dvs2 *(+sdot(mv,w1,wjz1) + sdot(mv,w2,wjz2)
     1                   +sdot(mv,w3,wjz3) + sdot(mv,w4,wjz4))
            xj2 = xj2 + 2.0d0 * v2i * (xxj - xj*xj)
            yj2 = yj2 + 2.0d0 * v2i * (yyj - yj*yj)
            zj2 = zj2 + 2.0d0 * v2i * (zzj - zj*zj)
            if (iv.gt.1) then
              do jv=1,iv-1                                ! sum over j < i only
                jwa = iwa-iv+jv
                v2j = v2(jwa)
                asq = v2j*(one-v2j)
                uviuvj = zero
                if (asq.gt.zero) uviuvj = uvi*sqrt(asq) 
                v2iv2j = v2j*v2i 
                call scopy (mq,a(1,jwa),w1)
                xj  = dvs2 *(-sdot(mv,w3,wjx1) + sdot(mv,w4,wjx2)
     1                       +sdot(mv,w1,wjx3) - sdot(mv,w2,wjx4))
                yj  = dvs2 *(-sdot(mv,w4,wjy1) - sdot(mv,w3,wjy2)
     1                       +sdot(mv,w2,wjy3) + sdot(mv,w1,wjy4))
                zj  = dvs2 *(+sdot(mv,w1,wjz1) + sdot(mv,w2,wjz2)
     1                       +sdot(mv,w3,wjz3) + sdot(mv,w4,wjz4))
                xj2 = xj2 - 4.0d0*xj*xj * (v2iv2j + uviuvj)
                yj2 = yj2 - 4.0d0*yj*yj * (v2iv2j + uviuvj)
                zj2 = zj2 - 4.0d0*zj*zj * (v2iv2j + uviuvj)

c               ..................................... Belyaev moment of inertia
c               NOTE: there is the square of the factor (ui vj - vi uj) missing
c                     in Eq. (6) of Bender, Bonche, Heenen, PRC70 (2004) 054304
c              (ui vj - vi uj)^2 = ui^2 vj^2 + vi^2 uj^2 - 2 ui vi uj vj
c                                = (1-vi^2) vj^2 + vi^2 (1-vj^2)- 2 ui vi uj vj
c                                = vj^2 + vi^2 - 2 vi^2 vj^2 - 2 ui vi uj vj
c               ...............................................................
                xeq = eqp(iwa) + eqp(jwa)
                fac = zero
                if (xeq.gt.zero) then 
                  fac = (v2i  + v2j - 2.0d0*(v2iv2j+uviuvj))/xeq
                endif
                fac = fac * 2.0d0       ! factor 2 from summing over i < j only
                fac = fac * 2.0d0       ! factor 2 from Belyaev formula
                xtheta = xtheta + fac * xj*xj
                ytheta = ytheta + fac * yj*yj
                ztheta = ztheta + fac * zj*zj
              enddo
            endif
          enddo
          djt = xj2 + yj2 + zj2
          if (ipridebug.eq.1) print 101,it,ib,xj2,yj2,zj2,djt
          xjis = xjis + xj2
          yjis = yjis + yj2
          zjis = zjis + zj2
        enddo
        djt = xjis + yjis + zjis
        if (ipridebug.eq.1) print 102,it,xjis,yjis,zjis,djt
        xjt = xjt + xjis
        yjt = yjt + yjis
        zjt = zjt + zjis
      enddo
      djt = xjt + yjt + zjt
      if (ipridebug.eq.1) print 103,xjt,yjt,zjt,djt
      if (ipridebug.eq.1) print 104

      erx = zero
      ery = zero
      erz = zero
      if (xtheta.gt.epst) erx = xjt/(2.0d0*xtheta)
      if (ytheta.gt.epst) ery = yjt/(2.0d0*ytheta)
      if (ztheta.gt.epst) erz = zjt/(2.0d0*ztheta)
      ert = erx + ery + erz
      print 111,xjt,yjt,zjt,djt
      print 112,xtheta,ytheta,ztheta
      
c     ......Calculation of the quantum mechanical rigid rotor moments of inertia
c      The formula is:
c       I_z = m \int d^3r ( x^2 + y^2) * rho
c      and analogous for the x and y directions.
c 
c      However, this number, when calculated by the code, is then in units 
c      of Mev * fm^2/c^2, while h^2/MeV is more customary. 
c      With some checking, the conversion coefficient is given by 
c                    (hbar c)^2
c    ..........................................................................
        qix = xm(2)*(y2p+z2p) + xm(1)*(y2n+z2n)                      
        qiy = xm(2)*(z2p+x2p) + xm(1)*(z2n+x2n)                      
        qiz = xm(2)*(x2p+y2p) + xm(1)*(x2n+y2n)
        qix = qix/(hbar*clum)**2
        qiy = qiy/(hbar*clum)**2
        qiz = qiz/(hbar*clum)**2
        print 120, qix,qiy,qiz


      return
      end subroutine jdeux
c______________________________________________________________________________      
      subroutine pipj

c..............................................................................
c     calculation of the full one- + two-body c.m. correction                 .
c     the subroutine 'fprtj' has to be called first.                          .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*20 afor

      parameter (zero=0.0d0,tt4=4.0d0,epsm5=1.0d-5)

      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu(3),ecm(2,3),ecmp(3)
      common /force/ tsk(16),hbar,hbm(2),xm(3),afor
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /nxyz / dx,dv
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)
      common /stor / a(mq,mw)
      common /stord/ da(3*mq,mw)
      common /wave / w1(mx,my,mz),w2(mx,my,mz),w3(mx,my,mz),w4(mx,my,mz)
     1              ,p1(mx,my,mz),p2(mx,my,mz),p3(mx,my,mz),p4(mx,my,mz)
      common /waved/ wx1(mx,my,mz),wx2(mx,my,mz),wx3(mx,my,mz)
     1              ,wx4(mx,my,mz),wy1(mx,my,mz),wy2(mx,my,mz)
     2              ,wy3(mx,my,mz),wy4(mx,my,mz),wz1(mx,my,mz)
     3              ,wz2(mx,my,mz),wz3(mx,my,mz),wz4(mx,my,mz)
      common /wcm2 / tx(mw,mw),ty(mw,mw),tz(mw,mw)
      dimension pp1(2),pp2(2),ep2(2),ecmh(3),ecm1(3),ecm2(3)

c..............................................................................
  200 format (  ' Centre-of-mass correction', 
     1        /,32x,'neutron',6x,'proton',7x,'total')
  201 format (  '  <P2>     (hbar/fm)^2   : ',3f12.3)
  202 format (  ' -<P2>/2M  (MeV)         : ',3f12.3)
  203 format (  ' CM correction used here : ',3f12.3)

c....................................... deduce 1-body part from kinetic energy
c     The kinetic energy has been already calculated in fprtj. Depending on the
c     option for cm correction, it may or may not contain the one-body part of
c     the cm correction.
c     ncm2 = -1, -2: ek(:) is the kinetic energy
c     ncm2 =  0,  1: ek(:) is the kinetic energy with the 1-body cm correction 
c     .........................................................................
      pp1 (:) = 0.0d0
      ecm1(:) = 0.0d0
      do it=1,2
        if ( ncm2.eq.-1 .or. ncm2 .eq.-2 ) then
          ecm1(it) = -(xm(it)/xm(3)) * ek(it)
        endif
        if ( ncm2.eq.0 .or. ncm2 .eq.1 ) then
          ecm1(it) = ecm(1,it)
        endif
        pp1 (it) = -ecm1(it) * (xm(3)/xm(it)) / hbm(it)
      enddo
      ecm1(3) = ecm1(1) + ecm1(2)

c     ........................................................... two-body part
      call nabxyz

      do it=1,2
        pp2(it) = zero
        ep2(it) = zero
        i1 = 1 + nwaven*(it-1)
        i2 =     nwaven        + nwavep*(it-1)
        do ibar = i1,i2
        do iket = i1,i2
c         if (ibar.ne.iket) then          ! old. unclear if it's more efficient
          if (kparz(ibar).ne.kparz(iket)) then
            u2ibar = 1.0d0 - v2(ibar)
            u2iket = 1.0d0 - v2(iket)
            v2v2 = v2(ibar) * v2(iket)
            uvuv = sqrt(v2(ibar) * u2ibar * v2(iket) * u2iket)
            txyz = tx(ibar,iket)**2+ty(ibar,iket)**2+tz(ibar,iket)**2
            pp2(it) = pp2(it) - (v2v2 + uvuv) * txyz
          endif
        enddo
        enddo
        ecm2(it) = -hbm(it)*(xm(it)/xm(3)) * pp2(it)  ! 2-body cm energy
        pp2 (it) = pp1(it) + pp2(it)                  ! 1 + 2-body part of <P2>
        ep2 (it) = -hbm(it)*(xm(it)/xm(3)) * pp2(it)  ! 1 + 2-body cm energy
      enddo
      ppt  = pp2(1) + pp2(2)
      ep2t = ep2(1) + ep2(2)
      ecm2(3) = ecm2(1) + ecm2(2)

c     ..................................................... diagnostic printing
      print 200
      print 201,pp2(1),pp2(2),ppt
      print 202,ep2(1),ep2(2),ep2t

c     .............................................. cm correction as used here
      ecmh(:) = 0.0d0
      if ( ncm2.eq.0 .or. ncm2.eq.1 .or.ncm2.eq.-2 ) then   ! add one-body part
        ecmh(:) = ecm1(:)
      endif
      if ( ncm2.eq.1 .or. ncm2.eq.-2 ) then                 ! add two-body part
        ecmh(:) = ecmh(:) + ecm2(:)
      endif
      print 203,ecmh(1:3)

      return
      end subroutine pipj

c______________________________________________________________________________
      subroutine writ8

c..............................................................................
c     write single-particle wave functions and other information on the       .
c     calculation to tape fort.13                                             .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      character*4 head
      character*20 afor

      common /champ/ qxxn,qyyn,qzzn,qrrn
     1              ,qxxp,qyyp,qzzp,qrrp
     2              ,qxxt,qyyt,qzzt,qrrt

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mv)
     1              ,imtd,imtg,icutq
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

      !Please note that the variable alpha contains the inverse of the 
      ! mathematical variable alpha!
      if (iver.ne.20) then
        write (13) npair,gn,gp,encut,epcut,dcut,xcut,alpha,alphap
     1            ,ambda,xlamb,0
      else
        write (13) npair,vn,vp,rangen,rangep,ambda,xlamb,0
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

      call densit
      call newpot

      return
      end subroutine writ8
c______________________________________________________________________________
      subroutine conver (iprint,iterm)

c..............................................................................
c  check the convergence of the imaginary-time evolution                      .
c  convergence criteria are                                                   .
c  - relative change of the binding energy in the last 7 iterations           .
c    (this avoids accidental termination in case of oscillating iteration)    .
c                                                                             .
c          (E(i) - E(i-n))/E(i) < epse    for n=1,...,7                       .
c                                                                             .
c  - sum of the dispersion of single-particle Hamiltonian weigthed with       .
c    the occupation number is small (divided by the mass number A)            .
c                                                                             .
c          sum_mu v_mu^2 [ (mu|h^2|mu) - (mu|h|mu)^2 ] / A < epsdh            .
c                                                                             .
c  - the change in Fermi energies is small in the last 7 iterations           .
c    for protons and neutrons                                                 .
c                                                                             .
c          lambda(i) - lambda(i-n) < epsl     for n=1,...,7                   .
c                                                                             .
c  - the constrained quadrupole moment approaches the desired value Qcon      .
c    in the last 7 iterations. For Qcon = 0, the absolute value is taken      .
c                                                                             .
c          (Q(i-n) - Qcon )        < epsq      for n=0,...,6                  .
c                                                                             .
c    otherwise its the relative change                                        .
c                                                                             .
c          (Q(i-n) - Qcon ) / Qcon < epsq     for n=0,...,6                   .
c                                                                             .
c    Note that the constraint is on the undampened quadrupole moment,         .
c    although the constraint potential is calculated with the damping!        .
c                                                                             .
c  When all of the above is fulfilled, iterm is set to +1                     .
c  otherwise it is 0                                                          .
c                                                                             .
c  - The occupation of the highest state with given quantum numbers is        .
c    smaller than epsm (this checks if the model space is large enough)       .
c                                                                             .
c  parameter handling through common block /condat/                           .
c                                                                             .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'
      parameter (epsm1=1.0d-1,epsm3=1.0d-3,epsm4=1.0d-4)

      logical let,lef,ldh,lq,lm,lterm

      common /condat/ epsq,epse,epsdh,epsf,epsm,iprintconver
      common /consto/ etold(7),efold(2,7),qxtold(7),qytold(7),qztold(7)
      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mx,my,mz),imtd
     1              ,imtg,icutq
      common /cstt / qxtc,qxcstt,qxfint,excstt,pentext
     1              ,qytc,qycstt,qyfint,eycstt,penteyt
     2              ,qztc,qzcstt,qzfint,ezcstt,pentezt
     3              ,qrtc,qrcstt,qrfint,ercstt,pentert
     4              ,q2fin,g2fin,q2cst,qtc0,gtc0
      common /evohe/ dt,nitert,nxmu,ndiag,itert,nprint,iverb
      common /fopt / nfunc,njmunu,ncm2,nmass,ndd,ncoex
      common /iwrit/ et0,ett,etd,et3a,et3b,esk,ecoul,eso(2),ek(3),ecoex
     1              ,ejmunu(3),ecm(2,3),ecmp(3)
      common /kfmom/ tfac,tfac1,tfac2,s3
      common /mmtc / x2p,y2p,z2p,x2n,y2n,z2n
      common /noyau/ nwaven,nwavep,nwave,npn,npp,npar(2,2)
      common /pair / ambda(2),xlamb(2),epair(3),eproj(3),disper(3)
      common /pairf/ gn,gp,delmax(2),dcut,encut,epcut,xcut,alpha,alphap
     1              ,npair,ntqp,ilqp,ifor
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)

      dimension sdh(3),dhmx(3),ocmnm(2),ocmnp(2)

c..............................................................................
   91 format(' summary',1i5,1f12.3,8es10.2,4es12.4)
   93 format(' Sum. it  :', ' Iter   =' , i5, 
     1     /,' Sum. E   :', ' Etot   =' , f10.3 ,' dE      =',e10.3 , 
     2     /,' Sum. Quad:', ' Qxx    =' , f10.3, ' Qy      =',f10.3,
     3        ' Qz      =',f10.3, 
     4     /,' Sum. Misc:', ' Sum D2H=',e10.3, ' dFermi N=',e10.3,
     5       ' dFermi P=', e10.3 )
   94 format (78('-'))

c...................................................change in quadrupole moment
      lq = .true.
      qxt = 0.0d0
      qyt = 0.0d0
      qzt = 0.0d0
      if (cq2.gt.0.0d0) then
        qxt = 2.0d0*x2p-y2p-z2p + 2.0d0*x2n-y2n-z2n
        qyt = 2.0d0*y2p-z2p-x2p + 2.0d0*y2n-z2n-x2n
        qzt = 2.0d0*z2p-x2p-y2p + 2.0d0*z2n-x2n-y2n

        do i=7,2,-1
          qxtold(i) = qxtold(i-1)
          qytold(i) = qytold(i-1)
          qztold(i) = qztold(i-1)
        enddo
        qxtold(1) = qxt
        qytold(1) = qyt
        qztold(1) = qzt
        lq = (1.eq.1)
        if ((abs(qxfint).lt.1.0d0).and.(abs(qyfint).lt.1.0d0)
     1                            .and.(abs(qzfint).lt.1.0d0)) then
          do i=1,7
            lq = (abs(qxtold(i)-qxfint).lt.epsm1.and.
     1            abs(qytold(i)-qyfint).lt.epsm1.and.
     2            abs(qztold(i)-qzfint).lt.epsm1.and.lq)
          enddo
        else
          do i=1,7
            lq = (abs((qxtold(i)-qxfint)/qxfint).lt.epsq.and.
     1            abs((qytold(i)-qyfint)/qyfint).lt.epsq.and.
     2            abs((qztold(i)-qzfint)/qzfint).lt.epsq.and.lq)
          enddo
        endif
      endif

c......................................................... total binding energy
c                                      in case of nfunc.eq.1, ejmunu(2) is zero
      call fprte(0)
      call fprtj
      call fprtz

      epa = 0.0d0
      if (npair.ne.0) epa = epair(1) + epair(2)
      ekin = ek (1) + ek (2)
      et   = esk + ecoul + ecoex + ekin
      et   = et  + epa   - (eproj(1)+eproj(2))
      if (ncm2.eq.-2) et = et + ecm(1,3) + ecm(2,3)
      if (ncm2.eq. 1) et = et            + ecm(2,3)

c............................................... change in total binding energy
      let = (1.eq.1)
      do i=1,7
        let = (abs((et-etold(i))/et).lt.epse.and.let)
      enddo
      etold(7) = etold(6)
      etold(6) = etold(5)
      etold(5) = etold(4)
      etold(4) = etold(3)
      etold(3) = etold(2)
      etold(2) = etold(1)
      etold(1) = et

c...................... calculate the sum of the occupation-weighted dispersion
c                            of the single-particle Hamiltonian and its maximum
c                                          and check if there are enough states
      do it=1,2
        if (it.eq.1) then
          na = 1
          nb = nwaven
        else
          na = nwaven + 1
          nb = nwaven + nwavep
        endif
        ocmnp(it) = 1.0d0
        ocmnm(it) = 1.0d0
        sdh  (it) = 0.0d0
        dhmx (it) = 0.0d0
        do i=na,nb
          esp4 = esp2(i) - esp1(i)*esp1(i)
          esp4 = abs(esp4)
!          if (esp4.lt.esp3(i)) esp4 = -esp4                      !         ????
          dh      = v2(i) * esp4
          sdh(it) = sdh(it) + dh
          if (dh.gt.dhmx(it)) dhmx(it) = dh
          if (kparz(i).eq.+1.and.v2(i).lt.ocmnp(it)) ocmnp(it)= v2(i)
          if (kparz(i).eq.-1.and.v2(i).lt.ocmnm(it)) ocmnm(it)= v2(i)
        enddo
      enddo
      dh    = (abs(sdh(1))+abs(sdh(2))) / (npn+npp)
      dhmax = sqrt(max(dhmx(1),dhmx(2)))
      ldh   = (dh.lt.epsdh)
      lm    = (ocmnp(1).lt.epsm3.and.ocmnp(2).lt.epsm3.and.
     1         ocmnm(1).lt.epsm3.and.ocmnm(2).lt.epsm3)

c..................................................... change in Fermi energies
      lef = (1.eq.1)
      do it=1,2
        do i=1,7
          lef = (abs(ambda(it)-efold(it,i)).lt.epsf.and.lef)
        enddo
      enddo
      do it=1,2
        efold(it,7) = efold(it,6)
        efold(it,6) = efold(it,5)
        efold(it,5) = efold(it,4)
        efold(it,4) = efold(it,3)
        efold(it,3) = efold(it,2)
        efold(it,2) = efold(it,1)
        efold(it,1) = ambda(it)
      enddo

c.............................................................. take a decision
      lterm = (let.and.lef.and.ldh.and.lq)
      iterm = 0
      if (lterm) iterm = 1
      if (iprint.gt.1) then
        print *,' logic',itert,let,ldh,lef,lq,lm,lterm
      endif

c.......................................................... diagnostic printout
      if (iprint.gt.-1 )then
        dee = (et-etold(2))/et
        dfn = efold(1,1)-efold(1,2)
        dfp = efold(2,1)-efold(2,2)
        df  = max(abs(dfn),abs(dfp))
        dqx = qxt - qxfint
        dqy = qyt - qyfint
        dqz = qzt - qzfint
        if (abs(qxfint).gt.1.0d0.and.
     1      abs(qyfint).gt.1.0d0.and.
     2      abs(qzfint).gt.1.0d0) then
          dqx = dqx/qxfint
          dqy = dqy/qyfint
          dqz = dqz/qzfint
        endif
        if( iprintconver.eq.1 ) then
        print 91,itert,et,dee,dh,dhmax,df,ocmnp(1),ocmnm(1),
     1                                    ocmnp(2),ocmnm(2),
     2           dqx,dqy,dqz
        else
          print 93,itert, et, dee, qxt, qyt, qzt, dhmax, 
     1            abs(ambda(:) - efold(:,2))
          print 94
        endif
      endif
      return
      end subroutine conver
      
c______________________________________________________________________________
      subroutine distance (func,sval)

c..............................................................................
c     cutoff of the constraint depending on the density distribution and      .
c     profile as proposed by K. Rutz et al, Nucl. Phys. A590 (1995) 690       .
c                                                                             .
c     this subroutine calculates the minimum distance of each point on the    .
c     mesh to the surface which is defined by having the value sval in the    .
c     mesh function func and then defines the cutoff function cutof2 from     .
c     this distance such that                                                 .
c                                                                             .
c                       1                                                     .
c          --------------------------                                         .
c           1 + exp[(dist-radd)/acut]                                         .
c                                                                             .
c     sval : switch value, recommended value is 1/10 of the maximum of func   .
c     information handling through common /cst/                               .
c       cutof2  : cutoff function                                             .
c       radd = rcut/     : additional distance from the surface               .
c       acut             : width of the cutoff                                .
c                                                                             .
c     the values radd = 4.0 fm and acut = 0.4 fm should be used.              .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (t1000=1000.d0,tt5=5.0d0)

      common /cst  / ral,epscst,cqr,cq2,rcut,acut,cutof2(mx,my,mz),imtd
     1               ,imtg, icutq
      common /nxyz / dx,dv
      common /mud  / xi (mx,my,mz),yi (mx,my,mz),zi (mx,my,mz)
     1              ,xii(mx,my,mz),yii(mx,my,mz),zii(mx,my,mz)
      common /wave / surf(mv,3),sig(mx,my,mz),reste(mv,4)

      dimension func(mx,my,mz),dist(mx,my,mz)

c........................................................... define the surface
c                         i[x,y,z]    are the coordinates of the current  point
c                         i[xx,yy,zz] are the coordinates of the adjacent point

      is = 0                               ! index of field with surface points
      do i=1,3*mv
        surf(i,1) = zero
      enddo
      
      surfx = 75.0d0   ! initialization - to be checked
      surfy = 75.0d0   ! sometimes, the values are not set below, which is not
      surfz = 75.0d0   ! a problem for the algorithm when they remain zero,
                       ! but might have some consequences on some platforms
      do iz=1,mz
        do iy=1,my
          do ix=1,mx
            x = xi(ix,iy,iz)
            y = yi(ix,iy,iz)
            z = zi(ix,iy,iz)
            fxyz = func(ix,iy,iz)
            if (fxyz.lt.sval) then
              sig(ix,iy,iz) =  one
            else
              sig(ix,iy,iz) = -one
            endif
            ixx = ix
            iyy = iy
            izz = iz

c           ......................................................  x direction
            if (ix.ne.1) ixx = ix - 1
            ffxyz = func(ixx,iyy,izz)
            if( (fxyz.lt.sval.and.ffxyz.lt.sval)
     1           .or.(fxyz.ge.sval.and.ffxyz.ge.sval)) then
              jx = ix                 ! ix? what is jx??? is this "do nothing"?
            else
              frac = (sval - fxyz) / (ffxyz - fxyz)
              if (ix.ne.ixx) surfx = x - frac*dx
              if (is.ge.mq) then
                call stp (' distance: is > mq!')
              else
                is = is + 1
                surf(is,1) = surfx
                surf(is,2) = y
                surf(is,3) = z
              endif
            endif
            ixx = ix

c           ......................................................  y direction
            if (iy.ne.1) iyy = iy - 1
            ffxyz = func(ixx,iyy,izz)
            if( (fxyz.lt.sval.and.ffxyz.lt.sval)
     1           .or.(fxyz.ge.sval.and.ffxyz.ge.sval)) then
              jx = ix                
            else
              frac = (sval - fxyz) / (ffxyz - fxyz)
              if (iy.ne.iyy) surfy = y - frac*dx
              if (is.ge.mq) then
                call stp (' distance: is > mq!')
              else
                is = is + 1
                surf(is,1) = x
                surf(is,2) = surfy
                surf(is,3) = z
              endif
            endif
            iyy = iy

c           ......................................................  z direction
            if (iz.ne.1) izz = iz - 1
            ffxyz = func(ixx,iyy,izz)
            if( (fxyz.lt.sval.and.ffxyz.lt.sval)
     1           .or.(fxyz.ge.sval.and.ffxyz.ge.sval)) then
              jx = ix                 
            else
              frac = (sval - fxyz) / (ffxyz - fxyz)
              if (iz.ne.izz) surfz = z - frac*dx
              if (is.ge.mq) then
                call stp (' distance: is > mq!')
              else
                is = is + 1
                surf(is,1) = x
                surf(is,2) = y
                surf(is,3) = surfz
              endif
            endif
            izz = iz
          enddo
        enddo
      enddo

      ismax = is

c     ...................................................... calculate distance
      do iz=1,mz
        do iy=1,my
          do ix=1,mx
            dmin = 1.d+12
            x = xi(ix,iy,iz)
            y = yi(ix,iy,iz)
            z = zi(ix,iy,iz)
            do is=1,ismax
              xx = surf(is,1)
              yy = surf(is,2)
              zz = surf(is,3)
              d = (x-xx)*(x-xx)+(y-yy)*(y-yy)+(z-zz)*(z-zz)
              dmin = min(dmin,d)
            enddo
            dist(ix,iy,iz) = sqrt(dmin) * sig(ix,iy,iz)
          enddo
        enddo
      enddo

c     ................................ calculate cutoff function for constraint
      do iz=1,mz
        do iy=1,my
          do ix=1,mx
            earg = (dist(ix,iy,iz) - rcut) / acut
            cutof2(ix,iy,iz) = 1.0d0 / (1.0d0 + exp(earg))
          enddo
        enddo
      enddo


      return
      end subroutine distance
      
c______________________________________________________________________________
      subroutine diagon (a,ndim,n,v,d,wd)

c..............................................................................
c     diagonalization of a real symmetric matrix a(i,j)                       .
c..............................................................................
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

c..............................................................................
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
        if (j.eq.jstop) then
            call stp (
     1      'Subroutine diagon did not succeed. ')
        endif
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
      end subroutine diagon

c______________________________________________________________________________
      subroutine matin1 (a,idim,n1,n2,index,mistak,determ)

c..............................................................................
c     matrix inversion with accompanying solution of linear equations
c..............................................................................
      implicit real*8 (a-h,o-z)

      parameter (zero=0.0d0,one=1.0d0)

      dimension a(idim),index(idim)

c..............................................................................
  100 format (//,' matin1 ..... the ',i4,
     1           'th. column of the matrix contains only zeros at the ',
     2         /,18x,i4,'. elemination step.')

c..............................................................................
      deter = one
      n     = n1
      iemat = n + n2
      nmin  = n - 1

      ipivc = 1 - idim

      do main=1,n
        pivot  = zero
        ipivc  = ipivc + idim
        ipivc1 = ipivc + main - 1
        ipivc2 = ipivc + nmin
        do i=ipivc1,ipivc2
          if (abs(a(i)).gt.abs(pivot)) then
            pivot = a(i)
            lpiv  = i
          endif
        enddo
        if (pivot.eq.zero) then
          mistak =-1
          determ = deter
          print 100,main,main
          return
        endif
        icol = lpiv - ipivc + 1
        index(main) = icol
        if (icol.gt.main) then
          deter =-deter
          icol  = icol - idim
          i3 = main - idim
          do i=1,iemat
            icol    = icol + idim
            i3      = i3   + idim
            swap    = a(i3)
            a(i3)   = a(icol)
            a(icol) = swap
          enddo
        endif
        deter = deter * pivot
        pivot = one   / pivot
        i3    = ipivc + nmin
        do i=ipivc,i3
          a(i) =-a(i) * pivot
        enddo
        a(ipivc1) = pivot
        i1   = main - idim
        icol =   1  - idim
        do i=1,iemat
          icol=icol+idim
          i1=i1+idim
          if (i.ne.main) then
            jcol = icol + nmin
            swap = a(i1)
            i3   = ipivc - 1
            do i2=icol,jcol
              i3    = i3 + 1
              a(i2) = a(i2) + swap * a(i3)
            enddo
            a(i1) = swap * pivot
          endif
        enddo
      enddo

      do i1=1,n
        main = n + 1 - i1
        lpiv = index(main)
        if (lpiv.ne.main) then
          icol  = (lpiv-1) * idim + 1
          jcol  = icol + nmin
          ipivc = (main-1) * idim + 1 - icol
          do i2=icol,jcol
            i3    = i2 + ipivc
            swap  = a(i2)
            a(i2) = a(i3)
            a(i3) = swap
          enddo
        endif
      enddo

      determ = deter
      mistak = 0

      return
      end subroutine matin1
c______________________________________________________________________________
      subroutine iniderlag

c..............................................................................
c     this routine initialises the coefficients of the Lagrange derivatives   .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,tt3=3.0d0,tt4=4.0d0)
      parameter (half=0.5d0,t0p08=0.08d0,eps2=1.0d-3)
      parameter (tt8=8.0d0,tt10=10.0d0,t1p75=1.75d0)

      common /der1 / drxp(mx,mx),drxm(mx,mx),
     1               dryp(my,my),drym(my,my),
     2               drzp(mz,mz),drzm(mz,mz)
      common /der2 / dlxp(mx,mx),dlxm(mx,mx),
     1               dlyp(my,my),dlym(my,my),
     2               dlzp(mz,mz),dlzm(mz,mz)
      common /kfcl / e2,e2eff,epscl,coexv,nnx,nny,nnz,iCoul,iprintcoul

      common /nxyz / dx,dv

c         ATTENTION, this works only if mw is the largest dimension !!!!!

      dimension s(mz+mz),r(mz+mz)

c..............................................................................
      pi   = tt4*atan2(one,one)
      dx2  = dx*dx
      dv   = two*two*two * dx*dx*dx

c     .................... initialisations of the coefficients for the first
c                                                  and second order derivatives
      drxp(:,:) = zero
      drxm(:,:) = zero
      dryp(:,:) = zero
      drym(:,:) = zero
      drzp(:,:) = zero
      drzm(:,:) = zero
      dlxp(:,:) = zero
      dlxm(:,:) = zero
      dlyp(:,:) = zero
      dlym(:,:) = zero
      dlzp(:,:) = zero
      dlzm(:,:) = zero

      nx2 = mx+mx-1
      pn  = pi/(nx2+one)
      t   = pn/dx
      pp  = t
      d2c =-pi*pi*(one-one/(tt4*mx*mx))/(tt3*dx2)
      do i=1,nx2
        t    =-t
        sa   = sin(pn*i)
        s(i) = t/sa
        r(i) =-two*cos(pn*i)*pp*t/sa**2
      enddo
      do i=1,mx
        drxp(i,i) =  s(i+i-1)
        dlxp(i,i) =  d2c+r(i+i-1)
        dlxm(i,i) =  d2c-r(i+i-1)
        drxm(i,i) = -s(i+i-1)
        do m=1,mx
          if (i.ne.m) then
            im = iabs(i-m)
            si = s(im)
            if (m.lt.i) si =-si
            drxp(m,i) = si+s(i+m-1)
            drxm(m,i) = si-s(i+m-1)
            dlxp(m,i) = r(im)+r(i+m-1)
            dlxm(m,i) = r(im)-r(i+m-1)
          endif
        enddo
      enddo

      ny2 = my+my-1
      pn  = pi/(ny2+one)
      t   = pn/dx
      pp  = t
      d2c =-pi*pi*(one-one/(tt4*my*my))/(tt3*dx2)
      do i=1,ny2
        t    =-t
        sa   = sin(pn*i)
        s(i) = t/sa
        r(i) =-two*t*pp*cos(pn*i)/sa**2
      enddo
      do i=1,my
        dryp(i,i) =  s(i+i-1)
        dlyp(i,i) =  d2c+r(i+i-1)
        dlym(i,i) =  d2c-r(i+i-1)
        drym(i,i) = -s(i+i-1)
        do m=1,my
          if (i.ne.m) then
            im = iabs(i-m)
            si = s(im)
            if (m.lt.i) si =-si
            dryp(i,m) = si+s(i+m-1)
            drym(i,m) = si-s(i+m-1)
            dlyp(i,m) = r(im)+r(i+m-1)
            dlym(i,m) = r(im)-r(i+m-1)
          endif
        enddo
      enddo

      nz2 = mz+mz-1
      pn  = pi/(nz2+one)
      t   = pn/dx
      pp  = t
      d2c =-pi*pi*(one-one/(tt4*mz*mz))/(tt3*dx2)
      do i=1,nz2
        t    =-t
        sa   = sin(pn*i)
        s(i) = t/sa
        r(i) =-two*t*pp*cos(pn*i)/sa**2
      enddo
      do i=1,mz
        drzp(i,i) =  s(i+i-1)
        dlzp(i,i) =  d2c+r(i+i-1)
        dlzm(i,i) =  d2c-r(i+i-1)
        drzm(i,i) = -s(i+i-1)
        do m=1,mz
          if (i.ne.m) then
            im = iabs(i-m)
            si = s(im)
            if (m.lt.i) si =-si
            drzp(i,m) = si+s(i+m-1)
            drzm(i,m) = si-s(i+m-1)
            dlzp(i,m) = r(im)+r(i+m-1)
            dlzm(i,m) = r(im)-r(i+m-1)
          endif
        enddo
      enddo

      return
      end subroutine iniderlag

c______________________________________________________________________________
      subroutine deriv (iz)

c..............................................................................
c     calculate first x, y and z derivatives of the four real functions       .
c     representing a spinor single-particle wave function of parity iz        .
c..............................................................................
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
      end subroutine deriv

c______________________________________________________________________________
      subroutine der (w,wx,wy,wz,ix,iy,iz)

c..............................................................................
c     first x, y and z derivatives of a function "w" on the mesh              .
c     with parities "ix", "iy" and "iz"                                       .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (ca=45.0d0,cb=9.0d0,cc=60.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),wx(mx,my,mz),wy(mx,my,mz),wz(mx,my,mz)

c................................................................. x derivative
      if (ix.gt.0) then
        do j=1,my*mz
          wx(1,j,1)    = ca*(w(2,j,1)   -w(1,j,1)   )
     1                  -cb*(w(3,j,1)   -w(2,j,1)   )
     2                     +(w(4,j,1)   -w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)   -w(1,j,1)   )
     1                  -cb*(w(4,j,1)   -w(1,j,1)   )
     2                     +(w(5,j,1)   -w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)   -w(2,j,1)   )
     1                  -cb*(w(5,j,1)   -w(1,j,1)   )
     2                     +(w(6,j,1)   -w(1,j,1)   )
          wx(mx-2,j,1) = ca*(w(mx-1,j,1)-w(mx-3,j,1))
     1                  -cb*(w(mx,j,1)  -w(mx-4,j,1))
     2                     +(           -w(mx-5,j,1))
          wx(mx-1,j,1) = ca*(w(mx,j,1)  -w(mx-2,j,1))
     1                  -cb*(           -w(mx-3,j,1))
     2                     +(           -w(mx-4,j,1))
          wx(mx,j,1)   = ca*(           -w(mx-1,j,1))
     1                  -cb*(           -w(mx-2,j,1))
     2                     +(           -w(mx-3,j,1))
        enddo
      else
        do j=1,my*mz
          wx(1,j,1)    = ca*(w(2,j,1)   +w(1,j,1)   )
     1                  -cb*(w(3,j,1)   +w(2,j,1)   )
     2                     +(w(4,j,1)   +w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)   -w(1,j,1)   )
     1                  -cb*(w(4,j,1)   +w(1,j,1)   )
     2                     +(w(5,j,1)   +w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)   -w(2,j,1)   )
     1                  -cb*(w(5,j,1)   -w(1,j,1)   )
     2                     +(w(6,j,1)   +w(1,j,1)   )
          wx(mx-2,j,1) = ca*(w(mx-1,j,1)-w(mx-3,j,1))
     1                  -cb*(w(mx,j,1)  -w(mx-4,j,1))
     2                     +(           -w(mx-5,j,1))
          wx(mx-1,j,1) = ca*(w(mx,j,1)  -w(mx-2,j,1))
     1                  -cb*(           -w(mx-3,j,1))
     2                     +(           -w(mx-4,j,1))
          wx(mx,j,1)   = ca*(           -w(mx-1,j,1))
     1                  -cb*(           -w(mx-2,j,1))
     2                     +(           -w(mx-3,j,1))
        enddo
      endif

      do i=4,mx-3
      do j=1,my*mz
        wx(i,j,1)    = ca*(w(i+1,j,1) -w(i-1,j,1) )
     1                -cb*(w(i+2,j,1) -w(i-2,j,1) )
     2                   +(w(i+3,j,1) -w(i-3,j,1) )
      enddo
      enddo

c     ............................................................ y derivative
      if (iy.gt.0) then
        do k=1,mz
        do i=1,mx
          wy(i,1,k)    = ca*(w(i,2,k)   -w(i,1,k)   )
     1                  -cb*(w(i,3,k)   -w(i,2,k)   )
     2                     +(w(i,4,k)   -w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)   -w(i,1,k)   )
     1                  -cb*(w(i,4,k)   -w(i,1,k)   )
     2                     +(w(i,5,k)   -w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)   -w(i,2,k)   )
     1                  -cb*(w(i,5,k)   -w(i,1,k)   )
     2                     +(w(i,6,k)   -w(i,1,k)   )
          wy(i,my-2,k) = ca*(w(i,my-1,k)-w(i,my-3,k))
     1                  -cb*(w(i,my,k)  -w(i,my-4,k))
     2                     +(           -w(i,my-5,k))
          wy(i,my-1,k) = ca*(w(i,my,k)  -w(i,my-2,k))
     1                  -cb*(           -w(i,my-3,k))
     2                     +(           -w(i,my-4,k))
          wy(i,my,k)   = ca*(           -w(i,my-1,k))
     1                  -cb*(           -w(i,my-2,k))
     2                     +(           -w(i,my-3,k))
        enddo
        enddo
      else
        do k=1,mz
        do i=1,mx
          wy(i,1,k)    = ca*(w(i,2,k)   +w(i,1,k)   )
     1                  -cb*(w(i,3,k)   +w(i,2,k)   )
     2                     +(w(i,4,k)   +w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)   -w(i,1,k)   )
     1                  -cb*(w(i,4,k)   +w(i,1,k)   )
     2                     +(w(i,5,k)   +w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)   -w(i,2,k)   )
     1                  -cb*(w(i,5,k)   -w(i,1,k)   )
     2                     +(w(i,6,k)   +w(i,1,k)   )
          wy(i,my-2,k) = ca*(w(i,my-1,k)-w(i,my-3,k))
     1                  -cb*(w(i,my,k)  -w(i,my-4,k))
     2                     +(           -w(i,my-5,k))
          wy(i,my-1,k) = ca*(w(i,my,k)  -w(i,my-2,k))
     1                  -cb*(           -w(i,my-3,k))
     2                     +(           -w(i,my-4,k))
          wy(i,my,k)   = ca*(           -w(i,my-1,k))
     1                  -cb*(           -w(i,my-2,k))
     2                     +(           -w(i,my-3,k))
        enddo
        enddo
      endif

      do j=4,my-3
      do i=1,mx
      do k=1,mz
        wy(i,j,k)    = ca*(w(i,j+1,k) -w(i,j-1,k) )
     1                -cb*(w(i,j+2,k) -w(i,j-2,k) )
     2                   +(w(i,j+3,k) -w(i,j-3,k) )
      enddo
      enddo
      enddo

c     ............................................................ z derivative
      if (iz.gt.0) then
        do i=1,mx*my
          wz(i,1,1)    = ca*(w(i,1,2)   -w(i,1,1)   )
     1                  -cb*(w(i,1,3)   -w(i,1,2)   )
     2                     +(w(i,1,4)   -w(i,1,3)   )
          wz(i,1,2)    = ca*(w(i,1,3)   -w(i,1,1)   )
     1                  -cb*(w(i,1,4)   -w(i,1,1)   )
     2                     +(w(i,1,5)   -w(i,1,2)   )
          wz(i,1,3)    = ca*(w(i,1,4)   -w(i,1,2)   )
     1                  -cb*(w(i,1,5)   -w(i,1,1)   )
     2                     +(w(i,1,6)   -w(i,1,1)   )
          wz(i,1,mz-2) = ca*(w(i,1,mz-1)-w(i,1,mz-3))
     1                  -cb*(w(i,1,mz)  -w(i,1,mz-4))
     2                     +(           -w(i,1,mz-5))
          wz(i,1,mz-1) = ca*(w(i,1,mz)  -w(i,1,mz-2))
     1                  -cb*(           -w(i,1,mz-3))
     1                     +(           -w(i,1,mz-4))
          wz(i,1,mz)   = ca*(           -w(i,1,mz-1))
     1                  -cb*(           -w(i,1,mz-2))
     2                     +(           -w(i,1,mz-3))
        enddo
      else
        do i=1,mx*my
          wz(i,1,1)    = ca*(w(i,1,2)   +w(i,1,1)   )
     1                  -cb*(w(i,1,3)   +w(i,1,2)   )
     2                     +(w(i,1,4)   +w(i,1,3)   )
          wz(i,1,2)    = ca*(w(i,1,3)   -w(i,1,1)   )
     1                  -cb*(w(i,1,4)   +w(i,1,1)   )
     2                     +(w(i,1,5)   +w(i,1,2)   )
          wz(i,1,3)    = ca*(w(i,1,4)   -w(i,1,2)   )
     1                  -cb*(w(i,1,5)   -w(i,1,1)   )
     2                     +(w(i,1,6)   +w(i,1,1)   )
          wz(i,1,mz-2) = ca*(w(i,1,mz-1)-w(i,1,mz-3))
     1                  -cb*(w(i,1,mz)  -w(i,1,mz-4))
     2                     +(           -w(i,1,mz-5))
          wz(i,1,mz-1) = ca*(w(i,1,mz)  -w(i,1,mz-2))
     1                  -cb*(           -w(i,1,mz-3))
     1                     +(           -w(i,1,mz-4))
          wz(i,1,mz)   = ca*(           -w(i,1,mz-1))
     1                  -cb*(           -w(i,1,mz-2))
     2                     +(           -w(i,1,mz-3))
        enddo
      endif

      do k=4,mz-3
      do i=1,mx*my
        wz(i,1,k)    = ca*(w(i,1,k+1) -w(i,1,k-1) )
     1                -cb*(w(i,1,k+2) -w(i,1,k-2) )
     2                   +(w(i,1,k+3) -w(i,1,k-3) )
      enddo
      enddo

      do i=1,mv
        wx(i,1,1) = wx(i,1,1) / (cc*dx)
        wy(i,1,1) = wy(i,1,1) / (cc*dx)
        wz(i,1,1) = wz(i,1,1) / (cc*dx)
      enddo

      return
      end subroutine der

c______________________________________________________________________________
      subroutine derx (w,wx,ix)

c..............................................................................
c     x derivative
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (ca=45.0d0,cb=9.0d0,cc=60.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),wx(mx,my,mz)

c..............................................................................
      if (ix.gt.0) then
        do j=1,my*mz
          wx(1,j,1)    = ca*(w(2,j,1)    - w(1,j,1)   )
     1                  -cb*(w(3,j,1)    - w(2,j,1)   )
     2                     +(w(4,j,1)    - w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)    - w(1,j,1)   )
     1                  -cb*(w(4,j,1)    - w(1,j,1)   )
     2                     +(w(5,j,1)    - w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)    - w(2,j,1)   )
     1                  -cb*(w(5,j,1)    - w(1,j,1)   )
     2                     +(w(6,j,1)    - w(1,j,1)   )
          wx(mx-2,j,1) = ca*(w(mx-1,j,1) - w(mx-3,j,1))
     1                  -cb*(w(mx,j,1)   - w(mx-4,j,1))
     2                     +(            - w(mx-5,j,1))
          wx(mx-1,j,1) = ca*(w(mx,j,1)   - w(mx-2,j,1))
     1                  -cb*(            - w(mx-3,j,1))
     2                     +(            - w(mx-4,j,1))
          wx(mx,j,1)   = ca*(            - w(mx-1,j,1))
     1                  -cb*(            - w(mx-2,j,1))
     2                     +(            - w(mx-3,j,1))
        enddo
      else
        do j=1,my*mz
          wx(1,j,1)    = ca*(w(2,j,1)    + w(1,j,1)   )
     1                  -cb*(w(3,j,1)    + w(2,j,1)   )
     2                     +(w(4,j,1)    + w(3,j,1)   )
          wx(2,j,1)    = ca*(w(3,j,1)    - w(1,j,1)   )
     1                  -cb*(w(4,j,1)    + w(1,j,1)   )
     2                     +(w(5,j,1)    + w(2,j,1)   )
          wx(3,j,1)    = ca*(w(4,j,1)    - w(2,j,1)   )
     1                  -cb*(w(5,j,1)    - w(1,j,1)   )
     2                     +(w(6,j,1)    + w(1,j,1)   )
          wx(mx-2,j,1) = ca*(w(mx-1,j,1) - w(mx-3,j,1))
     1                  -cb*(w(mx,j,1)   - w(mx-4,j,1))
     2                     +(            - w(mx-5,j,1))
          wx(mx-1,j,1) = ca*(w(mx,j,1)   - w(mx-2,j,1))
     1                  -cb*(            - w(mx-3,j,1))
     2                     +(            - w(mx-4,j,1))
          wx(mx,j,1)   = ca*(            - w(mx-1,j,1))
     1                  -cb*(            - w(mx-2,j,1))
     2                     +(            - w(mx-3,j,1))
        enddo
      endif

      do i=4,mx-3
      do j=1,my*mz
        wx(i,j,1)    = ca*(w(i+1,j,1) - w(i-1,j,1) )
     1                -cb*(w(i+2,j,1) - w(i-2,j,1) )
     2                   +(w(i+3,j,1) - w(i-3,j,1) )
      enddo
      enddo

      do i=1,mv
        wx(i,1,1) = wx(i,1,1) / (cc*dx)
      enddo

      return
      end subroutine derx

c______________________________________________________________________________
      subroutine dery (w,wy,iy)

c..............................................................................
c     y derivative
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (ca=45.0d0,cb=9.0d0,cc=60.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),wy(mx,my,mz)

c..............................................................................
      if (iy.gt.0) then
        do k=1,mz
        do i=1,mx
          wy(i,1,k)    = ca*(w(i,2,k)    - w(i,1,k)   )
     1                  -cb*(w(i,3,k)    - w(i,2,k)   )
     2                     +(w(i,4,k)    - w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)    - w(i,1,k)   )
     1                  -cb*(w(i,4,k)    - w(i,1,k)   )
     2                     +(w(i,5,k)    - w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)    - w(i,2,k)   )
     1                  -cb*(w(i,5,k)    - w(i,1,k)   )
     2                     +(w(i,6,k)    - w(i,1,k)   )
          wy(i,my-2,k) = ca*(w(i,my-1,k) - w(i,my-3,k))
     1                  -cb*(w(i,my,k)   - w(i,my-4,k))
     2                     +(            - w(i,my-5,k))
          wy(i,my-1,k) = ca*(w(i,my,k)   - w(i,my-2,k))
     1                  -cb*(            - w(i,my-3,k))
     2                     +(            - w(i,my-4,k))
          wy(i,my,k)   = ca*(            - w(i,my-1,k))
     1                  -cb*(            - w(i,my-2,k))
     2                     +(            - w(i,my-3,k))
        enddo
        enddo
      else
        do k=1,mz
        do i=1,mx
          wy(i,1,k)    = ca*(w(i,2,k)    + w(i,1,k)   )
     1                  -cb*(w(i,3,k)    + w(i,2,k)   )
     2                     +(w(i,4,k)    + w(i,3,k)   )
          wy(i,2,k)    = ca*(w(i,3,k)    - w(i,1,k)   )
     1                  -cb*(w(i,4,k)    + w(i,1,k)   )
     2                     +(w(i,5,k)    + w(i,2,k)   )
          wy(i,3,k)    = ca*(w(i,4,k)    - w(i,2,k)   )
     1                  -cb*(w(i,5,k)    - w(i,1,k)   )
     2                     +(w(i,6,k)    + w(i,1,k)   )
          wy(i,my-2,k) = ca*(w(i,my-1,k) - w(i,my-3,k))
     1                  -cb*(w(i,my,k)   - w(i,my-4,k))
     2                     +(            - w(i,my-5,k))
          wy(i,my-1,k) = ca*(w(i,my,k)   - w(i,my-2,k))
     1                  -cb*(            - w(i,my-3,k))
     2                     +(            - w(i,my-4,k))
          wy(i,my,k)   = ca*(            - w(i,my-1,k))
     1                  -cb*(            - w(i,my-2,k))
     2                     +(            - w(i,my-3,k))
        enddo
        enddo
      endif

      do k=1,mz
      do j=4,my-3
      do i=1,mx
        wy(i,j,k)    = ca*(w(i,j+1,k) - w(i,j-1,k) )
     1                -cb*(w(i,j+2,k) - w(i,j-2,k) )
     2                   +(w(i,j+3,k) - w(i,j-3,k) )
      enddo
      enddo
      enddo

      do i=1,mv
        wy(i,1,1) = wy(i,1,1) / (cc*dx)
      enddo

      return
      end subroutine dery

c______________________________________________________________________________
      subroutine derz (w,wz,iz)

c..............................................................................
c     z derivative
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (ca=45.0d0,cb=9.0d0,cc=60.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),wz(mx,my,mz)

c..............................................................................
      if (iz.gt.0) then
        do i=1,mx*my
          wz(i,1,1)    = ca*(w(i,1,2)    - w(i,1,1)   )
     1                  -cb*(w(i,1,3)    - w(i,1,2)   )
     2                     +(w(i,1,4)    - w(i,1,3)   )
          wz(i,1,2)    = ca*(w(i,1,3)    - w(i,1,1)   )
     1                  -cb*(w(i,1,4)    - w(i,1,1)   )
     2                     +(w(i,1,5)    - w(i,1,2)   )
          wz(i,1,3)    = ca*(w(i,1,4)    - w(i,1,2)   )
     1                  -cb*(w(i,1,5)    - w(i,1,1)   )
     2                     +(w(i,1,6)    - w(i,1,1)   )
          wz(i,1,mz-2) = ca*(w(i,1,mz-1) - w(i,1,mz-3))
     1                  -cb*(w(i,1,mz)   - w(i,1,mz-4))
     2                     +(            - w(i,1,mz-5))
          wz(i,1,mz-1) = ca*(w(i,1,mz)   - w(i,1,mz-2))
     1                  -cb*(            - w(i,1,mz-3))
     1                     +(            - w(i,1,mz-4))
          wz(i,1,mz)   = ca*(            - w(i,1,mz-1))
     1                  -cb*(            - w(i,1,mz-2))
     2                     +(            - w(i,1,mz-3))
        enddo
      else
        do i=1,mx*my
          wz(i,1,1)    = ca*(w(i,1,2)    + w(i,1,1)   )
     1                  -cb*(w(i,1,3)    + w(i,1,2)   )
     2                     +(w(i,1,4)    + w(i,1,3)   )
          wz(i,1,2)    = ca*(w(i,1,3)    - w(i,1,1)   )
     1                  -cb*(w(i,1,4)    + w(i,1,1)   )
     2                     +(w(i,1,5)    + w(i,1,2)   )
          wz(i,1,3)    = ca*(w(i,1,4)    - w(i,1,2)   )
     1                  -cb*(w(i,1,5)    - w(i,1,1)   )
     2                     +(w(i,1,6)    + w(i,1,1)   )
          wz(i,1,mz-2) = ca*(w(i,1,mz-1) - w(i,1,mz-3))
     1                  -cb*(w(i,1,mz)   - w(i,1,mz-4))
     2                     +(            - w(i,1,mz-5))
          wz(i,1,mz-1) = ca*(w(i,1,mz)   - w(i,1,mz-2))
     1                  -cb*(            - w(i,1,mz-3))
     1                     +(            - w(i,1,mz-4))
          wz(i,1,mz)   = ca*(            - w(i,1,mz-1))
     1                  -cb*(            - w(i,1,mz-2))
     2                     +(            - w(i,1,mz-3))
        enddo
      endif

      do k=4,mz-3
      do i=1,mx*my
        wz(i,1,k)    = ca*(w(i,1,k+1) - w(i,1,k-1) )
     1                -cb*(w(i,1,k+2) - w(i,1,k-2) )
     2                   +(w(i,1,k+3) - w(i,1,k-3) )
      enddo
      enddo

      do i=1,mv
        wz(i,1,1) = wz(i,1,1) / (cc*dx)
      enddo

      return
      end subroutine derz

c______________________________________________________________________________
      subroutine lapla (w,dw,xp,yp,zp)

c..............................................................................
c     calculation of the Laplacian of the four real functions representing a  .
c     spinor wave function of parity zp                                       .
c..............................................................................
      implicit real*8 (a-h,o-z)
      include 'param8.h'

      parameter (xxf=-8064.0d0,xxm=43050.0d0,xx2=1008.0d0)
      parameter (xx3=-128.0d0,xx4=9.0d0,xxd=5040.0d0)

      common /nxyz / dx,dv
      dimension w(mx,my,mz),dw(mx,my,mz)
      dimension w_t(my,mx,mz),d_T(my,mx,mz)

c..............................................................................
      xf = xxf                  ! = -8064.0d0
      xm = xxm / xf             ! = 43050.0d0 / xf
      x2 = xx2 / xf             ! =  1008.0d0 / xf
      x3 = xx3 / xf             ! =  -128.0d0 / xf
      x4 = xx4 / xf             ! =     9.0d0 / xf
      xd = xxd / xf * (dx*dx)   ! =  5040.0d0 / xf * (dx*dx)

c     .........................................................................
      do i=1,mv
        dw(i,1,1) = xm * w(i,1,1)
      enddo

c     ................................................................. x sweep
      do j=1,my*mz
        dw(1,j,1)    = dw(1,j,1)   +   w(1,j,1)*xp+x2*w(2,j,1)*xp
     1                             +x3*w(3,j,1)*xp+x4*w(4,j,1)*xp
        dw(2,j,1)    = dw(2,j,1)   +   w(1,j,1)   +x2*w(1,j,1)*xp
     1                             +x3*w(2,j,1)*xp+x4*w(3,j,1)*xp
        dw(3,j,1)    = dw(3,j,1)   +   w(2,j,1)   +x2*w(1,j,1)
     1                             +x3*w(1,j,1)*xp+x4*w(2,j,1)*xp
        dw(4,j,1)    = dw(4,j,1)   +   w(3,j,1)   +x2*w(2,j,1)
     1                             +x3*w(1,j,1)   +x4*w(1,j,1)*xp
        dw(mx-1,j,1) = dw(mx-1,j,1)+   w(mx,j,1)
        dw(mx-2,j,1) = dw(mx-2,j,1)+   w(mx-1,j,1)+x2*w(mx,j,1)
        dw(mx-3,j,1) = dw(mx-3,j,1)+   w(mx-2,j,1)+x2*w(mx-1,j,1)
     1                             +x3*w(mx,j,1)
      enddo
      do j=1,my*mz
      do i=1,mx-4
         dw(i,j,1) = dw(i,j,1)   +   w(i+1,j,1) +x2*w(i+2,j,1)
     1                           +x3*w(i+3,j,1) +x4*w(i+4,j,1)
      enddo
      enddo
      do j=1,my*mz
      do i=5,mx
         dw(i,j,1) = dw(i,j,1)   +   w(i-1,j,1) +x2*w(i-2,j,1)
     1                           +x3*w(i-3,j,1) +x4*w(i-4,j,1)
      enddo
      enddo

c     ................................................................. y sweep
c       transposition in y dimension in order to have vector lengh set to mx*mz
      do i=1,mx
      do j=1,my
      do k=1,mz
        w_t(j,i,k) =  w(i,j,k)
        d_t(j,i,k) = dw(i,j,k)
      enddo
      enddo
      enddo
      do k=1,mx*mz
        d_t(1,k,1)    = d_t(1,k,1)   +   w_t(1,k,1)*yp+x2*w_t(2,k,1)*yp
     1                               +x3*w_t(3,k,1)*yp+x4*w_t(4,k,1)*yp
        d_t(2,k,1)    = d_t(2,k,1)   +   w_t(1,k,1)   +x2*w_t(1,k,1)*yp
     1                               +x3*w_t(2,k,1)*yp+x4*w_t(3,k,1)*yp
        d_t(3,k,1)    = d_t(3,k,1)   +   w_t(2,k,1)   +x2*w_t(1,k,1)
     1                               +x3*w_t(1,k,1)*yp+x4*w_t(2,k,1)*yp
        d_t(4,k,1)    = d_t(4,k,1)   +   w_t(3,k,1)   +x2*w_t(2,k,1)
     1                               +x3*w_t(1,k,1)   +x4*w_t(1,k,1)*yp
!        d_t(mx-1,k,1) = d_t(mx-1,k,1)+   w_t(mx,k,1)
!        d_t(mx-2,k,1) = d_t(mx-2,k,1)+   w_t(mx-1,k,1)+x2*w_t(mx,k,1)
!        d_t(mx-3,k,1) = d_t(mx-3,k,1)+   w_t(mx-2,k,1)+x2*w_t(mx-1,k,1)
!     1                               +x3*w_t(mx,k,1)
        d_t(my-1,k,1) = d_t(my-1,k,1)+   w_t(my,k,1)
        d_t(my-2,k,1) = d_t(my-2,k,1)+   w_t(my-1,k,1)+x2*w_t(my,k,1)
        d_t(my-3,k,1) = d_t(my-3,k,1)+   w_t(my-2,k,1)+x2*w_t(my-1,k,1)
     1                               +x3*w_t(my,k,1)
      enddo
      do k=1,mx*mz
      do j=1,my-4
        d_t(j,k,1) = d_t(j,k,1)   +   w_t(j+1,k,1) +x2*w_t(j+2,k,1)
     1                            +x3*w_t(j+3,k,1) +x4*w_t(j+4,k,1)
      enddo
      enddo
      do k=1,mx*mz
      do j=5,my
        d_t(j,k,1) = d_t(j,k,1)   +   w_t(j-1,k,1) +x2*w_t(j-2,k,1)
     1                            +x3*w_t(j-3,k,1) +x4*w_t(j-4,k,1)
      enddo
      enddo
      do i=1,mx
      do j=1,my
      do k=1,mz
        dw(i,j,k) = d_t(j,i,k)
      enddo
      enddo
      enddo

c     ................................................................. z sweep
      do i=1,mx*my
        dw(i,1,1)    = dw(i,1,1)   +   w(i,1,1)*zp+x2*w(i,1,2)*zp
     1                             +x3*w(i,1,3)*zp+x4*w(i,1,4)*zp
        dw(i,1,2)    = dw(i,1,2)   +   w(i,1,1)   +x2*w(i,1,1)*zp
     1                             +x3*w(i,1,2)*zp+x4*w(i,1,3)*zp
        dw(i,1,3)    = dw(i,1,3)   +   w(i,1,2)   +x2*w(i,1,1)
     1                             +x3*w(i,1,1)*zp+x4*w(i,1,2)*zp
        dw(i,1,4)    = dw(i,1,4)   +   w(i,1,3)   +x2*w(i,1,2)
     1                             +x3*w(i,1,1)   +x4*w(i,1,1)*zp
        dw(i,1,mz-1) = dw(i,1,mz-1)+   w(i,1,mz)
        dw(i,1,mz-2) = dw(i,1,mz-2)+   w(i,1,mz-1)+x2*w(i,1,mz)
        dw(i,1,mz-3) = dw(i,1,mz-3)+   w(i,1,mz-2)+x2*w(i,1,mz-1)
     1                             +x3*w(i,1,mz)
      enddo
      do k=1,mz-4
      do i=1,mx*my
        dw(i,1,k) = dw(i,1,k)   +   w(i,1,k+1) +x2*w(i,1,k+2)
     1                          +x3*w(i,1,k+3) +x4*w(i,1,k+4)
      enddo
      enddo
      do k=5,mz
      do i=1,mx*my
        dw(i,1,k) = dw(i,1,k)   +   w(i,1,k-1) +x2*w(i,1,k-2)
     1                          +x3*w(i,1,k-3) +x4*w(i,1,k-4)
      enddo
      enddo

      do i=1,mv
        dw(i,1,1) = - dw(i,1,1) / xd
      enddo

      return
      end subroutine lapla

c______________________________________________________________________________
      subroutine derlagx (ipar,f,df)

c..............................................................................
c     calculation of the x derivative of an arbitrary function on the mesh
c     using lagrange mesh derivatives.
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      parameter (myz=my*mz)

      common /der1 / drxp(mx,mx),drxm(mx,mx),
     1               dryp(my,my),drym(my,my),
     2               drzp(mz,mz),drzm(mz,mz)

      dimension f(mv),df(mv)

c..............................................................................
      if (ipar.ne.1.and.ipar.ne.-1) call stp(' derx : ipar!')

      if (ipar.eq.1) then
        call mxm (drxp(1,1),mx,f(1),mx,df(1),myz)
      else
        call mxm (drxm(1,1),mx,f(1),mx,df(1),myz)
      endif

      return
      end subroutine derlagx

c______________________________________________________________________________
      subroutine derlagy (ipar,f,df)

c..............................................................................
c     calculation of the y derivative of an arbitrary function on the mesh    .
c     simplifies the handling of y derivatives as they require a              .
c     representation of the function to be derived as array of (x,y,z)        .
c     while x and z derivatives do not.
c     Using lagrange mesh derivatives!
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'

      common /der1 / drxp(mx,mx),drxm(mx,mx),
     1               dryp(my,my),drym(my,my),
     2               drzp(mz,mz),drzm(mz,mz)

      dimension f(mx,my,mz),df(mx,my,mz)

c..............................................................................
      if (ipar.ne.1.and.ipar.ne.-1) call stp(' dery : ipar!')

      if (ipar.eq.1) then
        do k=1,mz
          call mxm(f(1,1,k),mx,dryp(1,1),my,df(1,1,k),my)
        enddo
      else
        do k=1,mz
          call mxm(f(1,1,k),mx,drym(1,1),my,df(1,1,k),my)
        enddo
      endif

      return
      end subroutine derlagy

c______________________________________________________________________________
      subroutine derlagz (ipar,f,df)

c..............................................................................
c     calculation of the z derivative of an arbitrary function on the mesh    .
c     using lagrange mesh derivatives.
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      parameter (mxy=mx*my)

      common /der1 / drxp(mx,mx),drxm(mx,mx),
     1               dryp(my,my),drym(my,my),
     2               drzp(mz,mz),drzm(mz,mz)

      dimension f(mv),df(mv)

c..............................................................................
      if (ipar.ne.1.and.ipar.ne.-1) call stp(' derz : ipar!')

      if (ipar.eq.1) then
        call mxm (f(1),mxy,drzp(1,1),mz,df(1),mz)
      else
        call mxm (f(1),mxy,drzm(1,1),mz,df(1),mz)
      endif

      return
      end subroutine derlagz

c______________________________________________________________________________
      subroutine laplalag (w,d2w,ixp,iyp,izp)

c..............................................................................
c     calculation of the Laplacian of the function "w" with the Lagrange mesh .
c     derivatives.                                                            .
c         in  : w with parities ixp,iyp,izp                                   .
c         out : d2w                                                           .
c     This routine uses lagrange derivatives.                                 .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      parameter (one=1.0d0)
      parameter (mxy=mx*my,myz=my*mz)

      common /der2 / dlxp(mx,mx),dlxm(mx,mx),
     1               dlyp(my,my),dlym(my,my),
     2               dlzp(mz,mz),dlzm(mz,mz)
      dimension qlx(mx,my,mz),qly(mx,my,mz),qlz(mx,my,mz)
      dimension w(mx,my,mz),d2w(mx,my,mz)

c..............................................................................
      if (ixp.ne.1.and.ixp.ne.-1) call stp(' lapla xp !')
      if (iyp.ne.1.and.iyp.ne.-1) call stp(' lapla yp !')
      if (izp.ne.1.and.izp.ne.-1) call stp(' lapla zp !')

      if (ixp.eq.1) then
        call mxm(dlxp(1,1),mx,w(1,1,1),mx,qlx(1,1,1),myz)
      else
        call mxm(dlxm(1,1),mx,w(1,1,1),mx,qlx(1,1,1),myz)
      endif

      if (iyp.eq.1) then
        do k=1,mz
          call mxm(w(1,1,k),mx,dlyp(1,1),my,qly(1,1,k),my)
        enddo
      else
        do k=1,mz
          call mxm(w(1,1,k),mx,dlym(1,1),my,qly(1,1,k),my)
        enddo
      endif

      if (izp.eq.1) then
        call mxm(w(1,1,1),mxy,dlzp(1,1),mz,qlz(1,1,1),mz)
      else
        call mxm(w(1,1,1),mxy,dlzm(1,1),mz,qlz(1,1,1),mz)
      endif

      do i=1,mv
        d2w(i,1,1) = qlx(i,1,1) + qly(i,1,1) + qlz(i,1,1)
      enddo

      return
      end subroutine laplalag

c______________________________________________________________________________
      subroutine derivlag (ii)

c..............................................................................
c     calculation of the first derivatives of a single-particle wave function .
c     with the the exact Lagrange mesh derivatives.                           .
c     information handling through common blocks:                             .
c       in  : w1,...,w4                                                       .
c       out : wx1,...,wz4                                                     .
c     index iw of single-particle wave function needed to determine parity    .
c     the parity kparz (+ or -) determines the 3 plane symmetries:            .
c     (Note: all single-particle states have signature +1)                    .
c            1  (+,+, kparz)                                                  .
c            2  (-,-, kparz)                                                  .
c            3  (-,+,-kparz)                                                  .
c            4  (+,-,-kparz)                                                  .
c..............................................................................

      implicit real*8 (a-h,o-z)
      include 'param8.h'
      parameter (mxy=mx*my,myz=my*mz)

      common /der1 / drxp(mx,mx),drxm(mx,mx),
     1               dryp(my,my),drym(my,my),
     2               drzp(mz,mz),drzm(mz,mz)
      common /wave / w1(mx,my,mz),w2(mx,my,mz),
     1               w3(mx,my,mz),w4(mx,my,mz),
     2               psi1(mv),psi2(mv),psi3(mv),psi4(mv)
      common /waved/ wx1(mx,my,mz),wx2(mx,my,mz),
     1               wx3(mx,my,mz),wx4(mx,my,mz),
     2               wy1(mx,my,mz),wy2(mx,my,mz),
     3               wy3(mx,my,mz),wy4(mx,my,mz),
     4               wz1(mx,my,mz),wz2(mx,my,mz),
     5               wz3(mx,my,mz),wz4(mx,my,mz)
      common /spwf / esp1(mw),esp2(mw),esp3(mw),v2(mw),v22(mw),eqp(mw)
     1              ,delta(mw),ajzd(mw),aj2d(mw),kparz(mw),kiso(mw)

c..............................................................................
      do k=1,mz
        call mxm(w1(1,1,k),mx,dryp(1,1),my,wy1(1,1,k),my)
        call mxm(w2(1,1,k),mx,drym(1,1),my,wy2(1,1,k),my)
        call mxm(w3(1,1,k),mx,dryp(1,1),my,wy3(1,1,k),my)
        call mxm(w4(1,1,k),mx,drym(1,1),my,wy4(1,1,k),my)
      enddo
      call mxm(drxp(1,1),mx,w1(1,1,1),mx,wx1(1,1,1),myz)
      call mxm(drxm(1,1),mx,w2(1,1,1),mx,wx2(1,1,1),myz)
      call mxm(drxm(1,1),mx,w3(1,1,1),mx,wx3(1,1,1),myz)
      call mxm(drxp(1,1),mx,w4(1,1,1),mx,wx4(1,1,1),myz)
      if (kparz(ii).eq.1) then
        call mxm(w1(1,1,1),mxy,drzp(1,1),mz,wz1(1,1,1),mz)
        call mxm(w2(1,1,1),mxy,drzp(1,1),mz,wz2(1,1,1),mz)
        call mxm(w3(1,1,1),mxy,drzm(1,1),mz,wz3(1,1,1),mz)
        call mxm(w4(1,1,1),mxy,drzm(1,1),mz,wz4(1,1,1),mz)
      else
        call mxm(w1(1,1,1),mxy,drzm(1,1),mz,wz1(1,1,1),mz)
        call mxm(w2(1,1,1),mxy,drzm(1,1),mz,wz2(1,1,1),mz)
        call mxm(w3(1,1,1),mxy,drzp(1,1),mz,wz3(1,1,1),mz)
        call mxm(w4(1,1,1),mxy,drzp(1,1),mz,wz4(1,1,1),mz)
      endif

      return
      end subroutine derivlag

c______________________________________________________________________________
      function ssum (n,a)

      double precision zero,a(n),ssum

      parameter (zero=0.0d0)

      ssum = zero
      do i=1,n
        ssum = ssum + a(i)
      enddo

      return
      end function ssum

c______________________________________________________________________________
      function sdot (n,a,b)

      double precision zero,a(n),b(n),sdot

      parameter (zero=0.0d0)

      sdot = zero
      do i=1,n
        sdot = sdot + a(i) * b(i)
      enddo

      return
      end function sdot

c______________________________________________________________________________
      subroutine sswap (n,a,b)

      double precision a(n),b(n),tmp

      do i=1,n
        tmp  = b(i)
        b(i) = a(i)
        a(i) = tmp
      enddo

      return
      end subroutine sswap

c______________________________________________________________________________
      subroutine scopy (n,a,b)

      double precision a(n),b(n)

      do i=1, n
        b(i) = a(i)
      enddo

      return
      end subroutine scopy

c______________________________________________________________________________
      subroutine sscal (n,fac,a)

      double precision a(n),fac

      do i=1, n
        a(i) = fac * a(i)
      enddo

      return
      end subroutine sscal

c______________________________________________________________________________
      subroutine saxpy (n,fac,a,b)

      double precision a(n),b(n),fac

      do i=1, n
        b(i) = b(i) + fac * a(i)
      enddo

      return
      end subroutine saxpy

c______________________________________________________________________________
      subroutine stp (message)

      character (len = *) :: message

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

      subroutine mxm (a,n1,b,ks,c,n2)

      implicit none
      double precision a(n1,ks),b(ks,n2),c(n1,n2)
      integer n1,n2,i,i1,i2,k,ks

c..............................................................................
      do i=1,n1*n2
        c(i,1) = 0.0d0
      enddo

      do i2=1,n2
      do k=1,ks
      do i1=1,n1
        c(i1,i2) = c(i1,i2) + a(i1,k) * b(k,i2)
      enddo
      enddo
      enddo

      return
      end subroutine mxm

c______________________________________________________________________________
      subroutine to_upper (str, string)
      
        character(*)        :: str
        character(len(str)) :: string

        Integer :: ic, i
        !Ugly but effective and independent of platform and implementation.
        character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
        
        if(len(str) .ne. len(string)) then
          call stp('Strings of different length in to_upper')
        endif
        string = str
        do i = 1, len_trim(str)
          ic = INDEX(low, str(i:i)) !Note that ic = 0 when substring is not found
          if (ic > 0) then
            string(i:i) = cap(ic:ic)
          else
            string(i:i) = str(i:i)
          endif
        end do
      
      end subroutine to_upper
