PROGRAM RoughDiff

  USE mpi
  IMPLICIT NONE

  CHARACTER(200) :: wdir
  REAL*8 :: sigma0,sigYnow,Lo,alpha,eps ! Macroscopic load and scale
  REAL*8 :: Hu,qmax,qL,dz,mexp,nexp,Co,alphai,betai, &
       Anon,Apl,Aplpr,Apl2,Anonpr,dimfac,zetan,zetan2,fzetan,fzetan2,dzetan,dzetany,Al,Ael,eta,zprev
  REAL*8 :: dsig,dzeta,dmu,i2dsig,idsig,sigYmax,sigY0,Pszm,zetamax,fzmax,fzmin,dzetamax,dzetamin
  INTEGER :: nz,nsig,nsigmax,nspo,nsp,ierr, si, zi, myrank, NPU, nznpu,soindex
  
  INTEGER :: ntot,sni,sne,ni,ntotsub,nintz,next,Nred
  REAL*8  :: myni,ni2,dmuint,sfzint,rr2,Fel,Felall,Fpl,Ftot,ftilde1,f1
  REAL*8  :: s0sy, sns0, npsy, mysign,mypref,mysns0,mynpsy,nexpo,xx,zetarep,zetarep2
  REAL*8  :: ddsigy, idsigy, i2dsigy
  REAL*8, DIMENSION(:), ALLOCATABLE :: Pz,Pnon,Ppl
  
  INTEGER :: outstepsize,noutput,noutline
  REAL*8 :: zoutfac,zetaoprev,rr,bsigY,asigY,dsigY,sigYp,sigYp2,dApldz,dApldzold,sig,myend
  REAL*8 :: asigYm1, bsigYm1, csigYm1, dsigYm1
  REAL*8, DIMENSION(:), ALLOCATABLE :: Cq,zetas,mus,sigY,fz,sfz
  REAL*8, DIMENSION(:), ALLOCATABLE :: dpdzo,dpdzo2,Psz,Pszpr,PszB,PszE,Pszprall,Pzsum, &
       Pplsum,Pnonsum,cca,ccc,d
  INTEGER :: syindex,syrank,syrankm1,syrankall,syindexhome
  REAL*8 :: Ey, nu, pi, Eyp
  REAL*8 :: Pny, Pny1, Pny2, zoutmod
  INTEGER :: arraySize,nsize
  REAL*8 :: startTime,currTime,dTime,tfin
  
  PARAMETER(Ey = 1.d11, &  ! Young's modulus
       nu = 0.25)          ! Poisson's ratio
       
  PARAMETER (pi=3.141592653589793d0)
  
  ! ---------------------------------------------------
  ! Initialize MPI
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,nsize,ierr)

  startTime = MPI_Wtime()
  tfin = 23.d0*60.d0*60.d0 ! time to output backup of P distribution
!  tfin = 5.d0*60.d0 
  NPU = nsize ! Number of cores
  Nred = 2*NPU
  ! Discretization in stress space should be a multiple of NPU

  CALL readinput(wdir)
  CALL FLUSH(6)
  IF(myrank.EQ.0)THEN
     write(*,*) "opening files"
     write(*,*) trim(wdir)
     OPEN(100, FILE = trim(wdir)//'/1_contact')
     write(*,*) "opened contact"
     OPEN(101, FILE = trim(wdir)//'/1_Pzeta')
     OPEN(103, FILE = trim(wdir)//'/1_roughness')
     OPEN(109, FILE = trim(wdir)//'/distribution')
     OPEN(110, FILE = trim(wdir)//'/1_reloadDetails')
     write(*,*) "files opened"
  ENDIF
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  outstepsize = 1

  zoutfac = 5    ! Output every factor of 2 magnification 
  zoutmod = 0.1
  noutput = 100  ! output ever N magnification steps
  noutline = 1
  ! --------------------------------------------------
  ! Create input criteria for roughness power spectrum and diffusivity
  
  zetamax = 1.2d16 ! maximum magnification
  
  Lo = 1.d-3! 1.d0 ! Apparent area (1 m)
  qL = 2.d0*pi/ Lo

  Co = 1.d-6 !   ! Scaling for PSD (m^3)
  Co = 1.d-19   
  !Hu = 0.6d0   ! Hurst exponent
  !mexp = -1.d0 - 2.d0*Hu ! 1D PSD
  !mexp = -2.d0*(Hu+1.d0) ! 2D PSD

  ! Find representative magnification to calculate maximum diffusivity
  IF(Hu .LT. 0.5d0)THEN
     zetarep = zetamax
     zetarep2 = 1
  ELSE
     zetarep = 1
     zetarep2 = zetamax
  ENDIF
  ! -------------------------------------------------- 
  ! Create stress-space and structure I.C. and B.C.s
  ! B.C.s
  ! P(sig=0,zeta) = 0
  ! P(sig,zeta = 1) = delta(sig - sig0)
  ! Need to solve for B.C at sigma = sigma_{Y}(zeta)

  ! Yield stress
  sigY0 = 1.d10  ! Based on scaling from Chris's measurements
  sigY0 = 1.d8
!  sigY0 = 2.d8
  nexp = 0.2d0  ! Roughly based on Chris's measurements

  !Hu = -0.5d0*(mexp + 1.d0) ! Hurst exponent
  Hu = 0.8d0
  mexp = -2.d0*Hu - 1.d0
  
  ! Find representative magnification to calculate maximum diffusivity 
  IF(Hu .LT. 0.5d0)THEN
     zetarep = zetamax
     zetarep2 = 1
  ELSE
     zetarep = 1
     zetarep2 = zetamax
  ENDIF

  ! omega = 0 =>  n = m + 3 for 1D

  ! applied load
  sigma0 = 0.3d0 * sigY0
  
  ! Discrtization in stress-space
  dsig =0.1d-1*sigma0! sigY0!sigma0

  ! Initial diffusivity
  f1 = 0.125d0*((Ey/(1-nu*nu))**2) * (qL**3)*Co

  ! Non-dimensionalize stresses
  dimfac = dsig
  ftilde1 = f1/dimfac/dimfac
  sigY0 = sigY0 / dimfac
  sigma0 = sigma0 / dimfac
  dsig = dsig / dimfac
  
  ! calculate max yield stress
  CALL calcYieldBoundary(sigY0,sigYmax,zetamax,nexp)

  ! Total number of cells in stress-space
  CALL calcNSP(nsig,nsp,sigYmax,dsig,NPU)
  !dsig = sigYmax / dfloat(nsig-1)
  nsigmax = nsig
  CALL findYieldRank(dsig,sigY0,syrank,syindex,myrank,nsp,NPU)
  CALL MPI_Allreduce(syrank,syrankall,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  syrank = syrankall
  CALL MPI_Bcast(syindex,1,MPI_INTEGER,syrank,MPI_COMM_WORLD,ierr)
  IF(myrank.EQ.0)THEN
     WRITE(*,*) nsig, nsp, sigma0, sigY0, sigYmax
  ENDIF
  ALLOCATE(Psz(nsp),Pszpr(nsp),stat=ierr)
  ALLOCATE(d(nsp),cca(nsp),ccc(nsp),stat=ierr)
  IF(myrank.EQ.0)THEN
     ALLOCATE(Pszprall(nsigmax),stat=ierr)
  ENDIF
  
  DO si = 1,nsp
     Psz(si) = 0.d0
     Pszpr(si) = 0.d0
     d(si) = 0.d0
     cca(si) = 0.d0
     ccc(si) = 0.d0
  END DO
  nspo = nsp
  CALL MPI_Gather(Psz,nsp,MPI_REAL8,Pszprall,nsp,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

  ! --------------------------------------------------
  ! Create roughness power spectrum
  alpha = 0.125d0*0.25d0*0.25d0*0.25d0!*0.125d0*0.0625d0

  ! --------------------------------------------------
  ! Solve diffusion equation with magnification
  ALLOCATE(PszB(NPU),PszE(NPU),stat=ierr)

  ! Initial relative area of plastic contact is zero
  nz = 1
  Ael = 0.d0
  Al = 0.d0
  Fel = 0.d0
  Felall = 0.d0
  zetan = 1.d0
  Aplpr = 0.d0
  Apl = 0.d0
  Apl2 = 0.d0
  Anon = 0.d0
  Anonpr = 0.d0
  idsig = 1.d0/dsig
  i2dsig = 0.5d0*idsig*idsig ! Crank-Nicholson
  !i2dsig = idsig*idsig       ! Backward Euler
  zi = 1
  Fpl = 0.d0

  CALL calcYieldBoundary(sigY0,sigYnow,zetan,nexp)
  CALL calcNSP(nsig,nsp,sigYnow,dsig,NPU)

  ! Apply macroscopic load
  eps = dsig
  DO si = 1,nsp
     sig = dfloat(myrank*nsp + si-1)*dsig
     
     Psz(si) = diracDeltaGauss(sig,sigma0,eps)

!     IF((sig.GE.sigma0).AND.(sig.LT.(sigma0+dsig)))THEN
!        Psz(si) = 1.d0/dsig
!     ENDIF
     Al = Al + Psz(si)
     Fel = Fel + Psz(si)*sig
  END DO

  ! Need to save   Pszprall, zeta, Apl, Anon 

  Al = Al * dsig
  Fel = Fel * dsig

  CALL MPI_Reduce(Al,Ael,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Reduce(Fel,Felall,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Gather(Psz,nsp,MPI_REAL8,Pszprall,nsp,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  IF(myrank.EQ.syrank)THEN
     WRITE(*,*) "starting main loop"
     WRITE(*,*) zetan,nsig,nsp,syrank,syindex
  ENDIF

  ! Calculate diffusivity at current magnification, zeta_n
  CALL CalcFzMax(fzetan,zetan,ftilde1,mexp)
  
  ! Calculate change in yield stress at current magnification, sig_{Y}'_n
  sigYp = calcSigYp(sigY0,zetan,nexp)
  dzetan = alpha * dsig*dsig/fzetan                     ! calculate zeta step
  dzetamin = 1.d-7
!  dzetan = MAX(dzetan,dzetamin)

  Pny = Psz(syindex)
  syrankm1 = syrank-1
  IF(syrankm1.LT.0) syrankm1= 0
  CALL MPI_Bcast(Pny,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
  IF(syindex.GT.2)THEN
     Pny1 = Psz(syindex-1)
     Pny2 = Psz(syindex-2)
     CALL MPI_Bcast(Pny1,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     CALL MPI_Bcast(Pny2,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
  ELSEIF(syindex.EQ.2)THEN
     Pny1 = Psz(syindex-1)
     Pny2 = Psz(nsp)
     CALL MPI_Bcast(Pny1,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     CALL MPI_Bcast(Pny2,1,MPI_REAL8,syrankm1,MPI_COMM_WORLD,ierr)
  ELSE
     Pny1 = Psz(nsp)
     Pny2 = Psz(nsp-1)
     CALL MPI_Bcast(Pny1,1,MPI_REAL8,syrankm1,MPI_COMM_WORLD,ierr)
     CALL MPI_Bcast(Pny2,1,MPI_REAL8,syrankm1,MPI_COMM_WORLD,ierr)
  ENDIF
  
  ! Determine new magnification zeta_n+1
  zetan2 = zetan+dzetan                                 ! calculate new zeta, zeta_n+1
  CALL CalcFzMax(fzetan2,zetan2,ftilde1,mexp)           ! f_{n+1}
  sigYp2 = calcSigYp(sigY0,zetan2,nexp)                 ! sig_{Y}'_{n+1}
  dApldz = 0.d0
  dApldzold = 0.d0
  dsigY = (sigYp2/fzetan2)*(Apl + 0.5d0*dzetan*(sigYp*Pny - &   ! factor for yield B.C.
       fzetan*(idsigy*Pny -idsig*Pny1)))

  dsigYm1 =dzetan*fzetan*(Pny2*idsig*i2dsigy + Pny*i2dsigy*idsigy) + &
       (1.d0 -dzetan*fzetan*idsig*idsigy)*Pny1

  !Backward euler
  !dsigY = (sigYp2/fzetan2)*Apl
  !dsigYm1 = Pny1

  Felall = Felall*dimfac
  Fpl = Fpl *dimfac
  Ftot = Felall + Fpl
  IF(myrank.EQ.0)THEN
     WRITE(100,402) zetan, Ael, Anon, Apl, Felall, Fpl, Ftot   ! Output area of contact at every zeta step 
     CALL PRINTPSD(101,Pszprall/dimfac,zetan2,nsigmax,dsig*dimfac,outstepsize)
  ENDIF
  zetaoprev = zetan
  
  DO WHILE(zetan.LE.zetamax)
!  DO WHILE(nz.LT.100)
     ! calculate new yield boundary at zeta_{n+1} and adapt stress grid
     CALL calcYieldBoundary(sigY0,sigYnow,zetan2,nexp)
     CALL calcNSP(nsig,nsp,sigYnow,dsig,NPU)  

     ! calculate new grid and push pszprall to workers, Pszprall can be prepadded with zeros
     CALL MPI_Scatter(Pszprall,nsp,MPI_REAL8,Pszpr,nsp,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Allgather(Pszpr(nsp),1,MPI_REAL8,PszE,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
     CALL MPI_Allgather(Pszpr(1),1,MPI_REAL8,PszB,1,MPI_REAL8,MPI_COMM_WORLD,ierr)

     ! Find index and rank that holds the current yield boundary
     CALL findYieldRank(dsig,sigYnow,syrank,syindex,myrank,nsp,NPU)
     CALL MPI_Allreduce(syrank,syrankall,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr) 
     syrank = syrankall
     CALL MPI_Bcast(syindex,1,MPI_INTEGER,syrank,MPI_COMM_WORLD,ierr)

     ! We will use 2nd-order accurate Crank-Nicolson

     ! Want to update dsig for the last node spacing to reflect that the B.C. may be on the same node even if the spacing is increasing
     ddsigy = sigYnow - dfloat(syrank*nsp+syindex-1)*dsig    ! positive if sigYnow greater than grid point, negative if less
     idsigy = 1.d0/(dsig+ddsigy)
     i2dsigy = 1.d0/(2.d0*dsig+ddsigy)

     rr = dzetan*fzetan*i2dsig
     rr2 = dzetan*fzetan2*i2dsig
     !rr = 0.d0                 ! backward Euler 

     ! coefficients for 2nd-to-last node
     asigYm1 = -1.d0*dzetan*fzetan2*i2dsigy*idsig
     bsigYm1 = 1.d0 + dzetan*fzetan2*idsig*idsigy 
     csigYm1 = -1.d0*dzetan*fzetan2*idsigy*i2dsigy

     ! Backward Euler
     !asigYm1 = -2.d0*dzetan*fzetan2*i2dsigy*idsig
     !bsigYm1 = 1.d0 + 2.d0*dzetan*fzetan2*idsig*idsigy
     !csigYm1 = -2.d0*dzetan*fzetan2*idsigy*i2dsigy

     ! coefficient for yield boundary 
     asigY = -0.5d0*dzetan*sigYp2*idsigy+0.5d0*sigYp2*(dsig+ddsigy)/fzetan2
     bsigY = 1.d0 + 0.5d0*dzetan*sigYp2*(sigYp2/fzetan2 + idsigy)
     dsigY = dsigY + 0.5d0*sigYp2/fzetan2*(dsig+ddsigy)*Pny1

     ! Backward Euler
     !asigY = 0.5d0*sigYp2*(dsig+ddsigy)/fzetan2 - dzetan*sigYp2*idsigy
     !bsigY = 1.d0+ dzetan*sigYp2*(sigYp2/fzetan2 + idsigy)  

     CALL CalcLaplacian(Psz,Pszpr,PszB,PszE,nsp,rr,rr2,asigY,bsigY,dsigY, &
          asigYm1,bsigYm1,csigYm1,dsigYm1, &
          d,cca,ccc,myrank,syrank,syindex,NPU,Nred,ierr)

     ! Calculate area of elastic contact at zeta_{n+1}
     Al = 0.d0
     Fel = 0.d0
     DO si = 1, nsp
        Al = Al + Psz(si)
        Fel = Fel + Psz(si)*dfloat(myrank*nsp+si-1)*dsig
     END DO
     Al = Al * dsig
     Fel = Fel * dsig

     Pny = Psz(syindex)
     syrankm1 = syrank-1
     IF(syrankm1.LT.0) syrankm1 = 0


     CALL MPI_Bcast(Pny,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     IF(syindex.GT.2)THEN
        Pny1 = Psz(syindex-1)
        Pny2 = Psz(syindex-2)        
        CALL MPI_Bcast(Pny1,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
        CALL MPI_Bcast(Pny2,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     ELSEIF(syindex.EQ.2)THEN
        Pny1 = Psz(syindex-1)
        Pny2 = Psz(nsp)
        CALL MPI_Bcast(Pny1,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
        CALL MPI_Bcast(Pny2,1,MPI_REAL8,syrankm1,MPI_COMM_WORLD,ierr)
     ELSE
        Pny1 = Psz(nsp)
        Pny2 = Psz(nsp-1)
        CALL MPI_Bcast(Pny1,1,MPI_REAL8,syrankm1,MPI_COMM_WORLD,ierr)
        CALL MPI_Bcast(Pny2,1,MPI_REAL8,syrankm1,MPI_COMM_WORLD,ierr)
     END IF

     ! Correct yield boundary for elastic contact and force
     ! Solve for area of plastic and non-contact at zeta_{n+1}
     IF(myrank.EQ.syrank)THEN
        Al = Al + 0.5d0*(Pny + Pny1)*ddsigy - 0.5d0*Pny*dsig

        Fel = Fel + 0.5d0*Pny1*dfloat(myrank*nsp+syindex-2)*dsig*ddsigy - &
             Pny*dfloat(myrank*nsp+syindex-1)*dsig*dsig + &
             0.5d0*Pny*sigYnow*(dsig+ddsigy)

        dApldz = -1.d0*sigYp2*Pny - fzetan2*( &
             (3.d0*dsig+2.d0*ddsigy)*idsigy*i2dsigy* Pny - &
             (2.d0*dsig+ddsigy)*idsig*idsigy* Pny1 + &
             (dsig+ddsigy)*idsig*i2dsigy*Pny2)

        !Apl = Aplpr + dzetan*dApldz
        Apl = Aplpr + 0.5d0*dzetan*(dApldz+dApldzold)
        
        Apl2 = fzetan2/sigYp2
        Fpl = Apl*sigYnow*dimfac
     ENDIF
     Felall = 0.d0
     CALL MPI_Bcast(Apl,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     CALL MPI_Bcast(Fpl,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     CALL MPI_Bcast(dApldz,1,MPI_REAL8,syrank,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(Al,Ael,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     CALL MPI_Reduce(Fel,Felall,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     ! call mpi gather to collect Psz to pszprall to head node
     CALL MPI_Gather(Psz,nsp,MPI_REAL8,Pszprall,nsp,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

     Aplpr = Apl
     dApldzold = dApldz
     Felall = Felall*dimfac
     syindexhome = syrank*nsp+syindex 

     IF(myrank.EQ.0)THEN
        Anon = Anonpr + 0.5d0*idsig*dzetan*(fzetan* &
             (-0.5d0*Pszpr(3) + 2.d0*Pszpr(2))&
             +fzetan2* &
             (-0.5d0*Psz(3) + 2.d0*Psz(2)))
        !        Anon = Anonpr + idsig*dzetan*fzetan2* &
        !             (-0.5d0*Psz(3) + 2.d0*Psz(2))
        !arraySize = size(Pszprall)
       IF((zetan2.GE.(zoutfac*zetaoprev)).OR.(mod(nz,noutput).EQ.0))THEN
           Ftot = Felall + Fpl

           WRITE(100,402) zetan2, Ael, Anon, Apl, Felall, Fpl,Ftot                              ! Output relative area of contact
           WRITE(103,401) zetan2, fzetan2*dimfac, sigYp2*dimfac, sigYnow*dimfac                 ! Output diffusivity and yield stress
!           IF((zetan2-zprev).GE.zoutmod)THEN
!              CALL PRINTPSD(101,Pszprall/dimfac,zetan2,nsigmax,dsig*dimfac,outstepsize)
!              zprev = zetan2
!           ENDIF
           zetaoprev = zetan
        ENDIF
        Anonpr = Anon
     ENDIF

 
     ! Solve for new diffusivities and time step
     zetan = zetan2
     fzetan = fzetan2
     sigYp = sigYp2
     dzetan = alpha * dsig*dsig/fzetan  ! calculate zeta step
     zetan2 = zetan+dzetan
     CALL CalcFzMax(fzetan2,zetan2,ftilde1,mexp)
     sigYp2 = calcSigYp(sigY0,zetan2,nexp)
     !dsigY = (sigYp2/fzetan2)*(Apl+0.5d0*dzetan*dApldz)

     ! Information for calculating yield B.C at next zeta step
     dsigY = (sigYp2/fzetan2)*(Apl + 0.5d0*dzetan*(sigYp*Pny - & 
          fzetan*(idsigy*Pny -idsig*Pny1)))

     dsigYm1 =dzetan*fzetan*(Pny2*idsig*i2dsigy + Pny*i2dsigy*idsigy) + &
          (1.d0 -dzetan*fzetan*idsig*idsigy)*Pny1

     ! backward euler
     !dsigY = (sigYp2/fzetan2)*Apl
     !dsigYm1 = Pny1

     currTime = MPI_Wtime()
     dTime = currTime - startTime
     IF(dTime.GE.tfin)THEN
        IF(myrank.EQ.0)THEN
           WRITE(*,*) currTime,startTime,dTime
           arraySize = size(Pszprall) 
           DO si = 1,arraySize
              WRITE(109,404) Pszprall(si)
           ENDDO
           FLUSH(109)
           WRITE(110,403) zetan, zetan2, dzetan, dsigY, dsigYm1, &
                sigYp, sigYp2, fzetan, fzetan2, &
                Apl, Anon, dApldz, arraySize
           FLUSH(110)
        ENDIF
        startTime = currTime
     ENDIF

     nz = nz + 1
  END DO

  IF(myrank.EQ.0)THEN
     WRITE(*,*) "finished main loop ",nz
     CLOSE(101)
     CLOSE(100)
     CLOSE(103)
     CLOSE(109)
     CLOSE(110)
  ENDIF
400 format (3(1X,E15.7E3))  ! 1_contact
401 format (4(1X,E15.7E3))  ! analytic
402 format (1(1X,E15.7E3),6(1X,E15.7E3))
403 format (12(1X,E35.28),1X,I20) ! output details
404 format (1X,E35.28) 
  ! Deallocate arrays
  IF(ALLOCATED(Pszpr)) DEALLOCATE(Pszpr)
  IF(ALLOCATED(Psz)) DEALLOCATE(Psz)
  IF(ALLOCATED(Pszprall)) DEALLOCATE(Pszprall)
  IF(ALLOCATED(zetas)) DEALLOCATE(zetas)
  IF(ALLOCATED(sfz)) DEALLOCATE(sfz)
  IF(ALLOCATED(fz)) DEALLOCATE(fz)
  IF(ALLOCATED(Cq)) DEALLOCATE(Cq)
  IF(ALLOCATED(dpdzo)) DEALLOCATE(dpdzo)
  IF(ALLOCATED(dpdzo2)) DEALLOCATE(dpdzo2)
  IF(ALLOCATED(cca)) DEALLOCATE(cca,ccc,d)
  IF(ALLOCATED(PszB)) DEALLOCATE(PszB)
  IF(ALLOCATED(PszE)) DEALLOCATE(PszE)
  IF(ALLOCATED(Pz)) DEALLOCATE(Pz,Pnon,Ppl)
  IF(ALLOCATED(Pzsum)) DEALLOCATE(Pzsum,Pnonsum,Pplsum)
  
  CALL MPI_FINALIZE(ierr)
CONTAINS  
  !!--------------------------------------------------
  !!                  Sub-routines
  !! -------------------------------------------------

  !---------------------------------------------------
  ! Calculate max yield stress
  SUBROUTINE calcYieldBoundary(sigY0,sigYmax,zetamax,nexp)
    IMPLICIT NONE
    REAL*8 :: sigY0,sigYmax,zetamax,nexp
    sigYmax = sigY0 * zetamax**nexp
    RETURN
  END SUBROUTINE calcYieldBoundary


  !----------------------------------------------------
  ! calculate number of stress cells per worker
  SUBROUTINE calcNSP(nsig,nsp,sigYnow,dsig,NPU)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: dsig, sigYnow
    INTEGER, INTENT(IN) :: NPU
    INTEGER, INTENT(INOUT) :: nsig, nsp

    nsp = ceiling(sigYnow/dsig/dfloat(NPU))
    nsig = nsp * NPU
   
    RETURN
  END SUBROUTINE CALCNSP

  !--------------------------------------------------  
  ! Calculate initial distribution
  FUNCTION diracDeltaGauss(sig,sigma0,eps)
    IMPLICIT NONE
    REAL*8 :: sig, sigma0, eps, pi, diracDeltaGauss
    PARAMETER (pi=3.141592653589793d0)
    diracDeltaGauss = 1/(eps*dsqrt(pi))*dexp(-1.d0*(sig-sigma0)**2 / eps**2)
    RETURN
  END FUNCTION diracDeltaGauss

  !--------------------------------------------------
  ! Calculate sigma_{Y}'(\zeta)
  FUNCTION calcSigYp(sigY0,zeta,nexp)
    IMPLICIT NONE
    REAL*8 :: sigY0,zeta,nexp,calcSigYp
    calcSigYp = sigY0 * nexp * zeta**(nexp-1.d0)
    RETURN
  END FUNCTION calcSigYp

  !--------------------------------------------------
  ! Find rank and index hosting current yield stress
  SUBROUTINE findYieldRank(dsig,sigYnow,syrank,syindex,myrank,nsp,NPU)
    IMPLICIT NONE
    INTEGER:: i,NPU
    REAL*8 :: sig
    REAL*8, INTENT(IN) :: sigYnow, dsig
    INTEGER, INTENT(IN) :: myrank,nsp
    INTEGER, INTENT(INOUT) :: syrank,syindex
    syrank = -555
    syindex = -555

    DO i = 1,nsp
       sig = dfloat(myrank*nsp+i-1)*dsig
!       IF((sig.LE.sigYnow).AND.((sig+0.5d0*dsig).GE.sigYnow))THEN
!          syrank = myrank
!          syindex = i
!       ELSEIF((sigYnow.GT.(sig-0.5d0*dsig)).AND.(sigYnow.LT.sig))THEN
!          syrank = myrank
!          syindex = i
!       ENDIF
       IF((sigYnow.GE.(sig-0.5d0*dsig)).AND.(sigYnow.LT.sig+0.5d0*dsig))THEN
          syrank = myrank
          syindex = i
       ENDIF

!       ELSEIF(((sig+dsig).GT.sigYnow).AND.((sig+0.5d0*dsig).LT.sigYnow))THEN
!          IF(i.LT.nsp)THEN
!             syrank = myrank
!             syindex = i+1
!          ELSE
!             syrank = myrank+1
!             syindex = 1
!          ENDIF
!          
!       ENDIF
    END DO
    IF((myrank.EQ.(NPU-1)).AND.(sigYnow.GE.(dfloat(myrank*nsp+nsp-1)*dsig)))THEN
       syrank=myrank
       syindex=nsp
    ENDIF

    RETURN
  END SUBROUTINE findYieldRank

  !--------------------------------------------------
  ! Calculate the laplacian of P using modified Thomas Algorithm
  ! for parallelized solver for tridiagonal operation
  SUBROUTINE CalcLaplacian(Psz,Pszpr,PszB,PszE,nsp,rr,rr2,asigY,bsigY,dsigY, &
       asigYm1,bsigYm1,csigYm1,dsigYm1, &
       d,cca,ccc,myrank,syrank,syindex,NPU,Nred,ierr)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER :: i, nsp, myrank, ierr, NPU, Nred2, Nred, syrank, syindex
    REAL*8, DIMENSION(:) :: Psz,Pszpr
    REAL*8, DIMENSION(:) :: d,cca,ccc
    REAL*8 :: rr,bcen,r,rr2,bcen2,nrr2,asigY,bsigY,dsigY
    REAL*8 :: asigYm1,bsigYm1,csigYm1,dsigYm1
    REAL*8, DIMENSION(Nred) :: dred,ared,cred
    REAL*8, DIMENSION(2) :: dredsub,aredsub,credsub
    REAL*8, DIMENSION(NPU) :: PszE,PszB

    Nred2 = 2*(syrank+1)
    bcen2 = 1.d0 + 2.d0*rr2   ! b_ij for new time step
    nrr2  = -1.d0*rr2         ! a_ij, c_ij for new time step

    bcen  = 1.d0 - 2.d0*rr    ! b_ij for prev time step

    
    ! Problem is set up to use Crank-Nicolson implicit solver, 2nd-order acc in space and time
    ! Note that case where rr = 0 is backward Euler solution, 2nd-order acc in space, 1st-order in time

    IF(myrank.LT.syrank)THEN       
       IF(myrank.EQ.0)THEN
          ! Enforce P(0,zeta) = 0 B.C.
          d(1)   = 0.d0
          cca(1) = 0.d0
          ccc(1) = 0.d0
       ELSE
          d(1)   = rr*(PszE(myrank)+Pszpr(2)) + bcen*Pszpr(1)
          d(1)   = d(1) / bcen2
          ccc(1) = nrr2/bcen2
          cca(1) = nrr2/bcen2
       ENDIF
       
       d(2)   = rr*(Pszpr(1) + Pszpr(3)) + bcen*Pszpr(2)
       d(2)   = d(2) / bcen2
       ccc(2) = nrr2/bcen2
       cca(2) = nrr2/bcen2
       
       ! Modified Thomas algorithm
       ! We first want to transform every submatrix into a form
       ! where the diagonal will take b_{ij}' = 1 and the 
       ! equations for each row are a function of the first and last column
       ! while the lower and upper diagonal elements set to zero
       ! i.e. the equations for each line will relate the solution x_{i} to
       ! the first and last solution x_{1} and x_{nsp}
       
       ! remove lower diagonal element starting at third row
       DO i = 3,(nsp-1)
          d(i)   = rr*(Pszpr(i-1)+Pszpr(i+1)) + bcen*Pszpr(i)
          r      = 1.d0/(bcen2 + rr2*ccc(i-1))   
          d(i)   = r*(d(i) + rr2*d(i-1))
          ccc(i) = r*nrr2
          cca(i) = r*rr2*cca(i-1)
       END DO

       IF((syindex.EQ.1).AND.(myrank.EQ.(syrank-1)))THEN
          d(nsp) = dsigYm1
          r = 1.d0/(bsigYm1 - asigYm1*ccc(nsp-1))
          d(nsp) = r*(d(nsp) - asigYm1*d(nsp-1))
          ccc(nsp) = r*csigYm1
          cca(nsp) = -1.d0*r*asigYm1*cca(nsp-1)
       ELSE
          d(nsp) = rr*(Pszpr(nsp-1)+PszB(myrank+2)) + bcen*Pszpr(nsp)
          r = 1.d0/(bcen2 + rr2*ccc(nsp-1))
          d(nsp)   = r*(d(nsp) + rr2*d(nsp-1))
          ccc(nsp) = r*nrr2
          cca(nsp) = r*rr2*cca(nsp-1)
       ENDIF

       ! Back substitute to remove  upper diagonal element
       DO i = (nsp-2),2,-1
          d(i)   = d(i) - ccc(i)*d(i+1)
          ccc(i) = -1.d0*ccc(i)*ccc(i+1)         
          cca(i) = cca(i) - ccc(i)*cca(i+1)      
       END DO
       
       ! Note that b(i) = 1 now
       r = 1.d0/(1.d0-cca(2)*ccc(1))
       d(1) = r*(d(1)-ccc(1)*d(2))
       ccc(1) = -1.d0*r*ccc(1)*ccc(2)
       cca(1) = r*cca(1)
       
       aredsub(1) = cca(1)
       credsub(1) = ccc(1)
       dredsub(1) = d(1)
       aredsub(2) = cca(nsp)
       credsub(2) = ccc(nsp)
       dredsub(2) = d(nsp)
 
       ! Working on yield rank
    ELSEIF(myrank.EQ.syrank)THEN
       IF(syindex.EQ.1)THEN
          aredsub(1) = asigY/bsigY
          credsub(1) = 0.d0
          dredsub(1) = dsigY/bsigY
          aredsub(2) = 0.d0
          credsub(2) = 0.d0
          dredsub(2) = 0.d0
       ELSEIF(syindex.EQ.2)THEN
          aredsub(1) = asigYm1/bsigYm1
          credsub(1) = csigYm1/bsigYm1
          dredsub(1) = dsigYm1/bsigYm1
          
          dredsub(2) = dsigY/bsigY
          aredsub(2) = asigY/bsigY
          credsub(2) = 0.d0          
       ELSE
          IF(myrank.EQ.0)THEN
             d(1) = 0.d0
             cca(1) = 0.d0
             ccc(1) = 0.d0
          ELSE
             d(1) = rr*(PszE(myrank)+Pszpr(2)) + bcen*Pszpr(1)
             d(1)   = d(1) / bcen2
             ccc(1) = nrr2/bcen2
             cca(1) = nrr2/bcen2
          ENDIF
          d(2)   = rr*(Pszpr(1) + Pszpr(3)) + bcen*Pszpr(2)
          d(2)   = d(2) / bcen2
          ccc(2) = nrr2/bcen2
          cca(2) = nrr2/bcen2
          
          DO i = 3,(syindex-2)
             d(i)   = rr*(Pszpr(i-1)+Pszpr(i+1)) + bcen*Pszpr(i)
             r      = 1.d0/(bcen2 + rr2*ccc(i-1))
             d(i)   = r*(d(i) + rr2*d(i-1))
             ccc(i) = r*nrr2
             cca(i) = r*rr2*cca(i-1)
          END DO

          d(syindex-1) = dsigYm1
          r = 1.d0/(bsigYm1-asigYm1*ccc(syindex-2))
          d(syindex-1) = r*(d(syindex-1) - asigYm1*d(syindex-2))
          ccc(syindex-1) = r*csigYm1
          cca(syindex-1) = -1.d0*r*asigYm1*cca(syindex-2)
          
          ! Implement Sigma_{Y} B.C.
          d(syindex) = dsigY
          r = 1.d0/(bsigY - asigY*ccc(syindex-1))
          d(syindex) = r*(d(syindex) - asigY*d(syindex-1))
          ccc(syindex) = 0.d0
          cca(syindex) = -1.d0*r*cca(syindex-1)*asigY
          
          DO i = (syindex-2),2,-1
             d(i)   = d(i) - ccc(i)*d(i+1)
             ccc(i) = -1.d0*ccc(i)*ccc(i+1)
             cca(i) = cca(i) - ccc(i)*cca(i+1)
          END DO
          r = 1.d0/(1.d0-cca(2)*ccc(1))
          d(1) = r*(d(1)-ccc(1)*d(2))
          ccc(1) = -1.d0*r*ccc(1)*ccc(2)
          cca(1) = r*cca(1)
          
          aredsub(1) = cca(1)
          credsub(1) = ccc(1)
          dredsub(1) = d(1)
          aredsub(2) = cca(syindex)
          credsub(2) = ccc(syindex)
          dredsub(2) = d(syindex)
       ENDIF
    ELSE
       aredsub(1) = 0.d0
       credsub(1) = 0.d0
       dredsub(1) = 0.d0
       aredsub(2) = 0.d0
       credsub(2) = 0.d0
       dredsub(2) = 0.d0
    ENDIF

    ! Gather first and last entry and solver reduced tridiagonal problem with standard Thomas 
    CALL MPI_Gather(dredsub,2,MPI_REAL8,dred,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Gather(aredsub,2,MPI_REAL8,ared,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    CALL MPI_Gather(credsub,2,MPI_REAL8,cred,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    
    IF(myrank.EQ.0)THEN
       ! d1 = 0 and c1 = 0 for zero B.C. and b = 1 for all 
       DO i = 2, Nred2
          r = 1.d0/(1.d0 - ared(i)*cred(i-1))
          dred(i) = r*(dred(i) - ared(i)*dred(i-1))
          cred(i) = r*cred(i)
       END DO
       DO i = (Nred2-1),1,-1
          dred(i) = dred(i)-cred(i)*dred(i+1)
       END DO
    ENDIF

    ! Scatter to cores and solve for remaining unknowns
    CALL MPI_Scatter(dred,2,MPI_REAL8,dredsub,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    
    IF(myrank.LT.syrank)THEN
       Psz(1) = dredsub(1)
       Psz(nsp) = dredsub(2)
       DO i = 2, (nsp-1)
          Psz(i) = d(i) - cca(i)*Psz(1) - ccc(i)*Psz(nsp)
       END DO
    ELSEIF(myrank.EQ.syrank)THEN       
       Psz(syindex) = dredsub(2)
       Psz(1) = dredsub(1)
       DO i = 2,(syindex-1)
          Psz(i) = d(i) - cca(i)*Psz(1) - ccc(i)*Psz(syindex)
       END DO
       DO i = (syindex+1),nsp
          Psz(i) = 0.d0
       END DO
    ELSE
       DO i = 1,nsp
          Psz(i) = 0.d0
       END DO
    ENDIF
    
  END SUBROUTINE CalcLaplacian

  !--------------------------------------------------
  ! Calculate maximum diffusivity for zeta step
  !SUBROUTINE CalcFzMax(fzmax,zetarep,qL,Ey,nu,Co,mexp)
  SUBROUTINE CalcFzMax(fzmax,zetarep,ftilde1,mexp)
    IMPLICIT NONE
    REAL*8 :: ftilde1,mexp,zetarep,fzmax,zexp

    zexp = mexp + 2.d0
    fzmax = ftilde1 * zetarep**zexp
    
    RETURN
  END SUBROUTINE CalcFzMax
  
  !-------------------------------------------------
  SUBROUTINE PRINTPSD(nfile,Pszall,zetai,nsig,dsig,outstepsize)
    IMPLICIT NONE
    INTEGER nfile, nsig,outstepsize,si
    REAL*8 :: zetai,sig,dsig
    REAL*8, INTENT(IN) :: Pszall(nsig)
    
    DO si = 1,nsig,outstepsize
       sig = dfloat(si-1)*dsig
       WRITE(nfile,401) sig, Pszall(si)
    END DO
    WRITE(nfile,365) '# **',' above is for zeta ', zetai
    WRITE(nfile,'()')
    ENDFILE nfile
    BACKSPACE nfile

401 format (2(1X,E15.7E3))     ! 1_Pzeta
365 format (T1,A,A,E15.7)
    RETURN
  END SUBROUTINE PRINTPSD
  
  
  !--------------------------------------------------
  SUBROUTINE getdata(unit,line)
    IMPLICIT NONE
    INTEGER unit
    CHARACTER(200) :: line
    CHARACTER(1) :: char
    INTEGER i
    LOGICAL linecont

    char='#'
    linecont = .TRUE.
    DO WHILE(linecont)
       IF(char.EQ.'#')THEN
          read(unit,'(a)')line
          i=1
          char=line(1:1)
       ELSEIF(char.EQ.' ')THEN
          i=i+1
          char=line(i:i)
       ELSE
          linecont = .FALSE.
       ENDIF
    END DO
    
    RETURN
  END SUBROUTINE getdata
  
  !---------------------------------------------------
  SUBROUTINE readinput(wdir)
    INCLUDE 'mpif.h'
    INTEGER :: myrank, position, ierr
    CHARACTER(200) :: dataline, wdir
    INTEGER, PARAMETER :: psize = 512
    CHARACTER, DIMENSION(psize) ::  packed

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nsize,ierr)

    position = 0
    IF(0.EQ.myrank)THEN
       CALL getdata(5,dataline)
       READ(dataline,'(a)') wdir    
    ENDIF
    CALL MPI_Bcast(wdir,200,MPI_CHARACTER, &
         0,MPI_COMM_WORLD,ierr)
    RETURN
  END SUBROUTINE readinput
  
END PROGRAM RoughDiff

