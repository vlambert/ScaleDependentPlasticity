PROGRAM RoughDiff
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  CHARACTER(200) :: wdir
  REAL*8 :: Sigma0,Lo,alpha ! Macroscopic load and scale
  REAL*8 :: Hu,qmax,qL,dz,mexp,nexp,Co,alphai,betai, &
       Anon,Apl,Aplpr,Anonpr,dimfac,zetan,zetan2,fzetan,fzetan2,dzetan,Al,Ael,eta
  REAL*8 :: dsig,dzeta,dmu,i2dsig,idsig,sigYmax,sigY0,Pszm,zetamax,fzmax,fzmin,dzetamax
  INTEGER :: nz,nsig,nsp,ierr, si, zi, myrank, NPU, nznpu,size,soindex
  
  INTEGER :: ntot,sni,sne,ni,ntotsub,nintz,next,Nred
  REAL*8  :: myni,ni2,dmuint,sfzint,rr2
  REAL*8  :: s0sy, sns0, npsy, mysign,mypref,mysns0,mynpsy,nexpo,xx,zetarep,zetarep2
  REAL*8, DIMENSION(:), ALLOCATABLE :: Pz,Pnon,Ppl
  
  INTEGER :: outstepsize
  REAL*8 :: zoutfac,zetaoprev,rr
  REAL*8, DIMENSION(:), ALLOCATABLE :: Cq,zetas,mus,sigY,dsigY,fz,sfz,sig
  REAL*8, DIMENSION(:), ALLOCATABLE :: dpdzo,dpdzo2,Psz,Pszpr,PszB,PszE,Pszall,Pzsum, &
       Pplsum,Pnonsum
  INTEGER :: syindex,syrank,syrankall
  REAL*8 :: Ey, nu, pi, Eyp
  
  PARAMETER(Ey = 1.d11, &  ! Young's modulus
       nu = 0.25)          ! Poisson's ratio
       
  PARAMETER (pi=3.141592653589793d0)
  
  ! ---------------------------------------------------
  ! Initialize MPI
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
  NPU = size ! Number of cores
  Nred = 2*NPU
  ! Discretization in stress space should be a multiple of NPU

  CALL readinput(wdir)
  CALL FLUSH(6)
  IF(myrank.EQ.0)THEN
     write(*,*) "opening files"
     write(*,*) trim(wdir)
     OPEN(102, FILE = trim(wdir)//'/1_areaAnalytic')
     OPEN(103, FILE = trim(wdir)//'/1_roughness')
     write(*,*) "files opened"
  ENDIF
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  outstepsize = 10

  zoutfac = 2    ! Output every factor of 2 magnification 
  
  ! --------------------------------------------------
  ! Create input criteria for roughness power spectrum and diffusivity
  
  zetamax = 1.d4 ! maximum magnification
  
  Lo = 1.d0 ! Apparent area (1 m)
  qL = 2.d0*pi/ Lo

  Co = 1.d-5   !1.d-8 1.d-2  ! 1.d-12   ! Scaling for PSD (m^3)
  Hu = 0.6d0   ! Hurst exponent
  mexp = -1.d0 - 2.d0*Hu ! 1D PSD
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
  sigY0 = 1.d9
  sigY0 = 1.d10
  sigYmax = sigY0
  
  ! applied load
  sigma0 = 0.99d0 * sigY0
  
  ! Discrtization in stress-space
  !dsig = 5d-4 * sigma0 
  dsig = 1.d-3*sigma0
  ! sigma0 = 0.1, dsig = 0.1*sigma0

  ! Non-dimensionalize stresses
  dimfac = dsig
  sigma0 = sigma0 / dimfac
  sigY0 = sigY0 / dimfac
  dsig = dsig / dimfac

  ! Total number of cells in stress-space
  nsig = nint( sigY0 / dsig ) + 1
  nsp = nsig/NPU         ! find integer multiple of NPUs
  nsp = nsp + 1          ! number of stress discretization points per CPU
  nsig = nsp * NPU       ! total number of stress discretization points
  dsig = sigY0 / dfloat(nsig-1)

  IF(myrank.EQ.0)THEN
     WRITE(*,*) nsig, nsp, Sigma0, sigY0, sigYmax
  ENDIF
  ALLOCATE(sig(nsp),stat=ierr)
  DO si = 1,nsp
     sig(si) = (myrank*nsp + si-1)*dsig   
  END DO

  CALL MPI_Allreduce(syrank,syrankAll,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  syrank = syrankAll
  IF(myrank.EQ.syrank)THEN
     WRITE(*,*) syrank,syindex
  ENDIF

  Eyp = Ey / dimfac


  nz = 201
  nintz = nz*20
  dmu = log(zetamax)/dfloat(nz-1)
  dmuint = log(zetamax)/dfloat(nintz-1)
  ALLOCATE(sfz(nz),fz(nz),zetas(nz),Cq(nz),stat=ierr)
  sfz(1) = 0.d0
  next = 1
  zetan = 1.d0
  zetas(1) = zetan
  CALL calcFz(fzetan2,zetan,qL,Eyp,nu,Co,mexp)
  CALL calcCq(Cq(1),Co,mexp,zetan)
  fz(1) = fzetan2
  sfzint = 0.d0
  zetan2 = 1.d0
  DO zi = 2,nintz
     zetan = dexp(dmuint*dfloat(zi-1))
     CALL calcFz(fzetan,zetan,qL,Eyp,nu,Co,mexp)
     sfzint = sfzint + 0.5d0*dmuint*(zetan*fzetan+zetan2*fzetan2)
     IF(zetan.GE.dexp(dmu*dfloat(next)))THEN
        CALL calcCq(Cq(next+1),Co,mexp,zetan)
        fz(next+1) = fzetan
        sfz(next+1) = sfzint/sigY0/sigY0
        zetas(next+1) = zetan
        next = next + 1
     ENDIF
     fzetan2 = fzetan
     zetan2 = zetan
  END DO

  !! --------------------------------------------------------------
  !!                    Analytical solution
  !! --------------------------------------------------------------
  ALLOCATE(Pz(nz),Pnon(nz),Ppl(nz),STAT=ierr)
  IF(myrank.EQ.0)THEN
     ALLOCATE(Pzsum(nz),Pnonsum(nz),Pplsum(nz),STAT=ierr)
     DO zi = 1,nz
        Pzsum(zi) = 0.d0
        Pnonsum(zi) = 0.d0
        Pplsum(zi) = 0.d0
     END DO
     WRITE(*,*) "allocated"
  ENDIF

  DO zi = 1, nz
     Pz(zi) = 0.d0
     Pnon(zi) = 0.d0
     Ppl(zi) = 0.d0
  END DO
  

  xx = pi * sigma0/sigY0
  ntot = 100000

  eta = 1.d-3
  ntot = nint(-1.d0*log(eta)/(pi*pi*sfz(2)))
  ntotsub = ntot / NPU
  IF(ntotsub.LT.10) ntotsub = 10
  If(myrank.EQ.0)THEN
     WRITE(*,*) ntot,ntotsub,sfz(2)
  ENDIF
  sni = (myrank*ntotsub + 1)
  sne = (myrank+1)*ntotsub

  DO ni = sni, sne
     ni2  = dfloat(ni)
     sns0 = DSIN( ni2 * xx) / ni2
     npsy = (ni2 * pi) **2.d0
     mysign = (-1.d0)**(dfloat(ni)+1.d0)

     ! For Pz
     myni = dfloat(2*ni - 1)
     mysns0 = DSIN(myni * xx) / myni
     mynpsy = (myni * pi) **2
     
     DO zi = 1, nz
        nexpo =  dexp(-1.d0 * npsy * sfz(zi))
        
        Pz(zi) = Pz(zi) + mysns0 * dexp(-1.d0 * mynpsy * sfz(zi))
        Pnon(zi) = Pnon(zi) + sns0 * nexpo
        Ppl(zi) = Ppl(zi) + mysign * sns0 * nexpo
     END DO     
  END DO

  IF(myrank.EQ.0)THEN
     WRITE(*,*) "reducing"
  ENDIF
  CALL MPI_REDUCE(Pz,Pzsum,nz,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(Pnon,Pnonsum,nz,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(Ppl,Pplsum,nz,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  IF(myrank.EQ.0)THEN
     WRITE(*,*) "done counting"
     DO zi = 1, nz
        Pzsum(zi) = (4.d0/pi) * Pzsum(zi)
        Pplsum(zi) = sigma0/sigY0 - (2.d0/pi)*Pplsum(zi)
        Pnonsum(zi) = (1.d0 - sigma0/sigY0) - (2.d0/pi)*Pnonsum(zi)
        WRITE(102,401) zetas(zi),Pzsum(zi),Pnonsum(zi),Pplsum(zi)
        WRITE(103,400) zetas(zi),fz(zi),Cq(zi)
     END DO
  ! at each time step want to calculate Pzeta, Pnon, Ppl
  ! will also want to calculate P(sigma,zeta)
     WRITE(*,*) "done writing"
     CLOSE(102)
     CLOSE(103)
  ENDIF
400 format (3(1X,E15.7E3))  ! 1_contact
401 format (4(1X,E15.7E3))  ! analytic
  ! Deallocate arrays
  IF(ALLOCATED(Pszpr)) DEALLOCATE(Pszpr)
  IF(ALLOCATED(Psz)) DEALLOCATE(Psz)
  IF(ALLOCATED(Pszall)) DEALLOCATE(Pszall)
  IF(ALLOCATED(zetas)) DEALLOCATE(zetas)
  IF(ALLOCATED(sig)) DEALLOCATE(sig)
  IF(ALLOCATED(sfz)) DEALLOCATE(sfz)
  IF(ALLOCATED(fz)) DEALLOCATE(fz)
  IF(ALLOCATED(Cq)) DEALLOCATE(Cq)
  IF(ALLOCATED(dpdzo)) DEALLOCATE(dpdzo)
  IF(ALLOCATED(dpdzo2)) DEALLOCATE(dpdzo2)
  IF(ALLOCATED(PszB)) DEALLOCATE(PszB,PszE)
  IF(ALLOCATED(Pz)) DEALLOCATE(Pz,Pnon,Ppl)
  IF(ALLOCATED(Pzsum)) DEALLOCATE(Pzsum,Pnonsum,Pplsum)
  
  CALL MPI_FINALIZE(ierr)
CONTAINS  
  !!--------------------------------------------------
  !!                  Sub-routines
  !! -------------------------------------------------

  !---------------------------------------------------
  ! Calculate power spectral density for self-affine rough surface
  SUBROUTINE calcCq(Cq,Co,mexp,zeta)
    IMPLICIT NONE
    REAL*8 Cq,Co,mexp,zeta
    Cq = Co * zeta**mexp
    RETURN
  END SUBROUTINE calcCq

  !--------------------------------------------------
  ! Calculate diffusivity for zeta step
  SUBROUTINE calcFz(fzmax,zetarep,qL,Ey,nu,Co,mexp)
    IMPLICIT NONE
    REAL*8 :: Ey,nu,Co,mexp,qL,zetarep,fzmax,Cq,coeff,q

    q = zetarep*qL
    coeff = 0.125d0*((Ey/(1-nu*nu))**2) * qL*q*q
    Cq = Co* zetarep**mexp
    fzmax = coeff * Cq 
    
    RETURN
  END SUBROUTINE calcFz
  
  !-------------------------------------------------
  SUBROUTINE PRINTPSD(nfile,Pszall,zetai,nsig,outstepsize)
    IMPLICIT NONE
    INTEGER nfile, nsig,outstepsize,si
    REAL*8 :: zetai
    REAL*8 :: Pszall(nsig),zetas(nsig)
    
    DO si = 1,nsig,outstepsize
       WRITE(nfile,401) Pszall(si)
    END DO
    WRITE(nfile,365) '# **',' above is for zeta ', zetai
    WRITE(nfile,'()')
    ENDFILE nfile
    BACKSPACE nfile

401 format (1X,E15.7)     ! 1_Pzeta
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
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

    position = 0
    IF(0.EQ.myrank)THEN
       CALL getdata(5,dataline)
       READ(dataline,'(a)') wdir    
    ENDIF
    CALL MPI_BCAST(wdir,200,MPI_CHARACTER, &
         0,MPI_COMM_WORLD,ierr)
    RETURN
  END SUBROUTINE readinput
  
END PROGRAM RoughDiff

