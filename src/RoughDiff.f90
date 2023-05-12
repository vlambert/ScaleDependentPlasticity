PROGRAM RoughDiff
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  CHARACTER(200) :: wdir
  REAL*8 :: Sigma0,Lo ! Macroscopic load and scale
  REAL*8 :: Hu,qmax,qL,dz,mexp,nexp,Co,alphai,betai, &
       Anon,Apl,Aplpr,Anonpr
  REAL*8 :: dsig,dzeta,dmu,i2dsig,idsig,sigYmax,sigY0,Pszm,zetamax
  INTEGER :: nz,nsig,nsp,ierr, si, zi, myrank, NPU, nznpu,size,soindex

  INTEGER :: outstepsize
  REAL*8 :: zoutfac,zetaoprev 
  REAL*8, DIMENSION(:), ALLOCATABLE :: Cq,zetas,mus,sigY,dsigY,fz,sig
  REAL*8, DIMENSION(:), ALLOCATABLE :: dpdzo,Psz,Pszpr,PszB,PszE,Pszall
  INTEGER, DIMENSION(:), ALLOCATABLE :: syindex,syrank
  REAL*8 :: Ey, nu, pi
  
  PARAMETER(Ey = 1.d11, & ! Young's modulus
       nu = 0.3)          ! Poisson's ratio
       
  PARAMETER (pi=3.141592653589793d0)
  
  ! ---------------------------------------------------
  ! Initialize MPI
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
  NPU = size ! Number of cores
  ! Discretization in stress space should be a multiple of NPU

  CALL readinput(wdir)

  IF(myrank.EQ.0)THEN
     OPEN(100, FILE = trim(wdir)//'/1_contact')
     OPEN(101, FILE = trim(wdir)//'/1_Pzeta')
  ENDIF
  outstepsize = 10
  zoutfac = 2    ! Output every factor of 2 magnification 

  ! --------------------------------------------------
  ! Create or input a roughness distribution
  
  zetamax = 1.d4
  nz = 1000               ! minimum desired magnification increments
  nznpu = nz/NPU  ! find integer multiple of NPUs
  nznpu = nznpu + 1
  nz = nznpu*NPU         ! calculate # of magnification increments as multiple of NPU

  
  ALLOCATE(zetas(nz),Cq(nz),mus(nz),stat=ierr)
  ALLOCATE(syindex(nz),syrank(nz),stat=ierr)
  dzeta = (zetamax - 1.d0)/dfloat(nz-1)
  dmu = dlog(zetamax)/dfloat(nz-1)
  DO zi = 1,nz
     mus(zi) = (zi-1)*dmu
     zetas(zi) = dexp(mus(zi)) !dexp( (zi-1) )!1.d0 + (zi-1)*dz
     syindex(zi) = -555
     syrank(zi) = -555
  END DO


  IF(myrank.EQ.0)THEN
     WRITE(*,*) nz,NPU,mus(1),mus(nz),zetas(1),zetas(nz)
  ENDIF
  Lo = 1.d-2 ! Apparent area (1 cm)
  qL = 2.d0*pi/ Lo

  Co = 1.d-7 ! Scaling for PSD
  Hu = 0.6   ! Hurst exponent
  mexp = -2.d0*(Hu+1.d0)

  CALL createPSD(Cq,Co,mexp,zetas,nz)

  ! -------------------------------------------------- 
  ! Create stress-space and structure I.C. and B.C.s
  ! B.C.s
  ! P(sig=0,zeta) = 0
  ! P(sig,zeta = 1) = delta(sig - sig0)
  ! Need to solve for B.C at sigma = sigma_{Y}(zeta)
  Sigma0 = 1.d7
  dsig = Sigma0 / 100.d0

  sigY0 = 1.d9
  nexp = mexp + 4.d0 ! w = 0
  ALLOCATE(sigY(nz),dsigY(nz),stat=ierr)
  CALL createYieldBoundary(sigY,dsigY,zetas,sigY0,nexp,nz)
  sigYmax = sigY(nz)
  nsig = nint(sigYmax/dsig) + 1
  nsp = nsig/NPU         ! find integer multiple of NPUs
  nsp = nsp + 1          ! number of stress discretization points per core
  nsig = nsp * NPU       ! total number of stress discretization points
  dsig = sigYmax/dfloat(nsig-1)

  IF(myrank.EQ.0)THEN
     WRITE(*,*) nsig, nsp, Sigma0, sigY0, sigYmax
  ENDIF
  ALLOCATE(sig(nsp),stat=ierr)

  ALLOCATE(Psz(nsp),Pszpr(nsp),stat=ierr)
  IF(myrank.EQ.0)THEN
     ALLOCATE(Pszall(nsig),stat=ierr)
  ENDIF
  soindex = -555
  DO si = 1,nsp
     sig(si) = (myrank*nsp + si-1)*dsig   
     Psz(si) = 0.d0
     Pszpr(si) = 0.d0
     ! Find index for macroscopic load and apply I.C.
     IF((sig(si).GE.Sigma0).AND.(sig(si).LT.(Sigma0+dsig)))THEN
        Psz(si) = 1.d0/dsig
        Pszpr(si) = 1.d0/dsig
        soindex = si
     ENDIF
  END DO

  ! Find index for macroscopic load and apply I.C.  
  DO zi = 1,nz
     ! Find index for yield stress boundary
     DO si = 1,nsp
        IF((syindex(zi).LT.0).AND.(sig(si).GE.sigY(zi)))THEN
           syindex(zi) = si
           syrank(zi) = myrank
        ENDIF
     END DO
  END DO
  
  ! --------------------------------------------------
  ! Precalculate the effective diffusivity for varying
  ! magnifications f(zeta)
  ALLOCATE(fz(nz),stat=ierr)
  CALL CalcFz(fz,Cq,zetas,nz,qL,Ey,nu,dsig)

  ! --------------------------------------------------
  ! Solve diffusion equation with magnification
  ALLOCATE(dpdzo(nsp),stat=ierr)
  ALLOCATE(PszB(NPU),PszE(NPU),stat=ierr)

  ! Initial relative area of plastic contact is zero
  Aplpr = 0.d0
  Apl = 0.d0
  Anon = 0.d0
  idsig = 1.d0/dsig
  i2dsig = idsig*idsig
  zi = 1
  CALL MPI_GATHER(Psz,nsp,MPI_REAL8,Pszall,nsp,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  IF(myrank.EQ.0)THEN
     WRITE(100,400) zetas(zi), Apl, Anon, sigY(zi)  ! Output area of contact at every zeta step   
     CALL PRINTPSD(101,Pszall,zetas(zi),nsig,outstepsize)
  ENDIF
  zetaoprev = zetas(zi)
  
  ! We'll use 2nd-order accurate Heun's predictor-corrector scheme
  DO zi = 2,nz
     
     CALL MPI_Allgather(Pszpr(1),1,MPI_REAL8,PszB,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
     CALL MPI_Allgather(Pszpr(nsp),1,MPI_REAL8,PszE,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
     ! First prediction
     DO si = 1,nsp
        ! Diffusion equation valid between (0,sigmaY(zeta)
        IF((sig(si).GT.0.d0).AND.(sig(si).LT.sigY(zi)))THEN
           IF(si.EQ.1)THEN
              dpdzo(si) = Pszpr(2)-2.d0*Pszpr(si)+PszE(myrank)
           ELSEIF(si.EQ.nsp)THEN
              dpdzo(si) = Pszpr(nsp-1)-2.d0*Pszpr(si)+PszB(myrank+2)
           ELSE
              dpdzo(si) = Pszpr(si+1)-2.d0*Pszpr(si)+Pszpr(si-1)
           ENDIF
           Psz(si) = Pszpr(si) + dmu*zetas(zi-1)*i2dsig*fz(zi-1)*dpdzo(si)   ! hat(y)_n+1 = y_n + dt * f(y_n,t)
        ENDIF
     END DO
     
     ! If dsigma_{Y}/dzeta != 0 then need to solve for B.C. 
     IF(dsigY(zi).NE.0.d0)THEN
        ! Solve for yield boundary condition
        CALL MPI_Allgather(Psz(nsp),1,MPI_REAL8,PszE,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
        IF(myrank.EQ.syrank(zi))THEN
           IF(syindex(zi).EQ.1)THEN
              Pszm = PszE(myrank)
           ELSE
              Pszm = Pszpr(syindex(zi)-1)
           ENDIF
           ! First solve for update in relative area of plastic contact
           alphai = dsigY(zi)/fz(zi)
           betai = 1.d0/(dzeta*dsigY(zi))
           Apl = ( betai*Aplpr - 0.5d0*idsig*Pszm/alphai ) / (alphai + betai)
           
           Psz(syindex) = alphai*Apl
        ENDIF
     ENDIF
     
     CALL MPI_Allgather(Psz(1),1,MPI_REAL8,PszB,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
     CALL MPI_Allgather(Psz(nsp),1,MPI_REAL8,PszE,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
     ! Update prediction
     DO si = 1,nsp
        IF((sig(si).GT.0.d0).AND.(sig(si).LT.sigY(zi)))THEN
           IF(si.EQ.1)THEN
              dpdzo(si) = i2dsig*zetas(zi-1)*fz(zi-1)*dpdzo(si) + &
                   i2dsig * zetas(zi)*fz(zi)*(Psz(2)-2.d0*Psz(si)+PszE(myrank))
           ELSEIF(si.EQ.nsp)THEN
              dpdzo(si) = i2dsig*zetas(zi-1)*fz(zi-1)*dpdzo(si) + &
                   i2dsig *zetas(zi)* fz(zi)*(Psz(nsp-1)-2.d0*Psz(si)+PszB(myrank+2))
           ELSE
              dpdzo(si) = i2dsig * zetas(zi-1)*fz(zi-1)*dpdzo(si) + &
                   i2dsig * zetas(zi)*fz(zi)*(Psz(si+1)-2.d0*Psz(si)+Psz(si-1))
           ENDIF
           Psz(si) = Pszpr(si) + 0.5d0*dmu*dpdzo(si)  ! y_n+1 = y_n + (dt/2) * [ f(y_n,t) + f(hat(y)_n+1,t+1)
        ENDIF
     END DO
     
     ! Solve for updated yield boundary condition
     IF(dsigY(zi).NE.0.d0)THEN
        CALL MPI_Allgather(Psz(nsp),1,MPI_REAL8,PszE,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
        IF(myrank.EQ.syrank(zi))THEN
           IF(syindex(zi).EQ.1)THEN
              Pszm = PszE(myrank)
           ELSE
              Pszm = Psz(syindex(zi)-1)
           ENDIF
        ENDIF
        Apl = ( betai*Aplpr - 0.5d0*idsig*Pszm/alphai ) / (alphai + betai)
        Psz(syindex) = alphai*Apl
     ENDIF
     
     IF(myrank.EQ.0)THEN
        Anon = Anonpr + idsig*fz(zi)*Psz(2)
     ENDIF
     ! Save updated value as previous entry for next step
     DO si = 1,nsp
        Pszpr(si) = Psz(si)
     END DO
     Aplpr = Apl
     Anonpr = Anon
     
     CALL MPI_BCAST(Aplpr,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
     ! Need som routine to save data
     IF(myrank.EQ.0)THEN
        WRITE(100,400) zetas(zi), Apl, Anon, sigY(zi)  ! Output area of contact at every zeta step
     ENDIF
     
     IF(zetas(zi).GE.(zoutfac*zetaoprev))THEN
        CALL MPI_GATHER(Psz,nsp,MPI_REAL8,Pszall,nsp,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        IF(myrank.EQ.0)THEN
           CALL PRINTPSD(101,Pszall,zetas(zi),nsig,outstepsize)
        ENDIF
        zetaoprev = zetas(zi)
     ENDIF
  END DO
  
  CLOSE(101)
  CLOSE(100)
400 format (4(1X,E15.7))  ! 1_contact
  ! Deallocate arrays
  IF(ALLOCATED(Pszpr)) DEALLOCATE(Pszpr)
  IF(ALLOCATED(Psz)) DEALLOCATE(Psz)
  IF(ALLOCATED(dpdzo)) DEALLOCATE(dpdzo)
  IF(ALLOCATED(PszB)) DEALLOCATE(PszB,PszE)
  
  CALL MPI_FINALIZE(ierr)
CONTAINS  
  !!--------------------------------------------------
  !!                  Sub-routines
  !! -------------------------------------------------

  !---------------------------------------------------
  ! Create a self-affine power-spectrum
  SUBROUTINE createPSD(Cq,Co,mexp,zetas,nz)
    IMPLICIT NONE
    INTEGER :: nz,zi
    REAL*8 :: mexp,Co
    REAL*8 :: Cq(nz),zetas(nz)

    DO zi = 1, nz
       Cq(zi) = Co * zetas(zi)**mexp
    END DO

    RETURN
  END SUBROUTINE createPSD

  !--------------------------------------------------- 
  ! create the scale-dependent yield stress
  SUBROUTINE createYieldBoundary(sigY,dsigY,zetas,sigY0,nexp,nz)
    IMPLICIT NONE
    INTEGER :: nz,zi
    REAL*8 :: nexp, sigY0
    REAL*8 :: zetas(nz),sigY(nz),dsigY(nz)

    DO zi = 1,nz
       sigY(zi) = sigY0 * zetas(zi)**nexp
       dsigY(zi) = nexp * sigY0 * zetas(zi)**(nexp-1)
    END DO
    
    RETURN
  END SUBROUTINE createYieldBoundary
  
  !---------------------------------------------------
  ! Calculate the effective diffusivity as a function
  ! of magnification
  SUBROUTINE calcFz(fz,Cq,zetas,nz,qL,Ey,nu,dsig)
    IMPLICIT NONE
    INTEGER :: nz, zi
    REAL*8 :: Ey,nu,qL,coeff,dsig
    REAL*8 :: Cq(nz),fz(nz),dfz(nz),zetas(nz)

    coeff = 0.125d0*(Ey/(1-nu*nu))**2 * qL**3  ! Considering 1-D profile
    DO zi = 1,nz
       fz(zi) = coeff*zetas(zi)**2 * Cq(zi)
    END DO
    
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
