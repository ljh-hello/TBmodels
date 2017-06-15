program grapheneAF
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: Norb=2,Nspin=2,Nso=Norb*Nspin
  integer                                 :: Nk,Nktot,Lfreq,Nkpath
  real(8)                                 :: ts,tsp,phi,delta,xmu,beta,eps,wmax,Mh,Ndens(Nso)
  real(8)                                 :: a0,bklen
  integer                                 :: i,j,k,ik,ix,iy,iorb,jorb
  real(8)                                 :: kx,ky,Eshift
  real(8),dimension(:),allocatable        :: kxgrid,kygrid
  real(8),dimension(2)                    :: d1,d2,d3
  real(8),dimension(2)                    :: a1,a2,a3
  real(8),dimension(2)                    :: bk1,bk2,pointK,pointKp,kvec
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:,:),allocatable :: Gmats,Greal,Sfoo
  real(8),dimension(:),allocatable        :: wm,wr,Wtk
  real(8),dimension(:,:),allocatable      :: KPath
  complex(8)                              :: zeta
  complex(8),dimension(4,4)               :: GammaX,GammaY,GammaZ,Gamma0

  call parse_input_variable(Nk,"NK","inputGRAPHENE.conf",default=20)
  call parse_input_variable(nkpath,"NKPATH","inputGRAPHENE.conf",default=100)
  call parse_input_variable(ts,"TS","inputGRAPHENE.conf",default=1d0)
  call parse_input_variable(tsp,"TSP","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(mh,"MH","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(phi,"PHI","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(xmu,"XMU","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(wmax,"WMAX","inputGRAPHENE.conf",default=10d0)
  call parse_input_variable(eps,"EPS","inputGRAPHENE.conf",default=3.d-2)
  call parse_input_variable(beta,"BETA","inputGRAPHENE.conf",default=1000d0)
  call parse_input_variable(Lfreq,"LFREQ","inputGRAPHENE.conf",default=2000)
  call save_input_file("inputGRAPHENE.conf")
  phi=phi*pi

  Nktot=Nk*Nk
  write(*,"(A)")"Using Nk="//txtfy(Nktot)

  Gamma0 = kron_pauli(pauli_0,pauli_0)
  GammaX = kron_pauli(pauli_0,pauli_x)
  GammaY = kron_pauli(pauli_0,pauli_y)
  GammaZ = kron_pauli(pauli_0,pauli_z)


  !Lattice basis (see graphene.f90 for backupd version)
  !LATTICE BASIS:
  ! nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]
  !
  !
  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d2-d3                    !a*sqrt(3)[sqrt(3)/2,-1/2]
  a2 = d3-d1                    !a*sqrt(3)[-sqrt(3)/2,-1/2]
  a3 = d1-d2                    !a*sqrt(3)[0, 1]
  !
  !
  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/sqrt(3d0)
  bk1=bklen*[ sqrt(3d0)/2d0 ,  1d0/2d0 ]
  bk2=bklen*[ sqrt(3d0)/2d0 , -1d0/2d0 ]

  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


  Eshift=0d0

  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  allocate(kxgrid(Nk),kygrid(Nk))
  print*,"Build Hk(Nso,Nso) for the Graphene model"
  ik=0
  do iy=1,Nk
     ky = dble(iy-1)/Nk
     do ix=1,Nk
        ik=ik+1
        kx=dble(ix-1)/Nk
        kvec = kx*bk1 + ky*bk2
        kxgrid(ix) = kvec(1)
        kygrid(iy) = kvec(2)
        Hk(:,:,ik) = grapheneAF_model(kvec,Nso)
     enddo
  enddo
  Wtk = 1d0/Nktot


  print*,sum(Hk,3)/dble(Nktot)

  !THIS IS WRONG BECUASE THE GRIDS KXGRID,KYGRID
  !CAN ONLY BE DESCRIBED AS VECTORS GRID (i.e.
  ! you can not store the information with a single 
  ! real number) 
  ! call write_hk_w90("Hkrfile_grapheneAF.data",&
  !      No=Nso,&
  !      Nd=Norb,&
  !      Np=0,&
  !      Nineq=1,&
  !      Hk=Hk,&
  !      kxgrid=kxgrid,kygrid=kygrid,kzgrid=[0d0])


  Eshift=xmu

  allocate(Kpath(4,2))
  KPath(1,:)=[0,0]
  KPath(2,:)=pointK
  Kpath(3,:)=pointKp
  KPath(4,:)=[0d0,0d0]
  call TB_Solve_path(grapheneAF_model,Nso,KPath,Nkpath,&
       colors_name=[red1,green1,blue1,yellow1],&
       points_name=[character(len=10) :: "G","K","K`","G"],&
       file="Eigenbands.nint")



  ! allocate(wm(Lfreq),wr(Lfreq))
  allocate(Gmats(Nso,Nso,Lfreq),Greal(Nso,Nso,Lfreq),Sfoo(Nso,Nso,Lfreq))
  Gmats=zero
  Greal=zero
  Sfoo =zero
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=1)
  Sfoo=zero
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=1)


  do iorb=1,Nso
     Ndens(iorb)= fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  open(10,file="density.nint")
  write(10,"(10F20.12)")(Ndens(iorb),iorb=1,Nso),sum(Ndens)
  close(10)

  ! write(*,"(A,10F14.9)")"Occupations =",(Ndens(iorb),iorb=1,Nso),sum(Ndens)


contains


  function grapheneAF_model(kpoint,Nso) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: hk
    real(8)                       :: h0,hx,hy,hz
    real(8)                       :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = 2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = 2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh 
    hk = h0*Gamma0 + hx*GammaX + hy*GammaY + hz*GammaZ
  end function grapheneAF_model





end program grapheneAF
