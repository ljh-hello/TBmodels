program hm_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=1024,Norb=1,Nspin=1,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk

  real(8),dimension(:,:,:),allocatable    :: nkgrid
  real(8)                                 :: ts,xmu,beta,eps,wmax
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nso,Nso,L),Greal(Nso,Nso,L),Sfoo(Nso,Nso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,iener


  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(ts,"TS","inputHM.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputHM.conf",default="hkfile_bhz.in")
  call save_input_file("inputHM.conf")





  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  allocate(kxgrid(Nkx))
  write(*,*) "Using Nk_total="//txtfy(Nktot)


  Hloc  = zero
  kxgrid = kgrid(Nkx)
  Hk = TB_build_model(hk_model,Nso,kxgrid,kxgrid,[0d0])
  Wtk = 1d0/Nktot
  Hloc=sum(Hk(:,:,:),dim=3)/Nktot

  call write_hk_w90(trim(file),Nso,&
       Nd=Norb,&
       Np=1,   &
       Nineq=1,&
       hk=Hk,  &
       kxgrid=kxgrid,&
       kygrid=kxgrid,&
       kzgrid=[0d0])

  call write_Hloc(Hloc)

  Greal = zero
  Gmats = zero 
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  Sfoo =zero
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=1)
  Sfoo =zero
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=1)
  do iorb=1,Nso
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo

  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_M1
  kpath(3,:)=kpoint_X1
  kpath(4,:)=kpoint_Gamma
  call TB_Solve_path(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1],&
       points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")


  !plot observables
  open(10,file="observables.nint")
  write(10,"(10F20.12)")(n(iorb),iorb=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(n(iorb),iorb=1,Nso),sum(n)

contains




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -2d0*ts*(cos(kx)+cos(ky))
  end function hk_model




end program hm_2d


