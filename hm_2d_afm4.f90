! DESCRIPTION
!   Solve the non-interacting Hubbard model with AFM 2x2 basis 

program hm_afm_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                       :: L=2048,Norb=1,Nspin=1,Nlat=4,Nso=Nspin*Norb,Nlso=Nlat*Nso
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath,ktrims,ddk
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  real(8)                                 :: ts
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc,wmax
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nlso),h0(4,4)
  complex(8)                              :: w
  complex(8)                              :: Hloc(Nlso,Nlso)
  complex(8)                              :: Hlocsite(Nlat,Nso,Nso)
  complex(8)                              :: Gmats(Nlso,Nlso,L),Greal(Nlso,Nlso,L),Sfoo(Nlso,Nlso,L)
  character(len=32)                       :: file



  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(ts,"TS","inputHM.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputHM.conf",default="hkfile_bhz.in")
  call save_input_file("inputHM.conf")

  Nktot=Nkx*Nkx
  write(*,*) "Using Nk_total="//txtfy(Nktot)


  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(Wtk(Nktot))
  allocate(kxgrid(Nkx))
  Hloc   = zero
  kxgrid = kgrid(Nkx)
  Hk = TB_build_model(hk_model,Nlso,kxgrid,kxgrid,[0d0])

  Hloc=sum(Hk(:,:,:),dim=3)/Nktot
  Wtk = 1d0/Nktot


  !Build the local GF:
  Gmats=zero
  Greal=zero
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  Sfoo =zero
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=1)
  Sfoo =zero
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=1)

  do iorb=1,Nlso
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  !plot observables
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)



  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_X1
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_Gamma
  call solve_Hk_along_BZpath(Hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: '$\Gamma$', 'X', 'M', '$\Gamma$'],&
       file="Eigenbands_afm.nint")


contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = 2d0*kpoint(1)
    ky = 2d0*kpoint(2)
    !
    ! Hk =  | 0                  tx + tx^+.e^-ikx    0                   ty^+ + ty.e^iky |
    !       | tx^+ + tx.e^-ikx   0                   ty^+ + ty.e^-iky    0               |
    !       | 0                  ty + ty^+e^-iky     0                   tx^+ + tx.e^ikx |
    !       | ty + ty^+e^-iky    0                   tx + tx^+.e^-ikx    0               |
    !
    h0 = reshape([&
         0d0,1d0,0d0,1d0, &
         1d0,0d0,1d0,0d0, &
         0d0,1d0,0d0,1d0, &
         1d0,0d0,1d0,0d0],&
         [4,4])
    tk = reshape([ &
         zero, exp(-xi*kx), zero, exp(-xi*ky) , &
         exp(xi*kx), zero, exp(-xi*ky), zero,&
         zero, exp(xi*ky), zero, exp(xi*kx),&
         exp(xi*ky), zero, exp(-xi*kx), zero ],&
         [4,4])
    hk = ts*(h0 + tk)
  end function hk_model




end program hm_afm_2d







