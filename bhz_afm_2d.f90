! DESCRIPTION
!   Solve the non-interacting BHZ model with AFM 2x2 basis 
!   generate the hamiltoniana matrix H(k), 

program bhz_afm
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=2048,Norb=2,Nspin=2,Nlat=4,Nso=Nspin*Norb,Nlso=Nlat*Nso
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath,ktrims,ddk
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8)                                 :: mh,rh,lambda,delta
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                    :: wm,wr
  real(8)                                 :: n(Nlso)
  complex(8)                              :: w
  complex(8)                              :: Hloc(Nlso,Nlso)
  complex(8)                              :: Hlocsite(Nlat,Nso,Nso)
  complex(8)                              :: Gmats(Nlso,Nlso,L),Greal(Nlso,Nlso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag
  integer                                 :: icount,jcount,isite,jsite,ispin,jspin,row,col
  !Dirac matrices:
  complex(8),dimension(4,4)               :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz_afm.in")
  call save_input_file("inputBHZ.conf")

  if(Norb/=2)stop "Norb != 2"
  if(Nspin/=2)stop "Nspin != 2"
  if(Nso/=4)stop "Nso != 4"
  if(Nlso/=16)stop "Nlso != 16"

  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)

  Nktot=Nkx*Nkx
  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(kxgrid(Nkx))
  write(*,*) "Using Nk_total="//txtfy(Nktot)

  Greal = zero
  Gmats = zero 
  Hloc  = zero
  kxgrid = kgrid(Nkx)
  Hk = build_hk_model(hk_model,Nlso,kxgrid,kxgrid,[0d0])
  call write_hk_w90(trim(file),Nlso,&
       Nd=Nso,&
       Np=0,   &
       Nineq=4,&
       hk=Hk,  &
       kxgrid=kxgrid,&
       kygrid=kxgrid,&
       kzgrid=[0d0])
  Hloc=sum(Hk(:,:,:),dim=3)/Nktot
  !call write_Hloc(Hloc)



  !Build the local GF:
  wm = pi/beta*real(2*arange(1,L)-1,8)
  wr = linspace(-10.d0,10.d0,L)
  do ik=1,Nktot
     do i=1,L
        w = dcmplx(wr(i),eps)+xmu
        Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
        w = xi*wm(i)+xmu
        Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(w,Hk(:,:,ik))/Nktot
     enddo
  enddo
  do iorb=1,Nlso
     call splot("Gloc_lso"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_lso"//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
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
       colors_name=[character(len=10) :: 'red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue','red','blue'],&
       points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
       file="Eigenbands_afm.nint")



contains



  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    complex(8),dimension(Nso,Nso) :: M
    complex(8),dimension(Nso,Nso) :: tx,ty,thx,thy
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = 2d0*kpoint(1)
    ky = 2d0*kpoint(2)
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    ! H4 =  | m1                 tx + tx^+.e^-ikx    0                   ty^+ + ty.e^iky |
    !       | tx^+ + tx.e^-ikx   m2                  ty^+ + ty.e^-iky    0               |
    !       | 0                  ty + ty^+e^-iky     m3                  tx^+ + tx.e^ikx |
    !       | ty + ty^+e^-iky    0                   tx + tx^+.e^-ikx    m4              |
    !
    hk(1:4,1:4)    = M
    hk(1:4,5:8)    = tx  + thx*exp(-xi*kx)
    hk(1:4,9:12)   = zero
    hk(1:4,13:16)  = thy + ty*exp(xi*ky)
    !
    hk(5:8,1:4)    = thx + tx*exp(xi*kx)
    hk(5:8,5:8)    = M
    hk(5:8,9:12)   = thy + ty*exp(xi*ky) 
    hk(5:8,13:16)  = zero
    !
    hk(9:12,1:4)   = zero
    hk(9:12,5:8)   = ty  + thy*exp(-xi*ky)
    hk(9:12,9:12)  = M
    hk(9:12,13:16) = thx + tx*exp(xi*kx)
    !
    hk(13:16,1:4)  = ty  + thy*exp(-xi*ky)
    hk(13:16,5:8)  = zero
    hk(13:16,9:12) = tx  + thx*exp(-xi*kx)
    hk(13:16,13:16)= M
  end function hk_model



  function inverse_g0k(zeta,hk) result(g0k)
    complex(8)                      :: zeta
    complex(8),dimension(Nlso,Nlso) :: hk
    complex(8),dimension(Nlso,Nlso) :: g0k
    g0k = -hk
    forall(i=1:16)g0k(i,i)=zeta - hk(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k





  function iso2site(hktmp) result(hk)
    complex(8),dimension(16,16) :: hktmp,hk
    integer                     :: icount,isite,ispin,iband,iorb
    integer                     :: jcount,jsite,jspin,jband,jorb
    integer                     :: row,col
    icount=0
    do isite=1,4
       do ispin=1,2
          do iorb=1,2
             icount=icount+1
             row = (ispin-1)*2*4 + (isite-1)*2 + iorb
             jcount=0
             do jsite=1,4
                do jspin=1,2
                   do jorb=1,2
                      jcount=jcount+1
                      col = (jspin-1)*2*4 + (jsite-1)*2 + jorb
                      hk(icount,jcount) = hktmp(row,col)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function iso2site


  function site2iso(hk) result(hktmp)
    complex(8),dimension(16,16) :: hktmp,hk
    integer                     :: icount,isite,ispin,iband,iorb
    integer                     :: jcount,jsite,jspin,jband,jorb
    integer                     :: row,col
    icount=0
    do isite=1,4
       do ispin=1,2
          do iorb=1,2
             icount=icount+1
             row = (ispin-1)*2*4 + (isite-1)*2 + iorb
             jcount=0
             do jsite=1,4
                do jspin=1,2
                   do jorb=1,2
                      jcount=jcount+1
                      col = (jspin-1)*2*4 + (jsite-1)*2 + jorb
                      hktmp(row,col) = hk(icount,jcount)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function site2iso

end program bhz_afm







