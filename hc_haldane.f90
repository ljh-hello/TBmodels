program hc_haldane
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: Norb=2,Nso=Norb
  integer                                 :: Nk,Nktot,Lfreq,Nkpath
  real(8)                                 :: ts,tsp,phi,delta,xmu,beta,eps,wmax,Mh,Ndens(Nso)
  integer                                 :: i,j,k,ik,ix,iy,iorb,jorb
  real(8)                                 :: kx,ky,Eshift
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(2)                    :: a1,a2,a3
  real(8),dimension(2)                    :: b1,b2,b3
  complex(8),dimension(:,:,:),allocatable :: Hk
  complex(8),dimension(:,:,:),allocatable :: Gmats,Greal
  real(8),dimension(:),allocatable        :: wm,wr
  real(8),dimension(:,:),allocatable      :: KPath,berry_curvature
  complex(8)                              :: zeta,Eigvec(Norb,Norb)
  real(8)                                 :: chern,BZ_area,eigval(Norb)
  real(8),dimension(:,:,:),allocatable    :: nkgrid
  complex(8),dimension(:,:,:),allocatable :: State


  call parse_input_variable(Nk,"NK","inputHALDANE.conf",default=20)
  call parse_input_variable(nkpath,"NKPATH","inputHALDANE.conf",default=100)
  call parse_input_variable(ts,"TS","inputHALDANE.conf",default=1d0)
  call parse_input_variable(tsp,"TSP","inputHALDANE.conf",default=1.d0/3/sqrt(3d0))
  call parse_input_variable(mh,"MH","inputHALDANE.conf",default=0d0)
  call parse_input_variable(phi,"PHI","inputHALDANE.conf",default=0d0)
  call parse_input_variable(xmu,"XMU","inputHALDANE.conf",default=0d0)
  call parse_input_variable(wmax,"WMAX","inputHALDANE.conf",default=10d0)
  call parse_input_variable(eps,"EPS","inputHALDANE.conf",default=3.d-2)
  call parse_input_variable(beta,"BETA","inputHALDANE.conf",default=1000d0)
  call parse_input_variable(Lfreq,"LFREQ","inputHALDANE.conf",default=2000)
  call save_input_file("inputHALDANE.conf")
  phi=phi*pi

  Nktot=Nk*Nk
  write(*,"(A)")"Using Nk="//txtfy(Nktot)


  !Following:
  !http://www-personal.umich.edu/~sunkai/teaching/Fall_2012/chapter3_part7.pdf
  !Lattice basis (a=1) is:
  !\nu_1 = [ sqrt3  , 0  ]
  !\nu=2 = [-sqrt3/2,3/2 ]
  !nearest neighbor: A-->B, B-->A
  a1 = [0d0           ,  1d0    ]
  a2 = [-sqrt(3d0)/2d0, -1d0/2d0]
  a3 = [ sqrt(3d0)/2d0, -1d0/2d0]

  !next nearest-neighbor displacements: A-->A, B-->B \== \nu_1,\nu_2, \nu_3=\nu_1-\nu_2
  b1=[ sqrt(3d0), 0d0]
  b2=[-sqrt(3d0)/2d0, 3d0/2d0]
  b3=[-sqrt(3d0)/2d0,-3d0/2d0]

  eshift=0d0

  allocate(kxgrid(Nk))
  allocate(Hk(Nso,Nso,Nktot))

  print*,"Build Hk(Nso,Nso) for the Haldane model"
  kxgrid = kgrid(Nk)
  Hk = TB_build_model(haldane_model,Nso,kxgrid,kxgrid,[0d0])

  call write_hk_w90("Hkrfile_Haldane.data",&
       No=Nso,&
       Nd=Norb,&
       Np=0,&
       Nineq=1,&
       Hk=Hk,&
       kxgrid=kxgrid,kygrid=kxgrid,kzgrid=[0d0])


  Eshift=xmu

  allocate(Kpath(4,2))
  KPath(1,:)=[0,0]
  KPath(2,:)=[2d0/3d0/sqrt(3d0),0d0]  !K K`=-K
  Kpath(3,:)=[1d0/3d0/sqrt(3d0),1/3d0]
  KPath(4,:)=[0,0]
  KPath=Kpath*pi2
  call TB_Solve_path(haldane_model,Nso,KPath,Nkpath,&
       colors_name=[red1,blue1],&
       points_name=[character(len=10) :: "G","K","K`","G"],&
       file="Eigenbands.nint")




  allocate(wm(Lfreq),wr(Lfreq))
  allocate(Gmats(Nso,Nso,Lfreq),Greal(Nso,Nso,Lfreq))
  wm = pi/beta*(2*arange(1,Lfreq)-1)
  wr = linspace(-wmax,wmax,Lfreq)
  call start_timer
  do ik=1,Nktot
     do i=1,Lfreq
        zeta = dcmplx(wr(i),eps)+xmu
        Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(zeta,Hk(:,:,ik))/Nktot
        zeta = xi*wm(i)+xmu
        Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(zeta,Hk(:,:,ik))/Nktot
     enddo
     call eta(ik,Nktot)
  enddo
  call stop_timer

  do iorb=1,Nso
     Ndens(iorb)= fft_get_density(Gmats(iorb,iorb,:),beta)
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_iw.nint",wm,Gmats(iorb,iorb,:))
     call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(iorb))//"_realw.nint",wr,&
          -dimag(Greal(iorb,iorb,:))/pi,dreal(Greal(iorb,iorb,:)))
  enddo
  open(10,file="density.nint")
  write(10,"(10F20.12)")(Ndens(iorb),iorb=1,Nso),sum(Ndens)
  close(10)

  allocate(state(Norb,Nk,Nk))
  allocate(Berry_curvature(Nk,Nk))
  chern=0d0
  ik=0
  do i=1,Nk
     kx = kxgrid(i)
     do j=1,Nk
        ky = kxgrid(j)
        ik=ik+1
        Eigvec = Hk(:,:,ik)
        call eigh(Eigvec,Eigval)
        state(:,i,j) = Eigvec(:,1)
     enddo
  enddo

  ! BZ_area=(4*pi)**2/6/sqrt(3d0)
  ! chern=-simps2d(chern_nk,[-pi,pi-pi2/Nk],[-pi,pi-pi2/Nk])*BZ_area/pi2/pi2
  ! print*,chern
  call Get_Chern_number(state,chern,Berry_curvature,Nk/pi2*Nk/pi2)
  chern=chern/3
  call splot3d("Berry_Curvature.nint",kxgrid,kxgrid,Berry_Curvature)

  open(10,file="chern.nint")
  write(10,*)int(chern)
  close(10)

  write(*,"(A,10F14.9)")"Occupations =",(Ndens(iorb),iorb=1,Nso),sum(Ndens)
  write(*,"(A,I)")"Chern num.  =",nint(chern)


contains


  function haldane_model(kpoint,Nso) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: hk
    real(8)                       :: h0,hx,hy,hz
    real(8)                       :: kdota(3),kdotb(3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !(k.b_j)
    kdotb(1) = dot_product(kpoint,b1)
    kdotb(2) = dot_product(kpoint,b2)
    kdotb(3) = dot_product(kpoint,b3)
    !
    h0 = 2*tsp*cos(phi)*sum( cos(kdotb(:)) )
    hx =-ts*sum( cos(kdota(:)) )
    hy =-ts*sum( sin(kdota(:)) )
    hz = 2*tsp*sin(phi)*sum( sin(kdotb(:)) ) + Mh 
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
  end function haldane_model



  function nd_model(kpoint,M) result(dk)
    real(8),dimension(:),intent(in) :: kpoint
    integer                         :: M
    real(8),dimension(M)            :: dk
    real(8)                         :: kx,ky,norm
    real(8)                         :: h0,hx,hy,hz
    real(8)                         :: kdota(3),kdotb(3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !(k.b_j)
    kdotb(1) = dot_product(kpoint,b1)
    kdotb(2) = dot_product(kpoint,b2)
    kdotb(3) = dot_product(kpoint,b3)
    !
    h0 = 2*tsp*cos(phi)*sum( cos(kdotb(:)) )
    hx =-ts*sum( cos(kdota(:)) )
    hy =-ts*sum( sin(kdota(:)) )
    hz = 2*tsp*sin(phi)*sum( sin(kdotb(:)) ) + Mh 
    dk=[hx,hy,hz]
    norm = dot_product(dk,dk)
    dk = dk/sqrt(norm)
    where(abs(dk)<1.d-12)dk=0d0
  end function nd_model

  function chern_nk(kpoint) result(ck)
    real(8),dimension(:) :: kpoint
    real(8)              :: dk(3),dk_(3)
    real(8)              :: ddk(3,2)
    real(8)              :: ck,norm
    dk  = nd_model(kpoint,3)
    call djac_dk(kpoint,3,ddk)
    ck  = s3_product(dk,ddk(:,1),ddk(:,2))/4d0/pi
  end function chern_nk

  subroutine djac_dk(kpoint,M,ddk)
    real(8),dimension(:)            :: kpoint
    real(8),dimension(size(kpoint)) :: k_
    integer                         :: M
    real(8),dimension(M)            :: fvec,wa1
    real(8)                         :: ddk(M,size(kpoint))
    call djacobian(nd_model,kpoint,M,ddk)
  end subroutine djac_dk




  ! calcola il numero di chern di un generico stato dipendente da k con il metodo di Resta
  subroutine Get_Chern_number(State,Chern_number,Berry_curvatures,one_over_area)
    complex(8),intent(in),dimension(:,:,:)                     :: state !(nhilb, nk1, nk2)
    real(8),intent(out)                                        :: chern_number
    real(8),intent(out),dimension(size(state,2),size(state,3)) :: berry_curvatures
    real(8),intent(in)                                         :: one_over_area
    integer                                                    :: nhilb,nk1,nk2
    integer                                                    :: i1,i2,i3,ix,i1p,i1m,i2m,i2p,it
    complex(8)                                                 :: path_scalar_products(4)
    real(8)                                                    :: berry_phase
    !
    Nhilb= size(state,1)
    Nk1  = size(state,2)
    Nk2  = size(state,3)
    chern_number = zero
    do i1= 1, nk1
       i1p = modulo(i1,nk1) + 1
       i1m = modulo(i1-2,nk1) + 1
       do i2= 1, nk2           !faccio l'integrale sulla bz
          i2p = modulo(i2,nk2) + 1
          i2m = modulo(i2-2,nk2) + 1
          path_scalar_products(1) = dot_product(state(:,i1,i2),state(:,i1, i2p))
          path_scalar_products(2) = dot_product(state(:,i1,i2p),state(:,i1p, i2p))
          path_scalar_products(3) = dot_product(state(:,i1p,i2p),state(:,i1p, i2))
          path_scalar_products(4) = dot_product(state(:,i1p,i2),state(:,i1,i2))
          berry_phase = -dimag(zlog( product(path_scalar_products)  ))
          berry_curvatures(i1,i2) = berry_phase*one_over_area
          chern_number = chern_number + berry_phase
       enddo
    enddo
    chern_number = chern_number/pi2
  end subroutine get_chern_number






  function inverse_g0k(zeta,hk) result(g0k)
    complex(8)                    :: zeta
    complex(8),dimension(Nso,Nso) :: hk
    integer                       :: i
    complex(8),dimension(Nso,Nso) :: g0k
    complex(8)                    :: delta,ppi,vmix12,vmix21
    g0k=zero
    delta = zeta - hk(1,1)
    ppi   = zeta - hk(2,2)
    vmix12  = -hk(1,2)
    vmix21  = -hk(2,1)
    g0k(1,1) = 1d0/(delta - vmix12*vmix21/ppi)
    g0k(2,2) = 1d0/(ppi - vmix12*vmix21/delta)
    g0k(1,2) = -vmix12/(ppi*delta - vmix12*vmix21)
    g0k(2,1) = -vmix21/(ppi*delta - vmix12*vmix21)
  end function inverse_g0k



end program hc_haldane
