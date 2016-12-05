program hc_haldaneAF
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: Norb=1,Nspin=2,Nlat=2,Nso=Nspin*Norb,Nlso=Nlat*Nso
  integer                                 :: Nk,Nktot,Lfreq,Nkpath
  real(8)                                 :: ts,tsp,phi,delta,xmu,beta,eps,wmax,Mh,Ndens(Nlso)
  integer                                 :: i,j,k,ik,ix,iy,iorb,jorb
  real(8)                                 :: a,a0
  real(8)                                 :: kx,ky,Eshift,bklen
  real(8),dimension(:),allocatable        :: kxgrid,kygrid
  real(8),dimension(2)                    :: d1,d2,d3
  real(8),dimension(2)                    :: a1,a2,a3
  real(8),dimension(2)                    :: bk1,bk2,pointK,pointKp,kvec
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  complex(8),dimension(:,:,:),allocatable :: Gmats,Greal,Sfoo
  real(8),dimension(:,:),allocatable      :: KPath,berry_curvature
  complex(8)                              :: zeta,Eigvec(Nlso,Nlso)
  real(8)                                 :: chern,BZ_area,eigval(Nlso)
  real(8),dimension(:,:,:),allocatable    :: nkgrid
  complex(8),dimension(:,:,:),allocatable :: BlochStates
  complex(8),dimension(4,4)               :: GammaX,GammaY,GammaZ,Gamma0

  call parse_input_variable(Nk,"NK","inputHALDANE.conf",default=20)
  call parse_input_variable(nkpath,"NKPATH","inputHALDANE.conf",default=100)
  call parse_input_variable(ts,"TS","inputHALDANE.conf",default=1d0)
  call parse_input_variable(tsp,"TSP","inputHALDANE.conf",default=0d0)!1.d0/3/sqrt(3d0)
  call parse_input_variable(mh,"MH","inputHALDANE.conf",default=0d0)
  call parse_input_variable(phi,"PHI","inputHALDANE.conf",default=0d0,comment="In units of \pi")
  call parse_input_variable(xmu,"XMU","inputHALDANE.conf",default=0d0)
  call parse_input_variable(wmax,"WMAX","inputHALDANE.conf",default=10d0)
  call parse_input_variable(eps,"EPS","inputHALDANE.conf",default=3.d-2)
  call parse_input_variable(beta,"BETA","inputHALDANE.conf",default=1000d0)
  call parse_input_variable(Lfreq,"LFREQ","inputHALDANE.conf",default=2000)
  call save_input_file("inputHALDANE.conf")
  phi=phi*pi

  Nktot=Nk*Nk
  write(*,"(A)")"Using Nk="//txtfy(Nktot)

  Gamma0 = kron_pauli(pauli_0,pauli_0)
  GammaX = kron_pauli(pauli_0,pauli_x)
  GammaY = kron_pauli(pauli_0,pauli_y)
  GammaZ = kron_pauli(pauli_0,pauli_z)


  a=1d0
  a0=a*sqrt(3d0)

  !Lattice basis (a=1; a0=sqrt3*a) is:
  !\a_1 = a0 [ sqrt3/2 , 1/2 ]
  !\a_2 = a0 [ sqrt3/2 ,-1/2 ]
  !
  !
  !nearest neighbor: A-->B, B-->A
  d1= a*[  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= a*[  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= a*[ -1d0     , 0d0           ]
  !
  !next nearest-neighbor displacements: A-->A, B-->B \== \nu_1,\nu_2, \nu_3=\nu_1-\nu_2
  a1=a0*[ sqrt(3d0)/2d0, 1d0/2d0]
  a2=a0*[ sqrt(3d0)/2d0,-1d0/2d0]
  a3=a2-a1


  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/3d0
  bk1=bklen*[ 1d0/2d0 ,  sqrt(3d0)/2d0 ]
  bk2=bklen*[ 1d0/2d0 , -sqrt(3d0)/2d0 ]


  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


  Eshift=0d0
  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(Wtk(Nktot))
  allocate(kxgrid(Nk),kygrid(Nk))
  print*,"Build Hk(Nlso,Nlso) for the Graphene model"
  ik=0
  do iy=1,Nk
     ky = dble(iy-1)/Nk
     do ix=1,Nk
        ik=ik+1
        kx=dble(ix-1)/Nk
        kvec = kx*bk1 + ky*bk2
        kxgrid(ix)=kvec(1)
        kygrid(iy)=kvec(2) 
        Hk(:,:,ik) = hk_haldaneAF_model(kvec,Nlso)
     enddo
  enddo
  Wtk = 1d0/Nktot

  print*,sum(Hk,3)/dble(Nktot)

  !THIS IS WRONG BECUASE THE GRIDS KXGRID,KYGRID
  !CAN ONLY BE DESCRIBED AS VECTORS GRID (i.e.
  ! you can not store the information with a single 
  ! real number) 
  ! call write_hk_w90("Hkrfile_Haldane.data",&
  !      No=Nlso,&
  !      Nd=Nso,&
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
  call TB_Solve_path(hk_haldaneAF_model,Nlso,KPath,Nkpath,&
       colors_name=[red1,blue1],&
       points_name=[character(len=10) :: "G","K","K`","G"],&
       file="Eigenbands.nint")




  ! allocate(wm(Lfreq),wr(Lfreq))
  allocate(Gmats(Nlso,Nlso,Lfreq),Greal(Nlso,Nlso,Lfreq),Sfoo(Nlso,Nlso,Lfreq))
  Gmats=zero
  Greal=zero
  Sfoo =zero
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=1)
  Sfoo=zero
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=1)

  do iorb=1,Nlso
     Ndens(iorb)= fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  open(10,file="density.nint")
  write(10,"(10F20.12)")(Ndens(iorb),iorb=1,Nlso),sum(Ndens)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(Ndens(iorb),iorb=1,Nlso),sum(Ndens)


  call get_Chern_number_NEW(Hk,[Nk,Nk],2,Nk/pi2*Nk/pi2,Chern)

contains


  function hk_haldaneAF_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
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
    h0 = -2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = -2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh 
    !hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
    hk = h0*Gamma0 + hx*GammaX + hy*GammaY + hz*GammaZ
  end function hk_haldaneAF_model




  ! calcola il numero di chern di un generico stato dipendente da k con il metodo di Resta
  subroutine Get_Chern_number_NEW(Hk,Nkvec,Noccupied,one_over_area,Chern)
    complex(8),intent(in),dimension(:,:,:)    :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)           :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    integer,intent(in)                        :: Noccupied
    real(8),intent(in)                        :: one_over_area
    real(8),intent(out)                       :: Chern
    !
    integer                                   :: Nlso
    integer                                   :: Nktot
    integer                                   :: Nkx,Nky
    integer                                   :: ikx,iky
    integer                                   :: ikxP,ikyP
    integer                                   :: ik,iocc
    complex(8),dimension(:,:),allocatable     :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable          :: Eigval ![Nlso]
    complex(8),dimension(:,:),allocatable     :: Gmat
    complex(8),dimension(:,:,:,:),allocatable :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                   :: Ulink
    real(8),dimension(:,:),allocatable        :: BerryCurvature
    real(8)                                   :: berry_phase
    integer                                   :: unit
    !
    Nlso  = size(Hk,1)
    Nktot = size(Hk,3)
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    call assert_shape(Hk,[Nlso,Nlso,Nktot],"Get_Chern_NUmber_NEW","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number_NEW: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nlso,Nlso))
    allocate(Eigval(Nlso))
    allocate(BlochStates(Nkx,Nky,Noccupied,Nlso))
    allocate(BerryCurvature(Nkx,Nky))
    allocate(Gmat(Noccupied,Noccupied))
    ik=0
    do ikx=1,Nkx
       do iky=1,Nky
          ik=ik+1
          Eigvec = Hk(:,:,ik)
          call eigh(Eigvec,Eigval)
          do iocc=1,Noccupied
             BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
          enddo
       enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
       ikxP = modulo(ikx,Nkx) + 1
       !ikxM = modulo(ikx-2,Nkx) + 1
       do iky= 1, Nky
          ikyP = modulo(iky,Nky) + 1
          !ikyM = modulo(iky-2,Nky) + 1
          !
          if(Noccupied==1)then
             Ulink(1) = dot_product(BlochStates(ikx,iky,1,:)  , BlochStates(ikx,ikyP,1,:))
             Ulink(2) = dot_product(BlochStates(ikx,ikyP,1,:) , BlochStates(ikxP,ikyP,1,:))
             Ulink(3) = dot_product(BlochStates(ikxP,ikyP,1,:), BlochStates(ikxP,iky,1,:))
             Ulink(4) = dot_product(BlochStates(ikxP,iky,1,:) , BlochStates(ikx,iky,1,:))
             !
          else
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikx,iky,i,:)  , BlochStates(ikx,ikyP,j,:))
                enddo
             enddo
             Ulink(1) = det(gmat)
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikx,ikyP,i,:) , BlochStates(ikxP,ikyP,j,:))
                enddo
             enddo
             Ulink(2) = det(gmat)
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikxP,ikyP,i,:), BlochStates(ikxP,iky,j,:))
                enddo
             enddo
             Ulink(3) = det(gmat)
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikxP,iky,i,:) , BlochStates(ikx,iky,j,:))
                enddo
             enddo
             Ulink(4) = det(gmat)
             !
          endif
          !
          berry_phase = dimag(zlog( product(Ulink(:))  ))
          chern = chern + berry_phase
          BerryCurvature(ikx,iky) = berry_phase*one_over_area
          !
       enddo
    enddo
    !
    chern=chern/pi2
    !
    open(unit=free_unit(unit),file="Chern_Number.dat")
    write(unit,*)chern
    write(*,*)chern
    close(unit)
    !
    call splot3d("Berry_Curvature.dat",&
         linspace(0d0,pi2,Nkx,iend=.false.),&
         linspace(0d0,pi2,Nky,iend=.false.),&
         BerryCurvature(:Nkx,:Nky))
    !
  end subroutine Get_Chern_number_NEW



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
          berry_phase = dimag(zlog( product(path_scalar_products)  ))
          berry_curvatures(i1,i2) = berry_phase*one_over_area
          chern_number = chern_number + berry_phase
       enddo
    enddo
    chern_number = chern_number/pi2
  end subroutine get_chern_number










end program hc_haldaneAF
