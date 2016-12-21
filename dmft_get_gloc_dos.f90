MODULE DMFT_GET_GLOC_DOS
  USE SCIFOR
  USE DMFT_CTRL_VARS
  implicit none
  private


  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  interface select_block
     module procedure :: select_block_Nlso
     module procedure :: select_block_nnn
  end interface select_block
  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape
  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape
  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape
  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape


  real(8),dimension(:),allocatable          :: wm !Matsubara frequencies
  real(8),dimension(:),allocatable          :: wr !Real frequencies
  character(len=128)                        :: suffix
  character(len=128)                        :: gloc_suffix=".dat"
  integer                                   :: Lk,Nlso,Nlat,Nspin,Norb,Nso,Lreal,Lmats
  integer                                   :: i,j,ik,ilat,jlat,iorb,jorb,ispin,jspin,io,jo,is
  !

  public :: dmft_get_gloc_matsubara_normal_dos
  public :: dmft_get_gloc_realaxis_normal_dos

contains


  subroutine dmft_get_gloc_matsubara_normal_dos(Ebands,Dbands,Hloc,Gmats,Smats,iprint)
    real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)                  :: Smats     ![Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout)               :: Gmats     !as Smats
    integer,intent(in)                                          :: iprint
    !
    complex(8)                                                  :: gktmp
    complex(8),dimension(:,:,:),allocatable                     :: zeta_mats ![Nspin*Norb][Nspin*Norb][Lmats]
    !
    real(8)                                                     :: beta
    real(8)                                                     :: xmu,eps
    real(8)                                                     :: wini,wfin
    !
    !Retrieve parameters:
    call get_ctrl_var(beta,"BETA")
    call get_ctrl_var(xmu,"XMU")
    !
    Nspin = size(Smats,1)
    Norb  = size(Smats,3)
    Lmats = size(Smats,5)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_matsubara_normal_dos',"Ebands")
    call assert_shape(Smats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_dos',"Smats")
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,Lmats],'dmft_get_gloc_matsubara_normal_dos',"Gmats")
    !
    !Allocate and setup the Matsubara freq.
    allocate(zeta_mats(Nso,Nso,Lmats))
    if(allocated(wm))deallocate(wm);allocate(wm(Lmats))
    wm = pi/beta*dble(2*arange(1,Lmats)-1)
    !
    do i=1,Lmats
       zeta_mats(:,:,i)=(xi*wm(i)+xmu)*eye(Nso) - nn2so_reshape(Smats(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    call start_timer
    Gmats=zero
    do i=1,Lmats
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             do ik=1,Lk
                gktmp = Dbands(io,ik)/( zeta_mats(io,io,i)-Hloc(io)-Ebands(io,ik) )
                Gmats(ispin,ispin,iorb,iorb,i) = Gmats(ispin,ispin,iorb,iorb,i) + gktmp
             enddo
          enddo
       enddo
       call eta(i,Lmats)
    end do
    call stop_timer
    call dmft_gloc_print_matsubara(wm,Gmats,"Gloc",iprint)
  end subroutine dmft_get_gloc_matsubara_normal_dos




  subroutine dmft_get_gloc_realaxis_normal_dos(Ebands,Dbands,Hloc,Greal,Sreal,iprint)
    real(8),dimension(:,:),intent(in)                           :: Ebands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1),size(Ebands,2)),intent(in) :: Dbands    ![Nspin*Norb][Lk]
    real(8),dimension(size(Ebands,1)),intent(in)                :: Hloc      ![Nspin*Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)                  :: Sreal     ![Nspin][Nspin][Norb][Norb][Lreal]
    complex(8),dimension(:,:,:,:,:),intent(inout)               :: Greal     !as Sreal
    integer,intent(in)                                          :: iprint
    !
    complex(8)                                                  :: gktmp
    complex(8),dimension(:,:,:),allocatable                     :: zeta_real ![Nspin*Norb][Nspin*Norb][Lreal]
    !
    real(8)                                                     :: xmu,eps
    real(8)                                                     :: wini,wfin
    !
    !Retrieve parameters:
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(wini,"WINI")
    call get_ctrl_var(wfin,"WFIN")
    call get_ctrl_var(eps,"EPS")
    !
    Nspin = size(Sreal,1)
    Norb  = size(Sreal,3)
    Lreal = size(Sreal,5)
    Lk    = size(Ebands,2)
    Nso   = Nspin*Norb    
    !Testing part:
    call assert_shape(Ebands,[Nso,Lk],'dmft_get_gloc_realaxis_normal_dos',"Ebands")
    call assert_shape(Sreal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos',"Sreal")
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,Lreal],'dmft_get_gloc_realaxis_normal_dos',"Greal")
    !
    !Allocate and setup the Realaxis freq.
    allocate(zeta_real(Nso,Nso,Lreal))
    if(allocated(wr))deallocate(wr);allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    do i=1,Lreal
       zeta_real(:,:,i)=(wr(i)+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(:,:,:,:,i),Nspin,Norb)
    enddo
    !
    !invert (Z-Hk) for each k-point
    write(*,"(A)")"Get local Green's function (print mode:"//reg(txtfy(iprint))//")"
    call start_timer
    Greal=zero
    do i=1,Lreal
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb + (ispin-1)*Norb
             do ik=1,Lk
                gktmp = Dbands(io,ik)/( zeta_real(io,io,i)-Hloc(io)-Ebands(io,ik) )
                Greal(ispin,ispin,iorb,iorb,i) = Greal(ispin,ispin,iorb,iorb,i) + gktmp
             enddo
          enddo
       enddo
       call eta(i,Lreal)
    end do
    call stop_timer
    call dmft_gloc_print_realaxis(wr,Greal,"Gloc",iprint)
  end subroutine dmft_get_gloc_realaxis_normal_dos
















  !PRINT GLOC SINGLE SITE
  subroutine dmft_gloc_print_matsubara(w,Gmats,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Gmats
    real(8),dimension(size(Gmats,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    Nspin = size(Gmats,1)
    Norb  = size(Gmats,3)
    call assert_shape(Gmats,[Nspin,Nspin,Norb,Norb,size(Gmats,5)],"dmft_gloc_print_matsubara",reg(fname)//"_mats")
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Gloc matsubara: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,"(A)") "Gloc matsubara: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_iw"//reg(gloc_suffix)
             call splot(reg(suffix),w,Gmats(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,"(A)") "Gloc matsubara: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_iw"//reg(gloc_suffix)
                call splot(reg(suffix),w,Gmats(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default                  !print all off-diagonals
       write(*,"(A)") "Gloc matsubara: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_iw"//reg(gloc_suffix)
                   call splot(reg(suffix),w,Gmats(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gloc_print_matsubara

  subroutine dmft_gloc_print_realaxis(w,Greal,fname,iprint)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Greal
    real(8),dimension(size(Greal,5))           :: w
    character(len=*),intent(in)                :: fname
    integer,intent(in)                         :: iprint
    !
    Nspin = size(Greal,1)
    Norb  = size(Greal,3)
    call assert_shape(Greal,[Nspin,Nspin,Norb,Norb,size(Greal,5)],"dmft_gloc_print_realaxis",reg(fname)//"_real")
    !
    select case(iprint)
    case (0)
       write(*,"(A)") "Gloc real: not written to file."
       !
    case(1)                  !print only diagonal elements
       write(*,"(A)") "Gloc real: write spin-orbital diagonal elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix=reg(fname)//&
                  "_l"//reg(txtfy(iorb))//&
                  "_s"//reg(txtfy(ispin))//&
                  "_realw"//reg(gloc_suffix)
             call splot(reg(suffix),w,Greal(ispin,ispin,iorb,iorb,:))
          enddo
       enddo
       !
    case(2)                  !print spin-diagonal, all orbitals 
       write(*,"(A)") "Gloc real: write spin diagonal and all orbitals elements."
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                suffix=reg(fname)//&
                     "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//&
                     "_realw"//reg(gloc_suffix)
                call splot(reg(suffix),w,Greal(ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
       !
    case default                  !print all off-diagonals
       write(*,"(A)") "Gloc real: write all elements."
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix=reg(fname)//&
                        "_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
                        "_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//&
                        "_realw"//reg(gloc_suffix)
                   call splot(reg(suffix),w,Greal(ispin,jspin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine dmft_gloc_print_realaxis

















  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !--------------------------------------------------------------------!
  !PURPOSE:
  ! Bcast/Reduce a vector of Blocks [Nlat][Nso][Nso] onto a matrix [Nlat*Nso][Nlat*Nso]
  !--------------------------------------------------------------------!
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix

  function matrix_to_blocks(Matrix,Nlat,Nso) result(Vblocks)
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: i,j,ip
    Vblocks=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks






  !--------------------------------------------------------------------!
  !PURPOSE: select a single block of the diagonal from a large matrix.
  !--------------------------------------------------------------------!
  function select_block_Nlso(ip,Matrix,Nlat,Nso) result(Vblock)
    integer                                 :: ip
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nso,Nso)           :: Vblock
    integer                                 :: i,j
    Vblock=zero
    i = 1+(ip-1)*Nso
    j =       ip*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function select_block_nlso
  !
  function select_block_nnn(ip,Matrix,Nlat,Nspin,Norb) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    integer                                          :: Nlat,Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block_nnn








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,jlat,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso


END MODULE DMFT_GET_GLOC_DOS
