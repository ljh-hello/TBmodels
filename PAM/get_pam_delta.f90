program get_pam_teff
  implicit none
  real(8) :: tpd,delta,teff
  write(*,*)"tpd,teff:"
  read(*,*)tpd,teff
  teff=2.d0*teff
  delta = tpd**4/teff**2 - 4.d0*tpd**2
  delta=sqrt(delta)
  write(*,*)"delta=",delta
end program get_pam_teff
