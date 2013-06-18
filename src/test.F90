!-*-f90-*- 
program test

  use weak_rates
  implicit none
  include 'constants.inc'

  character*200 :: filename = "rates-ext.out"
  real*8 :: query_t9,query_lrYe,logECs,emissivity,eos_variables(14)
  integer reqnuc,rate,A,Z
  A = 65
  Z = 32
  query_t9 = 10.0d0
  query_lrYe = 1.0d0
  rate = 2
  eos_variables(1) = 2.0d0*10.0d0**10.0d0
  eos_variables(2) = 0.86
  eos_variables(3) = 0.5
  eos_variables(11) = 10.4

  call readrates_LMP(filename)
!  logECs =  weakrates(A,Z,query_t9,query_lrYe,rate)
!  write(*,*) logECs
!  call microphysical_electron_capture(emissivity)
  emissivity = emissivity_from_electron_capture_on_A(A,Z,eos_variables)
  write(*,*) "Emissivity: ", emissivity



  

end program test
