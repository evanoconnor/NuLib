!-*-f90-*- 
program test
  use weak_rates
  implicit none
  character*200 :: filename = "rates-ext.out"
  real*8 :: query_t9,query_lrYe,logECs
  integer reqnuc,rate,A,Z
  A = 65
  Z = 32
  query_t9 = 10.0d0
  query_lrYe = 1.0d0
!  reqnuc = 1
  rate = 2
  call readrates_LMP(filename)
  logECs =  weakrates(A,Z,query_t9,query_lrYe,rate)
  write(*,*) logECs

!  call microphysical_electron_capture()
  



  

end program test
