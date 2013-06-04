!-*-f90-*- 
program test
  use weak_rates
  implicit none
  character*200 :: filename = "rates-ext.out"
  DOUBLE PRECISION query_t9,query_lrYe,logECs
  integer reqnuc
  query_t9 = 10.0d0
  query_lrYe = 10.0d0
  reqnuc = 100
  call readrates_LMP(filename)
  call interpolant_2d(query_t9,query_lrYe,logECs)
  write(*,*) logECs

!  call microphysical_electron_capture()
  



  

end program test
