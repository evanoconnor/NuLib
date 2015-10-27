!-*-f90-*- 
program test
  use class_ratelibrary
  implicit none
  !NuLib parameters file (weak rates and EOS)
  character*200 :: parameters_filename = "../parameters"
  type(RateLibrary) :: weakrate_library

  
  call set_up_Hempel !set's up EOS for nuclear abundances
  weakrate_library = new_RateLibrary(parameters_filename)
  
  
end program test
