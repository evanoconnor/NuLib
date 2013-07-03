!-*-f90-*- 
program test

  use weak_rates
  use nulib
  implicit none

  character*200 :: filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/rates-ext.out"
  real*8 :: query_t9,query_lrYe,logECs,eos_variables(14)!,emissivity(number_groups)
  real*8, allocatable, dimension(:) :: emissivity
  real*8 dxfac,mindx
  integer reqnuc,rate,A,Z,i
  call initialize_nulib(1,6,24)
  allocate(emissivity(number_groups))

  !set up energies bins
  do_integrated_BB_and_emissivity = .false.
  add_nue_emission_weakinteraction_ecap = .true.
  mindx = 1.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 1.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 250.0d0
  
  call nulib_series2(number_groups-1,bin_bottom(2),bin_bottom(number_groups),mindx,dxfac)
  do i=4,number_groups
     bin_bottom(i) = bin_bottom(i-1)+(bin_bottom(i-1)-bin_bottom(i-2))*dxfac
  enddo
  
  !calculate bin widths & energies from the bottom of the bin & energy at top on bin
  do i=1,number_groups-1
     energies(i) = (bin_bottom(i)+bin_bottom(i+1))/2.0d0
     bin_widths(i) = bin_bottom(i+1)-bin_bottom(i)
     bin_top(i) = bin_bottom(i+1)
  enddo
  energies(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups-1)*dxfac/2.0d0
  bin_widths(number_groups) = 2.0d0*(energies(number_groups)-bin_bottom(number_groups))
  bin_top(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups)

  A = 65
  Z = 32
  query_t9 = 10.0d0
  query_lrYe = 1.0d0
  m_ref = m_amu !for SFHo_EOS (Hempel)
  eos_variables(1) = 2.0d0*10.0d0**10.0d0
  eos_variables(2) = 0.86
  eos_variables(3) = 0.5
  eos_variables(11) = 10.901772220438655

  call readrates(filename,table_bounds)
  call microphysical_electron_capture(1,eos_variables,emissivity)
  
  do i=1,number_groups
     write(*,*) energies(i),emissivity(i)
  end do

  write(*,*)  "Avg E from emissivity = ",Sum(emissivity(:)*bin_widths(:))/Sum(emissivity(:)*bin_widths(:)/energies(:))
  write(*,*)  "Summed rate from emissivity = ",Sum(emissivity(:)*bin_widths(:)/energies(:)*4.0d0*pi)
  write(*,*)  "Nu energy loss rate from emissivity = ",Sum(emissivity(:)*bin_widths(:)*4.0d0*pi)
  
  



end program test
