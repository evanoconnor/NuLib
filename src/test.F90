!-*-f90-*- 
program test

  use weak_rates
  use nulib

  implicit none
  
  ! real*8 :: rate
  ! real*8 :: temperature
  ! real*8 :: analytic_weakrates
  ! real*8 :: q_gs
  ! real*8 :: mue

  ! mue = 34.575d0  
  ! temperature = 7.50d9*kelvin_to_mev
  ! q_gs = -1.6240d0
  ! rate = analytic_weakrates(temperature,q_gs,mue)

  ! write(*,*) log10(rate), 7.124

  !Weak rate data (currently LMP rates only)
  character*200 :: weakrates_filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/rates-ext.out"
  real*8 :: t9
  real*8 :: lrhoye
  real*8 :: analytic_weakrates,rate
  integer i

  t9=(2.08d0/kelvin_to_mev)*10**(-9.0d0)
  lrhoye=11.607455d0
  open(11,file="analysis/qvrates.dat")
  open(22,file="analysis/qvrates_analytic.dat")

  call readrates(weakrates_filename,table_bounds)

  
  do i=1,100
     write(11,*) nuclear_species(i,1), 10**(weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),t9,lrhoye,2)),int(nuclear_species(i,2)),int(nuclear_species(i,3))
  end do

  rate=0.0d0
  do i=-200,200
     rate = analytic_weakrates((t9*10**9)*kelvin_to_mev,dble(i)/10.0d0,37.8d0)
     write(22,*) dble(i)/10.0d0,rate
  end do

  close(11)
  close(22)

end program test
