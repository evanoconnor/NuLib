!-*-f90-*- 
program test

  use weak_rates
  use nulib
  implicit none

  character*200 :: filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/parameters"
  real*8 :: query_t9,query_lrYe,logECs,eos_variables(14)!,emissivity(number_groups)
  real*8, allocatable, dimension(:) :: emissivity
  real*8 dxfac,mindx
  integer reqnuc,rate,A,Z,i,ilrhoye,itemp,j,k
  real*8, dimension(25) :: t9array


  t9array(1)=0.01
  t9array(2)=0.06
  t9array(3)=0.10
  t9array(4)=0.15
  t9array(5)=0.20
  t9array(6)=0.30
  t9array(7)=0.40
  t9array(8)=0.55
  t9array(9)=0.70
  t9array(10)=0.85
  t9array(11)=1.00
  t9array(12)=1.25
  t9array(13)=1.50
  t9array(14)=1.75
  t9array(15)=2.00
  t9array(16)=2.50
  t9array(17)=3.00
  t9array(18)=4.00
  t9array(19)=5.00
  t9array(20)=7.50
  t9array(21)=10.00
  t9array(22)=20.00
  t9array(23)=30.00
  t9array(24)=65.00
  t9array(25)=100.00

  call input_parser(filename)
  call initialize_nulib(1,6,24)
  call readrates(table_bounds)

  open(11,file='analysis/QvA3d.dat')
  
  do a=1,300
     do z=1,200
        if(nndc_mass_table(a,z).eq.0.0d0.or.nndc_mass_table(a,z-1).eq.0.0d0) cycle
        write(11,*) A-Z,"    ",Z,"    ",return_qec(a,z,z-1)
     end do
  end do
  
  close(11)

end program test
