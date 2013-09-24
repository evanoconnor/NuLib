!-*-f90-*- 
program test

  use weak_rates
  use nulib
  implicit none

  character*200 :: filename = "/projects/ceclub/gr1dnulib/GitHub/NuLib/src/extra_code_and_tables/rates-ext.out"
  real*8 :: query_t9,query_lrYe,logECs,eos_variables(14)!,emissivity(number_groups)
  real*8, allocatable, dimension(:) :: emissivity
  real*8 dxfac,mindx
  integer reqnuc,rate,A,Z,i,ilrhoye,itemp


  call initialize_nulib(1,6,24)
  call readrates(filename,table_bounds)
  write(*,*) "ListPlot[{"
  do ilrhoye=10,150
     do itemp=10,10
        write(*,*) "{",dble(ilrhoye)/10.0d0,",",weakrates(47,24,dble(itemp)/1.0d0,dble(ilrhoye)/10.0d0,2),"},"
     end do
  end do
  write(*,*) "}]"

end program test
