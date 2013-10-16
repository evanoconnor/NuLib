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

  open(11,file='lmsh-extended.dat')

  do i=1,nnuc
     write(11,"(A14,A3,I2,A3,I2,A3,I2,A3,F8.4)")"pos. daughter","z=",int(nuclear_species(i,3))-1,"n=",int(nuclear_species(i,2))-(int(nuclear_species(i,3)-1)),"a=",int(nuclear_species(i,2)),"Q=",nuclear_species(i,1)
     do j=2,25
        do k=1,25
           query_t9 = t9array(k)
           query_lrYe = dble(j)/2
           write(11,"(F7.2,F6.1,F8.3,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)") query_t9,query_lrYe,&
                0.0000d0,&
                weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),query_t9,query_lrYe,1),&
                weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),query_t9,query_lrYe,2),&
                weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),query_t9,query_lrYe,3),&
                weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),query_t9,query_lrYe,4),&
                weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),query_t9,query_lrYe,5),&
                weakrates(int(nuclear_species(i,2)),int(nuclear_species(i,3)),query_t9,query_lrYe,6)
        end do
     end do
     write(11,"(A9)") "end"
  end do

  close(11)

end program test
