!-*-f90-*-
! Psuedocode --  not everything declared

  real*8, allocatable, dimension(:,:,:,:,1:7) :: lambda  ! lambda[T,rhoYe,nucleus,rates+mu]
  real*8, allocatable, dimension(:,1,1,1) :: tao ! tao[nucleus, Q, A, Z]
  real*8, allocatable, dimension(:,:,1) :: eta ! eta(nucleus,energy group, emissivity)
  CHARACTER (LEN=*), PARAMETER :: file='rates-ext.out'


  call readfile_LMP(counter,lambda,tao,file)  ! counter indexes nuclei for which LMP data exists
  


  do i { for every, rho, T, Ye}

     do j=1,counter  ! counting over nuclei in LMP
        AvgE=lambda[i,j,6]/(lambda[i,j,5] + lambda[i,j,4])    ! Interpolation will be called. lpos?
        call nse_abundances(j,X)
        call AvgE_fit(AvgE,Q)
        do k=1,energy_groups
           eta(j,k,1)=(1/4*Pi)*AvgE*lambda(i,j,4)*rho(i)*X/M*nuSpectra(k)
        end do
     end do

  end do
