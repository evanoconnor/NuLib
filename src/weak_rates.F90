!-*-f90-*- 
module weak_rates
  implicit none
  real*8, allocatable,dimension(:,:,:,:) :: rates  ! rates[nuc,T,rhoYe,rates+uf(indexed by nrate)]
  real*8, allocatable,dimension(:,:) :: nucspec ! nucspec[nucleus, (Q, A, Z)]
  real*8, allocatable,dimension(:) :: t9dat
  real*8, allocatable,dimension(:) :: rhoYedat
  real*8, allocatable, dimension(:,:,:,:,:,:) :: C ! Matrix of spline coefficients (see desc. below)
  real*8, dimension(112,112) :: nucleus_index 
  integer nuc,nrho,nt9,nnuc,nrate


  contains

    subroutine readrates_LMP(filename)
 !      real*8, allocatable,dimension(:,:,:) :: eta ! eta(nucleus,energy group, emissivity
      character lindex
      character*200 :: filename,line
      integer i,dim,A,Z
      real*8 :: t9,lrho,uf,lbetap,leps,lnu,lbetam,lpos,lanu
      real*8 :: lrho_prior
      nuc = 0
      nrho = 0
      nt9 = 0
      dim = 1
      

      ! Count the dimension of the data in rhoYe and T9
      open(1,file=filename,status='old')
      do
         read(1,'(A)',end=10) line
         read(line,*) lindex
         if (lindex.eq.'n') then
            nuc = nuc + 1
            nrho = 0 
            nt9 = 0
         end if
         if (index('0123456789',lindex).ne.0) then
            lrho_prior = lrho
            read(line,*) t9,lrho,uf,lbetap,leps,lnu,lbetam,lpos,lanu            
            if (lrho.ne.lrho_prior) then
               nrho = nrho + 1
            end if
            nt9 = nt9 + 1
         end if
      end do
10    close(1)
      write(*,*) nuc,nrho,nt9/nrho

      ! allocate the array's based on dimension
      allocate(rates(nuc,nt9/nrho,nrho,7))
      allocate(nucspec(nuc,3))
      allocate(t9dat(nt9/nrho))
      allocate(rhoYedat(nrho))

      nuc = 0
      lrho = 0.0d0
      nrho = 0
      nt9 = 0
      ! fill the arrays
      open(1,file=filename,status='old')
      do
         read(1,'(A)',end=20) line
         read(line,*) lindex
         if (lindex.eq.'n') then
            nuc = nuc + 1
            read(line(index(line(1:),"Q= ")+3:index(line(1:),"Q= ")+10),*) nucspec(nuc,1)
            read(line(index(line(1:),"a=")+2:index(line(1:),"a=")+4),*) nucspec(nuc,2)
            A = nucspec(nuc,2)
            read(line(index(line(1:),"z=")+2:index(line(1:),"z=")+4),*) nucspec(nuc,3)
            Z = nucspec(nuc,3)
            nucleus_index(A,Z) = nuc
            nrho = 0 
            nt9 = 0
         end if
         if (index('0123456789',lindex).ne.0) then
            lrho_prior = lrho
            read(line,*) t9,lrho,uf,lbetap,leps,lnu,lbetam,lpos,lanu            
            if (lrho.ne.lrho_prior) then
               nrho = nrho + 1
               nt9 = 0
            end if
            nt9 = nt9 + 1
            rates(nuc,nt9,nrho,1)=lbetap
            rates(nuc,nt9,nrho,2)=leps
            rates(nuc,nt9,nrho,3)=lnu
            rates(nuc,nt9,nrho,4)=lbetam
            rates(nuc,nt9,nrho,5)=lpos
            rates(nuc,nt9,nrho,6)=lanu
            rates(nuc,nt9,nrho,7)=uf
            t9dat(nt9)=t9
            rhoYedat(nrho)=lrho
         end if
      end do
20    close(1)
      write(*,*) "Weak rate data loaded."

      ! build array of interpolating spline coefficients
      nnuc = nuc
      call monotonic_interp_2d(dim)
      write(*,*) "Interpolant functions built. Read-in is complete."
    end subroutine readrates_LMP


 ! #########################################################################
 ! #########################################################################
 ! #########################################################################

  ! This a 1-D Monotonic cubic spline interpolator. It takes 4 inputs and C 
  ! is global to the module:
       !       dim - Short for dimension, indicates which set of coeficients 
       !             to build. In the 1D case, dim = 1, and all interpolant
       !             coefficients are saved to C(1,:,:,:).
       !       spl - Short for spline function. This is an identifier for a
       !             particular interpolated function. In the 1D case, it is
       !             always equal to 1. In 2D, it will vary from 1 to N, where
       !             N is the number of data points.
       !         N - Number of data points to be interpolated over
       ! dataArray - Array, built by readfile, containing data in the
       !             form dataArray(n,1:2) where (x,y) :: (1,2)
       !         C - Array of spline coefficients. C = C(nucl.,rate,dim,
       !             splineseries,numdatapts,coef) where coef represents 
       !             coefficients a,b,c,d respectively used in the interpolation
       !   Reference:
       !       M. Steffen, "A simple method for monotonic interpolation in one dimension" 
       !       Astron. Astrophys. 239, 443-450 (1990)

       subroutine monotonic_interpolator(dim,spl,N,Data)
         INTEGER N      ! Number of data points (x,y) pairs
         real*8, dimension(N,2) :: Data
         real*8, dimension(N-1) :: S       ! Slope of secant passing through i+1 and i
         real*8, dimension(N-1) :: h       ! Size of interval, x_i+1 - x_i
         real*8, dimension(1:N) :: P     ! Slope of parabola passing through points i-1,i,i+1
         real*8, dimension(N) :: Dy         ! Value of 1st derivative at point i
         real*8 :: unity
         INTEGER i, j, dim, spl
         unity = 1.0d0


     ! Build the array of slopes
         do i=1,N-1
            S(i)=(Data(i+1,2)-Data(i,2))/(Data(i+1,1)-Data(i,1))
         end do
     ! Build the array of intervals
         do i=1,N-1
            h(i)=Data(i+1,1)-Data(i,1)
         end do
     ! Build the array of parabolic tangents
         P(1)=S(1)*(1+h(1)/(h(1)+h(2)))-S(2)*h(1)/(h(1)+h(2))
         do i=2,N-1
            P(i)=(S(i-1)*h(i)+S(i)*h(i-1))/(h(i-1)+h(i))
         end do
         P(N)=S(N-1)*(1+h(N-1)/(h(N-1)+h(N-2)))-S(N-2)*h(N-1)/(h(N-1)+h(N-2))
     ! Build the array of 1st derivatives
         do i=1,N
            if (i==1) then 
               if (P(1)*S(1) <= 0) then
                  Dy(1)=0.0d0
               else if (abs(P(1)) > 2.0d0*abs(S(1))) then
                  Dy(1)=2.0d0*S(1)
               else
                  Dy(1)=P(1)
               end if
            else if (i==N) then
               if (P(N)*S(N-1) <= 0.0d0) then
                  Dy(N)=0.0d0
               else if (abs(P(N)) > 2.0d0*abs(S(N-1))) then
                  Dy(N)=2.0d0*S(N-1)
               else
                  Dy(N)=P(N)
               end if
            else
               if (S(i-1)*S(i) <= 0.0d0) then
                  Dy(i)=0
               else if (abs(P(i)) > 2.0d0*abs(S(i-1))) then
                  Dy(i)=2.0d0*sign(unity,S(i-1))*min(abs(S(i-1)),abs(S(i)))
               else if (abs(P(i)) > 2.0d0*abs(S(i))) then
                  Dy(i)=2.0d0*sign(unity,S(i-1))*min(abs(S(i-1)),abs(S(i)))
               else
                  Dy(i)=P(i)
               end if
            end if
         end do

     ! Build the matrix of coefficients
          do i=1,N-1
            do j=1,4
               if (j==1) then
                  C(nuc,nrate,dim,spl,i,j)=(Dy(i)+Dy(i+1)-2.0d0*S(i))/(h(i)**2)
               else if (j==2) then
                  C(nuc,nrate,dim,spl,i,j)=(3.0d0*S(i)-2.0d0*Dy(i)-Dy(i+1))/h(i)
               else if (j==3) then
                  C(nuc,nrate,dim,spl,i,j)=Dy(i)
               else if (j==4) then
                  C(nuc,nrate,dim,spl,i,j)=Data(i,2)
               end if
            end do
         end do

         return

       end subroutine monotonic_interpolator

       subroutine interpolant(dim,spl,query)
         real*8 :: query,result
         INTEGER i,counter,dim,spl


         do i=1,nrho
            if (query <= rhoYedat(i)) then
               if (i==1) then
                  result = C(nuc,nrate,dim,spl,i,1)*(query-rhoYedat(i))**3 +&
                       C(nuc,nrate,dim,spl,i,2)*(query-rhoYedat(i))**2 +&
                       C(nuc,nrate,dim,spl,i,3)*(query-rhoYedat(i)) +&
                       C(nuc,nrate,dim,spl,i,4)
                  exit
               else
                  result = C(nuc,nrate,dim,spl,i-1,1)*(query-rhoYedat(i-1))**3 +&
                       C(nuc,nrate,dim,spl,i-1,2)*(query-rhoYedat(i-1))**2 +&
                       C(nuc,nrate,dim,spl,i-1,3)*(query-rhoYedat(i-1)) +&
                       C(nuc,nrate,dim,spl,i-1,4)
                  exit
               end if
            end if
         end do

         query = result
         return

       end subroutine interpolant

       subroutine monotonic_interp_2d(dim) 
         IMPLICIT NONE
         real*8, allocatable, dimension(:,:) :: Data
         INTEGER i,j,n,r,dim
         dim = 1

         allocate(Data(nt9,2))
         allocate(C(nuc,7,2,nrho,nt9,4))  ! 7 = 6 LMP weak rates + chemical potential

         nuc = 1
         do nuc=1,nnuc
            do i=1,nrho
               do nrate=1,7
                  do j=1,nt9
                     Data(j,1)=t9dat(j) ! (NEED TO ADD) prevent t9data being filled more than once
                     Data(j,2)=rates(nuc,j,i,nrate)
                  end do
                  call monotonic_interpolator(dim,i,nt9,Data)
               end do
            end do
         end do

         return
       end subroutine monotonic_interp_2d

       function weakrates(A,Z,query1,query2,rate_of_interest) result(interp_val) ! Efficiency consideration******

         real*8, allocatable, dimension(:,:) :: Data2d
         real*8 :: query1,query2,value,interp_val
         integer i,j,counter,dim,nucleus,rate_of_interest,A,Z
         dim = 1
         nrate = rate_of_interest
         allocate(Data2d(nrho,2))
         nucleus = nucleus_index(A,Z)

         do i=1,nrho
            do j=1,nt9 
                if (query1 <= t9dat(j)) then
                   if(j==1) then                ! Possible check to prevent extrapolation
                      value = C(nucleus,nrate,dim,i,j,1)*(query1-t9dat(j))**3 +&
                           C(nucleus,nrate,dim,i,j,2)*(query1-t9dat(j))**2 +&
                           C(nucleus,nrate,dim,i,j,3)*(query1-t9dat(j)) +&
                           C(nucleus,nrate,dim,i,j,4)
                      exit
                   else
                      value = C(nucleus,nrate,dim,i,j-1,1)*(query1-t9dat(j-1))**3 +&
                           C(nucleus,nrate,dim,i,j-1,2)*(query1-t9dat(j-1))**2 +&
                           C(nucleus,nrate,dim,i,j-1,3)*(query1-t9dat(j-1)) +&
                           C(nucleus,nrate,dim,i,j-1,4)
                      exit
                   end if
               end if
            end do
            Data2d(i,1) = rhoYedat(i)
            Data2d(i,2) = value
         end do

         interp_val = query2
         nuc = nucleus
         call monotonic_interpolator(2,1,nrho,Data2d)
         call interpolant(2,1,interp_val)

       end function weakrates

end module weak_rates

!use nulib, only : GLQ_n16_roots, GLQ_n16_weights
!call GaussLaguerreQuadrature_roots_and_weights(16,GLQ_n16_roots,GLQ_n16_weights)
!write(*,*) GLQ_n16_roots(1)
