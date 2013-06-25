 !-*-f90-*- 
 module weak_rates
   implicit none
   real*8, allocatable,dimension(:,:,:,:) :: rates  ! rates[nuc,T,rhoYe,rates+uf(indexed by nrate)]
   real*8, allocatable,dimension(:,:) :: nucspec ! nucspec[nucleus, (Q, A, Z)]
   real*8, allocatable,dimension(:) :: t9dat
   real*8, allocatable,dimension(:) :: rhoYedat
   real*8, allocatable, dimension(:,:,:,:,:,:) :: C ! Matrix of spline coefficients (see desc. below)
   integer, allocatable, dimension(:,:) :: nucleus_index 
   integer nuc,nrho,nt9,nnuc,nrate,nspecies


   contains

     subroutine readrates(filename)
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

       write(*,*) "Mark before hempel initialize"
       call set_up_Hempel ! set's up EOS for nuclear abundances
       write(*,*) "Mark before number_of_species"
       call get_Hempel_number_of_species(nspecies) ! returns the total number of nuclei
       write(*,*) "Mark after hempel in readrates, nspecies is ", nspecies
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
       allocate(nucleus_index(nspecies,nspecies))

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
       write(*,*) "out of loop"
                    
 20    close(1)
       write(*,*) "closed"
       write(*,*) "Weak rate data loaded."

       ! build array of interpolating spline coefficients
       nnuc = nuc
       call monotonic_interp_2d(dim)
       write(*,*) "Interpolant functions built. Read-in is complete."
     end subroutine readrates


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

        subroutine microphysical_electron_capture(ns,eos_variables,emissivity)
          
          use nulib, only : total_eos_variables, number_groups

          integer i,ns
          real*8, intent(in) :: eos_variables(total_eos_variables)
          integer, dimension(nspecies) :: nuclei_A
          integer, dimension(nspecies) :: nuclei_Z
          real*8, dimension(number_groups) :: emissivity
          real*8, dimension(nspecies) :: number_densities
          real*8, dimension(nspecies) :: mass_fractions

          if (ns == 1) then
             ! Hempel EOS and number of species are set up in readrates
             write(*,*) "Mark 1"
             call get_Hempel_As_and_Zs(nuclei_A,nuclei_Z)
             write(*,*) "Mark 2"
             call nuclei_distribution_Hempel(nspecies,nuclei_A(i),nuclei_Z(i),mass_fractions,number_densities,eos_variables)
             write(*,*) "Mark 3"
             emissivity = 0.0d0
             do i=1,nspecies
                emissivity = emissivity + emissivity_from_electron_capture_on_A(nuclei_A(i),nuclei_Z(i),number_densities(i),&
                     eos_variables)
                write(*,*) i
             enddo
          endif
        end subroutine microphysical_electron_capture


        function  emissivity_from_electron_capture_on_A(A,Z,number_density,eos_variables) result(emissivity)

          use nulib, only : total_eos_variables,energies,number_groups,do_integrated_BB_and_emissivity&
               ,mueindex,rhoindex,tempindex,yeindex,GLQ_n16_roots,GLQ_n16_weights
          include 'constants.inc'

          integer A,Z

          real*8, intent(in) :: eos_variables(total_eos_variables)
          real*8, intent(in) :: number_density
          real*8 :: emissivity(number_groups) !final answer in erg/cm^/s/srad/MeV
          real*8 :: avgenergy(2)
          real*8 :: qec_eff                   !effective Qec for the approximate neutrino spectra

          !spectrum integration
          real*8 :: normalization_constant    !nu spectra normalization, units of 1/MeV^5/s
          real*8 :: spectra
          real*8 :: t9
          real*8 :: lrhoYe
          real*8 :: nu_spectrum_eval
          integer i,ng

          ! if there is no data for a nucleus, this should prevent any further calculations for that species
          if (nucleus_index(A,Z) == 0) then
             do ng=1,number_groups
                emissivity(ng) = 0.0d0
             end do
             return
          endif

          !set local eos_variables for rate interpolation
          lrhoYe = log10(eos_variables(rhoindex)*eos_variables(yeindex))
          t9 = (eos_variables(tempindex)/kelvin_to_mev)/(10.0d0**9.0d0)  ! Conversion from MeV to GK

          !avgenergy(1) calculated from rates, avgenergy(2) calculated from the eff neutrino spectrum
          qec_eff = nucspec(nucleus_index(A,Z),1) !using Qgs from LMP table as seed
          avgenergy = average_energy(A,Z,qec_eff,eos_variables)
          qec_eff = qec_solver(avgenergy,qec_eff,eos_variables)

          !calculate normalization constant using effective neutrino spectra
          do i=1,16
             spectra = spectra + GLQ_n16_weights(i)*ec_neutrino_spectra(GLQ_n16_roots(i),qec_eff,eos_variables(mueindex),&
                  eos_variables(tempindex))
          end do
          spectra = (eos_variables(tempindex)**5)*spectra
          normalization_constant = (10.0d0**weakrates(A,Z,t9,lrhoYe,1)+10.0d0**weakrates(A,Z,t9,lrhoYe,2))&
               /spectra
          
          !
          if (do_integrated_BB_and_emissivity) then

          else
             do ng=1,number_groups
                nu_spectrum_eval =&
                     (eos_variables(tempindex)**4.0d0)*normalization_constant*&
                     ec_neutrino_spectra(energies(ng)/eos_variables(tempindex),qec_eff,eos_variables(mueindex),&
                     eos_variables(tempindex))
                emissivity(ng) = energies(ng)*number_density*nu_spectrum_eval/(4.0d0*pi) 
             end do
          endif

        end function emissivity_from_electron_capture_on_A

        function ec_neutrino_spectra(nu_energy_per_T,q,uf,T) result(nu_spectra)
          !definition: n(E) = (T^4)*N*ec_neutrino_spectra, where N is the normalization constant
          ![N] = #/MeV^5/s
          include 'constants.inc'

          !local variables
          real*8 :: T
          real*8 :: nu_energy_per_T
          real*8 :: q
          real*8 :: uf
          real*8 :: nu_spectra
          real*8 :: q_T
          real*8 :: uf_T

          !redefining q and the chemical potential to be dimensionless, T must be in MeV
          q_T = q/T
          uf_T = uf/T  

          !ec neutrino spectra (dimensionless)
          nu_spectra = (nu_energy_per_T**2.0d0)*((nu_energy_per_T-q_T)**2.0d0) &
               /(1.0d0+exp(nu_energy_per_T-q_T-uf_T))
          !note a factor of T^4 was removed to make the spectrum dimensionless, the returned result
          !should therefore be multiplied by T^4 (T in MeV)

        end function ec_neutrino_spectra


        function ec_neutrino_spectra_q_derivative(nu_energy_per_T,q,uf,T) result(nu_spectra_derivative)
          !definition: d/dq n(E) = (ec_neutrino_spectra_q_derivative/T) * n(E)

          include 'constants.inc'

          !local variables
          real*8 :: T !must be passed in with units of MeV
          real*8 :: nu_energy_per_T
          real*8 :: q
          real*8 :: uf
          real*8 :: nu_spectra_derivative
          real*8 :: q_T
          real*8 :: uf_T

          !redefining q and the chemical potential to be dimensionless, T must be in MeV
          q_T = q/T
          uf_T = uf/T

          !ec neutrino spectra derivative with respect to q (q=qec_eff) divided by the neutrino spectra
          nu_spectra_derivative = (1.0d0-1.0d0/(1.0d0+exp(nu_energy_per_T-q_T-uf_T))-2.0d0/(nu_energy_per_T-q_T))
          !note a factor of /T was removed to make the spectrum dimensionless, the returned result
          !should therefore be divided by T in MeV

        end function ec_neutrino_spectra_q_derivative

        function average_energy(A,Z,qec_eff,eos_variables) result(avgenergy)

          use nulib

          integer A,Z
          real*8, intent(in) :: eos_variables(total_eos_variables)

          !local rate variables
          real*8 :: lbetap
          real*8 :: leps
          real*8 :: lnu
          real*8 :: lrhoYe
          real*8 :: t9
          real*8 :: mu_e ! electron chem.pot.

          !local integration variables
          integer :: i
          real*8 :: avgenergy(2)
          real*8 :: qec_eff
          real*8 :: spectra
          real*8 :: number_density_integral
          real*8 :: energy_density_integral
          
          ! if there is no data for a nucleus, this should prevent any further calculations for that species
          if (nucleus_index(A,Z) == 0) then
             do i=1,2
                avgenergy(i) = 0.0d0
             end do
             return
          endif

          !setting local variables from eos_variables
          lrhoYe = log10(eos_variables(rhoindex)*eos_variables(yeindex))
          t9 = (eos_variables(tempindex)/kelvin_to_mev)/(10.0d0**9.0d0)  ! Conversion from MeV to GK
          mu_e = eos_variables(mueindex)

          !interpolating rates for given eos_variables and calculating average neutrino energy
          lbetap = weakrates(A,Z,t9,lrhoYe,1)
          leps = weakrates(A,Z,t9,lrhoYe,2)
          lnu = weakrates(A,Z,t9,lrhoYe,3)         
          avgenergy(1) = 10.0d0**lnu/(10.0d0**leps + 10.0d0**lbetap)

          !set weights and roots for quadrature integration, then calculate <E>
          avgenergy(2)=0.0d0
          energy_density_integral=0.0d0
          number_density_integral=0.0d0
          do i=1,16
             spectra = ec_neutrino_spectra(GLQ_n16_roots(i),qec_eff,mu_e,eos_variables(tempindex))
             energy_density_integral=energy_density_integral+GLQ_n16_weights(i)*GLQ_n16_roots(i)*spectra
             number_density_integral=number_density_integral+GLQ_n16_weights(i)*spectra
          end do
          avgenergy(2) = (eos_variables(tempindex))*(energy_density_integral/number_density_integral)

        end function average_energy

        function qec_solver(avgenergy,q,eos_variables) result(qec_eff)
          use nulib, only : GLQ_n16_roots, GLQ_n16_weights, total_eos_variables, mueindex, tempindex

          !starting point energies using a seed q = Qgs
          real*8, intent(in) :: avgenergy(2)
          real*8, intent(in) :: eos_variables(total_eos_variables)

          !local buffer variables for calculations
          real*8 :: avge_rates
          real*8 :: avge_spectra
          real*8 :: qec_eff
          real*8 :: q
          integer i

          !integration variables
          real*8 :: energy_density_integral
          real*8 :: number_density_integral
          real*8 :: dq_energy_density_integral
          real*8 :: dq_number_density_integral
          real*8 :: nu_spectra_deriv_coef
          real*8 :: nu_spectra

          avge_rates = avgenergy(1)
          avge_spectra = avgenergy(2)

          !Newton-Raphson technique to find zero (in q) of f(q) = <E>_rates - <E(q)>_spectra
          do while (abs(avge_rates - avge_spectra)/avge_rates > 1.0d-8) 
             energy_density_integral=0.0d0 
             number_density_integral=0.0d0
             dq_energy_density_integral=0.0d0
             dq_number_density_integral=0.0d0
             nu_spectra=0.0d0

             !Quadrature integration of n(E), n'(E), E*n(E), E*n'(E) such that E = kT*xi (xi is dimensionless)
             !q_n+1 = q_n - (<E>_rates - <E(q)>_spectra)/(-d/dq <E(q)>_spectra) 
             do i=1,16
                nu_spectra = ec_neutrino_spectra(GLQ_n16_roots(i),q,eos_variables(mueindex),&
                     eos_variables(tempindex))
                nu_spectra_deriv_coef = ec_neutrino_spectra_q_derivative(GLQ_n16_roots(i),q,&
                     eos_variables(mueindex),eos_variables(tempindex))
                energy_density_integral=energy_density_integral+GLQ_n16_weights(i)*GLQ_n16_roots(i)*&
                     nu_spectra
                number_density_integral=number_density_integral+GLQ_n16_weights(i)*nu_spectra
                dq_energy_density_integral=dq_energy_density_integral+GLQ_n16_weights(i)*GLQ_n16_roots(i)*&
                     nu_spectra*nu_spectra_deriv_coef
                dq_number_density_integral=dq_number_density_integral+GLQ_n16_weights(i)*nu_spectra*&
                    nu_spectra_deriv_coef
            end do 
            avge_spectra = eos_variables(tempindex)*(energy_density_integral/number_density_integral)
            q = q + (avge_rates - avge_spectra)/&
                 ((dq_energy_density_integral*number_density_integral-energy_density_integral*&
                 dq_number_density_integral)/(number_density_integral*number_density_integral))
         end do
         qec_eff = q
       end function qec_solver

end module weak_rates
