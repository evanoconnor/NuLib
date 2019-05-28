!-*-f90-*-
program bremsstrahlung_point_example

  use nulib
  use inputparser
  
#if NUCLEI_HEMPEL
  use nuclei_hempel, only : set_up_Hempel
#endif
  implicit none

  !many people use different number of species, this is to denote how they are devided up.
  ! mytable_neutrino_scheme = 1 (three output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino+antineutrino
  ! neutrino_scheme = 2 (four output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon+tau neutrino             #4 mu and tau antineutrino
  ! neutrino_scheme = 3 (six output species)
  ! species #1: electron neutrino             #2 electron antineutrino
  !         #3: muon neutrino                 #4 mu antineutrino
  !         #5: tau neutrino                  #6 tau antineutrino
  integer :: mypoint_neutrino_scheme = 1

  !number of species for nulib to calculate interactions for, must
  !be six currently, average is done via above parameter
  integer :: mypoint_number_species = 6

  !number of species for nulib to output interactions for, must
  !be commensurate with neutrino_scheme above
  integer :: mypoint_number_output_species = 3

  !number of energy groups
  integer :: mypoint_number_groups = 18 !18

  character*200 :: parameters = "./parameters"

  !local variables
  real*8, allocatable,dimension(:,:) :: local_emissivity
  real*8, allocatable,dimension(:,:) :: local_absopacity
  real*8, allocatable,dimension(:,:) :: local_scatopacity
  real*8, allocatable,dimension(:,:) :: local_delta
  real*8, allocatable,dimension(:,:) :: local_Phi0, local_Phi1
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_epannihil, local_Phi1_epannihil
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_bremsstrahlung
  real*8, allocatable,dimension(:,:,:) :: local_Phi0_bremsstrahlung2
  real*8, allocatable,dimension(:,:) :: blackbody_spectra
  real*8, allocatable,dimension(:) :: eos_variables
  real*8,allocatable, dimension(:) :: Phi0_brem,Phi0_brem2
  real*8,allocatable, dimension(:) :: Phi0_brem_ann,Phi0_brem2_ann
  real*8 :: matter_prs,matter_ent,matter_cs2,matter_dedt,matter_dpderho,matter_dpdrhoe
  real*8 :: single_neutrino_emissivity_from_NNBrem_given_energypoint
  real*8 :: single_neutrino_emissivity_from_epannihil_given_energypoint
  integer :: keytemp,keyerr
  real*8 :: precision = 1.0d-10
  real*8 :: xrho, xtemp, xye
  real*8 :: n_N,eta_star
  integer :: i,j,inde
  real*8 :: dxfac,mindx
  real*8 :: Q,Q_ann,Q2,Q2_ann,M,dE_dEdt,dE_dEdt2,dE_dEdt3,k
  real*8 :: Kuro_inv_lamb,Hann_inv_lamb,epan_inv_lamb
  real*8 :: dn_dEdt,dn_dEdt2,dn_dEdt3,ratio,ratio_2,dn_dEdt_burrows
  real*8 :: dE_dEdt_2,dE_dEdt2_2,dE_dEdt3_2,dE_dEdt_burrows
  real*8 :: Bremsstrahlung_Phi0_Hannestad
  real*8 :: Bremsstrahlung_Phi0_Kuroda
  real*8 :: epannihil_Phi_Bruenn
!~   real*8 :: find_s
  real*8,dimension(100) :: integral1,integral2,integral3,temp_array,dens_array,ener_array
  real*8:: integral
  real*8 :: J_1,J_1_bar,phi_a,phi_p,energy,energy_2
  real*8,dimension(100) :: M1_mom,M1_mom_inv
  
  
  !fermi function
  real*8 :: fermidirac,fermidirac_dimensionless
  real*8 :: x
  real*8 :: find_s,find_s2
  real*8 :: find_g
  real*8 :: nue_absorption_on_n
  real*8 :: anue_absorption_on_p
  

  real*8,parameter :: nulib_emissivity_gf = 5.59424238d-55/(6.77140812d-06**3*2.03001708d+05)
  real*8,parameter :: nulib_opacity_gf = 1.0d0/6.77140812d-06
  real*8,parameter :: nulib_energy_gf = 1.60217733d-6*5.59424238d-55
  real*8,parameter :: nulib_kernel_gf = 6.77140812d-06**3/2.03001708d+05
  
  !allocate the arrays for the point values
  allocate(local_emissivity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_absopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_scatopacity(mypoint_number_output_species,mypoint_number_groups))
  allocate(local_delta(mypoint_number_output_species,mypoint_number_groups))
  allocate(blackbody_spectra(mypoint_number_output_species,mypoint_number_groups))

  call input_parser(parameters)

  !this sets up many cooefficients and creates the energy grid (one
  !zone + log spacing) see nulib.F90:initialize_nulib
  call initialize_nulib(mypoint_neutrino_scheme,mypoint_number_species,mypoint_number_groups)
#if NUCLEI_HEMPEL
  call set_up_Hempel ! set's up EOS for nuclear abundances
#endif
  !read in EOS table & set reference mass
  call read_eos_table(eos_filename)
  m_ref = m_amu !for SFHo_EOS (Hempel)
  ! m_ref = m_n !for LS220

  !example point
  
  xrho =  4.0d14!366150704625021.44d0! 4.435d13 !g/cm^3
  xtemp =2.19d11*kelvin_to_mev!14.511804010728726d0! !MeV
  xye = 0.1d0!0.28092030837767307d0 !dimensionless
!~   xrho =  366150704625021.44d0! 4.435d13 !g/cm^3
!~   xtemp =14.511804010728726d0! !MeV
!~   xye = 0.28092030837767307d0 !dimensionless
	
  !set up energies bins
  do_integrated_BB_and_emissivity = .true.
  mindx = 2.0d0
  bin_bottom(1) = 0.0d0 !MeV
  bin_bottom(2) = 2.0d0 !MeV
  bin_bottom(3) = bin_bottom(2)+mindx
  bin_bottom(number_groups) = 250.0d0 ! MeV
  
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
  bin_widths(number_groups) = 2.0*(energies(number_groups)-bin_bottom(number_groups))
  bin_top(number_groups) = bin_bottom(number_groups)+bin_widths(number_groups)
  
  
  allocate(eos_variables(total_eos_variables))
  eos_variables(:) = 0.0d0
  eos_variables(rhoindex) = xrho
  eos_variables(tempindex) = xtemp
  eos_variables(yeindex) = xye
  
  !! EOS stuff
  call set_eos_variables(eos_variables)
  write(*,*) "Rho: ",eos_variables(rhoindex)," g/ccm"
  write(*,*) "T: ",eos_variables(tempindex)," MeV"
  write(*,*) "Ye: ",eos_variables(yeindex)
  write(*,*) "X_n: ",eos_variables(xnindex)
  write(*,*) "X_p: ",eos_variables(xpindex)
  write(*,*) "X_alpha: ",eos_variables(xaindex)
  

  
  if (add_nue_kernel_bremsstrahlung.or.add_anue_kernel_bremsstrahlung.or. &
       add_numu_kernel_bremsstrahlung.or.add_anumu_kernel_bremsstrahlung.or. &
       add_nutau_kernel_bremsstrahlung.or.add_anutau_kernel_bremsstrahlung) then
	 
     write(*,*)
     write(*,*)
     write(*,*) "brem call"
     
     
     allocate(local_Phi0_bremsstrahlung(mypoint_number_output_species,mypoint_number_groups,2))
     allocate(local_Phi0_bremsstrahlung2(mypoint_number_output_species,mypoint_number_groups,2))
     allocate(local_Phi0_epannihil(mypoint_number_output_species,mypoint_number_groups,2))
     allocate(local_Phi1_epannihil(mypoint_number_output_species,mypoint_number_groups,2))
     allocate(Phi0_brem(mypoint_number_groups))
     allocate(Phi0_brem2(mypoint_number_groups))
     allocate(Phi0_brem_ann(mypoint_number_groups))
     allocate(Phi0_brem2_ann(mypoint_number_groups))
     
     
     write(*,*) "brem emi", single_neutrino_emissivity_from_NNBrem_given_energypoint( &
     3,energies(mypoint_number_groups),eos_variables)
     
	write(*,*) "echem",eos_variables(mueindex),energies(11)
	write(*,*) "echem fermi",fermidirac_dimensionless( &
				energies(11),-eos_variables(mueindex))&
				,exp(-energies(11)-eos_variables(mueindex))
     

     open (unit = 2, file = "brem.txt")   
        
		! energy spectrum
     ratio_2 = 0.0d0   

     do i=1,mypoint_number_groups !neutrino energ

		
!~ 		dE_dEdt = 0.0d0
!~ 		dE_dEdt2 = 0.0d0
!~ 		dE_dEdt3 = 0.0d0
!~ 		dn_dEdt = 0.0d0
!~ 		dn_dEdt2 = 0.0d0
!~ 		dn_dEdt3 = 0.0d0
		Phi0_brem = 0.0d0
		Phi0_brem2 = 0.0d0
		Phi0_brem_ann = 0.0d0
		Phi0_brem2_ann = 0.0d0
		  call single_epannihil_kernel_point_return_all(i, &
				eos_variables(mueindex)/eos_variables(tempindex),&
!~ 				0.0d0/eos_variables(tempindex),&
	          eos_variables(tempindex),local_Phi0_epannihil &
	          ,local_Phi1_epannihil,mypoint_neutrino_scheme)  		

	    
	    
		do inde =1,1
			 local_Phi0_bremsstrahlung = 0.0d0
			 local_Phi0_bremsstrahlung2 = 0.0d0
			 n_N =0.0d0
			 
			 
			  ! neutron-neutron  n_n
		      if (inde .EQ. 1) then   
		      
				  n_N = eos_variables(rhoindex)&
						*eos_variables(xnindex)&
						/(m_n*mev_to_gram)
			  
			  
			  ! proton-proton n_n
			  else if (inde .EQ. 2 ) then 
			  
  				  n_N = eos_variables(rhoindex)&
						*eos_variables(xpindex)&
						/(m_p*mev_to_gram)
			  
			  ! neutron-proton  n_n
			  else if (inde .EQ. 3 ) then
  				  n_N = eos_variables(rhoindex)	&
						*sqrt(eos_variables(xnindex)*eos_variables(xpindex))&
  				       /(sqrt(m_p*m_n)*mev_to_gram)
			  else 
				write(*,*) "Wrong number of species for brem test"
			  endif


!~ 			write(*,*) "test",n_N

			  ! Hannestad version
	          call single_bremsstrahlung_kernel_point_return_all_Hannestad(i, &
	          n_N, &
	          eos_variables(tempindex) &
	          ,local_Phi0_bremsstrahlung,mypoint_neutrino_scheme) 
	          
	          ! Kuroda version
	          call single_bremsstrahlung_kernel_point_return_all_Kuroda(i, &
	          n_N, &
	          eos_variables(tempindex)&
	 	      ,local_Phi0_bremsstrahlung2,mypoint_neutrino_scheme) 


			  ! Hannestad production kernels
	          if (inde .LT. 3) then 
		          Phi0_brem = Phi0_brem + local_Phi0_bremsstrahlung(3,:,1)
		          
		          
		      else if ( inde .EQ. 3)  then      ! factor for neutron-proton process
				  Phi0_brem = Phi0_brem + &
							  28.0d0/3.0d0*local_Phi0_bremsstrahlung(3,:,1) 
		      endif
		      
		      ! Hannestad annihilation kernels
	          if (inde .LT. 3) then 
		          Phi0_brem_ann = Phi0_brem_ann + local_Phi0_bremsstrahlung(3,:,2)
		          
		        
		      else if ( inde .EQ. 3)  then 		! factor for neutron-proton process
				  Phi0_brem_ann = Phi0_brem_ann + &
							  28.0d0/3.0d0*local_Phi0_bremsstrahlung(3,:,2) 
		      endif
		      
		      
		      
		      ! Kuroda production kernels
	          if (inde .LT. 3) then 
		          Phi0_brem2 = Phi0_brem2 + local_Phi0_bremsstrahlung2(3,:,1)
		           
		      else if ( inde .EQ. 3)  then 		! factor for neutron-proton process
				  Phi0_brem2 = Phi0_brem2 + &
							  28.0d0/3.0d0*local_Phi0_bremsstrahlung2(3,:,1) 
		      endif
		      
		      ! Kuroda annihilation kernels
	          if (inde .LT. 3) then 
		          Phi0_brem2_ann = Phi0_brem2_ann + local_Phi0_bremsstrahlung2(3,:,2) 
		          
		      else if ( inde .EQ. 3)  then 		! factor for neutron-proton process
				  Phi0_brem2_ann = Phi0_brem2_ann + &
							  28.0d0/3.0d0*local_Phi0_bremsstrahlung2(3,:,2) 
		      endif
	     enddo 

		    
		    
			Q = 0.0d0
			Q_ann = 0.0d0
			Q2 = 0.0d0
			Q2_ann = 0.0d0
			M = 0.0d0
			ratio = 0.0d0
			ratio_2 = 0.0d0
			
			
			! Loop over antineutrino energies
	        do j=1,mypoint_number_groups 
	        
				! Kotake pair production rates
				
				! 		Hannestad 
				Q = Q + energies(j)**2	 /(2.0d0* pi * hbarc_mevcm)**3&
					 *2.0d0*pi * bin_widths(j) &
					 * Phi0_brem(j) 
!~ 					 * local_Phi0_bremsstrahlung(3,j,1)
						    
				Q_ann = Q_ann + energies(j)**2	 /(2.0d0* pi * hbarc_mevcm)**3&
					 *2.0d0*pi * bin_widths(j) &
					  * Phi0_brem_ann(j) 
						    
				!		Kuroda
				Q2 = Q2 + energies(j)**2  /(2.0d0* pi * hbarc_mevcm)**3 &
					 *2.0d0*pi* Phi0_brem2(j)* bin_widths(j) 
					 
	
				Q2_ann = Q2_ann + energies(j)**2  /(2.0d0* pi * hbarc_mevcm)**3 &
					 *2.0d0*pi* Phi0_brem2_ann(j)* bin_widths(j)
	
				! 		electron-positron 
				M = M + energies(j)**2  /(2.0d0* pi * hbarc_mevcm)**3&
					*2.0d0*pi* local_Phi0_epannihil(3,j,1) * bin_widths(j) 
					
	
	        enddo
			
			
			!!!number production spectra 
			
			!		Hannestad
	        dn_dEdt =1.0d0/(2.0d0 * pi *hbarc_mevcm)**3 * energies(i)**2 &
						* (4.0d0 * pi) * Q 
			!    	Kuroda
	        dn_dEdt2 =1.0d0/(2.0d0 * pi * hbarc_mevcm)**3 * energies(i)**2 &
						* 4.0d0 * pi * Q2 
			!		positron-electron pair
	        dn_dEdt3 =1.0d0/(2.0d0 * pi * hbarc_mevcm)**3 * energies(i)**2 &
						 * (4.0d0 * pi) * M
			!		Burrows
	        dn_dEdt_burrows =single_neutrino_emissivity_from_NNBrem_given_energypoint(3,&
							     energies(i),eos_variables)/mev_to_erg*(4.0d0*pi)/energies(i)
						 
						 
			!!!energy production spectra
			
			!		Hannestad
	        dE_dEdt =1.0d0/(2.0d0 * pi * hbarc_mevcm)**3 * energies(i)**3 &
						* (4.0d0 * pi) * (Q )
			
			!		Kuroda
	        dE_dEdt2 =1.0d0/(2.0d0 * pi * hbarc_mevcm)**3 * energies(i)**3 &
						 * 4.0d0 * pi * Q2
			
			!		positron-electron pair
	        dE_dEdt3 =1.0d0/(2.0d0 * pi *hbarc_mevcm)**3 * energies(i)**3 &
						 * (4.0d0 * pi) * M

			!		Energy production spectra Burrows
			dE_dEdt_burrows = single_neutrino_emissivity_from_NNBrem_given_energypoint(1,&
							     energies(i),eos_variables)/mev_to_erg*(4.0d0*pi)
							     
			!!! inverse mean free path for comparaison with Bruenn 2018
			
							     
			ratio = (dE_dEdt-dE_dEdt2) /dE_dEdt *100.0d0  ! Hannestad/Kuroda energy production spectra ratio 
!~ 			ratio = Q/Q2    ! Hannestad/Kuroda energy production spectra ratio 
			
			
			
			! Difference in emissivities from Burrows formula
			ratio_2 =  ratio_2 + dE_dEdt2/&
							single_neutrino_emissivity_from_NNBrem_given_energypoint(3,&
							     energies(i),eos_variables)*mev_to_erg/(4.0d0*pi)
							     
							 

							
        write(2,*) energies(i) &
					,dn_dEdt*1e-37,dn_dEdt2*1e-37,dn_dEdt3*1e-36 &
					,dn_dEdt_burrows*1e-37 &
					,dE_dEdt*1e-37 &!/n_N & !
					,dE_dEdt2*1e-37 &!/n_N & !
					,dE_dEdt3*1e-36 &!/n_N & !
					,dE_dEdt_burrows*1e-37 &!/n_N & !
					,dE_dEdt-dE_dEdt2,dE_dEdt-dE_dEdt_burrows &
					,ratio,dE_dEdt2-dE_dEdt_burrows
     enddo 
     

     close(2)


	 

     write(*,*) "find_s"
     !! reproduce Fig. 4.2  in Raffelt 1985 

	 open(unit=5,file="find_s.txt")
	 do  j = 10,200,10
	 
	 	eta_star = (hbarc_mevcm/clight  * (3.0d0 * pi **2 * n_N)**(1.0d0/3.0d0))**2 &
				/(2.0d0 *m_amu/clight**2 * (real(j)/10.0d0)) !adimensional
				
		 do inde = 1,mypoint_number_groups
			Q = 0.0d0
			Q2 = 0.0d0
				do i = 1,mypoint_number_groups
					Q = Q + find_s(eta_star,0.0d0,&
								(energies(i)+energies(inde))/(real(j)/10.0d0))&
								* bin_widths(i)
					Q2 = Q2 + find_s2(eta_star,0.0d0,&
								(energies(i)+energies(inde))/(real(j)/10.0d0))&
								* bin_widths(i)
				enddo
			write(5,*) (energies(inde)+energies(i))/(real(j)/10.0d0),&
				eta_star,find_s(eta_star,0.0d0,energies(inde)/(real(j)/10.0d0)),&
				ABS(find_s(eta_star,0.0d0,energies(inde)/(real(j)/10.0d0))&
				-find_s2(eta_star,0.0d0,energies(inde)/(real(j)/10.0d0)))&
				/find_s(eta_star,0.0d0,energies(inde)/(real(j)/10.0d0)),&
				real(j)/10.0d0,&
				(energies(inde)+energies(i))/(real(j)/10.0d0)*5.0d0/8.0d0&
				*exp(-(energies(inde)+energies(i))/(real(j)/10.0d0)) &
				*find_s(eta_star,0.0d0,(energies(inde)+energies(i))/(real(j)/10.0d0)),&
				(energies(inde)+energies(i))/(real(j)/10.0d0)*5.0d0/8.0d0&
				*exp(-(energies(inde)+energies(i))/(real(j)/10.0d0)) &
				*find_s2(eta_star,0.0d0,(energies(inde)+energies(i))/(real(j)/10.0d0))
		 enddo
	 enddo
	 close(5)
	 
	 
	 
	 
	 
     deallocate(Phi0_brem)
     deallocate(local_Phi0_bremsstrahlung)
     deallocate(local_Phi0_epannihil)
     deallocate(local_Phi1_epannihil)
  endif


end program bremsstrahlung_point_example
  
