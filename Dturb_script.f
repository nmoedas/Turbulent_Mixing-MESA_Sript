module run_star_extras
	! Turbulent diffusion   
      integer ::method_Dturb=0;						! Select the turbulent mixing method/prescription within the routine
      																		! 0 - Richer et al. (2000)
      																		! 1 - Proffit & Michaud (1991)
      integer :: D_turb_pms=0;						! Debugging during PMS when the star is fully convective
      																		! To avoid divergence (when using Proffit & Michaud; 1991)
      real(dp) :: W_turb=0;								! Value of the constant in the turbulent mixing
      real(dp) :: N_turb=0;								! Value of the in the exponent constant in the turbulent mixing
      logical :: use_mass_ref=.false.     ! Use mass envelope for reference depth for Richer et al. (2000)
      logical :: use_Temp_ref=.false.			! Use Temperature for reference depth for Richer et al. (2000)
      logical :: use_Pres_ref=.false.			! Use Pressure for reference depth for Richer et al. (2000)
      real(dp) :: Del_M_turb=0;						! Value of the reference mass enveloup
      real(dp) :: T_turb=0;								! Value of the reference Temperature
      real(dp) :: P_turb=0;								! Value of the reference  Pressure
      real(dp) :: m_sol = 1.9892d33				! solar mass (g)
      real(dp) :: r_sol = 6.9598d10				! solar radius (cm)
      real(dp) :: l_sol = 3.8418d33				! Luminosaty of the Sun (erg s^-1)
      integer(dp) :: counter = 0					! Counter for cycles
      integer(dp) :: d_len=0;							! Maximum value for the counter 
      real(dp) :: Dhe0=0;									! The diffusion value of helium at the reference depth Richer et al. (2000)
      real(dp) :: rho0=0;									! The density value at the reference depth Richer et al. (2000)
      real(dp) :: D_turb;									! Value of the calculate Turbulent Mixing
      logical :: selec_dturb_type=.false. ! This options will select witch turbulent mixing prescription the rotine uses depending of the size of the convective envelope at the ZAMS
      																		! M_CZ>10d-5 Msun Proffit & Michaud (1991); M_CZ<=10d-5 Msun Richer et al. (2000)
      logical :: om_n_select=.false. ;    ! This options select the constants depending of the prescription only uses when selec_dturb_type=.true.
																					! Values predifined in the routine

!####################################################################################################################################################################################################

	subroutine extras_controls(id, ierr)
	 	! Extra controls to take into account what was indicated in the list of models to select the recipe.
	 	! s%x_integer_ctrl(1)= 0 - no turbulent mixing
		! 										 1 - Richer et al. (2000)
		! 										 2 - Proffit & Michaud (1991)
		!											 3 - Selects depending on the mass of the convective envelope at ZAMS
		select case (s%x_integer_ctrl(1))
		  case(1) 
			s%use_other_D_mix = .true. 
			s% other_D_mix => turb_D_mix
			method_Dturb=0
			write(*,*) "------------Richer Prescreption------------"
			case(2) 
			s%use_other_D_mix = .true. 
			s% other_D_mix => turb_D_mix
			method_Dturb=1
			write(*,*) "------------Proffit & Michaud Prescreption------------"
		  case(3) 
			s%use_other_D_mix = .true. 
			s% other_D_mix => turb_D_mix
			selec_dturb_type = .true.
			om_n_select = .true.
			write(*,*) "------------Select  PRESCREPTION------------"
		case default
			write(*,*) "No Turbulent Mixing"
			s%use_other_D_mix = .false.
		end select 
	end subroutine extras_controls

!####################################################################################################################################################################################################

	subroutine turb_D_mix(id, ierr)
		integer, intent(in) :: id
		integer, intent(out) :: ierr
		real(dp) :: Xi										!Initial abundance of hydrogen in the star
		type (star_info), pointer :: s
		ierr = 0
		call star_ptr(id, s, ierr)
		if (ierr /= 0) return
		
		!Values of the constants added to the inlist of the model
		W_turb = s%x_ctrl(2)
		N_turb = s%x_ctrl(3)
		Del_M_turb= s%x_ctrl(1)

		if (selec_dturb_type) then !Selec prescription depending on the size of convective envelop depending on mass of the convective zone
															 !Is is predefined to use Proffit & Michaud only select Richer if convective zone smaller then the refference mass
															 !Only prepared for reference Mass
			method_Dturb = 1
			
			Xi=1-s% initial_z-s% initial_y
			if (abs((Xi-s%center_h1)/Xi)>=0.01) then !critery close to the ZAMS
				counter = 1
				do while (s%mixing_type(counter)==convective_mixing .and. (counter/=s%nz-1)) 
					counter=counter+1
				enddo
				 if ((s% m(1)/m_sol) -(s% m(counter)/m_sol) < Del_M_turb) then 
					 method_Dturb = 0
					 write(*,*) "---------CHANGE DTURB to Richer PRESCREPTION--------"
				 endif
				 write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
				 selec_dturb_type=.false.
			endif
		endif 

	counter = 1  !reset counter
	if (method_Dturb==0) then  !Richer prescriptio


		if (om_n_select) then !Values based in Vermas & Silva Aguirre (2019) calibration
				W_turb = 10000 
				N_turb = 4
		endif

		select case (s%x_integer_ctrl(2)) !select which parameter is used for the reference depth
			case(1) 
				use_mass_ref=.true.
				Del_M_turb= s%x_ctrl(1)
			case(2) 
				use_Temp_ref=.true.
				T_turb= s%x_ctrl(1)
			case(3) 
				use_Pres_ref=.true.
				P_turb= s%x_ctrl(1)
		end select 
		
		! Search of the refference depth
		if (use_mass_ref) then 
			do while ((s% m(1)/m_sol) -(s% m(counter)/m_sol) <= Del_M_turb) 
			counter=counter+1
			end do
		end if
		if (use_Temp_ref) then
			do while (s% T(counter)<= T_turb .and. counter <= s%nz) 
			counter=counter+1
			end do
		end if

		if (use_Pres_ref) then
			do while ((s% P(counter)) <= P_turb) 
			counter=counter+1
			end do
		end if

		!Calcule of the diffusion coeficient of helium and density at refference depth
		Dhe0=(3.3e-15*s% T(counter)**(2.5)) / (4*s% rho(counter) * log(1.0+1.125e-16* s% T(counter)**3/ s% rho(counter)))
		rho0=s% rho(counter)
		!Calculate turbulent mixing and add its value to diffusion in the model
		d_len=s%nz
		do counter=1,d_len

			if (s% rho(counter)<=0) then
				D_turb=0
			else 
				D_turb = W_turb*Dhe0 *(rho0/s% rho(counter))**N_turb
			end if
				s% D_mix(counter)=s% D_mix(counter)+D_turb
		end do

		!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

		elseif (method_Dturb==1) then  !Proffit & Michaud prescriptio
			if (om_n_select) then !Values based on solar calibration to reproduce the Lithium abundance.
					W_turb =  1615.4863
					N_turb =  1.3
			endif
			if (D_turb_pms>0) then  !avoid the divergences in the models when the star is fully convective in pms
					! verify if the stars is fully  convective
			 	  if (D_turb_pms==1 .and. s%mixing_type(s%nz)==convective_mixing .and. s%mixing_type(s%nz-1)==convective_mixing .and. buff_r_T<s% Teff) then
						D_turb_pms=2
					! verify if radiative zone appered in the star 
				  elseif (D_turb_pms==2 .and. s%mixing_type(s%nz)==no_mixing .and. s%mixing_type(s%nz-1)==no_mixing .and.  1.5>s% log_surface_luminosity) then
						D_turb_pms=0
				       
				  	end if
				    buff_r_T=s% Teff
			 	    return
			 end if
			
			!density at the bottom of the convective zone
			do while (s%mixing_type(counter)==convective_mixing .and. (counter/=s%nz)) 
				rho0=s% rho(counter)
				counter=counter+1
				if (counter>=s% nz-1) then 
					rho0=0
				end if
			end do
			d_len=s%nz
			
			!Calculate turbulent mixing and add its value to diffusion in the model
			do counter=1,d_len
				if (s% rho(counter)<=0) then
					D_turb=0
				else 
					D_turb = W_turb*(rho0/s% rho(counter))**N_turb
				end if
				s% D_mix(counter)=s% D_mix(counter)+D_turb
			end do   
		endif
	end subroutine turb_D_mix


end module run_star_extras
