module run_star_extras
	! Turbulent diffusion   
	integer ::method_Dturb=1;
	real(dp) :: W_turb=0;
	real(dp) :: N_turb=0;
	logical :: use_mass_ref=.false.
	logical :: use_Temp_ref=.false.
	logical :: use_Pres_ref=.false.
	real(dp) :: Del_M_turb=0;
	real(dp) :: T_turb=0;
	real(dp) :: P_turb=0;
	integer(dp) :: d_len=0;
	real(dp) :: m_sol = 1.9892d33  ! solar mass (g)  <<< gravitational mass, not baryonic
	real(dp) :: r_sol = 6.9598d10 ! solar radius (cm)
	real(dp) :: l_sol = 3.8418d33
	integer(dp) :: counter = 0
	real(dp) :: Dhe0=0;
	real(dp) :: rho0=0;
	real(dp) :: Dhe_in1=0;
	real(dp) :: Dhe_in2=0;
	real(dp) :: D_turb;


!####################################################################################################################################################################################################

	subroutine extras_controls(id, ierr)


		select case (s%x_integer_ctrl(1))
		case(1) 
		s%use_other_D_mix = .true. 
		s% other_D_mix => turb_D_mix
		method_Dturb=1
		write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!VERMA PRESCREPTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		case(2) 
		s%use_other_D_mix = .true. 
		s% other_D_mix => turb_D_mix
		method_Dturb=2
		write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Montréal PRESCREPTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		case default
		write(*,*) "EMPTY OPTION --> STANDARD MESA"
		s%use_other_D_mix = .false.
		end select 
			
	end subroutine extras_controls

!####################################################################################################################################################################################################

	subroutine turb_D_mix(id, ierr)
		integer, intent(in) :: id
		integer, intent(out) :: ierr
		integer :: nwu
		type (star_info), pointer :: s
		ierr = 0
		call star_ptr(id, s, ierr)
		if (ierr /= 0) return

		W_turb = s%x_ctrl(2) !Omega
		N_turb = s%x_ctrl(3) !N exponent
		counter = 1
		if (method_Dturb==1) then 
			select case (s%x_integer_ctrl(2))
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


			Dhe0=(3.3e-15*s% T(counter)**(2.5)) / (4*s% rho(counter) * log(1.0+1.125e-16* s% T(counter)**3/ s% rho(counter)))
			rho0=s% rho(counter)
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

		elseif (method_Dturb==2) then 
			!write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Montréal PRESCREPTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			do while (s%mixing_type(counter)==convective_mixing) 
				rho0=s% rho(counter)
				counter=counter+1
				if (counter>=s% nz) then 
					rho0=0
				end if
			end do
			d_len=s%nz

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

