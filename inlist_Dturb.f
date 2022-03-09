! Insert this in the MESA inlist use to deffine the stellar physics and MESA controls

&controls !

!-----------------------TURBULENT DIFFUSION---------------------------------

x_integer_ctrl(1)=1 ! Coment or 0 no Dturb, 1 refference parameter, 2 Montréal prescreption
! Verma prescription -> rho_0 is diffened at a reference parameter [Mass (based in Vermas prescription); or  Pressure; or Temperature] of the zoone 

x_integer_ctrl(2)=1 ! 1 --ref. mass; 2 -- ref. Temperature; 3 -- reff. Pressure 
x_ctrl(1) = 5d-4    ! Value od the reference parameter used. Defined in x_integer_ctrl(2)= 1 --ref mass; 2 -- ref Temperature; 3 -- reff Pressure

x_ctrl(2) = 10000   ! Omega Dturb (Verma and Silva Aguirre 2019)
x_ctrl(3) = 4       ! N Dturb     (Verma and Silva Aguirre 2019)

/ !END control

