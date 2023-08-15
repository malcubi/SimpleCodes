
advection: advection.f90
	gfortran -ffree-form advection.f90 -o advection

advimplicit: advimplicit.f90
	gfortran -ffree-form advimplicit.f90 -o advimplicit

wave: wave.f90
	gfortran -ffree-form wave.f90 -o wave

kleingordon: kleingordon.f90
	gfortran -ffree-form kleingordon.f90 -o kleingordon

burgers: burgers.f90
	gfortran -ffree-form burgers.f90 -o burgers

fluid: fluid.f90
	gfortran -ffree-form fluid.f90 -o fluid

riemann: riemann.f
	g77 riemann.f -o riemann

tunnel: tunnel.f90
	gfortran -ffree-form tunnel.f90 -o tunnel

tunnel2: tunnel2.f90
	gfortran -ffree-form tunnel2.f90 -o tunnel2

tortoise: tortoise.f90
	gfortran -ffree-form tortoise.f90 -o tortoise

newtongrav1: newtongrav1.f90
	gfortran -ffree-form newtongrav1.f90 -o newtongrav1

testtype: testtype.f90
	gfortran -ffree-form testtype.f90 -o testtype

testgrid: testgrid.f90
	gfortran -ffree-form testgrid.f90 -o testgrid

trumpet: trumpet.f90
	gfortran -ffree-form trumpet.f90 -o trumpet

toygammadriver: toygammadriver.f90
	gfortran -ffree-form toygammadriver.f90 -o toygammadriver

whitedwarf: whitedwarf.f90
	gfortran -ffree-form whitedwarf.f90 -o whitedwarf





