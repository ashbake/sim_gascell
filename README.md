Usage:

0 - import hapi/(maybe also some hapi fxns..) if importing code

1 - make location to store hitran files

    	 hit_path = './hitran/'
	 
2 - define species array

    species = np.array(['CH4', 'N2O', 'CO2'])

3 - define wavenumber bounadaries and environment variables

    v0, vf  = 2325.0, 5129.0			  # cm-1, 1.95-4.3 micron
    
    pres, temp, path_length = 0.33, 300, 10	# atm, Kelvin, centimeters

4 - # setup hapi, load tables once bc slow

    table   = setup_hapi(species)		  # this loads/downloads HITRAN data for molecules in species arrays and saves them to hit_path

5 - you can now generate a transmittance spectrum

    nu, spec = gen_transmission(mol,v0,vf,pres=pres, temp=temp, path_length=path_length)

6 - you can use hapi to convolve it like this where resolution is delta_nu in cm-1

     nu_,trans_,i1,i2,slit = hapi.convolveSpectrum(nu, spec,SlitFunction=hapi.SLIT_GAUSSIAN,Resolution=0.04,AF_wing=10.0)

* I also built in compare_sims whihc shows that the results are slightly different
* if you want multiple in one cell at total pressure P, evaluate each individual spectrum at P, but wth path_length equal to the fraction of the gas desired e.g. 1/3 of each gas at P=2, use P=2, with path_length = 2 * total_length/3 for each gas




