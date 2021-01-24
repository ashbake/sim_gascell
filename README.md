Usage:

0 - import hapi/(maybe also some hapi fxns..) if importing genspec.py instead of just running an edited version

1 - define paths - must create these folders if they dont already exist
	hit_path = './hitran/'              
	out_path = './HapiSimulations/'

2 - define gas cell params and wavenumber range
	species = np.array(['CH4', 'N2O', 'CO2'])
	v0, vf  = 2325.0, 5129.0           # cm-1, 1.95-4.3 micron
	pres, temp, path_length  = 0.3, 300, 10 # atm, K, cm

3 - run code - this saves it to out_path with option of plotting
	run(species,pres,temp,path_length,out_path,ploton=True)


* I also built in compare_sims which compares HAPI results to SpectralPlot online app and shows that the results are slightly different
* if you want multiple in one cell at total pressure P, evaluate each individual spectrum at P, but wth path_length equal to the fraction of the gas desired e.g. 1/3 of each gas at P=2, use P=2, with path_length = 2 * total_length/3 for each gas
* will probably want to save the hapi sims to file..



