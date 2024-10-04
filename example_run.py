# plot gas cell spectra
#https://docs.google.com/document/d/1Yy2HYR1b8KcXCd_9mAmnEvWcd6d2bwD3TzSvuuUzzoY/edit#
import numpy as np
import matplotlib.pylab as plt

from gascell import setup_hapi, degrade_spec, gen_transmission_all


## EDIT ME ###
# define paths
hit_path = './hitran/' # make sure this path exists, location to store hitran files from table query
out_path = './HapiSimulations/' # make sure this path exists, location to save hitran files

# Cell params
# https://www.wavelengthreferences.com/product/gas-cells/
# each gas cell has specific  length, pressure options
species = np.array([ 'CO2', 'C2H2','HCN','HF','CH4', 'H2O'])          # must name species as listed in hapi or wont work
lengths   =  np.array([80,    3,    5.5, 2.7,  5.5,   16.5]) # cm
pressures = np.array([500,   400,  100,  50,   100,    20]) * 0.00131579 # torr to mbar
#partial_pressures  = np.array([0.2])  # partial pressures of each gas, sums to total gas cell pressure
#other = 'NH3', 
#lengths   =       [80,   5.5,   16.5, 2.7, 5.5,   3/80]
t,l0 = 283, 1  # pressure total, temperature, cell length (do 1)

#instrument params
R = 100000 #resolving power 
v0, vf = 1e7/2000, 1e7/980 # cm-1, 1.95-4.3 micron , wavelength range to gen spectra over

##############


if __name__=='__main__':
	l0 = 1 # just always set l0=1 for generating spectra that will be used to scale to other cell lengths

	# make and save new spectra for each molecule for params defined above. maybe add check to see if files exist already
	table   = setup_hapi(species,v0,vf,hit_path=hit_path,rerun=False) # this seems to be globally defined
	specdic = gen_transmission_all(species, v0, vf, pressures, t, l0, [1]) # next fxn reloads outputs from the save file so dont need outputs here
	
	# plot each species
	plt.figure()
	for key in specdic.keys():
		# degrade spectrum version
		x,y = 1e7/specdic[key][1][0],specdic[key][1][1]
		lowres_spec = degrade_spec(x,y, R)
		#np.where()
		p = plt.plot(x,y,alpha=0.5,ls='--')
		plt.plot(x,lowres_spec,c=p[0].get_color(),label=key + ' %satm %scm'%(pressures[species==key],lengths[species==key]))
	
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmission')
	plt.legend()
	plt.show()
	plt.axhline(0.95, color='k', ls='--')

	species = np.array([ 'CO2'])         
	pressures = np.array([0.1,0.1])
	p = np.sum(pressures) 

	for specie in species: 
		#_, iso_nums, _ = hitran_ids(specie) #if want list of all iso_nums
		iso_nums = [1] # if just want the first, most common one, or can take subset. must be in list format though#
		_, _ = run_gen_spec(np.array([specie]),v0, vf,np.array([p]),t,l0,out_path,ploton=True,iso_nums=iso_nums) # next fxn reloads outputs from the save file so dont need outputs here

	# Combine into one cell, then degrades resolution,saves plot and text file with lam and tot_spec
	lam,tot_spec = run_combine_spec(pressures,t,10,R,species,iso_nums,plot=True,save=False,out_path=out_path)
