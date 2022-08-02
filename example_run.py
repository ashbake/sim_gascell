# plot gas cell spectra
#https://docs.google.com/document/d/1Yy2HYR1b8KcXCd_9mAmnEvWcd6d2bwD3TzSvuuUzzoY/edit#
import numpy as np
import matplotlib.pylab as plt

from gascell import setup_hapi, hitran_ids, run_gen_spec, run_combine_spec


## EDIT ME ###
# define paths
hit_path = './hitran/' # make sure this path exists, location to store hitran files
out_path = './HapiSimulations/' # make sure this path exists, location to save hitran files

# Cell params
species = np.array(['CH4', 'C2H2'])  # must name species as listed in hapi or wont work
partial_pressures  = np.array([0.15, 0.05])  # partial pressures of each gas, sums to total gas cell pressure
p,t,l = sum(partial_pressures), 300, 10  #pressure total, temperature, cell length

#instrument params
R = 100000 #resolving power
v0, vf = 2270, 5130 # cm-1, 1.95-4.3 micron , wavelength range to gen spectra over

##############


if __name__=='__main__':
	l0 = 1 # just always set l0=1 for generating spectra that will be used to scale to other cell lengths

	# make and save new spectra for each molecule for params defined above. maybe add check to see if files exist already
	table   = setup_hapi(species,v0,vf,hit_path=hit_path,rerun=False) # this seems to be globally defined
	for specie in species: 
		#_, iso_nums, _ = hitran_ids(specie) #if want list of all iso_nums
		iso_nums = [1] # if just want the first, most common one, or can take subset. must be in list format though
		_, _ = run_gen_spec(np.array([specie]),v0, vf,p,t,l0,out_path,ploton=True,iso_nums=iso_nums) # next fxn reloads outputs from the save file so dont need outputs here
	
	# Combine into one cell, then degrades resolution,saves plot and text file with lam and tot_spec
	lam,tot_spec = run_combine_spec(partial_pressures,t,l,R,species,iso_nums,plot=True,save=True,out_path='./HapiSimulations/')
