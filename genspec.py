import numpy as np
import matplotlib.pylab as plt

import hapi
from hapi import db_begin, absorptionCoefficient_Lorentz, select,transmittanceSpectrum, tableList, fetch_by_ids, getHelp

plt.ion()

def hitran_ids(molecule):
	"""
	given molecule name, return hitran molecule id and isotopologues list
	"""
	ids = hapi.ISO_ID

	mol_num, iso_nums, global_nums = [], [], []
	for key in ids.keys():
		if ids[key][-1] == molecule:
			mol_num = ids[key][0]
			iso_nums.append(ids[key][1])
			global_nums.append(key)

	return mol_num, iso_nums, global_nums

def gen_transmission(molecule,v0,vf,pres=1, temp=300, path_length=10):
	"""
	use hapi to  download hitran and generate transmission spectrum give pres, temp, pathlength from v0 to vf (cm-1)
	
	Inputs
	-------
	molecule - string of molecule to plot
	v0, vf - start and final wavenumbers (cm-1) to generate spectrum over
	pres - pressure in atm
	temp- temperature in kelvin
	path_length - length of tube in centimeters

	Outputs
	--------
	nu - wavenumbers of spectrum
	spec - transmission spectrum, 0 to 1

	Notes:
	-----
	** currently only plotting one isotopologue at a time, not sure how to do multiple at once
	** Lorentz is basically same as Voigt here and quicker
	** i'm cutting off lines weaker than 1e-28 even tho it doesnt seem to save time..

	# useful hapi functions
	db_begin()
	tableList()
	describeTable('H2O')
	select('tmp',ParameterNames=("gamma_self",'gamma_air'))   
	hapi.describeTable('tmp') 
	"""
	mol_id, iso_nums, _ = hitran_ids(molecule)

	# select subtable
	Cond = ('AND',('BETWEEN','nu',v0,vf),('>=','sw',1e-28)) 
	select(molecule,Conditions=Cond,DestinationTableName='tmp')

	# with table loaded, Loop through molecules again and save parameters 
	nu,coef  = absorptionCoefficient_Lorentz(Components=((mol_id,1),),
											SourceTables='tmp',HITRAN_units=False,
											Environment={'p':pres,'T':temp})
	nu,spec  = transmittanceSpectrum(nu,coef,Environment={'l':path_length,\
															'p':pres,'T':temp}) 

	return nu, spec

def setup_hapi(species):
	"""
	download HITRAN catalogs needed for each species

	returns table
	"""
	db_begin(hit_path)
	table = tableList()

	for molecule in species:
		if molecule in table:
			pass
		else:
			# Download so can work offline ..no idea where this stores shit and i htink it fails so just d/l it yourself
			_, iso_nums, global_nums = hitran_ids(molecule)
			print('Downloading Hitran data for ' + molecule)
			fetch_by_ids(molecule,global_nums,v0,vf,ParameterGroups=['160-char']) # first two isotopologues
			db_begin(hit_path)
			table = tableList()

	return table

def compare_sims():	
	"""
	Compare SpectralPlot simulatinos to HAPI ones

	ADB - HAPI seems to only take gamma-air which is the air broadening coefficient
	whereas spectralplot assumes one number for self broadening and air broadening
	Hapi's air broadening is more accurate (?), but doesnt do self broadening so IDK
	"""
	# load spectralplot sims
	f = np.loadtxt('SpectraPlotSimulations/CO2,x=.33,T=300K,P=1atm,L=10cm,simNum0.csv',delimiter=',')
	#gen corresponding hapi simulation
	pres, temp, path_length = 1, 300, 3.3
	nu, trans = gen_transmission('CO2',1e7/2010,1e7/2000,pres=pres, temp=temp, path_length=path_length)

	plt.figure()
	plt.plot(nu,trans,label='hapi')
	plt.plot(f[:,0],np.exp(-1*f[:,1]),label='SpectralPlot')
	plt.title('Pres %s Temp %s Length %s' %(pres, temp, path_length))
	plt.legend()

	f = np.loadtxt('SpectraPlotSimulations/CO2,x=1,T=300K,P=1atm,L=10cm,simNum3.csv',delimiter=',')
	#gen corresponding hapi simulation
	pres, temp, path_length = 1, 300, 10
	nu, trans = gen_transmission('CO2',1e7/2010,1e7/2000,pres=pres, temp=temp, path_length=path_length)

	plt.figure()
	plt.plot(nu,trans,label='hapi')
	plt.plot(f[:,0],np.exp(-1*f[:,1]),label='SpectralPlot')
	plt.title('Pres %s Temp %s Length %s' %(pres, temp, path_length))
	plt.legend()


	f = np.loadtxt('SpectraPlotSimulations/CO2,x=1,T=300K,P=.33atm,L=10cm,simNum4.csv',delimiter=',')
	#gen corresponding hapi simulation
	pres, temp, path_length = 0.33, 300, 10
	nu, trans = gen_transmission('CO2',1e7/2010,1e7/2000,pres=pres, temp=temp, path_length=path_length)

	plt.figure()
	plt.plot(nu,trans,label='hapi')
	plt.plot(f[:,0],np.exp(-1*f[:,1]),label='SpectralPlot')
	plt.title('Pres %s Temp %s Length %s' %(pres, temp, path_length))
	plt.legend()



if __name__=='__main__':
	hit_path = './hitran/' # make sure this path exists, location to store hitran files
	species = np.array(['CH4', 'N2O', 'CO2'])
	v0, vf  = 2325.0, 5129.0 # cm-1, 1.95-4.3 micron
	
	# setup hapi, load tables once bc slow
	table   = setup_hapi(species) # this seems to be globally defined

	# define environemnt
	pres, temp, path_length = 0.33, 300, 10

	# generate spec for each molecule "separately"
	plt.figure(1)
	specdic = {}
	for mol in species:
		specdic[mol] = gen_transmission(mol,v0,vf,pres=pres, temp=temp, path_length=path_length)
		#nu_,trans_,i1,i2,slit = hapi.convolveSpectrum(specdic[mol][0],specdic[mol][1],SlitFunction=hapi.SLIT_GAUSSIAN,Resolution=0.04,AF_wing=10.0)
		plt.plot(*specdic[mol],label=mol)

	plt.legend()
	plt.title('Pres %s Temp %s Length %s' %(pres, temp, path_length))

	# compare simulations hapi to spectralplot
	compare_sims()

	# if want multiple in one cell at total pressure P, evaluate each 
	# individual spectrum at P, but wth path_length equal to the fraction of the gas desired
	# e.g. 1/3 of each gas at P=2, use P=2, with path_length = 2 * total_length/3 for each gas


