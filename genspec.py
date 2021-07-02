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

def gen_transmission(molecule,v0,vf,pres=1, temp=300, path_length=10, iso_num=1):
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
	Cond = ('AND',('BETWEEN','nu',v0,vf),('>=','sw',1e-29)) 
	select(molecule,Conditions=Cond,DestinationTableName='tmp')

	# with table loaded, Loop through molecules again and save parameters 
	nu,coef  = absorptionCoefficient_Lorentz(Components=((mol_id,iso_num),),
											SourceTables='tmp',HITRAN_units=False,
											Environment={'p':pres,'T':temp})
	nu,spec  = transmittanceSpectrum(nu,coef,Environment={'l':path_length,\
															'p':pres,'T':temp}) 

	return nu, spec

def setup_hapi(species,rerun=True):
	"""
	download HITRAN catalogs needed for each species

	returns table
	"""
	db_begin(hit_path)
	table = tableList()

	for molecule in species:
		if molecule not in table or rerun==True:
			# Download so can work offline ..no idea where this stores shit and i htink it fails so just d/l it yourself
			_, iso_nums, global_nums = hitran_ids(molecule)
			print('Downloading Hitran data for ' + molecule)
			fetch_by_ids(molecule,global_nums,v0,vf,ParameterGroups=['160-char']) # first two isotopologues
			db_begin(hit_path)
			table = tableList()

	return table

def gen_gascell_spec(species, v0, vf, pres, temp, path_length, iso_nums,res=.11):
	"""
	save spectra to file
	"""
	specdic = {}
	specdic_lowres = {}
	nspecies=len(species)
	for mol in species:
		specdic[mol] = {}
		specdic_lowres[mol] = {}
		for iso_num in iso_nums:
			specdic[mol][iso_num] = gen_transmission(mol,v0,vf,pres=pres, 
													temp=temp, path_length=path_length/nspecies,
													iso_num=iso_num)
			nu_,trans_,i1,i2,slit = hapi.convolveSpectrum(specdic[mol][iso_num][0],specdic[mol][iso_num][1],SlitFunction=hapi.SLIT_GAUSSIAN,Resolution=res)
			specdic_lowres[mol][iso_num] = (nu_,trans_)

	return specdic, specdic_lowres

def save_gascell(specdic, savename):
	"""
	save specdics to file
	"""
	header = 'nu '
	# calculate total spectrum
	for i,mol in enumerate(specdic.keys()):
		for j, iso_num in enumerate(specdic[mol].keys()):
			if j==0: # if first isonum, start molspec total (spec total for that molecule)
				nu = specdic[mol][iso_num][0]
				molspec_total = specdic[mol][iso_num][1]
				if i==0: # if first molecule in cell, start spec total (all mol and isonums)
					spec_total = specdic[mol][iso_num][1]
					f  = np.zeros((len(specdic.keys()), len(nu))) # store everything here
				else:
					spec_total*= specdic[mol][iso_num][1]
			else:
				spec_total*= specdic[mol][iso_num][1]
				molspec_total *= specdic[mol][iso_num][1]
		f[i] = molspec_total
		header += mol  + ' '

	# write out to text file
	savedat = np.vstack((nu,f))
	np.savetxt(savename,savedat.T,header=header)

	return savedat

def plot_spectra(savedat, species, title, savename , fignum=-1):
	"""
	"""
	plt.figure(fignum,figsize=(12,6))
	specdic = {}
	for i,mol in enumerate(species):
		plt.plot(1e7/savedat[0],savedat[1+i],label=mol)

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')

	#plt.fill_between(1000*np.array([1.95,2.5]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	#plt.fill_between(1000*np.array([2.85,4.2]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)

	plt.title(title)
	plt.savefig(savename)

def run(species,v0, vf, pres,temp,path_length,out_path,ploton=True,iso_nums=[1],res=0.11):
	# setup hapi, load tables once bc slow
	# generate spec for each molecule "separately"
	specdic, specdic_lowres = gen_gascell_spec(species, v0, vf, pres, temp, path_length,iso_nums,res=res)

	# save it
	moltag = ''
	for mol in species:
		moltag += '_%s' %mol 

	savename = out_path + 'transpec_p_%s_t_%s_l_%s%s_%s' %(pres, temp, path_length, moltag,iso_nums)
	savedat  = save_gascell(specdic, savename + '.txt')
	savedat_lowres = save_gascell(specdic_lowres, savename + '_lowres.txt')

	title = 'Pres: %s Temp: %s Length: %s' %(pres, temp, path_length)
	#plot_spectra(savedat, title, savename + '.png', fignum=-1)
	if ploton:
		plot_spectra(savedat_lowres, species, title, savename + '_lowres.png', fignum=-1)

	return savename, specdic, specdic_lowres

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
	plt.savefig('SpectraPlotSimulations/comp_spec_p%s_t%s_l%s.png'%(pres,temp,path_length))

	f = np.loadtxt('SpectraPlotSimulations/CO2,x=1,T=300K,P=1atm,L=10cm,simNum3.csv',delimiter=',')
	#gen corresponding hapi simulation
	pres, temp, path_length = 1, 300, 10
	nu, trans = gen_transmission('CO2',1e7/2010,1e7/2000,pres=pres, temp=temp, path_length=path_length)

	plt.figure()
	plt.plot(nu,trans,label='hapi')
	plt.plot(f[:,0],np.exp(-1*f[:,1]),label='SpectralPlot')
	plt.title('Pres %s Temp %s Length %s' %(pres, temp, path_length))
	plt.legend()
	plt.savefig('SpectraPlotSimulations/comp_spec_p%s_t%s_l%s.png'%(pres,temp,path_length))


	f = np.loadtxt('SpectraPlotSimulations/CO2,x=1,T=300K,P=.33atm,L=10cm,simNum4.csv',delimiter=',')
	#gen corresponding hapi simulation
	pres, temp, path_length = 0.33, 300, 10
	nu, trans = gen_transmission('CO2',1e7/2010,1e7/2000,pres=pres, temp=temp, path_length=path_length)

	plt.figure()
	plt.plot(nu,trans,label='hapi')
	plt.plot(f[:,0],np.exp(-1*f[:,1]),label='SpectralPlot')
	plt.title('Pres %s Temp %s Length %s' %(pres, temp, path_length))
	plt.legend()
	plt.savefig('SpectraPlotSimulations/comp_spec_p%s_t%s_l%s.png'%(pres,temp,path_length))


def gen_spec_gascell1():
	"""
	wrapper to run molecules for first fabricated gas cell
	"""
	hit_path = './hitran/' # make sure this path exists, location to store hitran files
	out_path = './HapiSimulations/'

	v0, vf = 2270, 5130 # cm-1, 1.95-4.3 micron
	species = np.array(['CH4', 'N2O', 'CO2','C2H2', 'H2O'])

	torr_to_atm = 0.00131579 #atm/torr
	press    = torr_to_atm * np.array([23.3, 4.9, 15.2, 102.2, 2.2]) #Torr
	airpress = torr_to_atm *22.3 # torr, added air pressure

	press_tot = round(sum(press) + airpress ,3)
	temp = 300
	path_length = 1

	table   = setup_hapi(species) # this seems to be globally defined
	for specie in species:
		_, iso_nums, _ = hitran_ids(specie)
		savename, dic, diclow = run(np.array([specie]),v0, vf,press_tot,temp,path_length,out_path,ploton=True,iso_nums=[1])


def gen_spec_gascell_newrequest():
	"""
	wrapper to run molecules for first fabricated gas cell
	"""
	hit_path = './hitran/' # make sure this path exists, location to store hitran files
	out_path = './HapiSimulations/'

	v0, vf = 2270, 5130 # cm-1, 1.95-4.3 micron
	species = np.array(['CH4', 'N2O', 'CO2','C2H2'])

	press    = np.array([0.06, 0.07, 0.16, 0.06]) #atm

	press_tot = round(sum(press),3)
	temp = 300
	path_length = 1

	table   = setup_hapi(species) # this seems to be globally defined
	for specie in species:
		_, iso_nums, _ = hitran_ids(specie)
		savename, dic, diclow = run(np.array([specie]),v0, vf,press_tot,temp,path_length,out_path,ploton=True,iso_nums=[1])


if __name__=='__main__':
	# define paths
	hit_path = './hitran/' # make sure this path exists, location to store hitran files
	out_path = './HapiSimulations/'

	# define gas cell params
	species = np.array(['HCN']) #np.array(['CH4', 'N2O', 'CO2','HCN', 'H2O','C2H2'])
	v0, vf  = 6300, 6800 # cm-1, 1.95-4.3 micron
	pres, temp, path_length = 0.026, 300, 5.5 # atm, K, cm

	table   = setup_hapi(species,rerun=True) # this seems to be globally defined
	_, iso_nums, _ = hitran_ids(species[0])
	savename, dic, diclow = run(species,v0, vf, pres,temp,path_length,out_path,ploton=True,iso_nums=[1])

	# compare simulations hapi to spectralplot
	# compare_sims()




