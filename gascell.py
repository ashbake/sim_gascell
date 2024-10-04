import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate

import hapi,os
from hapi import db_begin, absorptionCoefficient_Lorentz, select,transmittanceSpectrum, tableList, fetch_by_ids, getHelp

plt.ion()

def hitran_ids(molecule):
	"""
	given molecule name, return hitran molecule id and isotopologues list
	
	inputs
	-------
	molecule - string of molecule to plot

	outputs
	--------	
	mol_num - hitran molecule number
	iso_nums - list of isotopologues
	global_nums - list of global ids
	"""
	ids = hapi.ISO_ID

	mol_num, iso_nums, global_nums = [], [], []
	for key in ids.keys():
		if ids[key][-1] == molecule:
			mol_num = ids[key][0]
			iso_nums.append(ids[key][1])
			global_nums.append(key)

	return mol_num, iso_nums, global_nums

def gen_transmission_one(molecule,v0,vf,pres=1, temp=300, path_length=10, iso_num=1,hit_path='./hitran/'):
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
	mol_id, _, _ = hitran_ids(molecule)
	if iso_num == 'all': _, iso_num, _ = hitran_ids(molecule)

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

def setup_hapi(species,v0,vf,rerun=True,hit_path='./hitran/'):
	"""
	download HITRAN catalogs needed for each species

	returns hapi specific table
	"""
	# make hitpath if doesnt exist
	if not os.path.exists(hit_path): os.makedirs(hit_path)

	db_begin(hit_path)
	table = tableList()

	for molecule in species:
		if molecule not in table or rerun:
			# Download so can work offline ..no idea where this stores shit and i htink it fails so just d/l it yourself
			_, iso_nums, global_nums = hitran_ids(molecule)
			print('Downloading Hitran data for ' + molecule)
			fetch_by_ids(molecule,global_nums,v0,vf,ParameterGroups=['160-char']) # first two isotopologues
			db_begin(hit_path)
			table = tableList()

	return table

def gen_transmission_all(species, v0, vf, pres, temp, path_length, iso_nums):#,res=.11):
	"""
		Run a bunch of species, each individual at one pressure and temperature and length.
		Parameters:
		-----------
		species : list
			List of species to run.
		v0 : float
			Initial frequency.
		vf : float
			Final frequency.
		pres : float
			Pressure.
		temp : float
			Temperature.
		path_length : float
			Path length.
		iso_nums : str or list
			Either 'all' or a list of isotope numbers.
		Returns:
		--------
		specdic
			A dictionary containing the transmission data for each species and isotope number
	"""
	specdic = {}
	for i,mol in enumerate(species):
		specdic[mol] = {}
		if iso_nums=='all':
			_, all_iso_nums, _ = hitran_ids(mol)
			specdic[mol]['all_iso'] = gen_transmission_one(mol,v0,vf,pres=pres[i], 
													temp=temp, path_length=path_length,#/nspecies,
													iso_num=all_iso_nums)	
		else:
			for iso_num in iso_nums:
				specdic[mol][iso_num] = gen_transmission_one(mol,v0,vf,pres=pres[i], 
													temp=temp, path_length=path_length,#/nspecies,
													iso_num=iso_num)
			# if want to degrade resolution through hapi tool:
			#nu_,trans_,i1,i2,slit = hapi.convolveSpectrum(specdic[mol][iso_num][0],specdic[mol][iso_num][1],SlitFunction=hapi.SLIT_GAUSSIAN,Resolution=res)

	return specdic


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

def run_gen_spec(species,v0, vf, pres,temp,path_length,out_path,ploton=True,iso_nums=[1],res=0.11):
	"""# setup hapi, load tables once bc slow
	# generate spec for each molecule "separately"

	inputs
	-------
	species - list of molecules to simulate
	v0, vf - start and final wavenumbers (cm-1) to generate spectrum over
	pres - pressure in atm
	temp- temperature in kelvin
	path_length - length of tube in centimeters
	out_path - path to save files
	ploton - plot the spectra
	iso_nums - list of isotopologues to simulate
	res - resolution of the spectra

	outputs
	--------
	savename - name of saved file
	specdic - dictionary of spectra
	"""
	if not os.path.exists(out_path): os.makedirs(out_path)
	
	specdic = gen_transmission_all(species, v0, vf, pres, temp, path_length,iso_nums)

	# save it
	moltag = ''
	for mol in species:
		moltag += '_%s' %mol 

	savename = out_path + 'transpec_p_%s_t_%s_l_%s%s_%s' %(pres, temp, path_length, moltag,iso_nums)
	savedat  = save_gascell(specdic, savename + '.txt')

	title = 'Pres: %s Temp: %s Length: %s' %(pres, temp, path_length)
	if ploton:
		plot_spectra(savedat, species, title, savename + '.png', fignum=-1)

	return savename, specdic

def load_sim(p,t,l,mol,iso_nums,lowres=True):
	f = np.loadtxt('HapiSimulations/' + 'transpec_p_%s_t_%s_l_%s_%s_%s.txt' %(p,t,l,mol,iso_nums))

	return f[:,0],f[:,1]

def define_lsf(v,res):
	# define gaussian in pixel elements to convolve resolved spectrum with to get rightish resolution
	dlam  = np.median(v)/res
	fwhm = dlam/np.abs(np.mean(np.diff(v))) # desired lambda spacing over current lambda spacing resolved to give sigma in array elements
	sigma = fwhm/2.634	
	x = np.arange(sigma*20) + 0.5
	if len(x)%2 == 1:
		x = np.arange(sigma*20  +sigma/2) + 0.5# if not symmetric then will shift everything
	gaussian = (1./sigma/np.sqrt(2*np.pi)) * np.exp(-0.5*( (x - 0.5*len(x))/sigma)**2 )

	if len(gaussian) < 20:
		raise(ValueError('Wavelength sampling too coarse for requestion resolving power - resample wavelength grid finer or lower resolution'))
	
	return gaussian

def degrade_spec(x,y,res):
	"""
	given wavelength, flux array, and resolving power R, return  spectrum at that R
	inputs
	-------
	x - wavelength array
	y - flux array
	res - resolving power

	outputs
	--------
	y_lowres - degraded spectrum
	"""
	lsf      = define_lsf(x,res=res)
	y_lowres = np.convolve(y,lsf,mode='same')

	return y_lowres

def run_combine_spec(partial_pressures,t,l,R,species,iso_nums,plot=True,save=True,out_path='./HapiSimulations/'):
	"""
	take the spectra made for given pres, temp and 
	"""
	if not os.path.exists(out_path): os.makedirs(out_path)

	p = sum(partial_pressures)
	ratios     = partial_pressures/p # ratio of pressures

	if plot: plt.figure()
	for i,mol in enumerate(species):
		xx,yy = load_sim(np.array([p]),t,1,mol,iso_nums) # set to high resolution spectrum, then convolve later
		lam,y = 1e7/xx[::-1], yy[::-1]**(l * ratios[i])
		if i==0: tot_spec = np.ones_like(y)
		if i==0: xtot_spec = lam*1.0

		if R < 300000:
			y_res      = degrade_spec(lam,y,R) 
		else:
			y_res = y

		if plot: plt.plot(lam,y_res,label=mol + ', %satm' %round(p * ratios[i],2),alpha=.7)
		# add to totl spectrum
		tck         = interpolate.splrep(lam,y_res, s=0)
		spec_interp = interpolate.splev(xtot_spec,tck,der=0,ext=1)
		spec_interp[np.where(spec_interp==0)] = 1

		tot_spec*=spec_interp

	if plot:
		plt.plot(xtot_spec,tot_spec,'k--',label='Combined',alpha=0.4)
		plt.xlabel('Wavelength (nm)')
		plt.ylabel('Transmittance')
	
		#plt.fill_between(1000*np.array([1.95,2.5]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
		#plt.fill_between(1000*np.array([2.85,4.2]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)

		plt.ylim(-0.01,1.1)
		title = 'Simulate Gas Cell, R=%s, T=%sK, $P_{tot}$=%satm'%(R,t,p)
		plt.title(title)

		plt.legend()
		if save: 
			savename = 'gascell_p_%s_t_%s_l_%s_R_%s_%s' %(p,t,l,R,species)
			plt.savefig(out_path + savename + '.pdf')
			#header = ??? add more info like partial pressures
			np.savetxt(out_path + savename + '.txt', np.vstack((xtot_spec,tot_spec)).T,
					header='lam transmittance') # maybe save each species spectrum too bc then can plot them separately

	return xtot_spec,tot_spec



if __name__=='__main__':
	# define paths
	hit_path = './hitran/' 
	out_path = './HapiSimulations/'

	# define gas cell params
	species = np.array(['CH4', 'NH3'])
	v0, vf = 1e7/2000, 1e7/980 # cm-1, 1.95-4.3 micron , wavelength range to gen spectra over
	pres, temp, path_length = 0.1, 273, 10 # atm, K, cm

	# make and save new spectra for each molecule for params defined above. maybe add check to see if files exist already
	table   = setup_hapi(species,v0,vf,hit_path=hit_path,rerun=False) # this seems to be globally defined
	specdic = gen_transmission_all(species, v0, vf, pres, temp, path_length, iso_nums) # next fxn reloads outputs from the save file so dont need outputs here
	
	# compare simulations hapi to spectralplot
	# compare_sims()




