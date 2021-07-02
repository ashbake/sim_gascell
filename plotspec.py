# plot gas cell spectra
#https://docs.google.com/document/d/1Yy2HYR1b8KcXCd_9mAmnEvWcd6d2bwD3TzSvuuUzzoY/edit#
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from hapi import convolveSpectrum, SLIT_GAUSSIAN

plt.ion()

def define_lsf(v,res):
	# define gaussian in pixel elements to convolve resolved spectrum with to get rightish resolution
	dlam  = np.median(v)/res
	fwhm = dlam/np.abs(np.mean(np.diff(v))) # desired lambda spacing over current lambda spacing resolved to give sigma in array elements
	sigma = fwhm/2.634	
	x = np.arange(sigma*20) + 0.5
	if len(x)%2 == 1:
		x = np.arange(sigma*20  +sigma/2) + 0.5# if not symmetric then will shift everything
	gaussian = (1./sigma/np.sqrt(2*np.pi)) * np.exp(-0.5*( (x - 0.5*len(x))/sigma)**2 )

	return gaussian

def degrade_spec(x,y,res):
	"""
	given wavelength, flux array, and resolving power R, return  spectrum at that R
	"""
	lsf      = define_lsf(x,res=res)
	y_lowres = np.convolve(y,lsf,mode='same')

	return y_lowres


def load_sim(p,t,l,mol,lowres=True):
	if lowres:
		f = np.loadtxt('HapiSimulations/' + 'transpec_p_%s_t_%s_l_%s_%s_lowres.txt' %(p,t,l,mol))
	else:
		f = np.loadtxt('HapiSimulations/' + 'transpec_p_%s_t_%s_l_%s_%s.txt' %(p,t,l,mol))

	return f[:,0],f[:,1]


def plot_spec_wl(species, p, t, l,ratios, lowres=True,ltot=10,res=30000):
	plt.figure(figsize=(12,6))
	res=res if lowres else 1e6

	for i,mol in enumerate(species):
		xx,yy = load_sim(p,t,l,mol,lowres=False) # set to high resolution spectrum, then convolve later
		lam,y = 1e7/xx[::-1], yy[::-1]**(ltot * ratios[i])
		if i==0: tot_spec = np.ones_like(y)
		if i==0: xtot_spec = lam*1.0

		if res < 1e6:
			#dnu = np.median(1e7/lam)/res #match genspec setting
			#test        = convolveSpectrum(1e7/lam[::-1],y[::-1],dnu,SlitFunction=SLIT_GAUSSIAN) 
			#tck         = interpolate.splrep(1e7/test[0][::-1], test[1][::-1], s=0)
			#y_res       = interpolate.splev(lam,tck,der=0,ext=0)
			y_res      = degrade_spec(lam,y,res) # this shifts data meh
		else:
			y_res = y

		plt.plot(lam,y_res,label=mol + ', %satm' %round(p * ratios[i],2),alpha=.7,zorder=100)
		# add to totl spectrum
		tck         = interpolate.splrep(lam,y_res, s=0)
		spec_interp = interpolate.splev(xtot_spec,tck,der=0,ext=1)
		spec_interp[np.where(spec_interp==0)] = 1

		tot_spec*=spec_interp

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	

	title = 'Pres: %satm Temp: %sK Length: %scm R:%s' %(p, t, ltot,res)
	plt.fill_between(1000*np.array([1.95,2.5]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	plt.fill_between(1000*np.array([2.85,4.2]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	plt.plot([1950,4300],[.95,0.95],'k--')

	plt.ylim(-0.01,1.1)
	plt.title(title)
	plt.savefig('HapiSimulations/plots/trn_spec_%s_p_%s_t_%s_l_%s_r_%s.pdf'%(species,p,t,ltot,res))

	return xtot_spec,tot_spec,title


def plot_spec_wl_old(species, p, t, l,ratios, lowres=True,ltot=10):
	plt.figure(-1,figsize=(12,6))
	res=30000 if lowres else 1e6

	for i,mol in enumerate(species):
		xx,yy = load_sim(p,t,l,mol,lowres=lowres)
		lam,y = 1e7/xx[::-1], yy[::-1]**(ltot * ratios[i])
		if i==0: tot_spec = np.ones_like(y)
		if i==0: xtot_spec = lam*1.0

		plt.plot(lam,y,label=mol + ', %satm' %round(p * ratios[i],2),alpha=.7,zorder=100)
		# add to totl spectrum
		tck         = interpolate.splrep(lam,y, s=0)
		spec_interp = interpolate.splev(xtot_spec,tck,der=0,ext=0)

		tot_spec*=spec_interp

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	
	title = 'Pres: %satm Temp: %sK Length: %scm R:%s' %(p, t, ltot,res)
	plt.fill_between(1000*np.array([1.95,2.5]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	plt.fill_between(1000*np.array([2.85,4.2]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	plt.plot([1950,4300],[.95,0.95],'k--')

	plt.ylim(-0.01,1.1)
	plt.title(title)
	plt.savefig('HapiSimulations/plots/trn_spec_%s_p_%s_t_%s_l_%s_r_%s.pdf'%(species,p,t,ltot,res))

	return xtot_spec,tot_spec,title


def overplot_online_spec(x,spec,p,t,l,lowres=False,ltot=10):
	"""
	overplot spectral plot spectrum for n2o, ch4, co2
	"""
	species = np.array(['CH4_[1]','N2O', 'CO2'])

	f, (ax1, ax2) = plt.subplots(2, 1, figsize=(12,6),sharex=True)	
	for i,mol in enumerate(species):
		x,y = load_sim(p,t,l,mol,lowres=lowres)
		if i==0: tot_spec = np.ones_like(y)
		if i==0: xtot_spec = 1e7/x[::-1]
		ax1.plot(1e7/x, y**(ltot/len(species)),label=mol + ', %satm' %round(p/len(species),2),alpha=.7,zorder=100)
		# add to totl spectrum
		tck         = interpolate.splrep(1e7/x[::-1], y[::-1]**(ltot/len(species)), s=0)
		spec_interp = interpolate.splev(xtot_spec,tck,der=0,ext=0)

		tot_spec*=spec_interp

	ax1.legend()
	ax2.set_xlabel('Wavelength (nm)')
	ax1.set_ylabel('Transmittance')
	ax2.set_ylabel('Transmittance')

	fco2 = np.loadtxt('SpectraPlotSimulations/CO2,x=.333,T=300K,P=.3atm,L=10cm,simNum0.csv',delimiter=',')
	fn2o = np.loadtxt('SpectraPlotSimulations/CH4,x=.333,T=300K,P=.3atm,L=10cm,simNum2.csv',delimiter=',')
	fch4 = np.loadtxt('SpectraPlotSimulations/N2O,x=.333,T=300K,P=.3atm,L=10cm,simNum1.csv',delimiter=',')

	ax2.plot(xtot_spec,tot_spec,'k-',label='hapi')
	ax2.plot(1e7/fco2[:,0],np.exp(-1*fch4[:,1])*np.exp(-1*fn2o[:,1])*np.exp(-1*fco2[:,1]),
			'r--',label='SpectralPlot')
	ax1.set_title('Pres %satm Temp %sK Length %scm' %(p, t, ltot))
	ax2.legend()
	ax2.set_xlim(np.min(1e7/fco2[:,0]),np.max(1e7/fco2[:,0]))
	ax2.set_ylim(0.1,1.1)

	plt.savefig('HapiSimulations/plots/comp_specplot_hapi_%s.pdf'%species)


def plot_totspec(x,spec,species,title,lowres):
	plt.figure(2,figsize=(12,6))
	plt.plot(x,spec)

	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	
	plt.fill_between(1000*np.array([1.95,2.5]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	plt.fill_between(1000*np.array([2.85,4.2]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)

	plt.ylim(-0.01,1.1)
	plt.title(title)
	plt.savefig('HapiSimulations/plots/trn_totspec_%s.pdf'%(species))

	# save x, spec to file
	if lowres:
		np.savetxt('HapiSimulations/final/trn_totspec_%s_lowres.txt'%species, np.vstack((x,spec)).T,
					header=title)


def plot_pres_change():
	"""
	"""
	x1,y1 = load_sim(0.1,300,10,'CH4_[1]',lowres=False)
	x3,y3 = load_sim(0.3,300,10,'CH4_[1]',lowres=False)
	x5,y5 = load_sim(0.5,300,10,'CH4_[1]',lowres=False)
	x7,y7 = load_sim(0.7,300,10,'CH4_[1]',lowres=False)


	plt.figure()
	plt.plot(x1,y1,'steelblue',ls='--',label='p=0.1atm')
	plt.plot(x3,y3,'r',ls='--',label='p=0.3atm')
	plt.plot(x5,y5,'g',ls='--',label='p=0.5atm')
	plt.plot(x7,y7,'orange',ls='--',label='p=0.7atm')


	x1,y1 = load_sim(0.1,300,10,'CH4_[1]',lowres=True)
	x3,y3 = load_sim(0.3,300,10,'CH4_[1]',lowres=True)
	x5,y5 = load_sim(0.5,300,10,'CH4_[1]',lowres=True)
	x7,y7 = load_sim(0.7,300,10,'CH4_[1]',lowres=True)	

	plt.plot(x1,y1,'steelblue', lw=2)
	plt.plot(x3,y3,'r', lw=2)
	plt.plot(x5,y5,'g', lw=2)
	plt.plot(x7,y7,'orange',lw=2)
	plt.plot(x1,y1-10,'k',lw=2,label='R~30k')

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	plt.title('Methane Pressure Dependence')

	plt.xlim(3896.8,3897.2)
	plt.ylim(0.2,1.1)

	# now plot low res versions were i've scaled the 0.1atm spectrum to be close

	# 
	plt.figure()
	x1,y1 = load_sim(0.1,300,10,'CH4_[1]',lowres=False)
	plt.plot(x1,y1,'steelblue',ls='--',label='p=0.1atm, 10cm')
	x1,y1 = load_sim(0.1,300,30,'CH4_[1]',lowres=False)
	plt.plot(x1,y1,'limegreen',ls='--',label='p=0.1atm, 30cm')
	plt.xlim(3894,3895.6)
	plt.ylim(0.2,1.1)
	x13,y13 = load_sim(0.1,300,30,'CH4_[1]',lowres=True)
	plt.plot(x13,y13,'steelblue',lw=2,label='p=0.1atm, 30cm, lowres')

	x3,y3 = load_sim(0.3,300,10,'CH4_[1]',lowres=True)
	plt.plot(x3,y3,'r', lw=2,label='p=0.3atm, 10cm, lowres')
	x3,y3 = load_sim(0.3,300,10,'CH4_[1]',lowres=False)
	plt.plot(x3,y3,'r',ls='--',label='p=0.3atm, 10cm')

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	plt.title('Methane Pressure Dependence')

	plt.xlim(3896.8,3897.2)
	plt.ylim(0.2,1.1)

	# again for p=0.5
	plt.figure()
	x1,y1 = load_sim(0.1,300,10,'CH4_[1]',lowres=False)
	plt.plot(x1,y1,'steelblue',ls='--',label='p=0.1atm, 10cm')
	x1,y1 = load_sim(0.1,300,50,'CH4_[1]',lowres=False)
	plt.plot(x1,y1,'pink',ls='--',label='p=0.1atm, 50cm')
	plt.xlim(3894,3895.6)
	plt.ylim(0.2,1.1)
	x13,y13 = load_sim(0.1,300,50,'CH4_[1]',lowres=True)
	plt.plot(x13,y13,'steelblue',lw=2,label='p=0.1atm, 50cm, lowres')

	x3,y3 = load_sim(0.5,300,10,'CH4_[1]',lowres=True)
	plt.plot(x3,y3,'g', lw=2,label='p=0.5atm, 10cm, lowres')
	x3,y3 = load_sim(0.5,300,10,'CH4_[1]',lowres=False)
	plt.plot(x3,y3,'g',ls='--',label='p=0.5atm, 10cm')

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	plt.title('Methane Pressure Dependence')

	plt.xlim(3896.8,3897.2)
	plt.ylim(0.2,1.1)


def load_gascell_scan(file='LabData/Spectra/B0183.20/B0183.1a.dpt'):
	"""
	load gas cell scan keeyoon made
	"""
	f = np.loadtxt(file,delimiter=',')

	x,y = f[:,0],f[:,1]

	return 1e7/x[::-1], y[::-1]
	#f = open(file,'r')
	#lines = f.readlines()
	#f.close()


def plot_scan_model():
	"""
	plot lab scan with corresponding model

	convolve with kpic resolution to see
	"""
	# parametrs from keeyoon's fits
	# EDIT THESE ONCE GET DATA
	# Using scanned gas cell properties from fit by keeyoon
	# gas cell
	species = np.array(['CH4_[1]','N2O_[1]', 'C2H2_[1]','CO2_[1]', 'H2O_[1]'])

	# define cell 1 pressures, convert to atm
	torr_to_atm = 0.00131579 #atm/torr
	press    = torr_to_atm * np.array([23.3, 4.9, 15.2, 102.2, 2.2]) #Torr
	airpress = torr_to_atm *22.3 # torr, added air pressure

	#define ratios
	p = round(sum(press) + airpress ,3)
	ratios = press/p
	lam,spec,title  = plot_spec_wl(species, p, 300, 1,ratios, lowres=False,ltot=10,res=500000)
	plt.plot([1950,4300],[0.95,0.95],'k--')
	plot_totspec(lam,spec,species,title,lowres)

	xx,yy = load_gascell_scan()
	plt.plot(1e7/xx,yy,'--',label='Cell 1 Scan Data',alpha=1)
	plt.xlim(1900,4300)
	plt.legend()

	# plot at low res
	lam,spec,title  = plot_spec_wl(species, p, 300, 1,ratios, lowres=True,ltot=10,res=30000)

def first_requested_gascell():
	species = np.array(['CH4_[1]','N2O_[1]', 'C2H2_[1]','CO2_[1]'])
	ratios  = np.array([0.21666, 0.1666, 0.1666, 0.45]) # ratios out of 1 of each gas
	# high res
	p,t,l,ltot,lowres = 0.3,300,1,10,False # **must** have l=1 for doing this
	x,spec,title  = plot_spec_wl(species, 0.3, 300, 1,ratios, lowres=lowres,ltot=10)
	plt.plot([1950,4300],[0.95,0.95],'k--')
	# low res
	p,t,l,ltot,lowres = 0.3,300,1,10,True # **must** have l=1 for doing this
	x,spec,title  = plot_spec_wl(species, 0.3, 300, 1, ratios, lowres=lowres,ltot=10)
	plt.plot([1950,4300],[0.95,0.95],'k--')



if __name__=='__main__':
	# new spectrum
	species = np.array(['CH4_[1]', 'C2H2_[1]','CO2_[1]'])
	pressures  = np.array([0.15, 0.05, 0.15])  # ratios out of 1 of each gas
	ratios     = pressures/np.sum(pressures)

	# high res
	#p,t,l,ltot,lowres = 0.3,300,1,10,False # **must** have l=1 for doing this
	#lam,spec,title  = plot_spec_wl(species, 0.3, 300, 1,ratios, lowres=lowres,ltot=10,res=1e6)
	#plt.plot([1950,4300],[0.95,0.95],'k--')

	# low res
	p,t,l,ltot,lowres = 0.35,300,1,10,True # **must** have l=1 for doing this
	lam,spec,title  = plot_spec_wl(species, p, 300, 1, ratios, lowres=lowres,ltot=10,res=30000)
	plt.plot([1950,4300],[0.95,0.95],'k--')

	plot_totspec(lam,spec,species,title,lowres)
	plt.plot([1950,4300],[0.95,0.95],'k--')





