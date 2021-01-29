# plot gas cell spectra
#https://docs.google.com/document/d/1Yy2HYR1b8KcXCd_9mAmnEvWcd6d2bwD3TzSvuuUzzoY/edit#
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate

plt.ion()

def load_sim(p,t,l,mol,lowres=True):
	if lowres:
		f = np.loadtxt('HapiSimulations/' + 'transpec_p_%s_t_%s_l_%s_%s_lowres.txt' %(p,t,l,mol))
	else:
		f = np.loadtxt('HapiSimulations/' + 'transpec_p_%s_t_%s_l_%s_%s.txt' %(p,t,l,mol))

	return f[:,0],f[:,1]


def plot_spec_wn(species, p, t, l):
	plt.figure(fignum,figsize=(12,6))
	for i,mol in enumerate(species):
		x,y = load_sim(p,t,l,mol)
		plt.plot(x, y,label=mol)

	plt.legend()
	plt.xlabel('Wavenumber')
	plt.ylabel('Transmittance')

	plt.title(title)
	plt.savefig(savename)


def plot_spec_wl(species, p, t, l,lowres=True,ltot=10):
	plt.figure(1,figsize=(12,6))

	for i,mol in enumerate(species):
		x,y = load_sim(p,t,l,mol,lowres=lowres)
		if i==0: tot_spec = np.ones_like(y)
		if i==0: xtot_spec = 1e7/x[::-1]
		plt.plot(1e7/x, y**(ltot/len(species)),label=mol + ', %satm' %round(p/len(species),2),alpha=.7,zorder=100)
		# add to totl spectrum
		tck         = interpolate.splrep(1e7/x[::-1], y[::-1]**(ltot/len(species)), s=0)
		spec_interp = interpolate.splev(xtot_spec,tck,der=0,ext=0)

		tot_spec*=spec_interp

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	
	res=30000 if lowres else 1e6
	title = 'Pres: %satm Temp: %sK Length: %scm R:%s' %(p, t, ltot,res)
	plt.fill_between(1000*np.array([1.95,2.5]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)
	plt.fill_between(1000*np.array([2.85,4.2]), -0.2, y2=1.2,zorder=-100,facecolor='gray',alpha=0.3)

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

if __name__=='__main__':
	species = np.array(['CH4_[1]','N2O', 'CO2'])
	p,t,l,ltot,lowres = 0.3,300,1,10,False
	x,spec,title  = plot_spec_wl(species, 0.3, 300, 1,lowres=True,ltot=10)

	plot_totspec(x,spec,species,title,True)

