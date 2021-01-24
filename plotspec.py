# plot gas cell spectra
#https://docs.google.com/document/d/1Yy2HYR1b8KcXCd_9mAmnEvWcd6d2bwD3TzSvuuUzzoY/edit#
import numpy as np
import matplotlib.pylab as plt

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


def plot_spec_wl(species, p, t, l,lowres=True):
	plt.figure(1,figsize=(12,6))
	
	for i,mol in enumerate(species):
		x,y = load_sim(p,t,l,mol,lowres=lowres)
		plt.plot(1e7/x, y**(10/len(species)),label=mol + ', %satm' %round(p/len(species),2),alpha=.7,zorder=100)

	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Transmittance')
	
	title = 'Pres: %satm Temp: %sK Length: %scm' %(p, t, 10)
	plt.fill_between(1000*np.array([1.95,2.5]), 0, y2=1,zorder=-100,facecolor='gray',alpha=0.3)
	plt.fill_between(1000*np.array([2.85,4.2]), 0, y2=1,zorder=-100,facecolor='gray',alpha=0.3)

	plt.title(title)
	#plt.savefig(savename)


if __name__=='__main__':
	species = np.array(['CH4_[1]','C2H2'])
	plot_spec_wl(species, 0.3, 300, 1,lowres=True)




