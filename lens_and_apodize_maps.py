#Code for lensing CMB maps for Joao's ML code
#The module performing the lenses is QuickLens.
#Joao requires unlensed Q & U maps, lensed E maps and kappa maps.
#To lens maps, QuickLens requires unlensed T maps and a real FFT of the phi map.

#Sourcing the required modules.
import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import quicklens as ql
import numpy
import matplotlib.pyplot as plt
import matplotlib
print("Finished loading modules")

#Defining map parameters
pixels = 192. #192 pixels on each side
nx = int(pixels)
side = 5 #5 degrees on each side
reso = side/pixels #resolution in degrees
reso_arcmin = reso*60 #resolution in arcminutes
dx = reso*numpy.pi/180.0 #resolution in radians
lmax = 180./reso #maximum l-mode achievable given these parameters
lstep = lmax*2/pixels #increase in l-mode from one pixel to the next
tfac = dx/pixels #converts from pixels to radians
pix = ql.maps.pix(nx, dx)
print("Finished setting parameters")

#Setting up data structure 
#cl_unl = ql.spec.get_camb_scalcl(lmax=int(lmax))
#teb_unl = ql.sims.tebfft(pix, cl_unl)
#tqu_unl = teb_unl.get_tqu()
#print("Finished setting up data structures")

#Loading in E & B maps 
temperature_map = numpy.load('all_unlensed_temperature_maps.npy')
print("Finished loading unlensed map")

#Loading in phi map FFTs
phi_map_ffts = numpy.load('all_phi_map_ffts.npy')
print("Finished loading in phi map FFTs")

#Setting up struture to save unlensed, lensed Q & U maps and kappa maps
kappa_maps = numpy.zeros((len(temperature_map),int(pixels),int(pixels)))
e_maps = numpy.zeros((len(temperature_map),int(pixels),int(pixels)))
q_maps = numpy.zeros((len(temperature_map),int(pixels),int(pixels)))
u_maps = numpy.zeros((len(temperature_map),int(pixels),int(pixels)))

#Loading in apodization mask
apod_mask = numpy.load("apod_highres.npy")
print("Finished loading apodization mask")

#Setting up a single 2D FFT (easiest way to calculate the l-modes)
emap_for_structure = ql.maps.rmap(nx, dx, map=temperature_map[-1,1])
efft_for_structure = emap_for_structure.get_rfft()

#Getting l-modes
lx, ly = efft_for_structure.get_lxly()
ell2D = efft_for_structure.get_ell()
print("Finished calculating l-modes")

#Calculating factor to convert from kappa to phi FFTs
#Setting the 0-mode to 1 so be don't get Runtime Errors later in the code
fac = (ell2D*(ell2D+1.0))/2

#Calculating angle
tpi  = 2.*numpy.arctan2(lx, -ly)

#Starting loop
print("Entering loop")
print("# of maps completed:")
for i, mapz in enumerate(temperature_map):

	if numpy.mod(i,100)==0:
		print i,

	#Setting E & B maps
	emap = ql.maps.rmap(nx, dx, map=mapz[1])
	bmap = ql.maps.rmap(nx, dx, map=mapz[2])
	
	#Calculating FFT of E & B maps
	efft = emap.get_rfft()
	bfft = bmap.get_rfft()

	#Calculating Q & U Maps
	qmap = numpy.fft.irfft2(numpy.cos(tpi)*efft.fft - numpy.sin(tpi)*bfft.fft) / tfac
	umap = numpy.fft.irfft2(numpy.sin(tpi)*efft.fft + numpy.cos(tpi)*bfft.fft) / tfac

	#Calculating unlensed apodized maps
	qmap_un = qmap * apod_mask
	umap_un = umap * apod_mask

	#Putting Q & U maps into QuickLens structure
	tqu_maps = ql.maps.tqumap(nx, dx, maps=[numpy.zeros((nx, nx)), qmap_un, umap_un])

	#Getting the T, E & B maps to get unlensed, apodized E map
	teb = tqu_maps.get_teb()

	#Getting apodized emap
	efft_apod = ql.maps.rfft(nx, dx, fft=teb.efft)
	emap_apod = efft_apod.get_rmap()

	#Putting phi map FFT into QuickLens structure
	phi_map_fft = ql.maps.rfft(nx, dx, fft=phi_map_ffts[i]*tfac)

	#Lensing the temperature/Q/U map using the phi map FFT and the QuickLens data structure
	lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_maps, phi_map_fft)

	#Applying apodization
	unlensed_qmap = apod_mask*qmap
	unlensed_umap = apod_mask*umap
	qmap_lensed_apod = apod_mask*lensed_tqu.qmap
	umap_lensed_apod = apod_mask*lensed_tqu.umap

	#Calculating the kappa map using QuickLens structure
	kappa_map_fft = ql.maps.rfft(nx, dx, fft=phi_map_fft.fft*tfac*fac)

	#Converting kappa map FFT to the kappa map
	kappa = kappa_map_fft.get_rmap()

	#Saving maps
	kappa_maps[i] = kappa.map*apod_mask
	q_maps[i] = qmap_lensed_apod
	u_maps[i] = umap_lensed_apod
	e_maps[i] = emap_apod.map

#Plotting! Lots of Plotting!
#Plotting the E map
plt.figure()
plt.title("Unlensed, Apodized E Map")
plt.imshow(emap_apod.map, cmap="plasma")
plt.colorbar()
plt.title("E Map")
plt.savefig("figures/e_map_final.png")

#Plotting Q map
plt.figure()
plt.title("Lensed, Apodized Q Map")
plt.imshow(qmap_lensed_apod, cmap="viridis")
plt.colorbar()
plt.savefig("figures/q_map_final.png")

#Plotting the U map
plt.figure()
plt.title("Lensed, Apodized U Map")
plt.imshow(umap_lensed_apod, cmap="viridis")
plt.colorbar()
plt.savefig("figures/u_map_final.png")

#Plotting Q map
plt.figure()
plt.title("Unlensed, Apodized Q Map")
plt.imshow(qmap_un, cmap="viridis")
plt.colorbar()
plt.savefig("figures/q_map_un_final.png")

#Plotting the U map
plt.figure()
plt.title("Unlensed, Apodized U Map")
plt.imshow(umap_un, cmap="viridis")
plt.colorbar()
plt.savefig("figures/u_map_un_final.png")

#Plotting the kappa map
plt.figure()
plt.title("Apodized Kappa Map")
plt.imshow(kappa.map*apod_mask, cmap="bwr")
plt.colorbar()
plt.savefig("figures/kappa_map_final.png")

print("Saving data")
numpy.save("q_maps", q_maps)
numpy.save("u_maps", u_maps)
numpy.save("e_maps", e_maps)
numpy.save("kappa_maps", kappa_maps)

exit()
