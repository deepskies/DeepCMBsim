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
cl_unl = ql.spec.get_camb_scalcl(lmax=int(lmax))
teb_unl = ql.sims.tebfft(pix, cl_unl)
tqu_unl = teb_unl.get_tqu()
print("Finished setting up data structures")

#Loading in T, Q & U maps 
temperature_map = numpy.load('test_temperature_map.npy')
tqu_unl.tmap = temperature_map[0]
tqu_unl.qmap = temperature_map[1]
tqu_unl.umap = temperature_map[2]
print("Finished loading unlensed map")

#Loading in the FFT of the phi map
phi_map_fft_raw = numpy.load("test_phi_map_fft.npy")
phi_map_fft = ql.maps.rfft(nx,dx, ny=nx, dy=dx, fft=phi_map_fft_raw*tfac)
print("Finished loading FFT of phi map")

#Generating the phi map from the FFT
phi_map = phi_map_fft.get_rmap()

#Calculating an array to represent the l-modes in one direction
ell = numpy.fft.fftfreq(192)*lstep

#Calculating the l-modes in a 2D-array for converting the phi map to kappa map
#(There is probably a more efficient way to do this but this was the way I thought of)
ell2D = numpy.zeros((192,97))
for m, i in enumerate(ell):
	for n, j in enumerate(ell[:97]):
		ell2D[m,n]= numpy.sqrt(i**2+j**2)

#Setting the 0-mode to 1 so be don't get Runtime Errors later in the code
ell2D = phi_map_fft.get_ell()
ell2D[0,0]=1

#Calculating the factor that relates the FFT of the phi map to the FFT of the kappa map
fac = (ell2D*(ell2D+1.0))/2

#Calculating the kappa map using QuickLens structure
kappa_map_fft = ql.maps.rfft(nx, dx, fft=phi_map_fft_raw*fac)
#try get_ell

#Converting kappa map FFT to the kappa map
kappa = kappa_map_fft.get_rmap()
print("Finished calculating kappa map")

#Lensing the temperature/Q/U map using the phi map FFT and the QuickLens data structure
lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_unl, phi_map_fft)
lensed_teb = lensed_tqu.get_teb()
print("Finished lensing map")

#Plotting! Lots of Plotting!
#Plotting the E map
plt.figure()
plt.title("Lensed E Map")
plt.imshow(numpy.fft.irfft2(lensed_teb.efft),norm=matplotlib.colors.LogNorm(), cmap="plasma")
plt.colorbar()
plt.savefig("figures/my_lensed_e_map.png")

numpy.save("E_map", numpy.fft.irfft2(lensed_teb.efft))

#Plotting the phi map
plt.figure()
plt.title("Phi Map")
plt.imshow(phi_map.map,norm=matplotlib.colors.LogNorm(), cmap="seismic")
plt.colorbar()
plt.savefig("figures/my_phi_map_fft.png")

#Plotting the FFT of the phi map
plt.figure()
plt.title("Phi Map FFT")
plt.imshow(abs(phi_map_fft_raw),norm=matplotlib.colors.LogNorm(), cmap="viridis")
plt.colorbar()
plt.savefig("figures/my_phi_map_fft.png")

#Plotting the FFT of the kappa map
plt.figure()
plt.title("Kappa Map FFT")
plt.imshow(abs(phi_map_fft_raw/fac/tfac), cmap="inferno", norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.savefig("figures/my_kappa_map_fft.png")

#Plotting the kappa map
plt.figure()
plt.title("Kappa Map")
plt.imshow(numpy.fft.irfft2(phi_map_fft_raw/fac/tfac), cmap="magma")
plt.colorbar()
plt.savefig("figures/my_192_kappa_map.png")

#Plotting the unlensed temperature map
plt.figure()
plt.title("Unlensed Temperature Map")
plt.imshow(tqu_unl.tmap, interpolation='None', cmap='viridis')
plt.colorbar()
#plt.show()
plt.savefig('figures/my_192_unlensed_map.png')

#Plotting the lensed temperature map
plt.figure()
plt.title("Lensed Temperature Map")
plt.imshow(lensed_tqu.tmap, interpolation='None', cmap='viridis')
plt.colorbar()
#plt.show()
plt.savefig('figures/my_192_lensed_t_map.png')

#Plotting the lensed Q map
plt.figure()
plt.title("Lensed Q Map")
plt.imshow(lensed_tqu.qmap, interpolation='None', cmap='viridis')
plt.colorbar()
#plt.show()
plt.savefig('figures/my_192_lensed_q_map.png')

#Plotting the lensed U map
plt.figure()
plt.title("Lensed U Map")
plt.imshow(lensed_tqu.umap, interpolation='None', cmap='viridis')
plt.colorbar()
#plt.show()
plt.savefig('figures/my_192_lensed_u_map.png')
print("Figures saved")

exit()
