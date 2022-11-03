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
tqu_unl = ql.maps.tqumap(nx, dx)
print("Finished setting up data structures")

#Loading in T, Q & U maps 
temperature_map = numpy.load('test_temperature_map.npy')
tqu_unl.tmap = temperature_map[0]
tqu_unl.qmap = temperature_map[1]
tqu_unl.umap = temperature_map[2]
print("Finished loading unlensed map")

#Getting unlensed E & B maps
teb_unlensed = tqu_unl.get_teb() #I believe this should be the same as teb_unl
unlensed_e_fft = ql.maps.rfft(nx, dx, fft=teb_unlensed.efft)
unlensed_b_fft = ql.maps.rfft(nx, dx, fft=teb_unlensed.bfft)
unlensed_e_map = unlensed_e_fft.get_rmap()
unlensed_b_map = unlensed_b_fft.get_rmap()
numpy.save("./unlensed_e_map", unlensed_e_map.map)
numpy.save("./unlensed_b_map", unlensed_b_map.map)

#Loading in the FFT of the phi map
phi_map_fft_raw = numpy.load("test_phi_map_fft.npy")
phi_map_fft = ql.maps.rfft(nx,dx, ny=nx, dy=dx, fft=phi_map_fft_raw*tfac)
print("Finished loading FFT of phi map")

#Generating the phi map from the FFT
phi_map = phi_map_fft.get_rmap()

#Getting l-modes for each pixel using built-in functionality
#Setting the 0-mode to 1 so be don't get Runtime Errors later in the code
ell2D = phi_map_fft.get_ell()
ell2D[0,0]=1

#Calculating the factor that relates the FFT of the phi map to the FFT of the kappa map
fac = (ell2D*(ell2D+1.0))/2

#Calculating the kappa map using QuickLens structure
kappa_map_fft = ql.maps.rfft(nx, dx, fft=phi_map_fft_raw*tfac*fac)

#Converting kappa map FFT to the kappa map
kappa = kappa_map_fft.get_rmap()
print("Finished calculating kappa map")

#Lensing the temperature/Q/U map using the phi map FFT and the QuickLens data structure
lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_unl, phi_map_fft)
lensed_teb = lensed_tqu.get_teb()
print("Finished lensing map")

#Use rfft structure to get out E map
e_fft = ql.maps.rfft(nx, dx, fft=lensed_teb.efft)
e_map = e_fft.get_rmap()
numpy.save("./lensed_e_map", e_map.map)

#Use rfft structure to get out B map
b_fft = ql.maps.rfft(nx, dx, fft=lensed_teb.bfft)
b_map = b_fft.get_rmap()
numpy.save("./lensed_b_map", b_map.map)

#Plotting! Lots of Plotting!
#Plotting the Ell map
plt.figure()
plt.title("L Map")
plt.imshow(ell2D,norm=matplotlib.colors.LogNorm(), cmap="plasma", vmin = 10**2)
plt.colorbar()
plt.savefig("figures/L-modes.png")

#Plotting the E map
plt.figure()
plt.title("E Map")
plt.imshow(e_map.map, cmap="bwr")
plt.colorbar()
plt.savefig("figures/my_lensed_e_map.png")

#Plotting the B map
plt.figure()
plt.title("B Map")
plt.imshow(b_map.map, cmap="bwr")
plt.colorbar()
plt.savefig("figures/my_lensed_b_map.png")

#Plotting the E map
plt.figure()
plt.title("Unlensed E Map")
plt.imshow(unlensed_e_map.map, cmap="bwr")
plt.colorbar()
plt.savefig("figures/my_unlensed_e_map.png")

#Plotting the B map
plt.figure()
plt.title("Unlensed B Map")
plt.imshow(unlensed_b_map.map, cmap="bwr")
plt.colorbar()
plt.savefig("figures/my_unlensed_b_map.png")

#Plotting the phi map
plt.figure()
plt.title("Phi Map")
plt.imshow(phi_map.map, cmap="seismic")
plt.colorbar()
plt.savefig("figures/my_phi_map.png")

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
plt.imshow(kappa.map, cmap="magma")
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
