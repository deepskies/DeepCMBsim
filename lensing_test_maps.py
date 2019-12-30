import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import quicklens as ql
import numpy
import matplotlib.pyplot as plt

print("Finished loading modules")
nx = 192
reso_arcmin = 5.0/nx*60
dx = reso_arcmin/60.0*numpy.pi/180.0 #reso_rad
pix        = ql.maps.pix(nx, dx)
lmax = 2500

tfac = numpy.sqrt((dx*dx) / (nx*nx))
print("Finished setting parameters")

cl_unl     = ql.spec.get_camb_scalcl(lmax=lmax)
teb_unl = ql.sims.tebfft(pix, cl_unl)
tqu_unl = teb_unl.get_tqu()
print("Finished setting up data structures")

#replace the tmap in tqu_unl by the checkboard 
checker = numpy.load('checkerboard_192.npy')
tqu_unl.tmap = checker
print("Finished loading unlensed map")

#loading fft of phi map
phi_map_fft_raw = numpy.load("phi_map_fft.npy")
phi_map_fft = ql.maps.rfft(nx,dx,fft=phi_map_fft_raw/tfac, ny=nx, dy=dx)
print("Finished loading FFT of phi map")

phi_map = phi_map_fft.get_rmap()

ell = numpy.arange(0, 191)
ell2D = numpy.zeros((192,97))
for m, i in enumerate(ell):
	for n, j in enumerate(ell[:97]):
		ell2D[m,n]= numpy.sqrt(i**2+j**2)

fac = 2.0/(ell2D*(ell2D+1.0))
kappamapfft = ql.maps.rfft(pix.nx, pix.dx, fft=phi_map_fft_raw/fac/tfac)
kappa = kappamapfft.get_rmap()
print("Finished calculating kappa map")

#lens
lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_unl, phi_map_fft)
print("Finished lensing map")

plt.figure()
plt.imshow(phi_map.map, interpolation='None', cmap='RdBu_r');
plt.colorbar()
#plt.show()
plt.savefig('192_phi_map.png')

plt.figure()
plt.imshow(kappa.map, interpolation='None', cmap='RdBu_r'); #colorbar() 
plt.colorbar()
#plt.show()
plt.savefig('192_kappa_map.png')

plt.figure()
plt.imshow(tqu_unl.tmap, interpolation='None', cmap='gray')
plt.colorbar()
#plt.show()
plt.savefig('192_unlensed_map.png')

plt.figure()
plt.imshow(lensed_tqu.tmap, interpolation='None', cmap='gray')
plt.colorbar()
#plt.show()
plt.savefig('192_lensed_map.png')
print("Figures saved")
