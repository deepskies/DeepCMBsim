import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import quicklens as ql
import numpy
import matplotlib.pyplot as plt
import matplotlib

print("Finished loading modules")
nx = 192
reso_arcmin = 5.0/nx*60
dx = reso_arcmin/60.0*numpy.pi/180.0 #reso_rad
pix        = ql.maps.pix(nx, dx)
lmax = numpy.pi / dx
lstep = lmax*2/nx
print(int(lstep))

tfac = numpy.sqrt((nx*nx) / (dx*dx))
print("Finished setting parameters")

cl_unl     = ql.spec.get_camb_scalcl(lmax=int(lmax))
teb_unl = ql.sims.tebfft(pix, cl_unl)
tqu_unl = teb_unl.get_tqu()
print("Finished setting up data structures")

#replace the tmap in tqu_unl by the checkboard 
checker = numpy.load('checkerboard_192.npy')
tqu_unl.tmap = checker
print("Finished loading unlensed map")

#loading fft of phi map
phi_map_fft_raw = numpy.load("phi_map_fft_0-mode.npy")
phi_map_fft = ql.maps.rfft(nx,dx, ny=nx, dy=dx, fft=phi_map_fft_raw/tfac)
print("Finished loading FFT of phi map")

phi_map = phi_map_fft.get_rmap()

#ell = numpy.arange(0, 192) * lstep
ell = numpy.fft.fftfreq(192)*lstep
ell2D = numpy.zeros((192,97))
for m, i in enumerate(ell):
	for n, j in enumerate(ell[:97]):
		ell2D[m,n]= numpy.sqrt(i**2+j**2)

fac = (ell2D*(ell2D+1.0))/2
kappamapfft = ql.maps.rfft(pix.nx, pix.dx, fft=phi_map_fft_raw*fac/tfac)
#try get_ell
kappa = kappamapfft.get_rmap()
print("Finished calculating kappa map")

#lens
lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_unl, phi_map_fft)
print("Finished lensing map")

plt.figure()
plt.title("Phi Map, Null 0-mode")
plt.imshow(phi_map.map, interpolation='None', cmap='RdBu_r');
plt.colorbar()
plt.xlabel("5 degrees")
plt.ylabel("5 degrees")
#plt.show()
plt.savefig('192_phi_map_0-mode.png')

plt.figure()
plt.title("Ell2D")
plt.imshow(ell2D, interpolation='None', cmap='magma');
plt.colorbar()
#plt.show()
plt.savefig('ell2D.png')

plt.figure()
plt.title("Ell Factor")
plt.imshow(fac, interpolation='None', cmap='inferno');
plt.colorbar()
#plt.show()
plt.savefig('ell_factor.png')

plt.figure()
plt.title("Phi Map FFT")
plt.imshow(abs(phi_map_fft_raw),norm=matplotlib.colors.LogNorm(), cmap="magma")
plt.colorbar()
plt.savefig("phi_map_fft.png")

plt.figure()
plt.title("Kappa Map FFT")
plt.imshow(abs(phi_map_fft_raw/fac/tfac), cmap="inferno", norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.savefig("kappa_map_fft.png")

plt.figure()
plt.title("Kappa Map, Null 0-mode")
plt.imshow(kappa.map, interpolation='None', cmap='RdBu_r'); #colorbar() 
plt.colorbar()
#plt.show()
plt.savefig('192_kappa_map_0-mode.png')


plt.figure()
plt.title("Unlensed Checkerboard Map")
plt.imshow(tqu_unl.tmap, interpolation='None', cmap='gray')
plt.colorbar()
#plt.show()
plt.savefig('192_unlensed_map.png')

plt.figure()
plt.title("Lensed Checkerboard")
plt.imshow(lensed_tqu.tmap, interpolation='None', cmap='gray')
plt.colorbar()
#plt.show()
plt.savefig('192_lensed_map.png')
print("Figures saved")
