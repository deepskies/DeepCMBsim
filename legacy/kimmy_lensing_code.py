import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import quicklens as ql
import numpy as np
import matplotlib.pyplot as plt


nx = 128
reso_arcmin = 5.0/nx*60
dx = reso_arcmin/60.0*np.pi/180.0 #reso_rad

pix        = ql.maps.pix(nx, dx)

lmax=2500


cl_unl     = ql.spec.get_camb_scalcl(lmax=lmax)

teb_unl = ql.sims.tebfft(pix, cl_unl)
tqu_unl = teb_unl.get_tqu()

#replace the tmap in tqu_unl by the checkboard 
checker = np.load('checkerboard.npy')
tqu_unl.tmap = checker

#kappa map load
kappa = np.load("kappa_100.npy")[0]

kapmap = ql.maps.rmap(pix.nx, pix.dx, map=kappa)
ell = np.arange(0, lmax+1)
fac = 2.0/(ell*(ell+1.0))
fac[0]=0.0
kappa_map_fft = kapmap.get_rfft()
phifft = kappa_map_fft*fac

#lens
lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_unl, phifft)

plt.figure()
plt.imshow(phifft.get_rmap().map, interpolation='None', cmap='RdBu_r'); #colorbar() 
plt.colorbar()
#plt.show()
plt.savefig('figures/kimmy_phi_map.png')

plt.figure()
plt.imshow(abs(kappa_map_fft.fft),interpolation='None', cmap='RdBu_r');
plt.colorbar()
plt.savefig("figures/kimmy_kappa_map_fft.png")

plt.figure()
plt.imshow(abs(phifft.fft), interpolation='None', cmap='viridis'); #colorbar() 
plt.colorbar()
#plt.show()
plt.savefig('figures/kimmy_phi_map_fft.png')

plt.figure()
plt.imshow(kappa, interpolation='None', cmap='RdBu_r'); #colorbar() 
plt.colorbar()
#plt.show()
plt.savefig('figures/kimmy_kappa_map.png')

plt.figure()
plt.imshow(tqu_unl.tmap, interpolation='None', cmap='gray')
#plt.show()
plt.savefig('figures/kimmy_checkerboard_unlensed.png')

plt.figure()
plt.imshow(lensed_tqu.tmap, interpolation='None', cmap='gray')
#plt.show()
plt.savefig('figures/kimmy_checkerboard_lensed.png')

