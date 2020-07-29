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

class tebmap(ql.maps.pix):

    def __init__(self, nx, dx, maps=None, ny=None, dy=None):
        """ class which contains temperature (T) and polarization (E, B) maps. """
        
        super( tebmap, self ).__init__(nx, dx, ny=ny, dy=dy)
        if maps is None:
            self.tmap = numpy.zeros( (self.ny, self.nx) )
            self.emap = numpy.zeros( (self.ny, self.nx) )
            self.bmap = numpy.zeros( (self.ny, self.nx) )
        else:
            [self.tmap, self.emap, self.bmap] = maps

        assert( (self.ny, self.nx) == self.tmap.shape )
        assert( (self.ny, self.nx) == self.emap.shape )
        assert( (self.ny, self.nx) == self.bmap.shape )

    def copy(self):
        return tebmap( self.nx, self.dx,
                       [self.tmap.copy(), self.emap.copy(), self.bmap.copy()],
                       ny = self.ny, dy = self.dy )

    def pad(self, nxp, nyp):
        """ make a new map with dimensions nxp (>nx), nyp (>ny) with this map at its center. """
        assert( nxp > self.nx )
        assert( nyp > self.ny )
        assert( numpy.mod( nxp - self.nx, 2 ) == 0 )
        assert( numpy.mod( nyp - self.ny, 2 ) == 0 )

        ret = tebmap( nx=nxp, dx=self.dx, ny=nyp, dy=self.dy )
        for this, that in [ [self.tmap, ret.tmap], [self.emap, ret.emap], [self.bmap, ret.bmap] ]:
            that[ (nyp-self.ny)/2:(nyp+self.ny)/2, (nxp-self.nx)/2:(nxp+self.nx)/2 ] = this
        return ret

    def threshold(self, vmin, vmax=None, vcut=0.):
        """ returns a new, thresholded version of the current map.
        threshold(v) -> set all pixels which don't satisfy (-|v| < val < |v|) equal to vcut.
        threshold(min,max) -> set all pixels which don't satisfy (vmin < val < vmax) equal to vcut.
        """
        if vmax is None:
            vmin = -numpy.abs(vmin)
            vmax = +numpy.abs(vmin)
        assert( vmin < vmax )

        ret = self.copy()
        for m in [ret.tmap, ret.emap, ret.bmap]:
            m[numpy.where(m < vmin)] = vcut
            m[numpy.where(m > vmax)] = vcut
        return ret

    def compatible(self, other):
        """ check whether this map can be added, subtracted, etc. to the map 'other'. """
        return ( hasattr(other, 'tmap') and
                 hasattr(other, 'emap') and
                 hasattr(other, 'bmap') and
                 super(tebmap, self).compatible(other) )

    def get_rfft(self):
        """ return real ffts of map objects """
        ret_t = ql.maps.rfft( self.nx, self.dx, ny = self.ny, dy = self.dy )
        ret_e = ql.maps.rfft( self.nx, self.dx, ny = self.ny, dy = self.dy )
        ret_b = ql.maps.rfft( self.nx, self.dx, ny = self.ny, dy = self.dy )
 
        #Calculating conversion factor
        tfac = numpy.sqrt((self.dx * self.dy) / (self.nx * self.ny))

        ret_t.fft[:] = numpy.fft.rfft2(self.tmap) * tfac
        ret_e.fft[:] = numpy.fft.rfft2(self.emap) * tfac
        ret_b.fft[:] = numpy.fft.rfft2(self.bmap) * tfac
        return ret_t, ret_e, ret_b

    def get_tqu(self):
        """ return a tqumap object containing the fourier transform of the T,Q,U maps. """

        #Calculating FFT of E & B maps
        tfft, efft, bfft = self.get_rfft()

	#Calculating l-modes in x- and y- directions
	lx, ly = efft.get_lxly()

	#Calculating angle
	tpi  = 2.*numpy.arctan2(lx, -ly)

	#Calculating FFT conversion factor
        tfac = numpy.sqrt((self.dx * self.dy) / (self.nx * self.ny))

        #Calculating Q & U Maps
        qmap = numpy.fft.irfft2(numpy.cos(tpi)*efft.fft - numpy.sin(tpi)*bfft.fft) / tfac
        umap = numpy.fft.irfft2(numpy.sin(tpi)*efft.fft + numpy.cos(tpi)*bfft.fft) / tfac

        ret = ql.maps.tqumap( self.nx, self.dx, maps=[self.tmap, qmap, umap])
        return ret

    def get_tebfft(self):
        """ return a QuickLens tebfft object """
        #Getting tebfft object
	ret = tebfft( self.nx, self.dx, ny = self.ny, dy = self.dy)

        #Calculating frequency-domain/spacial-domain conversion factor
        tfac = numpy.sqrt((self.dx * self.dy) / (self.nx * self.ny))

        #Calculating FFTs
        ret.tfft[:] = numpy.fft.rfft2(self.tmap) * tfac
        ret.efft[:] = numpy.fft.rfft2(self.emap) * tfac
        ret.bfft[:] = numpy.fft.rfft2(self.bmap) * tfac

        return ret


def plot_map(skymap, title, colorscheme="viridis", save_loc="figures/", filename=None):
     """Simple plotting function"""
     fig = plt.figure()
     plt.title(title)
     plt.imshow(skymap, cmap=colorscheme)
     plt.colorbar()
     if not filename:
          print("Saving figure "+title+".png"+" in "+save_loc)
          plt.savefig(save_loc+title+".png")
     else:
          print("Saving figure "+filename+" in "+save_loc)
          plt.savefig(save_loc+filename)
     return fig 

def get_tebmap(tqumap):
     """Function for getting tebmap object from tqumap object"""

     #Getting FFT object for TEB
     teb_fft = tqumap.get_teb()

     #Inverse FFT to get maps
     emap = numpy.fft.irfft2(teb_fft.efft)
     bmap = numpy.fft.irfft2(teb_fft.bfft)

     #Loading into tebmap structure
     teb_map = tebmap(tqumap.nx, tqumap.dx, [tqumap.tmap, emap, bmap])

     return teb_map

def lens_maps(cmb_maps, phi_fft, TEB_maps=False, TQU_maps=False, return_TQU=True, apodize_mask=None):
     """
     Function for lensing maps. Must be given QuickLens-style tebmap or tqumap
     object and QuickLens-style phi fft object
     """
     #Checking input maps
     if not (sum([TEB_maps,TQU_maps])==1):
           print("Do the maps that you're lensing have E & B modes or Q & U modes?")
           return

     #Prepping to-be-lensed-maps
     if TEB_maps:
           unlensed_maps = cmb_maps.get_tqu()

     elif TQU_maps:
           unlensed_maps = cmb_maps

     #This is the lens-y bit
     lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_maps, phi_fft)

     #Apodizing maps and/or converting to T, E and B maps
     if apodize_mask and return_TQU:
        qmap_lensed = apodize_mask*lensed_tqu.qmap
        umap_lensed = apodize_mask*lensed_tqu.umap
        lensed_maps = tqumap(cmb_maps.nx, cmb_maps.dx, [lensed_tqu.tmap, qmap_lensed, umap_lensed])

     elif apodize_mask and not return_TQU:
        qmap_lensed = apodize_mask*lensed_tqu.qmap
        umap_lensed = apodize_mask*lensed_tqu.umap
        lensed_tqu = tqumap(cmb_maps.nx, cmb_maps.dx, [lensed_tqu.tmap, qmap_lensed, umap_lensed])
        lensed_maps = get_tebmap(lensed_tqu)

     elif not apodize_mask and return_TQU:
        lensed_maps = lensed_tqu

     elif not apodize_mask and not return_TQU:
        lensed_maps = get_tebmap(lensed_tqu)

     else:
        print("How did I even get this far in the 'if' statements?")

     return lensed_maps
        

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

#Loading in E & B maps 
temperature_map = numpy.load('all_unlensed_temperature_maps.npy')
print("Finished loading unlensed map")

#Testing class structure
teb = tebmap(nx, dx, maps=temperature_map[0])

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
for i, mapz in enumerate(temperature_map[0:1]):

	if numpy.mod(i,100)==0:
		print i,

        #Trying out tebmap structure
	teb = tebmap(nx,dx,maps=[mapz[0],mapz[1],mapz[2]])

	#Putting Q & U maps into QuickLens structure
	tqu_maps = teb.get_tqu()

        #Old structure for comparison
	#Setting E & B maps
	emap = ql.maps.rmap(nx, dx, map=mapz[1])
	bmap = ql.maps.rmap(nx, dx, map=mapz[2])
	
	#Calculating FFT of E & B maps
	efft = emap.get_rfft()
	bfft = bmap.get_rfft()

	#Calculating Q & U Maps
	qmap = numpy.fft.irfft2(numpy.cos(tpi)*efft.fft - numpy.sin(tpi)*bfft.fft) / tfac
	umap = numpy.fft.irfft2(numpy.sin(tpi)*efft.fft + numpy.cos(tpi)*bfft.fft) / tfac

	#Putting Q & U maps into QuickLens structure
	tqu_maps_old = ql.maps.tqumap(nx, dx, maps=[numpy.zeros((nx, nx)), qmap, umap])

	#Putting phi map FFT into QuickLens structure
	phi_map_fft = ql.maps.rfft(nx, dx, fft=phi_map_ffts[i]*tfac)

	#Lensing the temperature/Q/U map using the phi map FFT and the QuickLens data structure
	lensed_tqu = ql.lens.make_lensed_map_flat_sky(tqu_maps, phi_map_fft)
        lensed_tqu_old = ql.lens.make_lensed_map_flat_sky(tqu_maps_old, phi_map_fft)

	#Applying apodization
	unlensed_qmap = apod_mask*tqu_maps.qmap
	unlensed_umap = apod_mask*tqu_maps.umap
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

#Plotting! Lots of Plotting!
#Plotting the E map
plot_map(tqu_maps.qmap, "New Code Unlensed Q Map", colorscheme="plasma") 
plot_map(tqu_maps_old.qmap, "Old Code Unlensed Q Map", colorscheme="plasma") 
plot_map(tqu_maps.qmap - tqu_maps_old.qmap, "Diff'd Unlensed Q Maps", colorscheme="plasma") 

exit()

#Plotting the E map
plt.figure()
plt.title("Unlensed, Apodized E Map")
plt.imshow(emap.map, cmap="plasma")
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

#Now with the class system
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

exit()

print("Saving data")
numpy.save("q_maps", q_maps)
numpy.save("u_maps", u_maps)
numpy.save("e_maps", e_maps)
numpy.save("kappa_maps", kappa_maps)

exit()
