#Code for lensing CMB maps
#The module performing the lensing is QuickLens.
#To lens maps, QuickLens requires unlensed TQU maps and a real FFT of the phi map.

#Sourcing the required modules.
import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import quicklens as ql
import numpy
import matplotlib.pyplot as plt
import matplotlib

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
     #Calculating FFT/spacial map conversion factor
     tfac = numpy.sqrt((tqumap.dx**2) /(tqumap.nx**2))

     #Getting FFT object for TEB
     teb_fft = tqumap.get_teb()

     #Inverse FFT to get maps
     emap = numpy.fft.irfft2(teb_fft.efft)/tfac
     bmap = numpy.fft.irfft2(teb_fft.bfft)/tfac

     #Loading into tebmap structure
     teb_map = tebmap(tqumap.nx, tqumap.dx, [tqumap.tmap, emap, bmap])

     return teb_map

def lens_maps(cmb_maps, phi_maps, dx, input_type="TQU", output_type="TQU", apodize_mask=[]):
    """
    Function for lensing multiple sets of maps. Takes set of maps (either T, Q
    & U maps or T, E, & B maps with a phi map) and lenses them, returning
    the maps as a large numpy array. Phi map is assumed to not be an FFT.
    Input maps are assumed to be .npy files.
    For a set of temperature maps, the expected dimension is (num_maps, 3, nx, nx)
    For a set of phi maps, the expected dimension is (num_maps, nx, nx)
    For just one map, the expected dimension is (3, nx, nx)
    For just one phi map, the expected dimension is (nx, nx)
    """
    #Checking for the right number of maps
    assert(cmb_maps.ndim == phi_maps.ndim+1), "You need to have the same number of CMB map sets and phi maps."

    assert(numpy.shape(cmb_maps)[-1] == numpy.shape(phi_maps)[-1]), "The CMB maps and phi maps should have the same dimensionality."

    if cmb_maps.ndim == 3:
       num_maps = 1

    elif cmb_maps.ndim == 4:
       num_maps = numpy.shape(cmb_maps)[0]

    nx = numpy.shape(cmb_maps)[-1]

    #Adding an extra dimension if only one map
    if num_maps == 1:
       cmb_maps = numpy.expand_dims(cmb_maps, axis=0)
       phi_maps = numpy.expand_dims(phi_maps, axis=0)

    #Creating structure for output maps
    output_maps = numpy.zeros((num_maps, 3, nx, nx))

    #Beginning loop for lensing
    print("Beginning lensing process. Number of map sets lensed:")
    for i in range(0,num_maps):

       if numpy.mod(i,100) == 0:
          print i, #print(i, end = ", ") #<--Python3 print statement

       #Loading phi map into QuickLens class structure
       phi_fft = load_phi(phi_maps[i], nx, dx)

       #Prepping to-be-lensed-maps
       if input_type=="TQU":
          unlensed_maps = ql.maps.tqumap(nx, dx, maps=[cmb_maps[i,0], cmb_maps[i,1], cmb_maps[i,2]])
       elif input_type=="TEB":
          teb_map_set = tebmap(nx, dx, maps=[cmb_maps[i,0], cmb_maps[i,1], cmb_maps[i,2]])
          unlensed_maps = teb_map_set.get_tqu()
       else:
          print("Are the input maps 'TQU' or 'TEB' maps?")
          return

       #This is the lens-y bit
       lensed_tqu = ql.lens.make_lensed_map_flat_sky(unlensed_maps, phi_fft)
   
       #Apodizing maps or not
       if numpy.shape(apodize_mask)[0]>0:
          tmap_lensed = apodize_mask*lensed_tqu.tmap
          qmap_lensed = apodize_mask*lensed_tqu.qmap
          umap_lensed = apodize_mask*lensed_tqu.umap
     
       else:
          tmap_lensed = lensed_tqu.tmap
          qmap_lensed = lensed_tqu.qmap
          umap_lensed = lensed_tqu.umap
     
       #Converting to TEB or leaving as TQU
       if output_type=="TQU":
          lensed_maps = ql.maps.tqumap(nx, dx, [tmap_lensed, qmap_lensed, umap_lensed])
          output_maps[i] = [lensed_maps.tmap, lensed_maps.qmap, lensed_maps.umap]
 
       elif output_type == "TEB":
          lensed_tqu = ql.maps.tqumap(nx, dx, [tmap_lensed, qmap_lensed, -umap_lensed])
          lensed_maps = get_tebmap(lensed_tqu)
          output_maps[i] = [lensed_maps.tmap, lensed_maps.emap, lensed_maps.bmap]

       else:
          print("Do you want the input maps to be 'TQU' maps or 'TEB' maps?")

    print("Finished lensing maps.")
    
    return output_maps

def load_phi(phi_map, nx, dx, is_fft=False):
    """
    Function for loading phi map or phi map fft into QuickLens structure.
    """

    if is_fft:
       #This assumes you've already multiplied by the conversion factor
       phi_map_fft = ql.maps.rfft(nx, dx, fft=phi_map)

    else:
       #Calculating FFT/spacial map conversion factor
       tfac = numpy.sqrt((dx**2) /(nx**2))

       #Calculating FFT
       phi_fft = numpy.fft.rfft2(phi_map)*tfac
       
       #Loading into QuickLens class structure
       phi_map_fft = ql.maps.rfft(nx, dx, fft=phi_fft)

    return phi_map_fft
