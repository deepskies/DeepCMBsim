#!/usr/bin/env python
# coding: utf-8

# # CMB Temperature Map Simulations
#Importing Necessary Modules

# This code creates unlensed Q & U maps of the CMB and phi map FFTs which can be used to lens them using QuickLens.
# It takes in a CMB spectrum generated by CAMB (https://camb.info/readme.html).
# #### Order matters here. The code breaks if you import pymaster after numpy or matplotlib.
# Import modules 
import pymaster as nmt
#import healpy
# If you prefer to work on spherical geometry you can use HEALPix/healpy -> conda install healpy -c conda-forge or pip install healpy
#Important! If you are using this in a notebook, import pymaster *before* any matplotlib/numpy to avoid OpenMP conflicts

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import h5py
import flatmaps as fm
from astropy.wcs import WCS
from scipy.interpolate import interp1d
import random
import time

def simple_visit_fn(a, b, d):
    try:
        _, d[a] = len(b), np.array(b)
    except Exception:
        d[a] = bool(b)

def hassign(file, outdict):
    with h5py.File(file, "r") as f:
        f.visititems(lambda x, y: simple_visit_fn(x, y, outdict))

#Defining functions
def build_checkerboard(w, h) :
      """
      Function for building a checkerboard pattern.
      This checkerboard can then be used for testing lensing functions.
      """
      re = np.r_[ w*[0,0,0,0,0,0,1,1,1,1,1,1] ]              # even-numbered rows
      ro = np.r_[ w*[1,1,1,1,1,1,0,0,0,0,0,0] ]              # odd-numbered rows
      return np.row_stack(h*(re,re,re,re,re,re,ro,ro,ro,ro,ro,ro))

def map_parameters(pixels, degrees, projection="ZEA"):
    """
    Function for storing map parameters for the map generation.
    """

    #Calculating some parameters
    reso = degrees/pixels 
    reso_arcmin = reso*60
    dx = reso*np.pi/180.0
    
    w = WCS(naxis=2)
    
    w.wcs.crpix = [pixels/2, pixels/2] #Center pixel X, Y
    w.wcs.cdelt = np.array([-reso, reso])
    w.wcs.crval = [0, 0] #Center coordinates RA, DEC at 0,0
    w.wcs.ctype = ["RA---"+projection, "DEC--"+projection]
    
    flat_map_info = fm.FlatMapInfo(w, nx=int(pixels), ny=int(pixels), lx=degrees, ly=degrees)
    
    return flat_map_info

def add_zeroed_modes(spectrum, l_array=False):
    """
        Function to concatenate [0,0] to spectrum to account for zeroed-out l-mode for 0 and 1 modes
    The simulation code, synfast_flat, interprets the first two numbers as the 0- and 1-mode.
    Since our spectrum starts at the second mode, I set the amplitude for the 0-mode and 1-mode to 0.
    
    """
    if l_array==False:
        concatenated_spectrum = np.concatenate(([0,0],spectrum))
    else:
        concatenated_spectrum = np.concatenate(([0,1],spectrum))
    return concatenated_spectrum

def load_cmb_spectra(location, fix_scaling=None, add_zeroed=True, lensedCls=False):
    """
    Loading in pregenerated spectra using CAMB The file has l-modes plus the
    following spectra are respectively: Temperature-Temperature, E mode-E mode,
    Temperature-E mode cross spectrum, Phi-Phi, Phi-Temperature cross
    """
    cmb_spectra = {}
    hassign(location, cmb_spectra)
    cmb_spectra['l'] = cmb_spectra['l'].astype(int)

    #CAMB scales spectra, so we fix that before feeding it into the map generation code
    if fix_scaling is not None:
        for keyword in cmb_spectra.keys():
            if keyword == "clPP" or keyword == "clPT":
                 cmb_spectra[keyword] = cmb_spectra[keyword] / (cmb_spectra["l modes"]**4 * fix_scaling)
            elif keyword == "l modes":
                 continue
            elif "c" in keyword:
                 cmb_spectra[keyword] = cmb_spectra[keyword] / ((cmb_spectra["l modes"]*(cmb_spectra["l modes"]+1)/(2*np.pi)))

    #We have to add 2 zeros in to account for CAMB starting at an l-mode of 2
    first_modes = 0
    for x in cmb_spectra:
        try:
            first_modes += np.sum(x[:2])
        except Exception:
            continue
    add_zeroed = (first_modes>1)
    if add_zeroed:
         for keyword in cmb_spectra.keys():
             if keyword == "l modes":
                  cmb_spectra[keyword] = add_zeroed_modes(cmb_spectra[keyword], l_array=True)
             else:
                  cmb_spectra[keyword] = add_zeroed_modes(cmb_spectra[keyword])

    return cmb_spectra

def combine_spectra(spectra_without_BB, spectra_with_BB):
    """
    Combining spectra so we can include BB spectra in our estimates
    """
    #The spectra may go to higher l-modes, even with the same .ini file
    new_length = min(len(spectra_without_BB["l modes"]), len(spectra_with_BB["l modes"]))

    new_spectra = {}
    new_spectra["l modes"] = spectra_without_BB["l modes"][:new_length]
    new_spectra["clTT"] = spectra_without_BB["clTT"][:new_length]
    new_spectra["clTE"] = spectra_without_BB["clTE"][:new_length]
    new_spectra["clEE"] = spectra_without_BB["clEE"][:new_length]
    new_spectra["clPP"] = spectra_without_BB["clPP"][:new_length]
    new_spectra["clPT"] = spectra_without_BB["clPT"][:new_length]
    new_spectra["clBB"] = spectra_with_BB["clBB"][:new_length]

    return new_spectra
 
def spectrum_plot(l_modes, spectrum, title):
    """
    Function for plotting spectrum
    """
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.plot(l_modes, spectrum, label=title)
    plt.xscale("log")
    plt.yscale("log")
    plt.tick_params(labelsize=15)
    plt.xlabel("$\mathit{l}$", fontsize=30)
    plt.ylabel("$\lambda K$", fontsize=30)
    plt.title("TT Spectrum", fontsize=24)
    ax.legend(bbox_to_anchor=(1.55, 1.0))
    
    return fig

def frequency_domain_plot(input_map_fft, title, scaling, colorbar_label):
    """
    #Function for plotting a 2D power spectrum or something in the frequency domain
    """
    fig = plt.figure()
    plt.imshow(np.log10(np.abs(input_map_fft)), vmin=scaling[0], vmax=scaling[1])
    plt.colorbar(label=colorbar_label)
    plt.title(title)
    
    return fig

def calculate_and_plot_fft(input_map, title, real=False, scaling=[0,6], colorbar_label="Log of Map Strength"):
    """
    Function for plotting 2D FFT
    This needs to be an rfft for the actual analysis and a complex fft for estimating spectra
    """
    if real:
        input_map_fft = np.fft.rfft2(input_map)
    else:
        input_map_fft = np.fft.fft2(input_map)

    fig = frequency_domain_plot(input_map_fft, title, scaling, colorbar_label)

    return fig

def calculate_map_power_2D(input_map, fmi, real=False):
    """
    Function for making 2D Power Map
    FFT needs to be real for the analysis and complex estimating spectra
    """
    if real:
        input_map_fft = np.fft.rfft2(input_map)
    else:
        input_map_fft = np.fft.fft2(input_map)
    
    #Calculate tfac (factor for converting from frequency to spacial domain and vice versa
    tfac = fmi.lx * np.pi/(180.*fmi.nx**2)

    #Taking the conjugate of the phi map fft
    input_map_fft_conj = np.conj(input_map_fft)
    
    #Calculating the 2D Map Power
    map_power_2D = input_map_fft*input_map_fft_conj*tfac**2
    #input tfac accounts for conversion between frequency-domain and real-domain
    #it's calculated using map and fft parameters: tfac = resolution_in radians /pixels

    return map_power_2D

def calculate_power_spectrum(input_map, fmi, bin_num=80, pixels=192, lstep=72, method_histogram=False):
    """
    Function for creating power spectrum plot. This function has two methods: one
    which involves summing over annuli to estimate the power in different l-modes
    and a second method (which I took from João Caldeira) that uses a histogram
    function to weight bins of l-modes
    """

    #Calculate tfac (factor for converting from frequency to spacial domain and vice versa
    tfac = fmi.lx * np.pi/(180.*fmi.nx**2)

    #Calculate the maximum achievable l-mode given map parameters
    peak_l = pixels/2*lstep
    
    #This is the first (and default) method of calculating the power spectrum:
    if method_histogram==False:
        map_power_2D = calculate_map_power_2D(input_map, fmi, real=False)

        #Dividing data into four quarters, rotating so lowest mode is in top left quarter
        halfpix = pixels//2

        quarter1 = map_power_2D[:halfpix, :halfpix]

        quarter2 = map_power_2D[:halfpix, halfpix:]
        quarter2 = np.rot90(quarter2)

        quarter3 = map_power_2D[halfpix:, :halfpix]
        quarter3 = np.rot90(quarter3, k=3)

        quarter4 = map_power_2D[halfpix:, halfpix:]
        quarter4 = np.rot90(quarter4, k=2)

        #Creating empty arrays:
        spectra_est = np.zeros((bin_num)) #spectum sum before averaging
        spectra_est_avg_bins = np.zeros((bin_num)) #spectrum estimate after averaging
        n = np.zeros((bin_num))#number of points going into each bin/annulus

        #Defining edges of bins
        bin_edges = np.linspace(0,np.sqrt(2)*peak_l, bin_num+1)
    
        #Defining width of 1 bin
        bin_step = bin_edges[1]

        #Calculating bin centers
        bin_centers = (bin_edges[1:]+bin_edges[:-1])/2

        #Creating array containing each quarter of data
        quarter_list = [quarter1, quarter2, quarter3, quarter4]

        #Looping over each quarter
        for quarter in quarter_list:
    
            #Looping over row
            for i, spec in enumerate(quarter):
                lx_ = 72.*i #the l-mode in the x-direction
    
                #Looping over each pixel in each row
                for j, pwr in enumerate(quarter):
            
                    ly_ = 72.*j #the l-mode in the y-direction
            
                    #Calculate overall l-mode for each pixel
                    point_on_circ = np.sqrt(lx_**2+ly_**2) #overall l-mode

                    #Determine which been l-mode belongs to by seeing which is the largest lefthand edge smaller than pixel's overall l-mode
                    left_bin_edge = np.where(bin_edges <= point_on_circ)[0][-1] #largest bin that's smaller than overall l-mode
            
                    #Add power from that pixel to spectrum estimate
                    spectra_est[left_bin_edge] += quarter[i,j] #adding power at coordinate (i,j) to bin
            
                    #Tally the pixel as being in this particular been (for averaging later)
                    n[left_bin_edge] += 1 #counting the number of points in the bin for averaging later
                
        #Divide the spectrum by the number of pixels that went into each bin
        spectra_est_avg_bins = spectra_est / n
        
    #This is the second (histogramming) method
    else:
        #Defining bin sizes
        lbins = np.linspace(0, peak_l, num=bin_num)
        bin_centers = 0.5*(lbins[0:-1] + lbins[1:])
        
        #Create grid to define l-mode in x- and y-directions and overall l as their geometric sum
        lx, ly = np.meshgrid(np.linspace(-peak_l+lstep, peak_l, int(pixels)),np.linspace(-peak_l+lstep, peak_l, int(pixels)))
        ell = np.sqrt(lx**2 + ly**2)
        
        #Histogramming l-modes into bins for normalization later 
        norm, bins = np.histogram(ell, bins=lbins)
        norm = norm.astype(float)
        norm[ np.where(norm != 0.0) ] = 1./norm[ np.where(norm != 0.0) ]
        
        #Create 2D power map
        temperature_map_real_fft = np.fft.fft2(input_map)
        temperature_map_real_fft_conj = np.conj(temperature_map_real_fft)
        real_map_power_2D = temperature_map_real_fft*temperature_map_real_fft_conj
        #"spectrum_2D" is the same as the powermap but recentered so the 0-mode is in the middle
        spectrum_2D = np.fft.fftshift(real_map_power_2D*tfac**2)
        
        #Histogramming l-modes from above, weighting by the 2D powerspectrum
        spectra_est_avg_bins, bins = np.histogram(ell, bins=lbins, weights=spectrum_2D)

        #Normalizing to account for different number of l-modes in eac bin
        spectra_est_avg_bins *= norm
    
    return bin_centers, spectra_est_avg_bins

def plot_spectra(l_modes, spectra, plot_title, plot_labels):
    """
    Function to plot a single or multiple spectra The spectra must be loaded with
    their own l-modes, you can not give a different number of l-modes and spectra
    """    
    fig = plt.figure()
    
    #If there are multiple sets of l-modes, we plot them and their respective spectra
    if (type(l_modes) == list) and len(l_modes)>1:
        for i, spectrum in enumerate(spectra):
            plt.plot(l_modes[i],spectra[i],label=plot_labels[i])
            
    elif type(l_modes == np.ndarray) and len(np.shape(l_modes))>1:
        for i, spectrum in enumerate(spectra):
            plt.plot(l_modes[i],spectra[i],label=plot_labels[i])
            
        
    #If there's only one l-mode given, it's assumed you only gave it one spectrum
    #If you gave it multiple spectra for one l-mode, you fucked up.
    else:
        plt.plot(l_modes,spectra,label=plot_labels)
    
    plt.legend(loc="lower left")
    plt.yscale("log")
    plt.xscale("log")
    plt.title(plot_title)
    plt.xlabel("$l$ mode")
    plt.ylabel("$C_{l}$")
    return fig

def generate_maps(spectra_dict, fmi, num_maps, pixels, temp_only=False, TQU_maps=False, TEB_maps=False, phi_map=False, give_seeds=None):
    """
    This function generates CMB maps.  The spectra are assumed to be from CAMB,
    with the numbers rescaled and zeroes added to the front to account for CAMB
    beginning its spectra at l=2.  This function can make temperature-only maps, OR
    temperature-only maps with Q and U modes, OR temperature-only maps with E and B
    modes, OR phi (gravitational deflection) maps only.  fmi is flat-map info. It
    is basically a class which contains information about themaps and is an
    imported module.
    """   
    #Estimating times for optimization purposes
    function_start_time = time.time()

    #Checking if given seeds are valid 
    if give_seeds is not None:
        #Are there the same number of seeds as maps?
        if not (len(give_seeds)==num_maps):
            print("You need to have the same number of seeds as maps.")
            return

        #Are all of the seeds integers?
        elif sum([np.mod(give_seeds[i],1) for i in range(0,num_maps)])>0:
            print("All seeds must be integers.")
            return

        #If all of these are good, set these as the input seeds
        else:
            input_seeds = give_seeds


    #Making sure we know what type of maps we're making
    if sum([temp_only, TQU_maps, TEB_maps, phi_map])>1:
        print("You can only pick only one of these options: temperature-map only, TQU maps, TEB maps or phi_map.")
        return None
   
    #This part makes only temperature maps
    elif temp_only:
    
        #Creating an array to store all of the maps
        all_maps = np.zeros([num_maps,int(pixels),int(pixels)])
        
        #Print statement for counting
        print("Number of maps generated: ")
 
        #Loop to create all maps
        for i in range(0,num_maps):

            #Creating map if seeds are not set            
            if give_seeds is None:
                #Creating each individual map
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad, [spectra_dict["clTT"]],[0], seed=random.randint(0,50000000))

            else:
                #Creating each individual map
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad, [spectra_dict["clTT"]],[0], seed=input_seeds[i])

            #Saving map in index i of array
            all_maps[i] = one_map_set

            #Printing how many maps have been finished (it can seem slow)
            if np.mod(i,100) == 0:
                print(i+1)
            
        
    #This part makes temperature maps with their respective Q and U modes
    elif TQU_maps:
        
        #Creating an array to store all of the maps
        all_maps = np.zeros([num_maps,3,int(pixels),int(pixels)])

        #Print statement for counting
        print("Number of maps generated: ")
        
        #Loop to create all maps
        for i in range(0,num_maps):

            #Creating map if seeds are not set            
            if give_seeds is None:
           
                #Creating *one set* of temperature maps plus its Q and U maps
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad,np.array([spectra_dict["clTT"], spectra_dict["clTE"], spectra_dict["clTB"], spectra_dict["clEE"], spectra_dict["clEB"], spectra_dict["clBB"]]),[0,2],seed=random.randint(0,50000000))

            else:
           
                #Creating *one set* of temperature maps plus its Q and U maps
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad,np.array([spectra_dict["clTT"], spectra_dict["clTE"], spectra_dict["clTB"], spectra_dict["clEE"], spectra_dict["clEB"], spectra_dict["clBB"]]),[0,2],seed=input_seeds[i])
            
            #Saving all three maps in index i of array
            all_maps[i] = one_map_set
            
            #Printing how many maps have been finished (it can seem slow)
            if np.mod(i,100) == 0:
                print(i+1)

    #This part makes temperature maps with their respective E and B modes
    elif TEB_maps:
        
        #Creating an array to store all the maps
        all_maps = np.zeros([num_maps,3,int(pixels),int(pixels)])

        #Print statement for counting
        print("Number of maps generated: ")
        
        #Loop to create all maps
        for i in range(0,num_maps):

            #Creating map if seeds are not set            
            if give_seeds is None:
 
                #Creating *one set* of temperature maps plus its E and B maps
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad,np.array([spectra_dict["clTT"], spectra_dict["clTE"], spectra_dict["clTB"], spectra_dict["clEE"], spectra_dict["clEB"], spectra_dict["clBB"]]),[0,0,0], seed=random.randint(0,50000000))
            
            else:

                #Creating *one set* of temperature maps plus its E and B maps
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad,np.array([spectra_dict["clTT"], spectra_dict["clTE"], spectra_dict["clTB"], spectra_dict["clEE"], spectra_dict["clEB"], spectra_dict["clBB"]]),[0,0,0], seed=input_seeds[i])

            #Saving all three maps in index i of aray
            all_maps[i] = one_map_set
            
            #Printing how many maps have been finished (it can seem slow)
            if np.mod(i,100) == 0:
                print(i+1)
           
    #This part makes phi maps only
    elif phi_map:
        
        #Creating an array to store all the maps
        all_maps = np.zeros([num_maps,int(pixels),int(pixels)])

        #Print statement for counting
        print("Number of maps generated: ")
        
        #Loop to create all maps
        #If these maps look odd, you may have forgotten to nomalize the phi spectrum from CAMB, which is done automatically using load_cmb_spectra.
        for i in range(0,num_maps):
            
            #Creating map if seeds are not set            
            if give_seeds is None:

                #Creating one phi (gravitational deflection) map
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad, [spectra_dict["clPP"]],[0], seed=random.randint(0,50000000))
            
            else:

                #Creating one phi (gravitational deflection) map
                one_map_set = nmt.synfast_flat(int(fmi.nx),int(fmi.ny),fmi.lx_rad,fmi.ly_rad, [spectra_dict["clPP"]],[0], seed=input_seeds[i])

            #Saving map to index i of array
            all_maps[i] = one_map_set
            
            #Printing how many maps have been finished (it can seem slow)
            if np.mod(i,100) == 0:
                print(i+1)
        
    else:
        print("What type of maps did you want?")
        return None
    
    maps_generation_time = time.time()

    return all_maps #, [function_start_time, maps_generation_time]
