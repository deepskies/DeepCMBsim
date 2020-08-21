#Code for lensing CMB maps
#The module performing the lensing is QuickLens.
#To lens maps, QuickLens requires unlensed TQU maps and a real FFT of the phi map.

#Sourcing the required modules.
import pymaster as nmt
import time
import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import flatmaps
import quicklens as ql
import numpy
import matplotlib.pyplot as plt
import matplotlib
import lens_and_apodize_maps as lens
import flat_cmb_sims as fcs

start_time = time.time()

#Defining map parameters
pixels = 192.
degrees = 5.
fmi = fcs.map_parameters(pixels, degrees, projection="AIR")
reso = float(degrees)/pixels #resolution in degrees
dx = reso*numpy.pi/180.0 #resolution in radians
num_maps = 10
output = "TQU"
spectra_loc = "./base_plikHM_TTTEEE_lowl_lowE_lensing_scalCls.dat"
timing = False
delete = True

print("Finished setting parameters")

parameters_time = time.time()

spectra_dict = fcs.load_cmb_spectra(spectra_loc, fix_scaling=True, add_zeroed=True)

spectra_dict["clTB"] = 0 * spectra_dict["clTT"]
spectra_dict["clBB"] = 0 * spectra_dict["clTT"]
spectra_dict["clEB"] = 0 * spectra_dict["clTT"]
print("Finished loading spectra")

spectra_time = time.time()

unlensed_maps_timing, tqu_maps = fcs.generate_maps(spectra_dict, fmi, num_maps, pixels, TQU_maps=True)
print("Finished generating unlensed maps")

unlensed_maps_time = time.time()

phi_maps_timing, phi_maps = fcs.generate_maps(spectra_dict, fmi, 1, pixels, phi_map=True)
print("Finished generating phi maps")

phi_maps_time = time.time()

#Setting up struture to save unlensed, lensed Q & U maps and kappa maps#Loading in apodization mask
apod_mask = numpy.load("apod_highres.npy")
print("Finished loading apodization mask")

load_apod_mask_time = time.time()

function_maps = lens.lens_maps(tqu_maps, phi_maps, dx, output_type=output, apodize_mask = apod_mask)

lens_maps_time = time.time()

if delete:
   del function_maps
   del tqu_maps
   del phi_maps

x = numpy.arange(6)
y = numpy.zeros((6))
y[0] = parameters_time - start_time
y[1] = spectra_time - parameters_time
y[2] = unlensed_maps_time - spectra_time
y[3] = phi_maps_time - unlensed_maps_time
y[4] = load_apod_mask_time - phi_maps_time
y[5] = lens_maps_time - load_apod_mask_time


if timing:
    plt.bar(x, y)
    plt.xlabel("Checkpoint #")
    plt.ylabel("Time (s)")
    plt.title("Timing for %d Lensed %s Map" %(num_maps,output))
    plt.savefig("figures/timing_%d_%s_map.png" %(num_maps,output))

    numpy.save("timing_data/timing_%d_%s_data.dat"%(num_maps,output), y)

exit()
