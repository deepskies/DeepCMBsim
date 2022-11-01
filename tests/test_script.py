#Code for lensing CMB maps
#The module performing the lensing is QuickLens.
#To lens maps, QuickLens requires unlensed TQU maps and a real FFT of the phi map.

#Sourcing the required modules.
import time
start_time = time.time()
import sys
sys.path.append("/home/samanthausman/quicklens")
sys.path.append("/home/samanthausman/quicklens/quicklens")
import numpy
import matplotlib.pyplot as plt
from simcmb import lens_and_apodize_maps as lens

imports_time = time.time()

#Defining map parameters
pixels = 192
degrees = 5 #5 degrees on each side
reso = float(degrees)/pixels #resolution in degrees
dx = reso*numpy.pi/180.0 #resolution in radians
tfac = numpy.sqrt((float(dx)**2) /(float(pixels)**2)) #converts from pixels to radians
print("Finished setting parameters")

parameters_time = time.time()

#Loading in E & B maps 
temperature_map = numpy.load('all_unlensed_temperature_maps.npy')
print("Finished loading unlensed map")

#Loading in phi map FFTs
phi_map_ffts = numpy.load('all_phi_map_ffts.npy')
print("Finished loading in phi map FFTs")

#Setting up struture to save unlensed, lensed Q & U maps and kappa maps#Loading in apodization mask
apod_mask = numpy.load("apod_highres.npy")
print("Finished loading apodization mask")

npy_load_time = time.time()

#Testing class structure
teb = lens.tebmap(pixels, dx, maps=temperature_map[0])

phi_maps = numpy.array([numpy.fft.irfft2(phi_map_ffts[0]), numpy.fft.irfft2(phi_map_ffts[2])])

function_maps = lens.lens_maps(temperature_map[0:2], phi_maps, dx, apodize_mask = apod_mask)

lens_time = time.time()

times = numpy.array([imports_time-start_time, parameters_time-imports_time, npy_load_time-parameters_time, lens_time-npy_load_time])
points = numpy.arange(1,len(times)+1)

#Plotting! Lots of Plotting!
fig, ax = plt.subplots()
ax.bar(points, times)
#ax.set_ticks(points)
#ax.set_xticklabels(['Start','Modules Imported','Parameters Set','.npy Files Loaded', 'Maps Lensed'], rotation=45)
plt.xlabel("Checkpoint #")
plt.ylabel("Seconds")
plt.savefig("timing.png")
exit()
