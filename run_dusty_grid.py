# morgan.fraser@ucd.ie
# Code to create grid of DUSTY models using MARCS spectra


from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import subprocess
from glob import glob
import pysynphot as S



# Function to smooth a spectrum
def  smooth(wave, flux, smoothing_factor=10):
    
    rescaled_spectum = []
    
    for i in range(0, len(wave), smoothing_factor):

        binned_flux = np.sum(flux[i:i+smoothing_factor])/smoothing_factor
        binned_wave = np.sum(wave[i:i+smoothing_factor])/smoothing_factor
        
        rescaled_spectum.append([binned_wave, binned_flux])

    rescaled_spectum = np.asarray(rescaled_spectum)
    rescaled_spectum_pd = pd.DataFrame({'Lambda': rescaled_spectum[:, 0], 'Flux': rescaled_spectum[:, 1]})
    
    rescaled_spectum_pd.drop(rescaled_spectum_pd.tail(1).index,inplace=True) 
    
    return (rescaled_spectum_pd)



bp_F658N = S.ObsBandpass('wfc3,uvis2,f658n')
bp_F555W = S.ObsBandpass('wfc3,uvis2,f555w')
bp_F814W = S.ObsBandpass('wfc3,uvis2,f814w')
bp_F160W = S.ObsBandpass('wfc3,ir,f160w')
bp_F350LP = S.ObsBandpass('wfc3,uvis2,f350lp')
#bp_Ch1 = S.FileBandpass('Spitzer_IRAC.I1.dat')
#bp_Ch2 = S.FileBandpass('Spitzer_IRAC.I2.dat')

#models = glob("models/s*flx")
models = glob("models/s*.flx")

wave = pd.read_csv("flx_wavelengths.vac", names=['lambda'])

counter = 0

results_array = []
    
for model in models:

    #  Just to keep track of progress
    counter = counter+1
    print ("[",counter,"/",len(models),"] ",model)

    # Get input spectrum details from filename. Assume MARCS syntax
    teff = model[8:12]
    logg = model[14:18]
    mass = model[20:23]
    
    # Read spectrum
    spectrum = pd.read_csv(model, delim_whitespace=True, names=['flux'])
    
    # Smooth spectrum, convert to microns and trim
    spectrum_smooth = smooth(wave['lambda']/1e4, s3800['flux'])
    spectrum_smooth = spectrum_smooth[spectrum_smooth['Lambda'] < 6]
    
    # Write this to a file in format ready for DUSTY
    spectrum_smooth.to_csv('smoothed_spectrum.dat', index=False, header=False, sep=' ')

    # Now loop over the range of thicknesses we want to try
    for thickness in [2, 10, 20, 50]:

    # Open a template for the DUSTY input file.
        fin = open('2016jbu_script_template.inp', "rt")
        fout = open("2016jbu_script_generated.inp", "wt")
        
        # Replace the thickness with the value we want
        for line in fin:
            fout.write(line.replace('[THICKNESS]', str(thickness)))

        fin.close()
        fout.close()
        
        # Run DUSTY
        dusty_process = subprocess.run(["/Applications/DUSTY/dusty", "2016jbu_script_generated.inp"])

        # Open the output from DUSTY
        fin = open('2016jbu_script_generated.out', "rt")

        model_info_lines = []
        for num, line in enumerate(fin):
            # See how many optical depths were calculated
            if "models with" in line:
                no_of_optical_depths = (line.split()[0])
            
            # Find what line numbers in the file corespond to the info on output spectra
            if "========================================================================================" in line:
                models_line_start = num     
            model_info_lines = list(range((int(models_line_start)+1),(int(models_line_start)+1+int(no_of_optical_depths))))
            
            if (num in model_info_lines):
                
                # Get index of output spectrum file, and value of tau_V
                spectrum_index = line.split()[0]
                tau_V = line.split()[1]

                # Read in output spectrum
                output_spectrum = pd.read_csv('2016jbu_script_generated.s00'+str(spectrum_index),
                                            delim_whitespace=True, comment='#', skiprows=6,
                        names=['lambda', 'fTot', 'xAtt', 'xDs', 'xDe', 'fInp', 'TauTot', 'albedo'])
                
                # Drop duplicate lines
                output_spectrum.drop_duplicates(subset=['lambda'], inplace=True)
                                       
                # Now get spectrum into format needed for Synphot
                w = (output_spectrum['lambda']*1e4).to_numpy()
                f = (output_spectrum['fTot']/output_spectrum['lambda']).to_numpy()
                sp = S.ArraySpectrum(w, f*1e-15, name='MySource',
                         fluxunits='flam', waveunits='angstrom')
                
                obs_f658n = S.Observation(sp, bp_F658N)
                obs_f350lp = S.Observation(sp, bp_F350LP)
                obs_f555w = S.Observation(sp, bp_F555W)
                obs_f814w = S.Observation(sp, bp_F814W)
                obs_f160w = S.Observation(sp, bp_F160W)
#                obs_ch1 = S.Observation(sp, bp_Ch1)
#                obs_ch2 = S.Observation(sp, bp_Ch2)

                mag_f350lp = obs_f350lp.effstim("vegamag")
                mag_f555w = obs_f555w.effstim("vegamag")
                mag_f814w = obs_f814w.effstim("vegamag")
                mag_f160w = obs_f160w.effstim("vegamag")
                mag_f658n = obs_f658n.effstim("vegamag")
#                mag_ch1    = obs_ch1.effstim("vegamag")
#                mag_ch2    = obs_ch2.effstim("vegamag")
                
 #               print (teff, logg, mass, thickness, tau_V, mag_f555w, mag_f814w, mag_f160w)
                results_array.append([teff, logg, mass, thickness, tau_V,
                                      mag_f350lp, mag_f555w, mag_f814w, mag_f160w,
                                      mag_f658n])  #, mag_ch1, mag_ch2])
                
        fin.close()

# Now write results to a csv. Use a pandas df to format nicely
df = pd.DataFrame(results_array,
                 columns=['teff', 'logg', 'mass', 'thickness', 'tau_V',
                          'mag_f350lp', 'mag_f555w', 'mag_f814w', 'mag_f160w',
                          'mag_f658n'])  # , 'mag_ch1', 'mag_ch2'])
df.to_csv("DUSTY_results.csv",  index=False)
