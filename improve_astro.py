#!/usr/bin/env pythonw

__version__ = "2.0"
__author__ = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import	argparse
from 	astropy.io import ascii, fits
import 	astropy.coordinates as coord
import 	astropy.units as u
import	fits_tools
from 	misc import bcolors
import	numpy as np
import	glob
import	os
import	subprocess
import 	sip_to_pv

# Prepare sextractor files

default_conv					= ("""CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
""")

default_nnw						= ("""NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:       9 for profile parameters + 1 for seeing.
# outputs:      ``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00 
 1.00000e+00
""")

default_param					= ("""NUMBER                   # Running object number
XWIN_IMAGE               # Windowed position estimate along x                        [pixel]
YWIN_IMAGE               # Windowed position estimate along y                        [pixel]
ERRX2WIN_IMAGE           # Variance of position along x                              [pixel**2]
ERRY2WIN_IMAGE           # Variance of position along y                              [pixel**2]
ERRXYWIN_IMAGE           # Covariance of position between x and y                    [pixel**2]
X2WIN_IMAGE              # Windowed variance along x                                 [pixel**2]
Y2WIN_IMAGE              # Windowed variance along y                                 [pixel**2]
XYWIN_IMAGE              # Windowed covariance between x and y                       [pixel**2]
ELONGATION               # A_IMAGE/B_IMAGE
ALPHAWIN_J2000           # Windowed right ascension (J2000)                          [deg]
DELTAWIN_J2000           # windowed declination (J2000)                              [deg]
FLAGS                    # Extraction flags
#IMAFLAGS_ISO             # FLAG-image flags ORed over the iso. profile
FWHM_IMAGE               # FWHM assuming a gaussian core                             [pixel]
CLASS_STAR               # S/G classifier output
FLUX_APER(3)             # Flux vector within fixed circular aperture(s)             [count]
FLUXERR_APER(3)          # RMS error vector for aperture flux(es)                    [count]
BACKGROUND               # Background at centroid position                           [count]
FLUX_MAX                 # Peak flux above background                                [count]
FLUX_AUTO                # Flux within a Kron-like elliptical aperture               [count]
FLUXERR_AUTO             # RMS error for AUTO flux                                   [count]
KRON_RADIUS              # Kron apertures in units of A or B
FLUX_ISO                 # Isophotal flux                                            [count]
FLUXERR_ISO              # RMS error for isophotal flux                              [count]
ISOAREA_IMAGE            # Isophotal area above Analysis threshold                   [pixel**2]
MU_MAX                   # Peak surface brightness above background                  [mag * arcsec**(-2)]
FLUX_RADIUS(2)           # Fraction-of-light radii                                   [pixel]
FLUX_PETRO               # Flux within a Petrosian-like elliptical aperture          [count]
FLUXERR_PETRO            # RMS error for PETROsian flux                              [count]
PETRO_RADIUS             # Petrosian apertures in units of A or B
SNR_WIN                  # Signal-to-noise ratio in a Gaussian window
""")

default_sex						= (
"""# Default configuration file for SExtractor 2.12.4
# EB 2010-10-10
#
 
#-------------------------------- Catalog ------------------------------------
 
CATALOG_NAME     test.cat       # name of the output catalog
CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,
                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC
PARAMETERS_NAME  default.param  # name of the file containing catalog contents
 
#------------------------------- Extraction ----------------------------------
 
DETECT_TYPE      CCD            # CCD (linear) or PHOTO (with gamma correction)
DETECT_MINAREA   3              # min. # of pixels above threshold

DETECT_THRESH    1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH  1.5            # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
 
FILTER           Y              # apply filter for detection (Y or N)?
FILTER_NAME      default.conv   # name of the file containing the filter
 
DEBLEND_NTHRESH  32             # Number of deblending sub-thresholds
DEBLEND_MINCONT  0.005          # Minimum contrast parameter for deblending
 
CLEAN            Y              # Clean spurious detections? (Y or N)?
CLEAN_PARAM      1.0            # Cleaning efficiency
 
#-------------------------------- WEIGHTing ----------------------------------

WEIGHT_TYPE      NONE           # type of WEIGHTing: NONE, BACKGROUND,
                                # MAP_RMS, MAP_VAR or MAP_WEIGHT
WEIGHT_IMAGE     weight.fits    # weight-map filename

#-------------------------------- FLAGging -----------------------------------

FLAG_IMAGE       flag.fits      # filename for an input FLAG-image
FLAG_TYPE        OR             # flag pixel combination: OR, AND, MIN, MAX
                                # or MOST

#------------------------------ Photometry -----------------------------------
 
PHOT_APERTURES   5              # MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS  2.5, 3.5       # MAG_AUTO parameters: <Kron_fact>,<min_radius>
PHOT_PETROPARAMS 2.0, 3.5       # MAG_PETRO parameters: <Petrosian_fact>,
                                # <min_radius>
PHOT_AUTOAPERS   0.0,0.0        # <estimation>,<measurement> minimum apertures
                                # for MAG_AUTO and MAG_PETRO

SATUR_LEVEL      50000.0        # level (in ADUs) at which arises saturation
SATUR_KEY        SATURATE       # keyword for saturation level (in ADUs)
 
MAG_ZEROPOINT    0.0            # magnitude zero-point
MAG_GAMMA        4.0            # gamma of emulsion (for photographic scans)
GAIN             0.0            # detector gain in e-/ADU
GAIN_KEY         GAIN           # keyword for detector gain in e-/ADU
PIXEL_SCALE      1.0            # size of pixel in arcsec (0=use FITS WCS info)
 
#------------------------- Star/Galaxy Separation ----------------------------
 
SEEING_FWHM      1.2            # stellar FWHM in arcsec
STARNNW_NAME     default.nnw    # Neural-Network_Weight table filename
 
#------------------------------ Background -----------------------------------
 
BACK_TYPE        AUTO           # AUTO or MANUAL
BACK_VALUE       0.0            # Default background value in MANUAL mode
BACK_SIZE        64             # Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE  3              # Background filter: <size> or <width>,<height>
 
#------------------------------ Check Image ----------------------------------
 
CHECKIMAGE_TYPE  NONE           # can be NONE, BACKGROUND, BACKGROUND_RMS,
                                # MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
                                # FILTERED, OBJECTS, -OBJECTS, SEGMENTATION,
                                # or APERTURES
CHECKIMAGE_NAME  check.fits     # Filename for the check-image
 
#--------------------- Memory (change with caution!) -------------------------
 
MEMORY_OBJSTACK  3000           # number of objects in stack
MEMORY_PIXSTACK  300000         # number of pixels in stack
MEMORY_BUFSIZE   1024           # number of lines in buffer
 
#------------------------------- ASSOCiation ---------------------------------

ASSOC_NAME       sky.list       # name of the ASCII file to ASSOCiate
ASSOC_DATA       2,3,4          # columns of the data to replicate (0=all)
ASSOC_PARAMS     2,3,4          # columns of xpos,ypos[,mag]
ASSOC_RADIUS     2.0            # cross-matching radius (pixels)
ASSOC_TYPE       NEAREST        # ASSOCiation method: FIRST, NEAREST, MEAN,
                                # MAG_MEAN, SUM, MAG_SUM, MIN or MAX
ASSOCSELEC_TYPE  MATCHED        # ASSOC selection type: ALL, MATCHED or -MATCHED

#----------------------------- Miscellaneous ---------------------------------
 
VERBOSE_TYPE     NORMAL         # can be QUIET, NORMAL or FULL
HEADER_SUFFIX    .head          # Filename extension for additional headers
WRITE_XML        N              # Write XML file (Y/N)?
XML_NAME         sex.xml        # Filename for XML output
XSL_URL          file:///usr/local/share/sextractor/sextractor.xsl
                                # Filename for XSL style-sheet


""")

# Get input arguments

parser							= argparse.ArgumentParser(description='Improving WCS calibration with astrometry.net and converting SIP terms to Sextractor format.')

# Object properties

parser.add_argument('--fits',			type	= str,
										help	= 'File name (Required)',
										required= True)

parser.add_argument('--ra',				type	= str,
										help	= 'RA(J2000) of the Object (HMS and DD system allowed). Note: if declination is negative, write " -12:20:20.2"). (Required)',
										required= True)

parser.add_argument('--dec', 			type	= str,
										help	= 'Dec (J2000) of the Object (HMS and DD system allowed). (Required)',
										required= True)

parser.add_argument('--radius', 		type	= float,
										help	= 'only search in indexes within \'radius\' of the field center (unit: deg; default: 0.125 deg)',
										default = 0.125)

parser.add_argument('--downsample',		type	= int,
										help	= 'downsample the image by factor <int> before running source extraction (default: 2)',
										default = 2)

parser.add_argument('--tweak-order',	type	= int,
										help	= 'Polynomial order of SIP WCS corrections (default: 2)',
										default = 2)

args							= parser.parse_args()

# Start

# Setup astrometry.net

output 							= open('default.conv', 'w')
output.write(default_conv)
output.close()

output 							= open('default.nnw', 'w')
output.write(default_nnw)
output.close()

output 							= open('default.param', 'w')
output.write(default_param)
output.close()

output 							= open('default.sex', 'w')
output.write(default_sex)
output.close()

sex_cfg							= 'default.sex'

log_astro						= open('astro.log', 'w')
log_sip2pv						= open('astro_sip2pv.log', 'w')

#for i in range(1):

ra_dd, dec_dd					= fits_tools.convert_hms_dd(args.ra, args.dec)

print(bcolors.HEADER + 'Process Image' + bcolors.ENDC)

# Check if fits file complies to standard
hdu								= fits.open(args.fits)
hdu_data						= hdu[0].data
hdu_header						= hdu[0].header
fits.writeto(args.fits, hdu_data, hdu_header, overwrite=True, output_verify='silentfix')

# Improve astro

sexcat							= args.fits.replace('.fits','.sexcat')

cmd								= ['solve-field', '--no-plots', \
								'--sextractor-config', sex_cfg, \
								'--x-column', 'XWIN_IMAGE', '--y-column', 'YWIN_IMAGE', \
								'--sort-column', 'FLUX_AUTO', \
								args.fits, \
								'--tweak-order', str(args.tweak_order), \
								'--uniformize', str(0), \
								'--cpulimit', str(60), \
								'--ra', str(ra_dd), '--dec', str(dec_dd), '--radius', str(args.radius), \
								'--new-fits', args.fits.replace('.fits', '_wcs.fits'), '--overwrite', '--downsample', str(args.downsample)]


print (' '.join(cmd))

process							= subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(stdoutstr,stderrstr) 			= process.communicate()

status							= process.returncode

if not os.path.isfile(args.fits.replace('.fits', '_wcs.fits')):
	# Log errors
	print(bcolors.FAIL + "ERROR: astrometry.net failed on " + str(args.fits) + bcolors.ENDC)
	log_astro.write(args.fits + '\n')
else:
	# Delete temp files
	cwd							= os.getcwd()
	if os.path.dirname(args.fits) != '':
		os.chdir(os.path.dirname(args.fits))
		os.system('rm *.axy *indx.xyls *.corr *.match *.rdls *.solved *.wcs')
		os.chdir(cwd)
	else:
		os.system('rm *.axy *indx.xyls *.corr *.match *.rdls *.solved *.wcs')
		
	# Transform SIP to PV to use distortion keywords in sextractor

	#print 
	status      = sip_to_pv.sip_to_pv(infile=args.fits.replace('.fits', '_wcs.fits'), outfile=args.fits.replace('.fits', '_wcs.fits'), tpv_format=True) #False gives warning line/breaks cloud in wcs.all_pix2world
	if status   == False:
		print(bcolors.FAIL + 'sip_to_pv failed. {} has only SIP keywords'.format(args.fits.replace('.fits', '_wcs.fits')))
		# log_sip2pv.write(args.fits.replace('.fits', '_wcs.fits') + '\n')
