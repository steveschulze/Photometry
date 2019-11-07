#!/usr/bin/env python

__version__ = "8.0"
__author__ = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import	argparse
from 	astropy import table, time, wcs
from	astropy.io import ascii, fits
import	fits_tools
import	logging
from 	matplotlib import pylab as plt
import	numpy as np
import	os
import	phot_routines
import	sys

# Get input arguments

parser								= argparse.ArgumentParser(description='Programme for aperture photometry of HST images.')

# Object properties

parser.add_argument('--ra',				type=str,
										help='RA(J2000) of the Object (HMS and DD system allowed). Note: if declination is negative, write " -12:20:20.2"). Required keyword.',
										required=True)

parser.add_argument('--dec', 			type=str,
										help='Dec (J2000) of the Object (HMS and DD system allowed). Required keyword.',
										required=True)

parser.add_argument('--fits',			type=str,
										help='Name of FITS file. Required keyword.',
										required=True)

# Photometry keywords

parser.add_argument('--ana-thresh',		type=float,
										help='Analysis threshold. Default: 3 sigma',
										default=5)

parser.add_argument('--det-thresh',		type=float,
										help='Detection threshold. Default: 3 sigma',
										default=5)

parser.add_argument('--ap-diam',		type=float,
										help='Aperture radii in arcsec (one or more value allowed). Default: [0.25, 0.5, 0.75, 1.00, 1.50, 2.00, 3.00]',
										default=[0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 2.50, 3.00], nargs='+')

parser.add_argument('--ap-inner-annulus',type=float,
										help='Diameter of the inner annulus of the background. Default: 1.25 x ap-diam',
										default=1.25)

parser.add_argument('--ap-outer-annulus',type=float,
										help='Diameter of the outer annulus of the background. Default: 3 x ap-diam',
										default=2.5)

# Other options

parser.add_argument('--auto',			action='store_true',
										help='Automatic mode? Default: False',
										default=False)

parser.add_argument('--bw',				action='store_true',
										help='Screen output in B/W? Default: False',
										default=False)

parser.add_argument('--centroid',		action='store_true',
										help='Centre on the most nearby object',
										default=False)

parser.add_argument('--keeptemp',		action	= 'store_true',
										help	= 'Keep temporary files',
										default	= False)

parser.add_argument('--loglevel',		type	= str,
										help	= 'Logger level (default: INFO, possible values: DEBUG, INFO, WARNING, ERROR, CRITICAL)',
										default	= 'INFO')

parser.add_argument('--outdir',			type=str,
										help='Output path. Default: \'results/\'',
										default='results/')

parser.add_argument('--sex-loglevel',	type	= str,
										help	= 'Sextractor logger level (default: WARNING, possible values: DEBUG, INFO, WARNING, ERROR, CRITICAL)',
										default	= 'WARNING')

parser.add_argument('--tol',			type=float,
										help='Tolerance of the cross-matching in arcsec. Default: 1 arcsec',
										default=2)

########

# Input parameters

args								= parser.parse_args()

# Set up logger

logger 								= logging.getLogger()
logger.setLevel(args.loglevel)

formatter 							= logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh 									= logging.FileHandler(args.fits.replace(args.fits.split('.')[-1], 'log'), mode='w')
fh.setLevel(args.loglevel)
fh.setFormatter(formatter)
logger.addHandler(fh)

# Initialise text output

if args.bw 							== "True":
	class bcolors:
	    HEADER 						= ''
	    OKBLUE 						= ''
	    OKGREEN 					= ''
	    WARNING 					= ''
	    FAIL 						= ''
	    ENDC 						= ''
	    BOLD 						= ''
	    UNDERLINE 					= ''
else:
	from misc import bcolors

# Create output directory

msg									= 'Step 1: Administration'
print(bcolors.HEADER + bcolors.BOLD +'\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

if args.outdir[-1] 					!= '/':
	args.outdir						+= '/'

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

logger.info('Output directory: %s' %args.outdir)

# Print command again

print(bcolors.OKGREEN + '\nCommand' + bcolors.ENDC)

cmd 								= 'photometry_hst.py '

for key in vars(args):

	if key in ['auto', 'bw', 'centroid']:
		value						= ''

		if vars(args)[key] 			== True:
			cmd += '--{key} {value} '.format(key=key.replace('_', '-'), value=value)

	else:
		if key == 'dec':
			dec 					= vars(args)[key] if (not ('-' in vars(args)[key]) and not('+' in vars(args)[key])) else '\"{value}\"'.format(value=vars(args)[key])
			value 					= dec

		elif key == 'ap_diam':
			value 					= ' '.join(str(x) for x in vars(args)[key])

		else:
			value 					= vars(args)[key]

		cmd 						+= '--{key} {value} '.format(key=key.replace('_', '-'), value=value)

print (cmd)

logger.info('Command: %s' %cmd)

# Consistency check

print(bcolors.OKGREEN + "\nIs the object in the image footprint?" + bcolors.ENDC)

try:
	x_exp, y_exp					= fits_tools.sky2xy (args.fits, RA=args.ra, DEC=args.dec)
	x_obs, y_obs					= x_exp, y_exp
	#print('Coordinates = %s, %s' %(x_exp, y_exp))
	catalog_cleaned_reg				= args.outdir+ os.path.basename(args.fits).split('.fits')[0]+'_xy.cat'
	ascii.write(np.array([x_obs, y_obs]).T, catalog_cleaned_reg, format='no_header', overwrite=True)

	print('Yes.')
	logger.info('Is the object in the image footprint: Passed')

except:
	msg								= 'Object (RA, DEC = {ra}, {dec}) is not in the image footprint. Check object coordinates and FITS file.'.format(
										ra=args.ra, dec=args.dec
										)
	print(bcolors.FAIL + msg + bcolors.ENDC)
	logger.error(msg)
	sys.exit()

# Convert from HMS to DD system

ra_dd, dec_dd						= fits_tools.convert_hms_dd(args.ra, args.dec)

# Generate a source catalogue with sextractor

msg									= 'Step 2: Generate source catalogue'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

sources								= phot_routines.sextractor_photometry(
																		ANALYSIS_THRESH	= args.ana_thresh,
																		DETECT_THRESH	= args.det_thresh,
																		FITS 			= args.fits, 
																		FLAG			= 'centroid',
																		GAIN			= 'CCDGAIN',
																		PATH			= args.outdir,
																		PHOT_APERTURES	= 1.0,
																		LOGGER=logger
																		)

# Processing header

hdulist								= fits.open(args.fits)
hdu_header							= hdulist[0].header

if len(hdulist) 					> 1:
	try: 
		hdulist[1].shape
		hdu_header					+= hdulist[1].header
		hdu_image					= hdulist[1].data
	except:
		hdu_header					+= hdulist[1].header
		hdu_image					= hdulist[0].data
else:
	hdu_image						= hdulist[0].data

pix2arcsec							= fits_tools.pix2arcsec(args.fits)

# Pick the science source

sources['DISTANCE']					= np.sqrt( (sources['XWIN_IMAGE'] - x_exp)**2 + (sources['YWIN_IMAGE'] - y_exp)**2 )
sources['DISTANCE']					*= pix2arcsec
sources['DISTANCE'].unit			= 'arcsec'
sources['DISTANCE'].format			= '.2f'

if len(sources) 					== 0:

	msg								= 'No objects detected in the image. Change sextractor thresholds'
	print(bcolors.FAIL + '\n{}\n'.format(msg) + bcolors.ENDC)
	logger.error(msg)

	sys.exit()

# If centroiding is chosen

if args.centroid:

	sources2						= sources[sources['DISTANCE'] <= args.tol]

	if len(sources2) 				== 0:
		msg							= 'No object detected {tol}\" from the source position. Change sextractor thresholds'.format(tol=args.tol)
		print(bcolors.FAIL + '\n{}\n'.format(msg) + bcolors.ENDC)
		logger.error(msg)
		sys.exit()

	else:
		sources2.sort('DISTANCE')
		x_obs, y_obs				= sources2['XWIN_IMAGE'][0], sources2['YWIN_IMAGE'][0]
		
		# Update object catalogue
		
		catalog_cleaned_reg			= args.outdir+ args.fits.split('.fits')[0]+'_xy.cat'
		ascii.write(np.array([x_obs, y_obs]).T, catalog_cleaned_reg, format='no_header', overwrite=True)

# Do photometry

msg									= 'Step 3: Background estimation (can take some time...)'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

image_rms							= phot_routines.background(hdu_image, SIGMA=5, SNR=5, NPIXEL=5, DILATE_SIZE=11)

msg									= 'Step 4: Perform aperture photometry'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

apertures							= np.array(args.ap_diam)
innerannulus						= args.ap_inner_annulus * apertures
outerannulus						= args.ap_outer_annulus * apertures

if args.ap_inner_annulus 			< 1:
	print(bcolors.WARNING + 'The sky annulus intersects with the source region. Check the \'ap_inner_annulus\' keyword' + bcolors.ENDC)

photometry							= phot_routines.hst_aperture_photometry(args.fits, np.array([x_obs, y_obs]), apertures / 2., innerannulus / 2., outerannulus / 2., pix2arcsec, image_rms)
ascii.write(photometry, args.outdir + args.fits.replace('fits', 'mag'), overwrite=True)

# Curve of growth

msg									= 'Step 5: Curve of growth analysis'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

phot_routines.hst_cog(args.fits, np.array([x_obs, y_obs]), args.ap_inner_annulus, args.ap_outer_annulus, pix2arcsec, image_rms, args.outdir)

# Make cutouts

msg									= 'Step 6: Make cutout'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

phot_routines.hst_make_cutout(args.fits, [x_obs, y_obs], [x_exp, y_exp], apertures / 2., innerannulus / 2., outerannulus / 2., pix2arcsec, args.outdir)

# Prepare output catalogue

msg									= 'Step 7: Prepare output catalogue'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

# Zeropoint table

zeropoint							= table.Table()
zeropoint['METHOD']					= ['MAG_APER_'+str(x) for x in range(len(apertures))]
zeropoint['ZP']						= ['{zp:.3f}'.format(zp=photometry['ZP_APER_' +  str(i)][0]) for i in range(len(apertures))]
zeropoint['ZP_ERRP']				= -99
zeropoint['ZP_ERRM']				= -99
zeropoint['NUMBER']					= -99
zeropoint['r(FWHM)']				= -99
zeropoint['d(px)']					= apertures / pix2arcsec
zeropoint['d(arcsec)']				= apertures
zeropoint['MAG_3UL_GLOB']			= [-2.5*np.log10(3*photometry['BKG_NOISE_' + str(i)][0]) + photometry['ZP_APER_' +  str(i)][0] for i in range(len(apertures))]

for key in ['d(px)', 'd(arcsec)']:
	zeropoint[key].format			= '.3f'

zeropoint['MAG_3UL_GLOB'].format	= '.3f'

# Science table

science								= table.Table(names=('PROPERTY', 'VALUE', 'ERROR+', 'ERROR-', 'COMMENT'), dtype=('S100', 'f', 'f', 'f', 'S100'))
science.add_row(['FILENAME', np.nan, np.nan, np.nan, args.fits])

for key in ['DATE-OBS', 'MJD', 'EXPTIME', 'NCOMBINE']:
	try:
		if key == 'DATE-OBS':
			science.add_row([key, np.nan, np.nan, np.nan, time.Time(hdu_header['DATE-OBS'], format='isot', scale='utc').isot])
		elif key == 'EXPTIME':
			science.add_row([key, np.nan, np.nan, np.nan, np.round(hdu_header[key], 2)])
		else:
			science.add_row([key, np.nan, np.nan, np.nan, hdu_header[key]])

	except:
		if key == 'NCOMBINE':
			science.add_row([key, np.nan, np.nan, np.nan, 1])
		elif key == 'MJD':
			science.add_row([key, np.nan, np.nan, np.nan, np.round(time.Time(hdu_header['DATE-OBS'], format='isot', scale='utc').mjd, 7)])
		else:
			science.add_row([key, np.nan, np.nan, np.nan, '...'])

science.add_row(['RA', ra_dd, np.nan, np.nan, 'degree'])
science.add_row(['DEC', dec_dd, np.nan, np.nan, 'degree'])

science.add_row(['X_IMAGE_EXP', x_exp, np.nan, np.nan, 'px'])
science.add_row(['Y_IMAGE_EXP', y_exp, np.nan, np.nan, 'px'])

if args.centroid and len(sources2) > 0:

	science.add_row(['X_IMAGE_OBS', np.round(x_obs, 3), np.nan, np.nan, 'px'])
	science.add_row(['Y_IMAGE_OBS', np.round(y_obs, 3), np.nan, np.nan, 'px'])

	science.add_row(['DISTANCE (px)', 	  np.round(sources2['DISTANCE'][0]/pix2arcsec, 3), 	np.nan, np.nan, 'px'])
	science.add_row(['DISTANCE (arcsec)', np.round(sources2['DISTANCE'][0], 3), 			np.nan, np.nan, 'arcsec'])

else:

	science.add_row(['X_IMAGE_OBS', np.round(x_exp, 3), np.nan, np.nan, 'px'])
	science.add_row(['Y_IMAGE_OBS', np.round(y_exp, 3), np.nan, np.nan, 'px'])

	science.add_row(['DISTANCE (px)', 	  0, np.nan, np.nan, 'px'])
	science.add_row(['DISTANCE (arcsec)', 0, np.nan, np.nan, 'arcsec'])

for i in range(len(apertures)):
	science.add_row(['FNU_APER_' + str(i), '{:.3e}'.format(float(photometry['FNU_APER_' + str(i)])), '{:.3e}'.format(float(photometry['FNUERR_APER_' + str(i)])), '{:.3e}'.format(float(photometry['FNUERR_APER_' + str(i)])), 'microJy'])
	if photometry['FNU_APER_' + str(i)] > 0:
		science.add_row(['MAG_APER_' + str(i), np.round(photometry['MAG_APER_' + str(i)], 3), np.round(photometry['MAGERRP_APER_' + str(i)], 3), np.round(photometry['MAGERRM_APER_' + str(i)], 3), 'mag'])
	elif photometry['FNU_APER_' + str(i)] <= 0:
		science.add_row(['MAG_APER_' + str(i), np.round(-2.5 * np.log10(3*photometry['FNUERR_APER_' + str(i)]) + 23.9, 3), np.nan, np.nan, 'mag'])
	science.add_row(['MAG_APER_' + str(i) + '_2sigma', np.round(-2.5 * np.log10(2*photometry['FNUERR_APER_' + str(i)]) + 23.9, 3), np.nan, np.nan, 'mag'])
	science.add_row(['MAG_APER_' + str(i) + '_3sigma', np.round(-2.5 * np.log10(3*photometry['FNUERR_APER_' + str(i)]) + 23.9, 3), np.nan, np.nan, 'mag'])

# Save to file

print(bcolors.OKGREEN + "\nZeropoint\n" + bcolors.ENDC)
print(zeropoint)

print(bcolors.OKGREEN + "\nScience \n" + bcolors.ENDC)
science.pprint(max_lines=-1)
print('\n')

msg									= 'Step 9: Write results to file'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

ascii.write(zeropoint,		args.outdir + args.fits.replace('.fits', '_zp.log'),	overwrite=True)
ascii.write(science,		args.outdir + args.fits.replace('.fits', '_phot.log'),	overwrite=True)

# Remove temporary files

if not args.keeptemp:
	msg 						= 'Step 10: Remove all temps'
	print(bcolors.HEADER + bcolors.BOLD + "\n{}\n".format(msg) + bcolors.ENDC)
	logger.info(msg)
	os.system('rm check_' + args.fits)
	os.system('rm ' + args.outdir + args.fits.replace('.fits', '_xy.cat'))
	os.system('rm ' + args.outdir + args.fits.replace('.fits', '_centroid.log'))
else:
	msg 							= 'Step 10: Keep all temps'
	print(bcolors.HEADER + bcolors.BOLD + "\n{}\n".format(msg) + bcolors.ENDC)
	logger.info(msg)

# Show plots

if args.auto:
	plt.show()
plt.close()
