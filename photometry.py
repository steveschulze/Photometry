#!/usr/bin/env pythonw

__version__ 						= "8.0"
__author__ 							= "Steve Schulze (steve.schulze@weizmann.ac.il)"

import	argparse
from 	astropy import table, time
from	astropy.io import ascii, fits
from 	astropy import units as u
import	cat_tools
import	fits_tools
import	logging
from 	matplotlib import pylab as plt
import	numpy as np
import	os
import	phot_routines
from	plotsettings_py36 import *
import	sys

# Get input arguments

parser								= argparse.ArgumentParser(description='Programme for aperture photometry.')

# Object properties

parser.add_argument('--ra',				type	= str,
										help	= 'RA(J2000) of the Object (HMS and DD system allowed). Note: if declination is negative, write " -12:20:20.2"). (required)',
										required= True)

parser.add_argument('--dec', 			type	= str,
										help	= 'Dec (J2000) of the Object (HMS and DD system allowed) (required)',
										required= True)

parser.add_argument('--fits',			type	= str,
										help	= 'Name of FITS file (required)',
										required= True)

parser.add_argument('--host-offset',	type	= float,
										help	= 'Host offset (default: 10 arcsec)',
										default	= 5)

# Properties of the reference catalogue

parser.add_argument('--ref-cat',		type	= str,
										help	= 'Which reference catalogue should be used (SDSS, 2MASS, PS1)?',
										default	= None)

parser.add_argument('--ref-filter',		type	= str,
										help	= 'Filter of the reference catalogue?',
										default	= None)

parser.add_argument('--ref-file',		type	= str,
										help	= 'Name of reference catalog file (if given overwrites \'--ref-cat\' option)',
										default	= None)

parser.add_argument('--ref-image',		type	= str,
										help	= 'Name of reference image to run sextractor in dual image mode?',
										default	= '')

parser.add_argument('--ref-radius',		type	= str,
										help	= 'Search radius in the reference catalogue query? (default: 10, unit: arcmin)',
										default	= 10)

# Photometry keywords

parser.add_argument('--ana-thresh',		type	= float,
										help	= 'Analysis threshold (default: 1 sigma)',
										default	= 1)

parser.add_argument('--ap-diam',		type	= float,
										help	= 'Aperture radii in FWHM (one or more) (default: [1.,2.,3.,4.])',
										default	= [1.,2.,3.,4.], nargs='+')

parser.add_argument('--ap-diam-ul',		type	= float,
										help	= 'Aperture radii in FWHM (one or more) (default: 2)',
										default	= 2,
										nargs	= '+')

parser.add_argument('--det-thresh',		type	= float,
										help	= 'Detection threshold (default: 1; unit: sigma)',
										default	= 1)

parser.add_argument('--gain',			type	= str,
										help	= 'Gain keyword in the header (default: Gain; unit: e-/ADU)',
										default	= None)

parser.add_argument('--back-size',		type	= int,
										help	= 'Background mesh: <size> (default: 64)',
										default	= 64)

parser.add_argument('--back-filtersize',type	= int,
										help	= 'Background filter: <size> (default: 3)',
										default	= 3)

parser.add_argument('--deblend-nthresh',type	= int,
										help	= 'Number of deblending sub-thresholds (default: 64)',
										default	= 64)

parser.add_argument('--deblend-mincont',type	= float,
										help	= 'Minimum contrast parameter for deblending (default: 0.00001)',
										default	= 0.00001)

# Zeropoint determination		
		
parser.add_argument('--mag-cut',		type	= float,
										help	= 'Remove all objects brighter than a given magnitude (default: 12)',
										default	= 12.)

parser.add_argument('--mag-stdfaint',  	type	= str, 
										help	= 'Lower magnitude cut for secondary standards (default: 0)',
										default	= 0)

parser.add_argument('--mag-stdbright', 	type	= str, 
										help	= 'Upper magnitude cut for secondary standards (default: 0)',
										default	= 0)

parser.add_argument('--maxstars', 		type	= int, 
										help	= 'Maximum number of stars used for building the local sequence (default: 200)',
										default	= 200)

# Other options

parser.add_argument('--auto',			action	= 'store_true',
										help	= 'Automatic mode? (default: False)',
										default	= False)

parser.add_argument('--bw',				action	= 'store_true',
										help	= 'Text output in color (default: False)',
										default	= False)

parser.add_argument('--keeptemp',		action	= 'store_true',
										help	= 'Keep temporary files',
										default	= False)

parser.add_argument('--loglevel',		type	= str,
										help	= 'Logger level (default: INFO, possible values: DEBUG, INFO, WARNING, ERROR, CRITICAL)',
										default	= 'INFO')

parser.add_argument('--outdir',			type	= str,
										help	= 'Output path efault: \'results/\'',
										default	= 'results/')

parser.add_argument('--sex-loglevel',	type	= str,
										help	= 'Sextractor logger level (default: WARNING, possible values: DEBUG, INFO, WARNING, ERROR, CRITICAL)',
										default	= 'WARNING')

parser.add_argument('--tol',			type	= float,
										help	= 'Tolerance of the cross-matching in arcsec (default: 1)',
										default	= 1)


########

# Input parameters
args								= parser.parse_args()

# Set up text output

if args.bw 							== True:
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

# Set up logger

logger 								= logging.getLogger()
logger.setLevel(args.loglevel)

formatter 							= logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh 									= logging.FileHandler(args.fits.replace(args.fits.split('.')[-1], 'log'), mode='w')
fh.setLevel(args.loglevel)
fh.setFormatter(formatter)
logger.addHandler(fh)

# 1) Start the programme

msg									= 'Step 1: Administration'
print(bcolors.HEADER + bcolors.BOLD +'\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

# Set up output dir

if args.outdir[-1] != '/':
	args.outdir		+= '/'

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

logger.info('Output directory: %s' %args.outdir)

print(bcolors.OKGREEN + '\nCommand' + bcolors.ENDC)

cmd 								= 'photometry.py '

for key in sorted(vars(args)):

	if key in ['auto', 'bw']:
		if vars(args)[key]:
			value					= ''
			cmd += '--{key} '.format(key=key)

	elif key == 'ref_image':
		if vars(args)[key] 			!= '':
			cmd += '--{key} {value} '.format(key=key.replace('_', '-'), value=vars(args)[key])

	elif key in ['ref_cat', 'ref_filter', 'ref_file']:

		if vars(args)['ref_file'] 	!= None:
			if not 'ref-file' in cmd:
				cmd 				+= '--{key} {value} '.format(key='ref-file', value=vars(args)['ref_file'])

		elif vars(args)['ref_file'] == None:
			if key 					!= 'ref_file':
				cmd 				+= '--{key} {value} '.format(key=key.replace('_', '-'), value=vars(args)[key])

	else:
		if key 						== 'dec':
			dec 					= vars(args)[key] if (not ('-' in vars(args)[key]) and not('+' in vars(args)[key])) else '\"{value}\"'.format(value=vars(args)[key])
			value 					= dec

		elif key 					== 'ap_diam':
			value 					= ' '.join(str(x) for x in vars(args)[key])

		else:
			value 					= vars(args)[key]

		cmd 						+= '--{key} {value} '.format(key=key.replace('_', '-'), value=value)

print (cmd)

logger.info('Command: %s' %cmd)

# Convert from HMS to DD system

ra_dd, dec_dd						= fits_tools.convert_hms_dd(args.ra, args.dec)

# Add information to the table 'object'

object_properties					= table.Table()
object_properties['OBJECT']			= [os.path.basename(args.fits).replace(args.fits.split('.')[-1], '')[:-1]]
object_properties['RA_HMS']			= args.ra
object_properties['DEC_DMS']		= args.dec
object_properties['RA']				= ra_dd
object_properties['DEC']			= dec_dd

# Consistency check

print(bcolors.OKGREEN + "\nIs the object in the image footprint?" + bcolors.ENDC)

try:
	x_exp, y_exp					= fits_tools.sky2xy (args.fits, RA=ra_dd, DEC=dec_dd)
	object_properties['X_EXP']		= x_exp
	object_properties['Y_EXP']		= y_exp
	filename_science_xy				= args.outdir + object_properties['OBJECT'][0] + '_xy.cat'
	ascii.write(np.array([x_exp, y_exp]).T, filename_science_xy, format='no_header', overwrite=True)

	# print('Coordinates = %s, %s' %(x_exp, y_exp))

	print('Yes.')
	logger.info('Is the object in the image footprint: Passed')

except:
	msg								= 'Object (RA, DEC = {ra}, {dec}) is not in the image footprint. Check object coordinates and FITS file.'.format(
										ra=args.ra, dec=args.dec
										)
	print(bcolors.FAIL + msg + bcolors.ENDC)
	logger.error(msg)
	sys.exit()

# If sextractor is run in dual image mode, process reference image first

'''
Get star catalogue for the reference image. Needed because the circular apertures are
defined according to the reference image
'''

'''
if args.ref_image != '':

	print(bcolors.HEADER + bcolors.BOLD + '\nProcess reference image' + bcolors.ENDC)

	# refimg_cat_name_orig		= args.fits.split('_')[0]+'_PS1_r.cat'
	# refimg_cat_name				= refimg_cat_name_orig.replace('.cat', '_refimg.cat')
	os.system('cp {input} {output}'.format(input=refimg_cat_name_orig, output=args.outdir+refimg_cat_name))

	# print('Reference image: %s' %args.ref_image)
	# print('Reference star cat: %s' %args.outdir+refimg_cat_name)

	# catalog_clean_ref			= fits_tools.sky2xy(args.ref_image, CAT=args.outdir+refimg_cat_name)
	# np.savetxt(args.outdir + refimg_cat_name.replace('.cat', '_xy.cat'), catalog_clean_ref)

	# ref_image					= phot_routines.sextractor_photometry(FITS=args.ref_image, PHOT_APERTURES="10", DETECT_THRESH=1, ANALYSIS_THRESH=1,\
	# 							ASSOC_NAME=args.outdir + refimg_cat_name.replace('.cat', '_xy.cat'), ASSOC_PARAMS="1,2", FLAG="refimg_stars", PATH=args.outdir)
'''

'''
Process photometric catalogue
'''

msg									= 'Step 2: Flux calibration'

print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

if args.ref_file != None:
	msg							= 'Use catalog file: ' + args.ref_file

	print(bcolors.OKGREEN + '\n' + msg + bcolors.ENDC)
	ref_cat						= ascii.read(args.ref_file, names=('RA', 'DEC', 'MAG', 'MAGERR'))

	logger.info(msg)

	object_properties['PHOTCAL']	= args.ref_file

	# Write catalogue to file

	filename_stars				= args.outdir + object_properties['OBJECT'][0] + '_' + os.path.basename(args.ref_file).replace(args.ref_file.split('.')[-1], '')[:-1] + '.cat'
	ascii.write(ref_cat, filename_stars.replace('.cat', '_refcat.cat'),
									overwrite=True,
									format='no_header')

	logger.info('Copy of the catalogue is stored in %s' %filename_stars)

else:
	# Retrieve catalogue from the VIZIER database

	msg								= 'Retrieve catalogue from VIZIER'
	print(bcolors.OKGREEN + "\n" + msg + bcolors.ENDC)
	
	logger.info(msg)

	object_properties['PHOTCAL']	= args.ref_cat + '/'+ args.ref_filter
	object_properties['CAT']		= args.ref_cat
	object_properties['FILTER']		= args.ref_filter

	# Unaltered input catalogue will be stored in the following file
	
	filename_stars					= args.outdir + object_properties['OBJECT'][0] + '_' + args.ref_cat + '_' + args.ref_filter + '.cat'
	
	cat_tools.retrieve_photcat(object_properties, args.ref_cat, phot_routines.catalog_prop,
									FILENAME= filename_stars,
									RADIUS	= args.ref_radius * u.arcmin,
									OUTDIR	= args.outdir)

	# Make copy of the catalogue
	os.system('cp {file1} {file2}'.format(file1=filename_stars, file2=filename_stars.replace('.cat', '_refcat.cat')))

	logger.info('Copy of the catalogue is stored in %s' %filename_stars)

# 2) Building the local sequence

msg									= 'Building the local sequence'

ref_cat								= ascii.read(filename_stars.replace('.cat', '_refcat.cat'), names=('RA', 'DEC', 'MAG', 'MAGERR'))

print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
logger.info(msg)

# Postprocess local sequence
# Remove objects that are brighter than a certain magnitude

ref_cat 							= ref_cat[ref_cat['MAG'] > args.mag_cut]
logger.info('Remove all stars that are brighter than {magcut:.2f} mag'.format(magcut=args.mag_cut))
ascii.write(ref_cat, filename_stars.replace('.cat', '_refcat.cat'), overwrite=True, format='no_header')

# Remove all objects that lie outside of the image footprint

print(filename_stars.replace('.cat', '_refcat.cat'))
catalog_clean 						= fits_tools.sky2xy(args.fits, CAT=filename_stars.replace('.cat', '_refcat.cat'))

if len(catalog_clean) == 0:
	msg								= 'No stars in the catalogue'
	print(bcolors.FAIL + msg + bcolors.ENDC)
	logger.error(msg)
	sys.exit()

# Write catalogue to file. Needed for forced photometry.

ascii.write(catalog_clean, filename_stars.replace('.cat', '_refcat_xy.cat'), overwrite=True, format='no_header')

print(bcolors.OKGREEN + '\nSelecting stars for local sequence' + bcolors.ENDC)

ref_stars							= phot_routines.sextractor_photometry(
																		ANALYSIS_THRESH	= 3,
																		ASSOC_NAME		= filename_stars.replace('.cat', '_refcat_xy.cat'),
																		ASSOC_PARAMS	= "1,2",
																		ASSOC_RADIUS	= args.tol/fits_tools.pix2arcsec(args.fits),
																		DETECT_THRESH	= 3,
																		FITS 			= args.fits,
																		FLAG			= "loc_seq",
																		GAIN			= args.gain,
																		LOGGER			= logger,
																		LOGLEVEL		= args.sex_loglevel,
																		PATH			= args.outdir,
																		PHOT_APERTURES	= "10",
																		REF_FILE		= args.ref_image)

logger.info('{numstars} stars identified to build local sequence'.format(numstars=len(ref_stars)))

# Apply quality cuts to select the good stars to build the local sequence

if args.mag_stdbright == 0 and args.mag_stdfaint == 0:

	# Mask stars
	# 1) Take stars that are not blended (FLAG == 0)
	# 2) Remove objects with an eccentricity of > 20%
	# 3) Take stars with PSFs within a certain range
	# 4) Apply photometry cuts
	# 5) Only consider stars with measurement errors between 0.001 and 0.2 mag

	logger.info('Pruning:')
	logger.info('1) Ellipticity of < 20%')
	logger.info('2) Objects with sextractor flag == 0')
	logger.info('3) FHWM(PSF) < MEDIAN + 1.5 x IQR')
	logger.info('4) Limiting the number to {numstars}'.format(numstars=args.maxstars))
	logger.info('5) Only consider stars with measurements errors between 0.001 and 0.2 mag')


	# Flag clipping
	ref_stars						= ref_stars[ref_stars['FLAGS'] == 0]
	
	# Remove elongated objects (eccentricity < 0.2)
	ref_stars						= ref_stars[np.sqrt(1 - (ref_stars['A_IMAGE'] / ref_stars['A_IMAGE'])**2) < 0.2 ]

	# PSF clipping
	# print (ref_stars)
	
	fwhm_p25						= np.percentile(ref_stars['FWHM_IMAGE'], 25)
	fwhm_p50						= np.percentile(ref_stars['FWHM_IMAGE'], 50)
	fwhm_p75						= np.percentile(ref_stars['FWHM_IMAGE'], 75)
		
	fwhm_mask_good					= np.where(ref_stars['FWHM_IMAGE'] < fwhm_p75 + 1.5 * (fwhm_p75 - fwhm_p25))[0]
	ref_stars 						= ref_stars[fwhm_mask_good]

	# Magnitude cut

	cut_high						= 0.20
	cut_low							= 0.001

	ref_stars						= ref_stars[ref_stars['MAGERR_APER'] <= cut_high]
	ref_stars						= ref_stars[ref_stars['MAGERR_APER'] >= cut_low]

	# print (len(ref_stars))

	# sigma clipping

	if len(ref_stars) > args.maxstars:
		print(bcolors.WARNING + 'List contains more than {maxstars} stars. Truncate the faint end.'.format(maxstars=args.maxstars) + bcolors.ENDC)

		ref_stars.sort('MAGERR_AUTO')
		ref_stars 					= ref_stars[:args.maxstars]

		# print (len(ref_stars))

	logger.info('{numstars} stars remain after pruning'.format(numstars=len(ref_stars)))

	# Crossmatch the sextractor catalog with the reference catalog

	msg								= 'Cross-match catalogues'
	print(bcolors.OKGREEN + '\n'+ msg + bcolors.ENDC)
	logger.info(msg)

	# matched_standard				= phot_routines.merge_table(ref_stars, ref_cat, 'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'RA', 'DEC', args.tol/3600.)['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE', 'MAG_APER', 'MAGERR_APER', 'MAG', 'MAGERR', 'FWHM_IMAGE', 'FLUX_RADIUS']
		
	ref_stars_keys					= ['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE', 'MAG_AUTO', 'MAGERR_AUTO',
									'MAG_PETRO', 'MAGERR_PETRO', 'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_PETRO', 'FLUXERR_PETRO',
									'FWHM_IMAGE', 'FWHM_WORLD', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'KRON_RADIUS', 'FLUX_RADIUS',
									'FLAGS', 'VECTOR_ASSOC', 'NUMBER_ASSOC', 'MAG_APER', 'MAGERR_APER', 'FLUX_APER', 'FLUXERR_APER']
		
	ref_cat_keys					= ref_cat.keys()
	
	ref_stars 						= np.asarray(ref_stars[ref_stars_keys]).view((float, len(ref_stars.dtype.names)))
	ref_cat 						= np.asarray(ref_cat[ref_cat_keys]).view((  float, len(ref_cat.dtype.names)))
		
	matched_standard				= cat_tools.wrapper_crossmatch(ref_stars, ref_cat, args.tol)
	matched_standard				= table.Table(matched_standard, names=ref_cat_keys+ref_stars_keys+['DIST'])
	matched_standard['DIST']		= matched_standard['DIST']*3600
	matched_standard['DIST'].format	= '.3f'

	# merge_table(ref_stars, ref_cat, 'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'RA', 'DEC', args.tol/3600.)['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE', 'MAG_APER', 'MAGERR_APER', 'MAG', 'MAGERR', 'FWHM_IMAGE', 'FLUX_RADIUS']

	matched_standard['MAG_APER'].name 		= 'MAG_INS'
	matched_standard['MAGERR_APER'].name 	= 'MAGERR_INS'
	matched_standard['MAG'].name 			= 'MAG_CAT'
	matched_standard['MAGERR'].name 		= 'MAGERR_CAT'

	# Define local sequence

	local_sequence					= phot_routines.local_sequence(matched_standard,
																		AUTO			= args.auto,
																		FILENAME		= filename_stars.replace('.cat', '_loc_cleaned.cat'),
																		FITS			= args.fits,
																		LOGGER			= logger,
																		LOWER			= 5,
																		PATH			= args.outdir,
																		UPPER			= 90)
#																		BW				= args.bw)

	matched_standard				= local_sequence['CAT']

	logger.info('Final number of stars in the local sequenc: {num_stars}'.format(num_stars=local_sequence['NUMSTARS']))
	print(bcolors.OKGREEN + '\nFinal number of stars in the local sequence: ' + str(local_sequence['NUMSTARS']) + bcolors.ENDC)

else:

	# If magnitude boundaries are defined == auto mode

	logger.info('Pruning:')
	logger.info('1) Ellipticity of < 20%')
	logger.info('2) Objects with sextractor flag == 0')
	logger.info('3) FHWM(PSF) < MEDIAN + 1.5 x IQR')
	logger.info('4) Limiting the number to {numstars}'.format(numstars=args.maxstars))
	logger.info('5) Only consider stars with measurements errors between 0.001 and 0.2 mag')

	ref_stars						= ref_stars[ref_stars['MAG_APER'] >= float(args.mag_stdbright)]
	ref_stars						= ref_stars[ref_stars['MAG_APER'] <= float(args.mag_stdfaint)]

	logger.info('Magnitude range set to: {magbright:.2f}, {magfaint:.2f}'.format(magbright=args.mag_stdbright, magfaint=args.mag_stdfaint))

	# Mask stars
	# 1) Take stars that are not blended (FLAG == 0)
	# 2) Remove objects with an eccentricity of > 20%
	# 3) Take stars with PSFs within a certain range
	# 4) Apply photometry cuts
	# 5) Only consider stars with measurement errors between 0.001 and 0.2 mag

	# Flag clipping
	ref_stars						= ref_stars[ref_stars['FLAGS'] == 0]
	
	# Remove elongated objects (eccentricity < 0.2)
	ref_stars						= ref_stars[np.sqrt(1 - (ref_stars['A_IMAGE'] / ref_stars['A_IMAGE'])**2) < 0.2 ]
			
	# PSF clipping			
			
	fwhm_p25						= np.percentile(ref_stars['FWHM_IMAGE'], 25)
	fwhm_p50						= np.percentile(ref_stars['FWHM_IMAGE'], 50)
	fwhm_p75						= np.percentile(ref_stars['FWHM_IMAGE'], 75)
			
	fwhm_mask_good					= np.where(ref_stars['FWHM_IMAGE'] < fwhm_p75 + 1.5 * (fwhm_p75 - fwhm_p25))[0]
	ref_stars 						= ref_stars[fwhm_mask_good]

	# sigma clipping

	if len(ref_stars) > args.maxstars:
		print(bcolors.WARNING + 'List contains more than {maxstars} stars. Truncate the faint end.'.format(maxstars=args.maxstars) + bcolors.ENDC)

		cut_high					= 0.20
		cut_low						= 0.001

		ref_stars					= ref_stars[ref_stars['MAGERR_APER'] <= cut_high]
		ref_stars					= ref_stars[ref_stars['MAGERR_APER'] >= cut_low]

		ref_stars.sort('MAGERR_AUTO')
		ref_stars 					= ref_stars[:args.maxstars]

	# Match the sextractor catalog with the reference catalog

	msg								= 'Cross-match catalogues'
	print(bcolors.OKGREEN + '\n'+ msg + bcolors.ENDC)
	logger.info(msg)

	ref_stars_keys					= ['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'XWIN_IMAGE', 'YWIN_IMAGE', 'MAG_AUTO', 'MAGERR_AUTO',
									'MAG_PETRO', 'MAGERR_PETRO', 'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_PETRO', 'FLUXERR_PETRO',
									'FWHM_IMAGE', 'FWHM_WORLD', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'KRON_RADIUS', 'FLUX_RADIUS',
									'FLAGS', 'VECTOR_ASSOC', 'NUMBER_ASSOC', 'MAG_APER', 'MAGERR_APER', 'FLUX_APER', 'FLUXERR_APER']
		
	ref_cat_keys					= ref_cat.keys()
	
	ref_stars 						= np.asarray(ref_stars[ref_stars_keys]).view((float, len(ref_stars.dtype.names)))
	ref_cat 						= np.asarray(ref_cat  [ref_cat_keys]).view((  float, len(ref_cat.dtype.names)))
		
	matched_standard				= cat_tools.wrapper_crossmatch(ref_stars, ref_cat, args.tol)
	matched_standard				= table.Table(matched_standard, names=ref_cat_keys+ref_stars_keys+['DIST'])
	matched_standard['DIST']		= matched_standard['DIST']*3600
	matched_standard['DIST'].format	= '.3f'

	matched_standard['MAG_APER'].name 		= 'MAG_INS'
	matched_standard['MAGERR_APER'].name 	= 'MAGERR_INS'
	matched_standard['MAG'].name 			= 'MAG_CAT'
	matched_standard['MAGERR'].name 		= 'MAGERR_CAT'

	# Define local sequence

	local_sequence					= phot_routines.local_sequence(matched_standard,
																		AUTO	= args.auto,
																		FILENAME= filename_stars.replace('.cat', '_loc_cleaned.cat'),
																		FITS	= args.fits,
																		LOGGER	= logger,
																		LOWER	= 10,
																		PATH	= args.outdir,
																		UPPER	= 90)
#																		BW		= args.bw)

	matched_standard				= local_sequence['CAT']

	logger.info('Final number of stars in the local sequenc: {num_stars}'.format(num_stars=local_sequence['NUMSTARS']))
	print(bcolors.OKGREEN + '\nFinal number of stars in the local sequence: ' + str(local_sequence['NUMSTARS']) + bcolors.ENDC)

# Image properties

if args.ref_image == '':
	if all(matched_standard['FWHM_IMAGE'] == 0):
		# Sometimes FWHM_IMAGE only return zeros. The reasons are not clear:
		# https://www.astromatic.net/forum/showthread.php?tid=615
		# http://www.astromatic.net/forum/showthread.php?tid=458
		# https://www.astromatic.net/forum/showthread.php?tid=947

		print(bcolors.WARNING + '\nAll FWHM_IMAGE values are identical to zero. FLUX_RADIUS is used instead.' + bcolors.ENDC)

		FWHM_median					= np.median(matched_standard['FLUX_RADIUS'])*2/1.1	
	else:		
		FWHM_median					= np.median(matched_standard['FWHM_IMAGE'])
else:		
	FWHM_median						= fwhm_p50

apertures							= FWHM_median * np.array(args.ap_diam)

# Step 3: Compute zeropoint

msg									= 'Step 3: Zeropoint calculation'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

msg									= 'Run Sextractor'
print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
logger.info(msg)


phot_stars							= phot_routines.sextractor_photometry(
																		ANALYSIS_THRESH	= 3,
																		ASSOC_NAME		= filename_stars.replace('.cat', '_loc_cleaned.cat'),
																		ASSOC_PARAMS	= "1,2",
																		ASSOC_RADIUS	= args.tol/fits_tools.pix2arcsec(args.fits),
																		DETECT_THRESH	= 3,
																		FITS			= args.fits,
																		FLAG			= 'ref_star',
																		GAIN			= args.gain,
																		LOGGER			= logger,
																		LOGLEVEL		= args.sex_loglevel,
																		PATH			= args.outdir,
																		PHOT_APERTURES	= 2*np.array(args.ap_diam)*FWHM_median,
																		REF_FILE		= args.ref_image)

msg									= 'Compute zeropoint'
print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
logger.info(msg)

# print(matched_standard, phot_stars, 1000,
# 									args.fits,
# 									logger,
# 									args.outdir,
# 									args.tol)

summary_zeropoint					= phot_routines.zeropoint(matched_standard, phot_stars,
																		NITER			= 1000,
																		FITS			= args.fits,
																		LOGGER			= logger,
																		PATH			= args.outdir,
																		TOLERANCE		= args.tol)
summary_zeropoint['r(FWHM)']		= np.nan

msg									= 'Step 4: Aperture photometry'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

msg									= 'Run sextractor'
print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
logger.info(msg)

phot_all							= phot_routines.sextractor_photometry(
																		ANALYSIS_THRESH	= 1,
																		BACK_SIZE		= args.back_size,
																		BACK_FILTERSIZE	= args.back_filtersize,
																		DEBLEND_NTHRESH	= args.deblend_nthresh,
																		DEBLEND_MINCONT	= args.deblend_mincont, 
																		DETECT_THRESH	= 1,
																		FLAG			= 'all',
																		FITS			= args.fits,
																		GAIN			= args.gain,
																		LOGGER			= logger,
																		LOGLEVEL		= args.sex_loglevel,
																		PATH			= args.outdir,
																		PHOT_APERTURES	= 2*np.array(args.ap_diam)*FWHM_median,
																		REF_FILE		= args.ref_image)

phot_science						= phot_routines.sextractor_photometry(
																		ANALYSIS_THRESH	= args.ana_thresh,
																		ASSOC_NAME		= filename_science_xy,
																		ASSOC_PARAMS	= "1,2",
																		ASSOC_RADIUS	= args.host_offset/fits_tools.pix2arcsec(args.fits),
																		BACK_SIZE		= args.back_size,
																		BACK_FILTERSIZE	= args.back_filtersize,
																		DEBLEND_NTHRESH	= args.deblend_nthresh,
																		DEBLEND_MINCONT	= args.deblend_mincont, 
																		DETECT_THRESH	= args.det_thresh,
																		FITS			= args.fits,
																		FLAG			= 'science',
																		GAIN			= args.gain,
																		LOGGER			= logger,
																		LOGLEVEL		= args.sex_loglevel,
																		PATH			= args.outdir,
																		PHOT_APERTURES	= 2*np.array(args.ap_diam)*FWHM_median,
																		REF_FILE		= args.ref_image)

msg									= 'Post-process Sextractor output'
print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
logger.info(msg)

# Postprocess Sextractor output
# 1) Set negavitve fluxes to NaN
# 2) Recompute magnitude errors (becomes asymmetric). Done to be consistent with forced photometry.

msg									= '1) Set negavitve fluxes to NaN'
logger.info(msg)
msg									= '2) Recompute magnitude errors (becomes asymmetric). Done to be consistent with forced photometry.'
logger.info(msg)

phot_all							= phot_routines.sextractor_postprocess(phot_all)
phot_science						= phot_routines.sextractor_postprocess(phot_science, PRINTHELP=False)

if len(phot_science)				== 0:

	msg								= 'No object detected at RA, DEC = {ra}, {dec} within {radius:.2f} arcsec.'.format(ra=object_properties['RA_HMS'][0], dec=object_properties['DEC_DMS'][0], radius=args.host_offset)
	print(bcolors.FAIL + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
	logger.info(msg)

	msg								= 'Forced photometry at the object coordinates.'
	print(bcolors.FAIL + bcolors.BOLD + msg + bcolors.ENDC)
	logger.info(msg)

	# Open FITS file

	hdu								= fits.open(args.fits)
	hdu_data						= hdu[0].data
	hdu_header						= hdu[0].header

	# Do forced photometry for all apertures
	# Instrumental magnitudes

	forced_phot						= phot_routines.aperture_photometry(hdu_data,
																		[object_properties['X_EXP'][0], object_properties['Y_EXP'][0]],
																		apertures,
																		1.2 * apertures,
																		2.0 * apertures,
																		1,
																		GAIN		= 10,
																		ZEROPOINT	= np.array(summary_zeropoint['ZP'][2:]),
																		FA			= 1)

	# Calibrated magnitudes

	for i in range(len(apertures)):

		if i == 0:
			forced_phot['MAG_APER_' + str(i)] 		+= 0#summary_zeropoint['ZP'][summary_zeropoint['METHOD'] == 'MAG_APER'][0]
			forced_phot['MAGERRP_APER_' + str(i)] 	= np.sqrt(forced_phot['MAGERRP_APER_' + str(i)][0]**2 + summary_zeropoint['ZP_ERRP'][summary_zeropoint['METHOD'] == 'MAG_APER'][0]**2)
			forced_phot['MAGERRM_APER_' + str(i)] 	= np.sqrt(forced_phot['MAGERRM_APER_' + str(i)][0]**2 + summary_zeropoint['ZP_ERRM'][summary_zeropoint['METHOD'] == 'MAG_APER'][0]**2)

		else:

			forced_phot['MAG_APER_' + str(i)]		+= 0#summary_zeropoint['ZP'][summary_zeropoint['METHOD'] == 'MAG_APER_'+str(i)][0]
			forced_phot['MAGERRP_APER_' + str(i)] 	= np.sqrt(forced_phot['MAGERRP_APER_' + str(i)][0]**2 + summary_zeropoint['ZP_ERRP'][summary_zeropoint['METHOD'] == 'MAG_APER_'+str(i)][0]**2)
			forced_phot['MAGERRM_APER_' + str(i)] 	= np.sqrt(forced_phot['MAGERRM_APER_' + str(i)][0]**2 + summary_zeropoint['ZP_ERRM'][summary_zeropoint['METHOD'] == 'MAG_APER_'+str(i)][0]**2)

else:
	forced_phot						= None

msg									= 'Convert instrumental to calibrated magnitudes'
print(bcolors.OKGREEN + bcolors.BOLD + '\n' + msg + bcolors.ENDC)
logger.info(msg)

for key in [x for x in phot_science.keys() if 'MAG_' in x]:

	# Calibrate magnitudes

	if forced_phot == None:

		phot_science[key]									= phot_science[key] 	+ summary_zeropoint['ZP'][summary_zeropoint['METHOD'] == key]
		phot_science[key.replace('MAG_', 'MAGERRP_')]		= np.sqrt(phot_science[key.replace('MAG_', 'MAGERRP_')]**2 + summary_zeropoint['ZP_ERRP'][summary_zeropoint['METHOD'] == key]**2)
		phot_science[key.replace('MAG_', 'MAGERRM_')]		= np.sqrt(phot_science[key.replace('MAG_', 'MAGERRM_')]**2 + summary_zeropoint['ZP_ERRM'][summary_zeropoint['METHOD'] == key]**2)

		phot_science[key].format 							= '7.3f'
		phot_science[key.replace('MAG_', 'MAGERRM_')].format= '7.3f'
		phot_science[key.replace('MAG_', 'MAGERRP_')].format= '7.3f'

	phot_all[key]										= phot_all[key] 		+ summary_zeropoint['ZP'][summary_zeropoint['METHOD'] == key]
	phot_all[key.replace('MAG_', 'MAGERRP_')]			= np.sqrt(phot_all[key.replace('MAG_', 'MAGERRP_')]**2 + summary_zeropoint['ZP_ERRP'][summary_zeropoint['METHOD'] == key]**2)
	phot_all[key.replace('MAG_', 'MAGERRM_')]			= np.sqrt(phot_all[key.replace('MAG_', 'MAGERRM_')]**2 + summary_zeropoint['ZP_ERRM'][summary_zeropoint['METHOD'] == key]**2)

	phot_all[key].format 								= '7.3f'
	phot_all[key.replace('MAG_', 'MAGERRM_')].format 	= '7.3f'
	phot_all[key.replace('MAG_', 'MAGERRP_')].format 	= '7.3f'

msg									= 'Step 5: Summaries'
print(bcolors.HEADER + bcolors.BOLD + '\n{}\n'.format(msg) + bcolors.ENDC)
logger.info(msg)

# Zeropoint summary table

j									= 0

for i in range(len(summary_zeropoint['METHOD'])):
	if 'APER' in summary_zeropoint['METHOD'][i]:
		summary_zeropoint['r(FWHM)'][i]		= args.ap_diam[j]
		j							+= 1

summary_zeropoint['diam(px)']				= 2 * summary_zeropoint['r(FWHM)'] * FWHM_median

summary_zeropoint['MAG_3UL_GLOB']			= np.nan*len(summary_zeropoint)
summary_zeropoint['AP_cor']					= summary_zeropoint['ZP'][-1] - summary_zeropoint['ZP'] 

for key in summary_zeropoint['METHOD']:

	temp_measurements 				= phot_all[key]# + zeropoint['ZP'][zeropoint['METHOD'] == key]
	mask_bad						= np.where( (phot_all['MAGERR_'+key.split('_')[1]] >= 0.3) & (phot_all['MAGERR_'+key.split('_')[1]] <= 0.34))[0]

	if len(temp_measurements[mask_bad]) > 0:
		summary_zeropoint['MAG_3UL_GLOB'][summary_zeropoint['METHOD'] == key] 	= np.min(temp_measurements[mask_bad])
	else:
		summary_zeropoint['MAG_3UL_GLOB'][summary_zeropoint['METHOD'] == key] 	= np.nan

for key in ['ZP', 'ZP_ERRP', 'ZP_ERRM', 'MAG_3UL_GLOB', 'AP_cor', 'diam(px)']:
	summary_zeropoint[key].format 			= '.03f'

# Print summary

print('\n')
print(bcolors.OKGREEN + "\nZeropoint \n" + bcolors.ENDC)

logger.info('Zeropoint summary')
logger.info(summary_zeropoint)
print(summary_zeropoint)

# Build output catalogue

print(bcolors.OKGREEN + "\nScience \n" + bcolors.ENDC)

summary_science 				= phot_routines.make_scicat (args.fits, object_properties, phot_science, forced_phot, summary_zeropoint, args.host_offset, logger)
logger.info('Photometry summary')
logger.info(summary_science)

# Print summary
summary_science.pprint(max_lines=-1)

# Make cutouts

msg								= 'Step 6: Make poststamp'
print(bcolors.HEADER + bcolors.BOLD + "\n{}\n".format(msg) + bcolors.ENDC)
logger.info(msg)

if 'DISTANCE (arcsec)' in phot_science.keys() and phot_science['DISTANCE (arcsec)'][0] <= args.host_offset:
	phot_routines.make_poststamp(args.fits, [x_exp, y_exp], [summary_science['VALUE'][summary_science['PROPERTY'] == 'XWIN_IMAGE_OBS'][0], summary_science['VALUE'][summary_science['PROPERTY'] == 'YWIN_IMAGE_OBS'][0]], PATH=args.outdir)
else:
	phot_routines.make_poststamp(args.fits, [x_exp, y_exp], [0, 0], PATH=args.outdir)

# Save results to file

msg								= 'Step 7: Save to file'
print(bcolors.HEADER + bcolors.BOLD + "\n{}\n".format(msg) + bcolors.ENDC)
logger.info(msg)

ascii.write(summary_science, 	args.outdir + args.fits.replace('.fits', '_phot.log'), 							overwrite=True)
ascii.write(summary_zeropoint, 	args.outdir + args.fits.replace('.fits', '_zp.log'),   							overwrite=True)

ascii.write(phot_all,			args.outdir + args.fits.replace('.fits', '_all_abs_cal.phot'),     				overwrite=True)

# if forced_phot 					!= None:
# 	ascii.write(forced_phot,	args.outdir + args.fits.replace('.fits', '_science_forcedphot_abs_cal.phot'),	overwrite=True)

# Remove temporary files

if not args.keeptemp:
	msg 						= 'Step 8: Remove all temps'
	print(bcolors.HEADER + bcolors.BOLD + "\n{}\n".format(msg) + bcolors.ENDC)
	logger.info(msg)
	os.system('rm check_' + args.fits)
	os.system('rm {}*refreg*'.format(args.outdir + args.fits.replace(args.fits.split('.')[-1], '')[:-1]))
	os.system('rm {}*refcat*'.format(args.outdir + args.fits.replace(args.fits.split('.')[-1], '')[:-1]))
	os.system('rm {}*_loc_*'.format(args.outdir + args.fits.replace(args.fits.split('.')[-1], '')[:-1]))
	os.system('rm ' + args.outdir + args.fits.replace('.fits', '_ref_star.*'))
	os.system('rm ' + args.outdir + args.fits.replace('.fits', '_xy.cat'))
	os.system('rm ' + args.outdir + args.fits.replace('.fits', '_all.phot'))	
	os.system('rm ' + args.outdir + '*science*')

else:
	msg 							= 'Step 8: Keep all temps'
	print(bcolors.HEADER + bcolors.BOLD + "\n{}\n".format(msg) + bcolors.ENDC)
	logger.info(msg)

# Show plots

if not(args.auto) and not(float(args.mag_stdbright) != 0.) and not(float(args.mag_stdfaint) != 0.):
	plt.show()

plt.close()
