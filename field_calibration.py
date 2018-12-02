#!/usr/bin/env pythonw

import	argparse
from 	astroquery.vizier import Vizier
import 	astropy.units as u
import 	astropy.coordinates as coord
from 	astropy.io import ascii
from 	astropy import table
import	cat_tools
import	fits_tools
from	misc import bcolors
import	numpy as np
import	os
import	phot_routines

##########################

parser				 = argparse.ArgumentParser(description='Retrieve photometric catalogues. 2MASS, PS1 and SDSS are the input sources. Bessel catalogue is generated through colour equations. PS1 and SDSS cats are in the AB system, whereas Bessel and 2MASS cats are in the Vega system.')

parser.add_argument('--ra',				type	= str,
										help	= 'RA(J2000) of the Object (HMS and DD system allowed). Note: if declination is negative, write " -12:20:20.2"). (Required)',
										required= True)

parser.add_argument('--dec', 			type	= str,
										help	= 'Dec (J2000) of the Object (HMS and DD system allowed). (Required)',
										required= True)

parser.add_argument('--radius', 		type	= float,
										help	= 'only search in indexes within \'radius\' of the field center (unit: arcmin; default: 20)',
										default = 10)

parser.add_argument('--outdir',			type	= str,
										help	= 'Output directory (default: results/)',
										default	= 'results/')

parser.add_argument('--type',			type	= str,
										help	= 'Generate catalogues for optical/nir/all (default: all)',
										default	= 'all')

args				= parser.parse_args()

if args.outdir[-1] != '/':
	args.outdir		+= '/'

if not os.path.exists(args.outdir):
	os.mkdir(args.outdir)

# Convert coordinates

ra_dd, dec_dd		= fits_tools.convert_hms_dd(args.ra, args.dec)

# Set up output dir

if not(os.path.isdir(args.outdir)):
	os.mkdir(args.outdir)

# Get SDSS catalogue
if args.type == 'all' or args.type == 'optical':

	print(bcolors.OKGREEN + 'Generate optical catalogues' + bcolors.ENDC)


	print(bcolors.HEADER + 'SDSS catalogues' + bcolors.ENDC)

	try:
		# Retrieve catalogue

		v 				= Vizier(row_limit=10000)
		result 			= v.query_region(coord.SkyCoord(ra_dd, dec_dd, unit=u.deg), radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['SDSS']['CATID'])
		result 			= result[cat_tools.catalog_prop['SDSS']['CATID_OUT']]
		result			= result[cat_tools.catalog_prop['SDSS']['KEYWORDS']]

		result			= result[(result['class'] == 6) & (result['mode'] == 1)]

		# Relabeling

		for filter in ['u', 'g', 'r', 'i', 'z']:
			result.rename_column(filter + 'mag', filter + '_SDSS')
			result.rename_column('e_' + filter + 'mag', filter + '_SDSS_ERR')

		# Formatting

		for key in [x for x in result.keys() if 'SDSS' in x]:
			result[key].format= '.4f'

		# Write SDSS catalogues

		for filter in ['u', 'g', 'r', 'i', 'z']:			

			mask_good	= np.where((result[filter + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (result[filter + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]

			filename	= args.outdir + 'SDSS_SDSS_' + filter + '.cat'
			ascii.write(result[['RA_ICRS', 'DE_ICRS', filter+'_SDSS', filter+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

		# Convert to Bessel system

		result 			= cat_tools.SDSS_to_Bessel(result)

		# Formatting

		for key in [x for x in result.keys() if 'BESSEL' in x]:
			result[key].format= '.4f'

		# Write Bessel catalogues

		for filter in ['B', 'V', 'R', 'I']:

			mask_good	= np.where((result[filter + '_BESSEL_ERR'] > 0.) & (result[filter + '_BESSEL_ERR'] < 0.3))[0]

			filename	= args.outdir + 'SDSS_BESSEL_' + filter + '.cat'
			ascii.write(result[['RA_ICRS', 'DE_ICRS', filter+'_BESSEL', filter+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

	except:
		print(bcolors.FAIL + 'Field not covered by SDSS or VizieR query failed.' + bcolors.ENDC)
		pass

# Get PS1 catalogue

if args.type == 'all' or args.type == 'optical':

	print(bcolors.HEADER + 'PS1 catalogues' + bcolors.ENDC)

	try:

		v 				= Vizier(row_limit=10000, columns=['all'])
		result 			= v.query_region(coord.SkyCoord(ra_dd, dec_dd, unit=u.deg), radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['PanSTARRS']['CATID'])
		result 			= result[cat_tools.catalog_prop['PanSTARRS']['CATID_OUT']]
		result			= result[cat_tools.catalog_prop['PanSTARRS']['KEYWORDS']]

		# print(result.keys())

		result 			= result[
						(result['o_gmag'] >= 0.85) & (result['o_rmag'] >= 0.85) &
						(result['o_imag'] >= 0.85) & (result['o_zmag'] >= 0.85) &
						(result['o_ymag'] >= 0.85) &
	 					(np.round(result['gPSFf'], 0) >= 1) & (np.round(result['rPSFf'], 0) >= 1) &
	 					(np.round(result['iPSFf'], 0) >= 1) & (np.round(result['zPSFf'], 0) >= 1) &
	 					(np.round(result['yPSFf'], 0) >= 1) &
						(abs(result['imag'] - result['iKmag']) < 0.05) &
						(result['imag'] > 14) & (result['imag'] < 21)
						]

		for filter in ['g', 'r', 'i', 'z', 'y']:
			result.rename_column(filter + 'mag', filter + '_PS1')
			result.rename_column('e_' + filter + 'mag', filter + '_PS1_ERR')

		# Formatting

		for key in [x for x in result.keys() if 'PS1' in x]:
			result[key].format= '.4f'

		# Write PS1 catalogues

		for filter in ['g', 'r', 'i', 'z', 'y']:

			mask_good	= np.where((result[filter + '_PS1_ERR'] > cat_tools.catalog_prop['PanSTARRS']['SIGMA_HIGH']) & (result[filter + '_PS1_ERR'] < cat_tools.catalog_prop['PanSTARRS']['SIGMA_LOW']))[0]

			filename	= args.outdir + 'PS1_PS1_' + filter + '.cat'
			ascii.write(result[['RAJ2000', 'DEJ2000', filter+'_PS1', filter+'_PS1_ERR']][mask_good], filename, overwrite=True, format='no_header')

		# Write SDSS catalogues

		result 			= cat_tools.PS1_to_SDSS(result)

		for key in [x for x in result.keys() if 'SDSS' in x]:
			result[key].format= '.4f'

		for filter in ['u', 'g', 'r', 'i', 'z']:

			mask_good	= np.where((result[filter + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (result[filter + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]

			filename	= args.outdir + 'PS1_SDSS_' + filter + '.cat'
			ascii.write(result[['RAJ2000', 'DEJ2000', filter+'_SDSS', filter+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

		# Convert to Bessel system

		result 			= cat_tools.SDSS_to_Bessel(result)

		# Formatting

		for key in [x for x in result.keys() if 'BESSEL' in x]:
			result[key].format= '.4f'

		# Write Bessel catalogues

		for filter in ['B', 'V', 'R', 'I']:

			mask_good	= np.where((result[filter + '_BESSEL_ERR'] > 0.) & (result[filter + '_BESSEL_ERR'] < 0.3))[0]

			filename	= args.outdir + 'PS1_BESSEL_' + filter + '.cat'
			ascii.write(result[['RAJ2000', 'DEJ2000', filter+'_BESSEL', filter+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

	except:
		print(bcolors.WARNING + 'Field not covered by PS1' + bcolors.ENDC)
		pass


# Get 2MASS catalogue

if args.type == 'all' or args.type == 'nir':

	print(bcolors.OKGREEN + 'Generate NIR catalogues' + bcolors.ENDC)
	print(bcolors.HEADER + '2MASS catalogues' + bcolors.ENDC)

	try:
		v 		= Vizier(row_limit=10000)

		result 	= v.query_region(coord.SkyCoord(ra_dd, dec_dd, unit=u.deg), radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['2MASS']['CATID'])#
		result	= result[cat_tools.catalog_prop['2MASS']['CATID_OUT']]
		result	= result[cat_tools.catalog_prop['2MASS']['KEYWORDS']]

		for key in [x for x in result.keys() if 'mag' in x]:

			result[key].format= '.4f'

		for filter in ['J', 'H', 'K']:			

			mask_good	= np.where((result['e_'+filter+'mag'] > cat_tools.catalog_prop['2MASS']['SIGMA_HIGH']) & (result['e_'+filter+'mag'] < cat_tools.catalog_prop['2MASS']['SIGMA_LOW']))[0]

			filename	= args.outdir + '2MASS_' + filter + '.cat'
			ascii.write(result[['RAJ2000', 'DEJ2000', filter+'mag', 'e_'+filter+'mag']][mask_good], filename, overwrite=True, format='no_header')

	except:
		print(bcolors.FAIL + 'VizieR query failed.' + bcolors.ENDC)
		pass
