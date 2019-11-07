#!/usr/bin/env python

__version__ 						= "2019-11-04"
__author__ 							= "Steve Schulze (steve.schulze@weizmann.ac.il)"

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
import	urllib

##########################

parser				 = argparse.ArgumentParser(description='Retrieve photometric catalogues. 2MASS, PS1, SDSS and SkyMapper are the input catalogues. Bessel catalogue is generated through colour equations. PS1, SDSS and SkyMapper cats are in the AB system, whereas Bessel and 2MASS cats are in the Vega system.')

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
coordinates			= coord.SkyCoord(ra_dd, dec_dd, unit=u.deg)

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
		result 			= v.query_region(coordinates, radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['SDSS']['CATID'])
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
# Based on crossmatching Gaia and PS1
# Objects are considered to be stars if either the parallax or one of the proper motion measurements has a significance of >3 sigma.

if args.type == 'all' or args.type == 'optical':

	print(bcolors.HEADER + 'PS1 catalogues' + bcolors.ENDC)

	try:

		# Gaia

		v					= Vizier(columns = ['all'], row_limit = -1)
		result_gaia			= v.query_region(coordinates, radius = args.radius * u.arcmin, catalog=cat_tools.catalog_prop['GAIA']['CATID'])
		result_gaia			= result_gaia[cat_tools.catalog_prop['GAIA']['CATID_OUT']]
		result_gaia			= result_gaia[cat_tools.catalog_prop['GAIA']['KEYWORDS']]
		result_gaia			= result_gaia[(result_gaia['Plx'] > 0) & (result_gaia['e_pmRA'] > 0) & (result_gaia['e_pmDE'] > 0)]
		result_gaia			= result_gaia[(result_gaia['Plx'] / result_gaia['e_Plx'] > 3) | (abs(result_gaia['pmRA'] / result_gaia['e_pmRA']) > 3) | (abs(result_gaia['pmDE'] / result_gaia['e_pmDE']) > 3)]

		# PS1

		v					= Vizier(columns = ['all'], row_limit = -1)
		result_panstarrs	= v.query_region(coordinates, radius = args.radius * u.arcmin, catalog=cat_tools.catalog_prop['PanSTARRS']['CATID'])
		result_panstarrs	= result_panstarrs[cat_tools.catalog_prop['PanSTARRS']['CATID_OUT']]
		result_panstarrs	= result_panstarrs[cat_tools.catalog_prop['PanSTARRS']['KEYWORDS']]

		result_panstarrs	= result_panstarrs[
							(result_panstarrs['o_gmag'] >= 0.85) & (result_panstarrs['o_rmag'] >= 0.85) &
							(result_panstarrs['o_imag'] >= 0.85) & (result_panstarrs['o_zmag'] >= 0.85) &
							(result_panstarrs['o_ymag'] >= 0.85) &
							(np.round(result_panstarrs['gPSFf'], 0) >= 1) & (np.round(result_panstarrs['rPSFf'], 0) >= 1) &
							(np.round(result_panstarrs['iPSFf'], 0) >= 1) & (np.round(result_panstarrs['zPSFf'], 0) >= 1) &
							(np.round(result_panstarrs['yPSFf'], 0) >= 1)
							]

		result_panstarrs_keys=['RAJ2000', 'DEJ2000',
								'gmag', 'e_gmag',
								'rmag', 'e_rmag',
								'imag', 'e_imag',
								'zmag', 'e_zmag',
								'ymag', 'e_ymag']

		result_panstarrs	= result_panstarrs[result_panstarrs_keys]

		# Cross-match catalogues

		ref_stars 						= np.array(np.asarray(result_gaia).tolist())#.view((  float, len(result_gaia.dtype.names)))#.reshape((-1,len(result_gaia)))#, len(result_gaia)))
		ref_ps1							= np.array(np.asarray(result_panstarrs).tolist())

		matched_standard				= cat_tools.wrapper_crossmatch(ref_stars, ref_ps1, 0.5)
		matched_standard				= table.Table(matched_standard, names=result_panstarrs_keys+['GAIA_'+x for x in cat_tools.catalog_prop['GAIA']['KEYWORDS']]+['DIST'])
		matched_standard['DIST']		*= 3600
		matched_standard				= matched_standard[result_panstarrs_keys]

		for filter in ['g', 'r', 'i', 'z', 'y']:
			matched_standard.rename_column(filter + 'mag', filter + '_PS1')
			matched_standard.rename_column('e_' + filter + 'mag', filter + '_PS1_ERR')

		# Formatting

		for key in [x for x in matched_standard.keys() if 'PS1' in x]:
			matched_standard[key].format= '.4f'

		result							= matched_standard

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

# Get SkyMapper

if args.type == 'all' or args.type == 'optical':

	# Consistency check

	if coordinates.dec.deg <0:

		print(bcolors.HEADER + 'SkyMapper catalogue' + bcolors.ENDC)

		try:

			# Retrieve catalogue

			get_skymapper					= urllib.request.urlopen("http://skymapper.anu.edu.au/sm-cone/public/query?RA={ra:.3f}&DEC={dec:.3f}&SR={radius:.1f}&RESPONSEFORMAT=CSV".format(
																	ra		= coordinates.ra.deg,
																	dec		= coordinates.dec.deg,
																	radius	= args.radius/60.
																	)).read()

			# Write to file

			file							= open(args.outdir + 'skymapper.csv', 'w')
			file.write(get_skymapper.decode())
			file.close()

			cat_skymapper					= ascii.read(args.outdir + 'skymapper.csv', format='csv')
			os.system('rm ' + args.outdir + 'skymapper.csv')

			# Filter

			cat_skymapper					= cat_skymapper[(cat_skymapper['class_star'] > 0.95) & (cat_skymapper['flags'] == 0)]
			cat_skymapper					= cat_skymapper[(cat_skymapper['g_psf'] > 0) & (cat_skymapper['r_psf'] > 0) & (cat_skymapper['i_psf'] > 0) & (cat_skymapper['z_psf'] > 0)]

			keywords						= ['raj2000','dej2000','u_psf','e_u_psf','v_psf','e_v_psf','g_psf','e_g_psf','r_psf','e_r_psf','i_psf','e_i_psf','z_psf','e_z_psf']

			cat_skymapper					= cat_skymapper[keywords]

			cat_skymapper.rename_column('raj2000', 'RAJ2000')
			cat_skymapper.rename_column('dej2000', 'DEJ2000')

			for filter in ['u', 'g', 'r', 'i', 'z']:
				cat_skymapper.rename_column(filter + '_psf', filter + '_SDSS')
				cat_skymapper.rename_column('e_' + filter + '_psf', filter + '_SDSS_ERR')

			# Write SDSS catalogues

			for key in [x for x in cat_skymapper.keys() if 'SDSS' in x]:
				cat_skymapper[key].format= '.4f'

			for filter in ['u', 'g', 'r', 'i', 'z']:
				mask_good	= np.where((cat_skymapper[filter + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (cat_skymapper[filter + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]
				filename	= args.outdir + 'SkyMapper_SDSS_' + filter + '.cat'
				ascii.write(cat_skymapper[['RAJ2000', 'DEJ2000', filter+'_SDSS', filter+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

			# Convert to Bessel system

			cat_skymapper 			= cat_tools.SDSS_to_Bessel(cat_skymapper)

			# Formatting

			for key in [x for x in cat_skymapper.keys() if 'BESSEL' in x]:
				cat_skymapper[key].format= '.4f'

			# Write Bessel catalogues

			for filter in ['B', 'V', 'R', 'I']:
				mask_good	= np.where((cat_skymapper[filter + '_BESSEL_ERR'] > 0.) & (cat_skymapper[filter + '_BESSEL_ERR'] < 0.3))[0]
				filename	= args.outdir + 'SkyMapper_BESSEL_' + filter + '.cat'
				ascii.write(cat_skymapper[['RAJ2000', 'DEJ2000', filter+'_BESSEL', filter+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

		except:
			print(bcolors.WARNING + 'SkyMapper query failed' + bcolors.ENDC)
			pass

	else:
		cmd							= 'Declination is above 0 deg. No entries in SkyMapper catalogue.'
		print(bcolors.FAIL + cmd + bcolors.ENDC)

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
