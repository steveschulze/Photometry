#!/usr/bin/env python

__version__                         = "2021-07-20"
__author__                          = "Steve Schulze (steve.schulze@weizmann.ac.il)"

import    argparse
from      astroquery.vizier import Vizier
import    astropy.units as u
import    astropy.coordinates as coord
from      astropy.io import ascii
from      astropy import table
import    cat_tools
import    copy
import    fits_tools
from      misc import bcolors
import    numpy as np
import    os
import    routines_noir
import    sys
import    urllib

def get_parser():

    parser                 = argparse.ArgumentParser(description='Retrieve photometric catalogues. 2MASS, PS1, SDSS and SkyMapper are the input catalogues. Bessel catalogue is generated through colour equations. PS1, SDSS and SkyMapper cats are in the AB system, whereas Bessel and 2MASS cats are in the Vega system.')

    parser.add_argument('--ra',             
                        type    = str,
                        help    = 'RA(J2000) of the Object (HMS and DD system allowed). Note: if declination is negative, write " -12:20:20.2")',
                        required= True)

    parser.add_argument('--dec',            
                        type    = str,
                        help    = 'Dec (J2000) of the Object (HMS and DD system allowed)',
                        required= True)

    parser.add_argument('--radius',         
                        type    = float,
                        help    = 'only search in indexes within \'radius\' of the field center (unit: arcmin; default: 20)',
                        default = 10)

    parser.add_argument('--outdir',         
                        type    = str,
                        help    = 'Output directory (default: results/)',
                        default = 'results/')

    parser.add_argument('--type',           
                        type    = str,
                        help    = 'Generate catalogues for optical/nir/all (default: all)',
                        default = 'all')

    return parser

def main(args):

    parser             = get_parser()
    args               = parser.parse_args(args)

    if args.outdir[-1] != '/':
        args.outdir        += '/'

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)

    # Convert coordinates

    ra_dd, dec_dd        = fits_tools.convert_hms_dd(args.ra, args.dec.replace(' ', '').replace('"', ''))
    coordinates          = coord.SkyCoord(ra_dd, dec_dd, unit=u.deg)

    # Set up output dir

    if not(os.path.isdir(args.outdir)):
        os.mkdir(args.outdir)

    # Get SDSS catalogue

    print(bcolors.OKBLUE + '\nSDSS catalogues' + bcolors.ENDC)

    try:
        # Retrieve catalogue

        v                 = Vizier(row_limit=100000)
        result             = v.query_region(coordinates, radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['SDSS']['CATID'])
        result             = result[cat_tools.catalog_prop['SDSS']['CATID_OUT']]
        result            = result[cat_tools.catalog_prop['SDSS']['KEYWORDS']]

        result            = result[(result['class'] == 6) & (result['mode'] == 1)]

        # Relabeling

        for f in ['u', 'g', 'r', 'i', 'z']:
            result.rename_column(f + 'mag', f + '_SDSS')
            result.rename_column('e_' + f + 'mag', f + '_SDSS_ERR')

        # Formatting

        for key in [x for x in result.keys() if 'SDSS' in x]:
            result[key].format= '.4f'

        result_sdss                        = copy.deepcopy(result)

        # Write SDSS catalogues

        for f in ['u', 'g', 'r', 'i', 'z']:            

            mask_good    = np.where((result[f + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (result[f + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]

            filename    = args.outdir + 'SDSS_SDSS_' + f + '.cat'
            ascii.write(result[['RA_ICRS', 'DE_ICRS', f+'_SDSS', f+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to Bessel system

        result             = cat_tools.SDSS_to_Bessel(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'BESSEL' in x]:
            result[key].format= '.4f'

        # Write Bessel catalogues

        for f in ['B', 'V', 'R', 'I']:

            mask_good    = np.where((result[f + '_BESSEL_ERR'] > 0.) & (result[f + '_BESSEL_ERR'] < 0.3))[0]

            filename    = args.outdir + 'SDSS_BESSEL_' + f + '.cat'
            ascii.write(result[['RA_ICRS', 'DE_ICRS', f+'_BESSEL', f+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to GROND system

        print(bcolors.WARNING + 'Convert SDSS -> GROND' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_GROND(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'GROND' in x]:
            result[key].format= '.4f'

        # Write GROND catalogues

        for f in ['g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_GROND_ERR'] > 0.) & (result[f + '_GROND_ERR'] < 0.3))[0]
            filename    = args.outdir + 'SDSS_GROND_' + f + '.cat'
            ascii.write(result[['RA_ICRS', 'DE_ICRS', f+'_GROND', f+'_GROND_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to DES system
        
        print(bcolors.WARNING + 'Convert SDSS -> GROND' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_DES(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'DES' in x]:
            result[key].format= '.4f'

        # Write DES catalogues

        for f in ['g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_DES_ERR'] > 0.) & (result[f + '_DES_ERR'] < 0.3))[0]
            filename    = args.outdir + 'SDSS_DES_' + f + '.cat'
            ascii.write(result[['RA_ICRS', 'DE_ICRS', f + '_DES', f + '_DES_ERR']][mask_good], filename, overwrite=True, format='no_header')

    except:
        print(bcolors.FAIL + 'Field not covered by SDSS or VizieR query failed.' + bcolors.ENDC)
        pass

    # Get PS1 catalogue
    # Based on crossmatching Gaia and PS1
    # Objects are considered to be stars if either the parallax or one of the proper motion measurements has a significance of >3 sigma.

    print(bcolors.OKBLUE + '\nPS1 catalogues' + bcolors.ENDC)

    try:
    #if True:

        # Gaia

        v                    = Vizier(columns = ['all'], row_limit = -1)
        result_gaia            = v.query_region(coordinates, radius = args.radius * u.arcmin, catalog=cat_tools.catalog_prop['GAIA']['CATID'])
        result_gaia            = result_gaia[cat_tools.catalog_prop['GAIA']['CATID_OUT']]
        result_gaia            = result_gaia[cat_tools.catalog_prop['GAIA']['KEYWORDS']]
        result_gaia            = result_gaia[(result_gaia['Plx'] > 0) & (result_gaia['e_pmRA'] > 0) & (result_gaia['e_pmDE'] > 0)]
        result_gaia            = result_gaia[(result_gaia['Plx'] / result_gaia['e_Plx'] > 3) | (abs(result_gaia['pmRA'] / result_gaia['e_pmRA']) > 3) | (abs(result_gaia['pmDE'] / result_gaia['e_pmDE']) > 3)]

        # PS1

        v                    = Vizier(columns = ['all'], row_limit = -1)
        result_panstarrs    = v.query_region(coordinates, radius = args.radius * u.arcmin, catalog=cat_tools.catalog_prop['PanSTARRS']['CATID'])
        result_panstarrs    = result_panstarrs[cat_tools.catalog_prop['PanSTARRS']['CATID_OUT']]
        result_panstarrs    = result_panstarrs[cat_tools.catalog_prop['PanSTARRS']['KEYWORDS']]

        result_panstarrs_keys=['RAJ2000', 'DEJ2000',
                                'gmag', 'e_gmag',
                                'rmag', 'e_rmag',
                                'imag', 'e_imag',
                                'zmag', 'e_zmag',
                                'ymag', 'e_ymag']

        result_panstarrs    = result_panstarrs[result_panstarrs_keys]

        # Cross-match catalogues

        ref_stars                         = np.array(np.asarray(result_gaia).tolist())#.view((  float, len(result_gaia.dtype.names)))#.reshape((-1,len(result_gaia)))#, len(result_gaia)))
        ref_ps1                            = np.array(np.asarray(result_panstarrs).tolist())

        matched_standard                = cat_tools.wrapper_crossmatch(ref_stars, ref_ps1, 0.5)
        matched_standard                = table.Table(matched_standard, names=result_panstarrs_keys+['GAIA_'+x for x in cat_tools.catalog_prop['GAIA']['KEYWORDS']]+['DIST'])
        matched_standard['DIST']        *= 3600
        matched_standard                = matched_standard[result_panstarrs_keys]

        for f in ['g', 'r', 'i', 'z', 'y']:
            matched_standard.rename_column(f + 'mag', f + '_PS1')
            matched_standard.rename_column('e_' + f + 'mag', f + '_PS1_ERR')

        # Formatting

        for key in [x for x in matched_standard.keys() if 'PS1' in x]:
            matched_standard[key].format= '.4f'

        result                            = matched_standard
        result_ps1                        = copy.deepcopy(result)

        # Write PS1 catalogues

        for f in ['g', 'r', 'i', 'z', 'y']:

            mask_good    = np.where((result[f + '_PS1_ERR'] > cat_tools.catalog_prop['PanSTARRS']['SIGMA_HIGH']) & (result[f + '_PS1_ERR'] < cat_tools.catalog_prop['PanSTARRS']['SIGMA_LOW']))[0]

            filename    = args.outdir + 'PS1_PS1_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f + '_PS1', f + '_PS1_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert SDSS catalogues

        print(bcolors.WARNING + 'Convert PS1 -> SDSS' + bcolors.ENDC)

        result             = cat_tools.PS1_to_SDSS(result_ps1)

        for key in [x for x in result.keys() if 'SDSS' in x]:
            result[key].format= '.4f'

        result_sdss        = copy.deepcopy(result)

        for f in ['u', 'g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (result[f + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]

            filename    = args.outdir + 'PS1_SDSS_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_SDSS', f+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to Bessel system

        print(bcolors.WARNING + 'Convert PS1 -> SDSS -> BESSEL' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_Bessel(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'BESSEL' in x]:
            result[key].format= '.4f'

        # Write Bessel catalogues

        for f in ['B', 'V', 'R', 'I']:

            mask_good    = np.where((result[f + '_BESSEL_ERR'] > 0.) & (result[f + '_BESSEL_ERR'] < 0.3))[0]

            filename    = args.outdir + 'PS1_BESSEL_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_BESSEL', f+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to ZTF system

        print(bcolors.WARNING + 'Convert PS1 -> ZTF' + bcolors.ENDC)

        result             = cat_tools.PS1_to_ZTF(result_ps1)

        # Formatting

        for key in [x for x in result.keys() if 'ZTF' in x]:
            result[key].format= '.4f'

        # Write ZTF catalogues

        for f in ['g', 'r']:

            mask_good    = np.where((result[f + '_ZTF_ERR'] > 0.) & (result[f + '_ZTF_ERR'] < 0.3))[0]
            filename    = args.outdir + 'PS1_ZTF_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_ZTF', f+'_ZTF_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to GROND system

        print(bcolors.WARNING + 'Convert PS1 -> SDSS -> GROND' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_GROND(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'GROND' in x]:
            result[key].format= '.4f'

        # Write ZTF catalogues

        for f in ['g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_GROND_ERR'] > 0.) & (result[f + '_GROND_ERR'] < 0.3))[0]
            filename    = args.outdir + 'PS1_GROND_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_GROND', f+'_GROND_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to DES system

        print(bcolors.WARNING + 'Convert PS1 -> SDSS -> DES' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_DES(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'DES' in x]:
            result[key].format= '.4f'

        # Write ZTF catalogues

        for f in ['g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_DES_ERR'] > 0.) & (result[f + '_DES_ERR'] < 0.3))[0]
            filename    = args.outdir + 'PS1_DES_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_DES', f+'_DES_ERR']][mask_good], filename, overwrite=True, format='no_header')

    except:
        print(bcolors.FAIL + 'Field not covered by PS1' + bcolors.ENDC)
        pass

    print(bcolors.OKBLUE + '\nDESI Legacy Imaging Surveys DR9' + bcolors.ENDC)

    try:
    #if True:

        result         = routines_noir.query_noir_datalab_stars(coordinates.ra.value, coordinates.dec.value, args.radius)

        if len(result) > 0:

            for f in ['g', 'r', 'z']:

                temp_result = copy.deepcopy(result)
                temp_result = temp_result['ra', 'dec', 'mag_' + f, 'snr_' + f, 'allmask_' + f]
                temp_result['ra'].format='.6f'
                temp_result['dec'].format='.6f'

                # Select stars with a S/N >= 3 and good data quality (i.e., no saturation etc.)

                mask        = np.where( (temp_result['snr_' + f] >= 3) & (temp_result['allmask_' + f] == 0))[0]
                temp_result = temp_result[mask]

                if len(temp_result) >= 0:

                    # Calculate photometric error
                    temp_result['err_' + f] = 2.5 * np.log10( 1 + 1 / temp_result['snr_' + f])

                    # Prettify

                    for key in ['mag_' + f, 'err_' + f]:
                        temp_result[key].format= '.4f'

                    # Write to file
                    filename    = args.outdir + 'LS_LS_' + f + '.cat'
                    ascii.write(temp_result['ra', 'dec', 'mag_' + f, 'err_' + f], filename, overwrite=True, format='no_header')

            print(bcolors.WARNING + 'Query successful.' + bcolors.ENDC)

            #####

            mask = (result['allmask_g'] == 0) & (result['snr_g'] >= 3) & (result['allmask_r'] == 0) & (result['snr_r'] >= 3) & (result['allmask_z'] == 0) & (result['snr_z'] >= 3)
            result = result[mask]
            
            for f in ['g', 'r', 'z']:
                result['err_' + f] = 2.5 * np.log10( 1 + 1 / result['snr_' + f])

            # Relabeling

            for f in ['g', 'r', 'z']:
                result.rename_column('mag_' + f, f + '_DES')
                result.rename_column('err_' + f, f + '_DES_ERR')

            # Formatting

            for key in [x for x in result.keys() if 'DES' in x]:
                result[key].format= '.4f'

            # Convert SDSS catalogues

            print(bcolors.WARNING + 'Convert LS -> SDSS' + bcolors.ENDC)

            result['RAJ2000'] = result['ra']
            result['DEJ2000'] = result['dec']
            result['i_DES'] = result['z_DES']
            result['i_DES_ERR'] = result['z_DES_ERR']

            result             = cat_tools.DES_to_SDSS(result)

            for key in [x for x in result.keys() if 'SDSS' in x]:
                result[key].format= '.4f'

            result_sdss        = copy.deepcopy(result)

            for f in ['g', 'r']:

                mask_good    = np.where((result[f + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (result[f + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]

                filename    = args.outdir + 'LS_SDSS_' + f + '.cat'
                ascii.write(result[['RAJ2000', 'DEJ2000', f+'_SDSS', f+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

            # Convert to Bessel system

            print(bcolors.WARNING + 'Convert LS -> SDSS -> BESSEL' + bcolors.ENDC)

            result             = cat_tools.SDSS_to_Bessel(result_sdss)

            # Formatting

            for key in [x for x in result.keys() if 'BESSEL' in x]:
                result[key].format= '.4f'

            # Write Bessel catalogues

            for f in ['B', 'V', 'R']:

                mask_good    = np.where((result[f + '_BESSEL_ERR'] > 0.) & (result[f + '_BESSEL_ERR'] < 0.3))[0]

                filename    = args.outdir + 'LS_BESSEL_' + f + '.cat'
                ascii.write(result[['RAJ2000', 'DEJ2000', f+'_BESSEL', f+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

        else:
            print(bcolors.FAIL + 'Field not covered.' + bcolors.ENDC)
    
    except:
        print(bcolors.FAIL + 'NOIR appears to be offline.' + bcolors.ENDC)
        pass
    

    print(bcolors.OKBLUE + '\nDES catalogues' + bcolors.ENDC)

    #if True:
    try:

        v                = Vizier(columns = ['all'], row_limit = -1)
        result             = v.query_region(coordinates, radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['DES']['CATID'])
        result             = result[cat_tools.catalog_prop['DES']['CATID_OUT']]
        result            = result[cat_tools.catalog_prop['DES']['KEYWORDS']]
        result            = result[(result['S_Gg'] > 0.85) & (result['S_Gr'] > 0.85) & (result['S_Gi'] > 0.85)  & (result['gFlag'] <= 4) & (result['rFlag'] <= 4) & (result['iFlag'] <= 4) & (result['zFlag'] <= 4)]

        # Relabeling

        for f in ['g', 'r', 'i', 'z', 'Y']:
            result.rename_column(f + 'mag', f + '_DES')
            result.rename_column('e_' + f + 'mag', f + '_DES_ERR')

        # Formatting

        for key in [x for x in result.keys() if 'DES' in x]:
            result[key].format= '.4f'

        # Write DES catalogues

        for f in ['g', 'r', 'i', 'z', 'Y']:

            mask_good    = np.where((result[f + '_DES_ERR'] > cat_tools.catalog_prop['DES']['SIGMA_HIGH']) & (result[f + '_DES_ERR'] < cat_tools.catalog_prop['DES']['SIGMA_LOW']))[0]

            filename    = args.outdir + 'DES_DES_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_DES', f+'_DES_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert SDSS catalogues

        print(bcolors.WARNING + 'Convert DES -> SDSS' + bcolors.ENDC)

        result             = cat_tools.DES_to_SDSS(result)

        for key in [x for x in result.keys() if 'SDSS' in x]:
            result[key].format= '.4f'

        result_sdss        = copy.deepcopy(result)

        for f in ['g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (result[f + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]

            filename    = args.outdir + 'DES_SDSS_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_SDSS', f+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to Bessel system

        print(bcolors.WARNING + 'Convert DES -> SDSS -> BESSEL' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_Bessel(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'BESSEL' in x]:
            result[key].format= '.4f'

        # Write Bessel catalogues

        for f in ['B', 'V', 'R', 'I']:

            mask_good    = np.where((result[f + '_BESSEL_ERR'] > 0.) & (result[f + '_BESSEL_ERR'] < 0.3))[0]

            filename    = args.outdir + 'DES_BESSEL_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_BESSEL', f+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

        # Convert to GROND system

        print(bcolors.WARNING + 'Convert DES -> SDSS -> GROND' + bcolors.ENDC)

        result             = cat_tools.SDSS_to_GROND(result_sdss)

        # Formatting

        for key in [x for x in result.keys() if 'GROND' in x]:
            result[key].format= '.4f'

        # Write GROND catalogues

        for f in ['g', 'r', 'i', 'z']:

            mask_good    = np.where((result[f + '_GROND_ERR'] > 0.) & (result[f + '_GROND_ERR'] < 0.3))[0]
            filename    = args.outdir + 'DES_GROND_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'_GROND', f+'_GROND_ERR']][mask_good], filename, overwrite=True, format='no_header')

    except:
        print(bcolors.FAIL + 'Field not covered by DES' + bcolors.ENDC)
        pass

    # Get SkyMapper

    # Consistency check
    
    print(bcolors.OKBLUE + '\nSkyMapper catalogue' + bcolors.ENDC)

    if coordinates.dec.deg <0:

        try:

            # Retrieve catalogue

            get_skymapper                    = urllib.request.urlopen("http://skymapper.anu.edu.au/sm-cone/public/query?RA={ra:.3f}&DEC={dec:.3f}&SR={radius:.1f}&RESPONSEFORMAT=CSV".format(
                                                                    ra        = coordinates.ra.deg,
                                                                    dec        = coordinates.dec.deg,
                                                                    radius    = args.radius/60.
                                                                    )).read()

            # Write to file

            file                            = open(args.outdir + 'skymapper.csv', 'w')
            file.write(get_skymapper.decode())
            file.close()

            cat_skymapper                    = ascii.read(args.outdir + 'skymapper.csv', format='csv')
            os.system('rm ' + args.outdir + 'skymapper.csv')

            # Filter

            cat_skymapper                    = cat_skymapper[(cat_skymapper['class_star'] > 0.95) & (cat_skymapper['flags'] == 0)]
            cat_skymapper                    = cat_skymapper[(cat_skymapper['g_psf'] > 0) & (cat_skymapper['r_psf'] > 0) & (cat_skymapper['i_psf'] > 0) & (cat_skymapper['z_psf'] > 0)]

            keywords                        = ['raj2000','dej2000','u_psf','e_u_psf','v_psf','e_v_psf','g_psf','e_g_psf','r_psf','e_r_psf','i_psf','e_i_psf','z_psf','e_z_psf']

            cat_skymapper                    = cat_skymapper[keywords]

            cat_skymapper.rename_column('raj2000', 'RAJ2000')
            cat_skymapper.rename_column('dej2000', 'DEJ2000')

            for f in ['u', 'g', 'r', 'i', 'z']:
                cat_skymapper.rename_column(f + '_psf', f + '_SDSS')
                cat_skymapper.rename_column('e_' + f + '_psf', f + '_SDSS_ERR')

            # Write SDSS catalogues

            for key in [x for x in cat_skymapper.keys() if 'SDSS' in x]:
                cat_skymapper[key].format= '.4f'

            for f in ['u', 'g', 'r', 'i', 'z']:
                mask_good    = np.where((cat_skymapper[f + '_SDSS_ERR'] > cat_tools.catalog_prop['SDSS']['SIGMA_HIGH']) & (cat_skymapper[f + '_SDSS_ERR'] < cat_tools.catalog_prop['SDSS']['SIGMA_LOW']))[0]
                filename    = args.outdir + 'SkyMapper_SDSS_' + f + '.cat'
                ascii.write(cat_skymapper[['RAJ2000', 'DEJ2000', f +'_SDSS', f+'_SDSS_ERR']][mask_good], filename, overwrite=True, format='no_header')

            # Convert to Bessel system

            cat_skymapper             = cat_tools.SDSS_to_Bessel(cat_skymapper)

            # Formatting

            for key in [x for x in cat_skymapper.keys() if 'BESSEL' in x]:
                cat_skymapper[key].format= '.4f'

            # Write Bessel catalogues

            for f in ['B', 'V', 'R', 'I']:
                mask_good    = np.where((cat_skymapper[f + '_BESSEL_ERR'] > 0.) & (cat_skymapper[f + '_BESSEL_ERR'] < 0.3))[0]
                filename    = args.outdir + 'SkyMapper_BESSEL_' + f + '.cat'
                ascii.write(cat_skymapper[['RAJ2000', 'DEJ2000', f+'_BESSEL', f+'_BESSEL_ERR']][mask_good], filename, overwrite=True, format='no_header')

        except:
            print(bcolors.WARNING + 'SkyMapper query failed' + bcolors.ENDC)
            pass

    else:
        cmd                            = 'Declination is above 0 deg. No entries in SkyMapper catalogue.'
        print(bcolors.FAIL + cmd + bcolors.ENDC)

    # Get 2MASS catalogue

    print(bcolors.OKBLUE + '\n2MASS catalogues' + bcolors.ENDC)

    try:
        v         = Vizier(row_limit=10000)

        result     = v.query_region(coord.SkyCoord(ra_dd, dec_dd, unit=u.deg), radius=args.radius/60.*u.deg, catalog=cat_tools.catalog_prop['2MASS']['CATID'])#
        result    = result[cat_tools.catalog_prop['2MASS']['CATID_OUT']]
        result    = result[cat_tools.catalog_prop['2MASS']['KEYWORDS']]

        for key in [x for x in result.keys() if 'mag' in x]:

            result[key].format= '.4f'

        for f in ['J', 'H', 'K']:            

            mask_good    = np.where((result['e_'+f+'mag'] > cat_tools.catalog_prop['2MASS']['SIGMA_HIGH']) & (result['e_'+f+'mag'] < cat_tools.catalog_prop['2MASS']['SIGMA_LOW']))[0]

            filename    = args.outdir + '2MASS_' + f + '.cat'
            ascii.write(result[['RAJ2000', 'DEJ2000', f+'mag', 'e_'+f+'mag']][mask_good], filename, overwrite=True, format='no_header')

    except:
        print(bcolors.FAIL + 'VizieR query failed.' + bcolors.ENDC)
        pass

    # Get WISE catalogue

    print (bcolors.OKBLUE + '\nALLWISE' + bcolors.ENDC)

    try:
        v                    = Vizier(columns = ['all'], catalog = cat_tools.catalog_prop['WISE']['CATID'], row_limit = -1)

        data_wise            = v.query_region(coordinates, radius = 50*u.arcmin)[cat_tools.catalog_prop['WISE']['CATID_OUT']]
        data_wise            = data_wise[cat_tools.catalog_prop['WISE']['KEYWORDS']]

        mask                 = np.where( (data_wise['e_W1mag'] < 0.3) & (data_wise['e_W1mag'] > 0.0) & \
                                        (data_wise['e_W2mag'] < 0.3) & (data_wise['e_W2mag'] > 0.0)
                                        )[0]

        data_wise             = data_wise[mask]

        ascii.write(data_wise[['RAJ2000', 'DEJ2000', 'W1mag', 'e_W1mag']], 'photcat/unWISE_W1.cat', format='no_header', overwrite=True)
        ascii.write(data_wise[['RAJ2000', 'DEJ2000', 'W2mag', 'e_W2mag']], 'photcat/unWISE_W2.cat', format='no_header', overwrite=True)
        ascii.write(data_wise, 'photcat/unWISE.cat', overwrite=True)

    except:
        print(bcolors.FAIL + 'VizieR query failed.' + bcolors.ENDC)
        pass


    print('\n')


##############################################

if __name__ == "__main__":
    
    main(sys.argv[1:])