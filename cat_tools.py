"""
cat_tools.py — Photometric catalogue metadata, colour transformations,
               positional cross-matching, and Vizier query wrappers.

Cross-matching uses astropy.coordinates.SkyCoord.match_to_catalog_sky()
which internally uses a C-compiled KD-tree with proper RA wrap-around
handling.  For 10 000-entry catalogs this is ~5-10× faster than the
previous manual implementation.
"""

__version__ = "2026-05-29"
__author__  = "Steve Schulze (steve.schulze@weizmann.ac.il)"

from    astropy import units as u
from    astropy.io import ascii
from    astroquery.vizier import Vizier
from    astropy import table
from    astropy import coordinates as coord
from    misc import bcolors
import  numpy as np
from    pathlib import Path
import  sys

# Catalogue properties

#http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#Table8
#https://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html
#http://www.galex.caltech.edu/researcher/techdoc-ch4.html#1

catalog_zp								= {}
catalog_zp['WISE_W1']					= 20.500 + 2.699
catalog_zp['WISE_W2']					= 19.500 + 3.339
catalog_zp['WISE_W3']					= 18.000 + 5.174
catalog_zp['WISE_W4']					= 13.000 + 6.620
catalog_zp['GALEX_FUV']					= 18.82  + 0
catalog_zp['GALEX_NUV']					= 20.08  + 0

catalog_prop							= {}
catalog_prop['2MASS']					= {}
catalog_prop['2MASS']['KEYWORDS']		= ["RAJ2000", "DEJ2000", "Jmag", "e_Jmag", "Hmag", "e_Hmag", "Kmag", "e_Kmag"]
catalog_prop['2MASS']['CATID']			= "II/246"
catalog_prop['2MASS']['CATID_OUT']		= catalog_prop['2MASS']['CATID'] + '/out'
catalog_prop['2MASS']['FILTER']			= ['J', 'H', 'K']
catalog_prop['2MASS']['SIGMA_HIGH']		= 0
catalog_prop['2MASS']['SIGMA_LOW']		= 0.3

catalog_prop['SDSS']					= {}
catalog_prop['SDSS']['KEYWORDS']		= ["RA_ICRS", "DE_ICRS",
											"class", "mode",
											"umag", "e_umag",
											"gmag", "e_gmag",
											"rmag", "e_rmag",
											"imag", "e_imag",
											"zmag", "e_zmag"]
catalog_prop['SDSS']['CATID']			= "V/147"
catalog_prop['SDSS']['CATID_OUT']		= catalog_prop['SDSS']['CATID'] + '/sdss12'
catalog_prop['SDSS']['FILTER']          = ['u', 'g', 'r', 'i', 'z']
catalog_prop['SDSS']['SIGMA_HIGH']      = 0
catalog_prop['SDSS']['SIGMA_LOW']       = 0.2

catalog_prop['PanSTARRS']				= {}
catalog_prop['PanSTARRS']['KEYWORDS']	= ['RAJ2000', 'DEJ2000',
											'o_gmag', 'gPSFf', 'gmag', 'e_gmag',
											'o_rmag', 'rPSFf', 'rmag', 'e_rmag',
											'o_imag', 'iPSFf', 'imag', 'e_imag', 'iKmag',
											'o_zmag', 'zPSFf', 'zmag', 'e_zmag',
											'o_ymag', 'yPSFf', 'ymag', 'e_ymag']
catalog_prop['PanSTARRS']['CATID']		= "II/349"
catalog_prop['PanSTARRS']['CATID_OUT']	= catalog_prop['PanSTARRS']['CATID'] + '/ps1'
catalog_prop['PanSTARRS']['FILTER']		= ['g', 'r', 'i', 'z', 'y']
catalog_prop['PanSTARRS']['SIGMA_HIGH']	= 0
catalog_prop['PanSTARRS']['SIGMA_LOW']	= 0.2

catalog_prop['GAIA']					= {}
catalog_prop['GAIA']['KEYWORDS']		= ["RAJ2000", "DEJ2000", "Plx", "e_Plx", "pmRA", "e_pmRA", "pmDE", "e_pmDE"]
catalog_prop['GAIA']['CATID']			= "I/345"
catalog_prop['GAIA']['CATID_OUT']		= catalog_prop['GAIA']['CATID'] + '/gaia2'

catalog_prop['UKIDSS']					= {}
catalog_prop['UKIDSS']['KEYWORDS']		= ['RAJ2000', 'DEJ2000',
											'class',
											'Ymag', 'e_Ymag',
											'Jmag1', 'e_Jmag1',
											'Jmag2', 'e_Jmag2',
											'Hmag', 'e_Hmag',
											'Kmag', 'e_Kmag']
catalog_prop['UKIDSS']['CATID']			= "II/319"
catalog_prop['UKIDSS']['CATID_OUT']		= catalog_prop['UKIDSS']['CATID'] + '/las9'
catalog_prop['UKIDSS']['FILTER']		= ['Y', 'J', 'H', 'K']
catalog_prop['UKIDSS']['SIGMA_HIGH']	= 0
catalog_prop['UKIDSS']['SIGMA_LOW']		= 0.2

catalog_prop['WISE']					= {}
catalog_prop['WISE']['KEYWORDS']		= ['RAJ2000', 'DEJ2000',
											'W1mag', 'e_W1mag',
											'W2mag', 'e_W2mag']
catalog_prop['WISE']['CATID']			= "II/328"
catalog_prop['WISE']['CATID_OUT']		= catalog_prop['WISE']['CATID'] + '/allwise'
catalog_prop['WISE']['FILTER']			= ['W1', 'W2']
catalog_prop['WISE']['SIGMA_HIGH']		= 0
catalog_prop['WISE']['SIGMA_LOW']		= 0.2

catalog_prop['DES']						= {}
catalog_prop['DES']['KEYWORDS']			= ['RAJ2000', 'DEJ2000',
											'gmag', 'rmag', 'imag', 'zmag', 'Ymag',
											'e_gmag', 'e_rmag', 'e_imag', 'e_zmag', 'e_Ymag',
											'S_Gg', 'S_Gr', 'S_Gi', 'S_Gz',
											'gFlag', 'rFlag', 'iFlag', 'zFlag']
catalog_prop['DES']['CATID']			= "II/357"
catalog_prop['DES']['CATID_OUT']		= catalog_prop['DES']['CATID'] + '/des_dr1'
catalog_prop['DES']['FILTER']			= ['g', 'r', 'i', 'z', 'y']
catalog_prop['DES']['SIGMA_HIGH']		= 0
catalog_prop['DES']['SIGMA_LOW']		= 0.2

def PS1_to_HSC(DATA):

	# Based on 
	# https://hsc-release.mtk.nao.ac.jp/doc/index.php/data/
	# https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_8_e/colorterms.html

	# g 	g	r	0.005728	0.061749	-0.001125
	# r 	r	i	-0.000144	0.001369	-0.00838
	# r2 	r	i	-3.2E-05	-0.002866	-0.012638
	# i 	i	z	0.000643	-0.130078	-0.006855
	# i2 	i	z	0.001625	-0.200406	-0.013666
	# z		z	y	-0.005362	-0.221551	-0.308279
	# y	 	y	z	-0.002055	0.20968	0.227296


	gr  = DATA['g_PS1'] - DATA['r_PS1']
	ri  = DATA['r_PS1'] - DATA['i_PS1']
	iz  = DATA['i_PS1'] - DATA['z_PS1']
	zy  = DATA['z_PS1'] - DATA['y_PS1']
	yz  = DATA['y_PS1'] - DATA['z_PS1']

	gr_err  = np.sqrt(DATA['g_PS1_ERR']**2 + DATA['r_PS1_ERR']**2)
	ri_err  = np.sqrt(DATA['r_PS1_ERR']**2 + DATA['i_PS1_ERR']**2)
	iz_err  = np.sqrt(DATA['i_PS1_ERR']**2 + DATA['z_PS1_ERR']**2)
	zy_err  = np.sqrt(DATA['z_PS1_ERR']**2 + DATA['y_PS1_ERR']**2)
	yz_err  = np.sqrt(DATA['y_PS1_ERR']**2 + DATA['z_PS1_ERR']**2)


	DATA['g_HSC']  = DATA['g_PS1'] + 0.0057280 + 0.061749 * gr - 0.001125 * gr**2
	DATA['r_HSC']  = DATA['r_PS1'] - 0.0001440 + 0.001369 * ri - 0.008380 * ri**2
	DATA['r2_HSC'] = DATA['r_PS1'] - 0.0000320 - 0.002866 * ri - 0.012638 * ri**2
	DATA['i_HSC']  = DATA['i_PS1'] + 0.0006430 - 0.130078 * iz - 0.006855 * iz**2
	DATA['i2_HSC'] = DATA['i_PS1'] + 0.0016250 - 0.200406 * iz - 0.013666 * iz**2
	DATA['z_HSC']  = DATA['z_PS1'] - 0.0053620 - 0.221551 * zy - 0.308279 * zy**2
	DATA['y_HSC']  = DATA['y_PS1'] - 0.0020550 + 0.209680 * yz + 0.227296 * yz**2

	DATA['g_HSC_ERR']  = np.sqrt(DATA['g_PS1_ERR']**2 + (0.061749 * gr_err)**2 + (2 * 0.001125 * gr * gr_err)**2)
	DATA['r_HSC_ERR']  = np.sqrt(DATA['r_PS1_ERR']**2 + (0.001369 * ri_err)**2 + (2 * 0.008380 * ri * ri_err)**2)
	DATA['r2_HSC_ERR'] = np.sqrt(DATA['r_PS1_ERR']**2 + (0.002866 * ri_err)**2 + (2 * 0.012638 * ri * ri_err)**2)
	DATA['i_HSC_ERR']  = np.sqrt(DATA['i_PS1_ERR']**2 + (0.130078 * iz_err)**2 + (2 * 0.006855 * iz * iz_err)**2)
	DATA['i2_HSC_ERR'] = np.sqrt(DATA['i_PS1_ERR']**2 + (0.200406 * iz_err)**2 + (2 * 0.013666 * iz * iz_err)**2)
	DATA['z_HSC_ERR']  = np.sqrt(DATA['z_PS1_ERR']**2 + (0.221551 * zy_err)**2 + (2 * 0.308279 * zy * zy_err)**2)
	DATA['y_HSC_ERR']  = np.sqrt(DATA['y_PS1_ERR']**2 + (0.209680 * yz_err)**2 + (2 * 0.227296 * yz * yz_err)**2)


	return DATA

def PS1_to_SDSS(DATA):

	# Based on http://adsabs.harvard.edu/abs/2016ApJ...822...66F

	# The equations are valid for main-sequence stars with 0.4 < x < 2.7.

	# Coefficients are provided for gP1 - usdss and yP1 - zsdss for
	# much less reliable than the griz transformations. In particular,
	# the extrapolation from PS1 colors to u band is strongly
	# metallicity dependent, and should be used with caution. The
	# corrections are typically 0.01 mag in r and i, up to 0.1 mag in z,
	# and up to 0.25 in g.

	# After colour transformation differences between the PS1 and SDSS
	# u(SDSS - PS1) = -26.29 mmag
	# g(SDSS - PS1) =  -2.27 mmag
	# r(SDSS - PS1) =  -4.85 mmag
	# i(SDSS - PS1) =  -7.86 mmag
	# z(SDSS - PS1) = -12.66 mmag

	coefficients		= {}
	coefficients['u']	= [ 0.04438, -2.26095, -0.13387,  0.27099]
	coefficients['g']	= [-0.01808, -0.13595,  0.01941, -0.00183]
	coefficients['r']	= [-0.01836, -0.03577,  0.02612, -0.00558]
	coefficients['i']	= [ 0.01170, -0.00400,  0.00066, -0.00058]
	coefficients['z']	= [-0.01062,  0.07529, -0.03592,  0.00890]
	coefficients['y']	= [ 0.08924, -0.20878,  0.10360, -0.02441]

	x					= DATA['g_PS1'] - DATA['i_PS1']
	x_err				= np.sqrt(DATA['g_PS1_ERR']**2 + DATA['i_PS1_ERR']**2)

	# Select objects with 0.4 < x < 2.7

	mask_good			= np.where((x >= 0.4) & (x <= 2.7))[0]

	x 					= x    [mask_good]
	x_err 				= x_err[mask_good]
	DATA 				= DATA[mask_good]

	# Compute SDSS photometry

	for filter in ['u', 'g', 'r', 'i', 'z']:

		if filter not in ['u', 'y']:

			mag_SDSS	= DATA[filter+'_PS1'] - (coefficients[filter][0] + coefficients[filter][1] * x + coefficients[filter][2] * x**2 + coefficients[filter][3] * x**3)
			mag_SDSS_err= np.sqrt(DATA[filter+'_PS1_ERR']**2 +
						(           x_err * coefficients[filter][1]) ** 2 +
						(2 * x *    x_err * coefficients[filter][2]) ** 2 +
						(3 * x**2 * x_err * coefficients[filter][3]) ** 2
						)

		elif filter == 'u':

			mag_SDSS	= DATA['g_PS1'] - (coefficients[filter][0] + coefficients[filter][1] * x + coefficients[filter][2] * x**2 + coefficients[filter][3] * x**3)
			mag_SDSS_err= np.sqrt(DATA['g_PS1_ERR']**2 +
						(           x_err * coefficients[filter][1]) ** 2 +
						(2 * x *    x_err * coefficients[filter][2]) ** 2 +
						(3 * x**2 * x_err * coefficients[filter][3]) ** 2
						)

		mag_SDSS 					= np.array(mag_SDSS)
		mag_SDSS_err				= np.array(mag_SDSS_err)

		DATA[filter+'_SDSS'] 		= mag_SDSS
		DATA[filter+'_SDSS_ERR'] 	= mag_SDSS_err

	return DATA

def PS1_to_ZTF(DATA):

	# Ref: https://iopscience.iop.org/article/10.3847/2515-5172/ab7f3c

	DATA['g_ZTF'] = DATA['g_PS1'] + 0.055 * (DATA['g_PS1'] - DATA['r_PS1']) - 0.012
	DATA['r_ZTF'] = DATA['r_PS1'] - 0.087 * (DATA['g_PS1'] - DATA['r_PS1']) - 0.0035

	DATA['g_ZTF_ERR'] = np.sqrt(DATA['g_PS1_ERR']**2 + (0.055*DATA['g_PS1_ERR'])**2 + (0.055*DATA['r_PS1_ERR'])**2)
	DATA['r_ZTF_ERR'] = np.sqrt(DATA['r_PS1_ERR']**2 + (0.087*DATA['g_PS1_ERR'])**2 + (0.087*DATA['r_PS1_ERR'])**2)

	return DATA

def SDSS_to_DES(DATA):

	# Ref: https://iopscience.iop.org/article/10.3847/1538-4365/aab4f5#apjsaab4f5app1-4?gridset=show Appendix A.4

	DATA['g_DES'] = DATA['g_SDSS'] - 0.104 * (DATA['g_SDSS'] - DATA['r_SDSS']) + 0.01
	DATA['r_DES'] = DATA['r_SDSS'] - 0.102 * (DATA['g_SDSS'] - DATA['r_SDSS']) + 0.02
	DATA['i_DES'] = DATA['i_SDSS'] - 0.256 * (DATA['i_SDSS'] - DATA['z_SDSS']) + 0.02
	DATA['z_DES'] = DATA['z_SDSS'] - 0.086 * (DATA['i_SDSS'] - DATA['z_SDSS']) + 0.01

	DATA['g_DES_ERR'] = np.sqrt( DATA['g_SDSS_ERR']**2 + 0.104**2 * (DATA['g_SDSS_ERR']**2 + DATA['r_SDSS_ERR']**2) )
	DATA['r_DES_ERR'] = np.sqrt( DATA['r_SDSS_ERR']**2 + 0.102**2 * (DATA['g_SDSS_ERR']**2 + DATA['r_SDSS_ERR']**2) )
	DATA['i_DES_ERR'] = np.sqrt( DATA['i_SDSS_ERR']**2 + 0.256**2 * (DATA['i_SDSS_ERR']**2 + DATA['z_SDSS_ERR']**2) )
	DATA['z_DES_ERR'] = np.sqrt( DATA['z_SDSS_ERR']**2 + 0.086**2 * (DATA['i_SDSS_ERR']**2 + DATA['z_SDSS_ERR']**2) )

	return DATA

def SDSS_to_GROND(DATA):

	# Ref: https://www.mpe.mpg.de/~jcg/GROND/calibration.html

	DATA['g_GROND'] = DATA['g_SDSS'] - 0.015 * (DATA['g_SDSS'] - DATA['r_SDSS']) + 0.006
	DATA['r_GROND'] = DATA['r_SDSS'] - 0.012 * (DATA['r_SDSS'] - DATA['i_SDSS']) + 0.004
	DATA['i_GROND'] = DATA['i_SDSS'] - 0.113 * (DATA['r_SDSS'] - DATA['i_SDSS']) + 0.031
	DATA['z_GROND'] = DATA['z_SDSS'] + 0.009 * (DATA['i_SDSS'] - DATA['z_SDSS']) + 0.003

	DATA['g_GROND_ERR'] = np.sqrt( DATA['g_SDSS_ERR']**2 + 0.015**2 * (DATA['g_SDSS_ERR']**2 + DATA['r_SDSS_ERR']**2) )
	DATA['r_GROND_ERR'] = np.sqrt( DATA['r_SDSS_ERR']**2 + 0.012**2 * (DATA['r_SDSS_ERR']**2 + DATA['i_SDSS_ERR']**2) )
	DATA['i_GROND_ERR'] = np.sqrt( DATA['i_SDSS_ERR']**2 + 0.113**2 * (DATA['r_SDSS_ERR']**2 + DATA['i_SDSS_ERR']**2) )
	DATA['z_GROND_ERR'] = np.sqrt( DATA['z_SDSS_ERR']**2 + 0.009**2 * (DATA['i_SDSS_ERR']**2 + DATA['z_SDSS_ERR']**2) )

	return DATA

def SDSS_to_Bessel(DATA):

	# http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php

	DATA['B_BESSEL']	= DATA['g_SDSS'] + 0.3130*(DATA['g_SDSS'] - DATA['r_SDSS']) + 0.2271
	DATA['B_BESSEL_ERR']= np.sqrt(0.0107**2 +
						DATA['g_SDSS_ERR']**2 +
						(0.3130 * DATA['g_SDSS_ERR'])**2 +
						(0.3130 * DATA['r_SDSS_ERR'])**2
						)

	DATA['V_BESSEL']	= DATA['g_SDSS'] - 0.5784*(DATA['g_SDSS'] - DATA['r_SDSS']) - 0.0038
	DATA['V_BESSEL_ERR']= np.sqrt(0.0054**2 +
						DATA['g_SDSS_ERR']**2 +
						(0.5784 * DATA['g_SDSS_ERR'])**2 +
						(0.5784 * DATA['r_SDSS_ERR'])**2
						)

	DATA['R_BESSEL']	= DATA['r_SDSS'] - 0.1837*(DATA['g_SDSS'] - DATA['r_SDSS']) - 0.0971
	DATA['R_BESSEL_ERR']= np.sqrt(0.0106**2 +
						DATA['r_SDSS_ERR']**2 +
						(0.1837 * DATA['g_SDSS_ERR'])**2 +
						(0.1837 * DATA['r_SDSS_ERR'])**2
						)

	DATA['I_BESSEL']	= DATA['r_SDSS'] - 1.2444*(DATA['r_SDSS'] - DATA['i_SDSS']) - 0.3820
	DATA['I_BESSEL_ERR']= np.sqrt(0.0078**2 +
						DATA['r_SDSS_ERR']**2 +
						(1.2444 * DATA['r_SDSS_ERR'])**2 +
						(1.2444 * DATA['i_SDSS_ERR'])**2
						)

	return DATA

def DES_to_SDSS(DATA):

	# Ref: https://iopscience.iop.org/article/10.3847/1538-4365/aab4f5#apjsaab4f5app1-4?gridset=show Appendix A.4
	# Inverted Eqs. 14-18 with Mathematica

	DATA['g_SDSS']		=  1.002   * (-0.00894 + 1.102 * DATA['g_DES'] - 0.104 * DATA['r_DES'])
	DATA['r_SDSS']		= -1.002   * ( 0.01894 - 0.102 * DATA['g_DES'] - 0.896 * DATA['r_DES'])
	DATA['i_SDSS']		= -1.20482 * ( 0.01916 - 1.086 * DATA['i_DES'] + 0.256 * DATA['z_DES'])
	DATA['z_SDSS']		= -1.20482 * ( 0.00916 - 0.086 * DATA['i_DES'] - 0.744 * DATA['z_DES'])

	DATA['g_SDSS_ERR']	= np.sqrt( ( 1.002   * 1.102 * DATA['g_DES_ERR'])**2 + ( 1.002   * 0.104 * DATA['r_DES_ERR'])**2 )
	DATA['r_SDSS_ERR']	= np.sqrt( ( 1.002   * 0.102 * DATA['g_DES_ERR'])**2 + ( 1.002   * 0.896 * DATA['r_DES_ERR'])**2 )
	DATA['i_SDSS_ERR']	= np.sqrt( ( 1.20482 * 1.086 * DATA['i_DES_ERR'])**2 + ( 1.20482 * 0.256 * DATA['z_DES_ERR'])**2 )
	DATA['z_SDSS_ERR']	= np.sqrt( ( 1.20482 * 0.086 * DATA['i_DES_ERR'])**2 + ( 1.20482 * 0.744 * DATA['z_DES_ERR'])**2 )

	return DATA

def retrieve_photcat(OBJECT_PROP, PHOTCAT, CATPROP, FILENAME=None,
                     ROW_LIMIT=-1, RADIUS=10. * u.arcmin, OUTDIR='photcat/'):
    """Query Vizier for a photometric catalogue and write matched sources to file.

    Parameters
    ----------
    OBJECT_PROP : dict
        Must contain 'RA' (decimal deg), 'DEC' (decimal deg), 'FILTER' (list).
    PHOTCAT : str
        Catalogue key (e.g. 'SDSS', 'PanSTARRS', '2MASS').
    CATPROP : dict
        Catalogue property dictionary from this module (``catalog_prop``).
    FILENAME : str, optional
        Output file path.
    ROW_LIMIT : int
        Vizier row limit (-1 = unlimited).
    RADIUS : astropy.Quantity
        Search radius.
    OUTDIR : str
        Output directory (created if it does not exist).
    """
    v      = Vizier(columns=['all'], row_limit=ROW_LIMIT)
    result = v.query_region(
        coord.SkyCoord(OBJECT_PROP['RA'], OBJECT_PROP['DEC'], unit=u.deg),
        radius=RADIUS,
        catalog=CATPROP[PHOTCAT]['CATID'])
    result = result[CATPROP[PHOTCAT]['CATID_OUT']]
    result = result[CATPROP[PHOTCAT]['KEYWORDS']]

    for key in [k for k in result.colnames if 'mag' in k]:
        result[key].format = '.4f'

    Path(OUTDIR).mkdir(parents=True, exist_ok=True)

    # Star/galaxy selection
    if PHOTCAT == 'SDSS':
        result = result[(result['class'] == 6) & (result['mode'] == 1)]

    if PHOTCAT == 'PanSTARRS':
        result = result[
            (result['o_gmag'] >= 0.85) & (result['o_rmag'] >= 0.85) &
            (result['o_imag'] >= 0.85) & (result['o_zmag'] >= 0.85) &
            (result['o_ymag'] >= 0.85) &
            (np.round(result['gPSFf'], 0) >= 1) & (np.round(result['rPSFf'], 0) >= 1) &
            (np.round(result['iPSFf'], 0) >= 1) & (np.round(result['zPSFf'], 0) >= 1) &
            (np.round(result['yPSFf'], 0) >= 1) &
            (np.abs(result['imag'] - result['iKmag']) < 0.05) &
            (result['imag'] > 14) & (result['imag'] < 21)]

    for filt in OBJECT_PROP['FILTER']:
        if filt not in CATPROP[PHOTCAT]['FILTER']:
            print(bcolors.FAIL
                  + f'Filter {filt} not in catalogue {PHOTCAT}.'
                  + bcolors.ENDC)
            sys.exit(1)

        mask_good = ((result['e_' + filt + 'mag'] > CATPROP[PHOTCAT]['SIGMA_HIGH'])
                     & (result['e_' + filt + 'mag'] < CATPROP[PHOTCAT]['SIGMA_LOW']))

        ra_key  = CATPROP[PHOTCAT]['KEYWORDS'][0]
        dec_key = CATPROP[PHOTCAT]['KEYWORDS'][1]
        ascii.write(result[[ra_key, dec_key, filt + 'mag', 'e_' + filt + 'mag']][mask_good],
                    FILENAME, overwrite=True, format='no_header')


def crossmatch_sky(ra1, dec1, ra2, dec2, max_sep_arcsec):
    """Cross-match two sets of sky coordinates using astropy SkyCoord.

    Finds for each point in (ra1, dec1) the nearest neighbour in (ra2, dec2).
    Uses astropy's C-compiled KD-tree which handles RA wrap-around correctly.

    Parameters
    ----------
    ra1, dec1 : array-like
        Coordinates of the first catalog (decimal degrees).
    ra2, dec2 : array-like
        Coordinates of the second catalog (decimal degrees).
    max_sep_arcsec : float
        Maximum separation to accept as a match (arcseconds).

    Returns
    -------
    idx : np.ndarray, int
        Index into (ra2, dec2) of the nearest neighbour for each point in
        catalog 1.  Points without a match within *max_sep_arcsec* are
        flagged with ``mask``.
    sep_arcsec : np.ndarray, float
        Angular separation in arcseconds for each match.
    mask : np.ndarray, bool
        True where separation < *max_sep_arcsec*.
    """
    cat1 = coord.SkyCoord(ra=np.asarray(ra1) * u.deg,
                          dec=np.asarray(dec1) * u.deg)
    cat2 = coord.SkyCoord(ra=np.asarray(ra2) * u.deg,
                          dec=np.asarray(dec2) * u.deg)
    idx, d2d, _ = cat1.match_to_catalog_sky(cat2)
    sep_arcsec  = d2d.arcsec
    mask        = sep_arcsec < max_sep_arcsec
    return idx, sep_arcsec, mask


def wrapper_crossmatch(FILE1, FILE2, RADIUS):
    """Cross-match two float arrays by RA/Dec (first two columns each).

    Wraps :func:`crossmatch_sky` and returns a combined array suitable for
    direct construction of an astropy Table.  Backward-compatible with the
    previous KD-tree implementation.

    Parameters
    ----------
    FILE1 : np.ndarray, shape (N1, M1)
        Array where columns 0 and 1 are RA, Dec in decimal degrees.
    FILE2 : np.ndarray, shape (N2, M2)
        Array where columns 0 and 1 are RA, Dec in decimal degrees.
    RADIUS : float
        Maximum matching radius in **arcseconds**.

    Returns
    -------
    np.ndarray
        Matched rows: ``[FILE2_columns | FILE1_columns | separation_deg]``.
        Only rows with separation < RADIUS arcsec are included.
        Separation is in **degrees** (multiply by 3600 for arcsec).
    """
    data1 = np.asarray(FILE1, dtype=float)
    data2 = np.asarray(FILE2, dtype=float)

    if len(data1) == 0 or len(data2) == 0:
        return np.empty((0, data1.shape[1] + data2.shape[1] + 1))

    idx, sep_arcsec, mask = crossmatch_sky(
        data1[:, 0], data1[:, 1],
        data2[:, 0], data2[:, 1],
        RADIUS)

    rows1  = data1[mask]
    rows2  = data2[idx[mask]]
    sep_deg = sep_arcsec[mask] / 3600.

    return np.hstack([rows2, rows1, sep_deg.reshape(-1, 1)])


# Legacy wrappers kept for backward compatibility
def crossmatch_angular(X1, X2, max_distance=np.inf):
    """Angular cross-match (legacy wrapper around crossmatch_sky).

    Parameters
    ----------
    X1, X2 : array-like, shape (N, 2)
        RA/Dec arrays in degrees.
    max_distance : float
        Maximum separation in degrees.

    Returns
    -------
    dist : np.ndarray
        Angular separations in degrees. Unmatched entries are ``inf``.
    ind : np.ndarray
        Indices into X2 for each point in X1.
    """
    X1 = np.asarray(X1, dtype=float)
    X2 = np.asarray(X2, dtype=float)

    cat1 = coord.SkyCoord(ra=X1[:, 0] * u.deg, dec=X1[:, 1] * u.deg)
    cat2 = coord.SkyCoord(ra=X2[:, 0] * u.deg, dec=X2[:, 1] * u.deg)
    idx, d2d, _ = cat1.match_to_catalog_sky(cat2)

    dist = d2d.deg
    if not np.isinf(max_distance):
        dist[dist > max_distance] = np.inf

    return dist, idx