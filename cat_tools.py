from    astropy import units as u
from    astropy.io import ascii
from    astroquery.vizier import Vizier
from    astropy import table
from    astropy import coordinates as coord
from    astropy import units as u
import  numpy as np
import  os
from    scipy.spatial import cKDTree

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

catalog_prop['UKIDSS']					= {}
catalog_prop['UKIDSS']['KEYWORDS']		= ['RAJ2000', 'DEJ2000',
											'Ymag', 'e_Ymag',
											'Jmag1', 'e_Jmag1',
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
    
def retrieve_photcat(OBJECT_PROP, PHOTCAT, CATPROP, FILENAME=None, ROW_LIMIT=-1, RADIUS=10. * u.arcmin, OUTDIR='photcat/'):

	# Query VIZIER
	v					=  Vizier(columns = ['all'], row_limit = ROW_LIMIT)
	#result				= v.query_region(coord.SkyCoord(OBJECT_PROP['RA'], OBJECT_PROP['DEC'], unit=(u.hour, u.deg)), radius = RADIUS, catalog=CATPROP[PHOTCAT]['CATID'])
	result				= v.query_region(coord.SkyCoord(OBJECT_PROP['RA'], OBJECT_PROP['DEC'], unit=u.deg), radius = RADIUS, catalog=CATPROP[PHOTCAT]['CATID'])
	result				= result[CATPROP[PHOTCAT]['CATID_OUT']]
	result				= result[CATPROP[PHOTCAT]['KEYWORDS']]

	# Some formatting

	for key in [x for x in result.keys() if 'mag' in x]:
		result[key].format= '.4f'

	# Check if output dir exists

	if not os.path.isdir(OUTDIR) and OUTDIR != './':
		os.system('mkdir %s' %OUTDIR)

	# Filter output

	if PHOTCAT == 'SDSS':
		result          = result[(result['class'] == 6) & (result['mode'] == 1)]

	if PHOTCAT == 'PanSTARRS':
		result			= result[
						(result['o_gmag'] >= 0.85) & (result['o_rmag'] >= 0.85) &
						(result['o_imag'] >= 0.85) & (result['o_zmag'] >= 0.85) &
						(result['o_ymag'] >= 0.85) &
						(np.round(result['gPSFf'], 0) >= 1) & (np.round(result['rPSFf'], 0) >= 1) &
						(np.round(result['iPSFf'], 0) >= 1) & (np.round(result['zPSFf'], 0) >= 1) &
						(np.round(result['yPSFf'], 0) >= 1) &
						(abs(result['imag'] - result['iKmag']) < 0.05) &
						(result['imag'] > 14) & (result['imag'] < 21)
						]

	# Write result to file

	for filter in OBJECT_PROP['FILTER']:
		if filter in CATPROP[PHOTCAT]['FILTER']:

			mask_good	= np.where((result['e_'+filter+'mag'] > CATPROP[PHOTCAT]['SIGMA_HIGH']) &
							(result['e_'+filter+'mag'] < CATPROP[PHOTCAT]['SIGMA_LOW']))[0]

			filename	= FILENAME
			ascii.write(result[[CATPROP[PHOTCAT]['KEYWORDS'][0], CATPROP[PHOTCAT]['KEYWORDS'][1], filter+'mag', 'e_'+filter+'mag']][mask_good], filename, overwrite=True, format='no_header')

		else:
			print(bcolors.ERROR + 'Filter {filter} not in catalogue {catalog}.'.format(filter=filter, catalog=PHOTCAT) + bcolors.ENDC)
			sys.exit()


	return None

def crossmatch(X1, X2, max_distance=np.inf):
	"""Cross-match the values between X1 and X2
	
	By default, this uses a KD Tree for speed.
	
	Parameters
	----------
	X1 : array_like
		first dataset, shape(N1, D)
	X2 : array_like
		second dataset, shape(N2, D)
	max_distance : float (optional)
		maximum radius of search.  If no point is within the given radius,
		then inf will be returned.
	
	Returns
	-------
	dist, ind: ndarrays
		The distance and index of the closest point in X2 to each point in X1
		Both arrays are length N1.
		Locations with no match are indicated by
		dist[i] = inf, ind[i] = N2
	
	Taken from astroML. Add multi-processing capabilities
	
	"""
	X1			= np.asarray(X1, dtype=float)
	X2	 		= np.asarray(X2, dtype=float)
	
	N1, D		= X1.shape
	N2, D2 		= X2.shape
	
	if D != D2:
		raise ValueError('Arrays must have the same second dimension')
	
	kdt			= cKDTree(X2)
	
	dist, ind	= kdt.query(X1, k=1, distance_upper_bound=max_distance, n_jobs=-1)
	
	return dist, ind
 
def crossmatch_angular(X1, X2, max_distance=np.inf):
	"""Cross-match angular values between X1 and X2

	by default, this uses a KD Tree for speed.  Because the
	KD Tree only handles cartesian distances, the angles
	are projected onto a 3D sphere.

	Parameters
	----------
	X1 : array_like
		first dataset, shape(N1, 2). X1[:, 0] is the RA, X1[:, 1] is the DEC,
		both measured in degrees
	X2 : array_like
		second dataset, shape(N2, 2). X2[:, 0] is the RA, X2[:, 1] is the DEC,
		both measured in degrees
	max_distance : float (optional)
		maximum radius of search, measured in degrees.
		If no point is within the given radius, then inf will be returned.

	Returns
	-------
	dist, ind: ndarrays
		The angular distance and index of the closest point in X2 to
		each point in X1.  Both arrays are length N1.
		Locations with no match are indicated by
		dist[i] = inf, ind[i] = N2

	Taken from astroML.
	"""

	X1 = X1 * (np.pi / 180.)
	X2 = X2 * (np.pi / 180.)
	max_distance = max_distance * (np.pi / 180.)

	# Convert 2D RA/DEC to 3D cartesian coordinates
	Y1 = np.transpose(np.vstack([np.cos(X1[:, 0]) * np.cos(X1[:, 1]),
								np.sin(X1[:, 0]) * np.cos(X1[:, 1]),
								np.sin(X1[:, 1])]))
	Y2 = np.transpose(np.vstack([np.cos(X2[:, 0]) * np.cos(X2[:, 1]),
								np.sin(X2[:, 0]) * np.cos(X2[:, 1]),
								np.sin(X2[:, 1])]))

	# law of cosines to compute 3D distance
	max_y = np.sqrt(2 - 2 * np.cos(max_distance))
	dist, ind = crossmatch(Y1, Y2, max_y)

	# convert distances back to angles using the law of tangents
	not_inf = ~np.isinf(dist)
	x = 0.5 * dist[not_inf]
	dist[not_inf] = (180. / np.pi * 2 * np.arctan2(x,
								np.sqrt(np.maximum(0, 1 - x ** 2))))

	return dist, ind
 
def wrapper_crossmatch(FILE1, FILE2, RADIUS):
 
    # # Load data file
 
    data_1              = FILE1#np.loadtxt(FILE1)
    data_2              = FILE2#np.loadtxt(FILE2)
 
    # Make matricies of coordinates
 
    dist, idx           = crossmatch_angular(data_1[:,:2], data_2[:,:2], 1)
    result              = np.hstack([data_2[idx,:], data_1])
    result              = np.hstack([result, dist.reshape(-1,1)])
 
    # Filter data
 
    result              = result[result[:,-1] < RADIUS/3600.,:]
 
    return result#table.Table(result, names=('RA_1', 'DEC_1', 'MAG_1', 'MAGERR_1', 'RA_2', 'DEC_2', 'MAG_2', 'MAGERR_2', 'DIST'))