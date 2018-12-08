import 	numpy as np

def statNclip(A, NITER=1000):

	"""
	This function performs statistics on a numpy array A.
	The array can contain NaN

	"""

	# MC and bootstrap from observations

	i = 0

	random_array	= []

	p25				= np.percentile(A[:,0], 25)
	p50				= np.percentile(A[:,0], 50)
	p75				= np.percentile(A[:,0], 75)
	iqr				= abs(p75 - p25)

	mask			= np.where((A[:,0] > p25 - 1.5 * iqr) & (A[:,0] < p75 + 1.5 * iqr))[0]

	while(i < NITER):
		random_MC	= np.random.normal(A[mask,0], A[mask,1])
		if len(A[:,0]) > 10:
			random_boot	= np.random.choice(random_MC, len(random_MC))
			random_array.append(np.median(random_boot))
		else:
			random_array.append(np.median(random_MC))
		i			= i+1

	zp_med			= np.percentile(random_array,50)
	zp_inf 			= np.percentile(random_array,50-34.1)
	zp_sup 			= np.percentile(random_array,50+34.1)

	return np.hstack([zp_med, zp_sup - zp_med, zp_med - zp_inf, len(A[mask,0])])


# def statNclip(A, NITER=1000):

# 	"""
# 	This function performs statistics on a numpy array A.
# 	The array can contain NaN

# 	"""

# 	# MC and bootstrap from observations

# 	i = 0

# 	random_array		= []

# 	p25				= np.percentile(A['MAG_DIFF'], 25)
# 	p50				= np.percentile(A['MAG_DIFF'], 50)
# 	p75				= np.percentile(A['MAG_DIFF'], 75)
# 	iqr				= abs(p75 - p25)

# 	mask			= np.where((A['MAG_DIFF'] > p25 - 1.5 * iqr) & (A['MAG_DIFF'] < p75 + 1.5 * iqr))[0]

# 	while(i < NITER):
# 		# random_MC	= np.random.normal(A['MAG_DIFF'], A['MAGERR_DIFF'])
# 		random_MC	= np.random.normal(A['MAG_DIFF'][mask], A['MAGERR_DIFF'][mask])
# 		random_boot	= np.random.choice(random_MC, len(random_MC))
# 		random_array.append(np.median(random_boot))
# 		i	= i+1

# 	zp_med		= np.percentile(random_array,50)
# 	zp_inf 		= np.percentile(random_array,50-34.1)
# 	zp_sup 		= np.percentile(random_array,50+34.1)

# 	return np.array([zp_med, zp_sup - zp_med, zp_med - zp_inf, len(A)])
