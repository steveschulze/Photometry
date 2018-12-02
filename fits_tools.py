from 	astropy import coordinates as coord
from 	astropy import wcs
from 	astropy.io import fits
from 	astropy import units as u
from 	misc import bcolors
import	numpy as np
import	os

def convert_hms_dd(RA, DEC):

	'''
	Convert HMS to DD system
	'''

	if (':' in RA) and (':' in DEC):
		Coord_dd	= coord.SkyCoord(RA, DEC, unit=(u.hour,u.degree), frame='icrs')
		RA_dd		= Coord_dd.ra.deg
		Dec_dd		= Coord_dd.dec.deg

	elif (not (':' in RA) and not (':' in DEC)) and (('.' in RA) and ('.' in DEC)):
		RA_dd, Dec_dd	= float(RA), float(DEC)

	else:
		print(bcolors.FAIL + 'Coordinates have wrong format.' + bcolors.ENDC)
		sys.exit()

	return RA_dd, Dec_dd

def get_header(FILE, KEYWORD):

	'''
	Get keyword from fits file
	'''

	header	= fits.getheader(FILE)
	return header[KEYWORD]

def pix2arcsec(FITS):

	'''
	Get pixel scale
	'''

	hdu		= fits.open(FITS)
	if len(hdu) > 1:
		header	= fits.getheader(FITS, 0)
		header	+= fits.getheader(FITS, 1)
	else:
		header	= fits.getheader(FITS)

	hdu_wcs			= wcs.WCS(header)
	return np.median(wcs.utils.proj_plane_pixel_scales(hdu_wcs)) * 3600
	

def sky2xy (FITS, RA=False, DEC=False, CAT=None):

	'''
	Coordinate transformation: sky -> xy
	'''

	if CAT == None:
		if RA != False and  DEC != False:
			cmd=('sky2xy %s %s %s | grep -v off' %(FITS, RA, DEC))
			program_call = os.popen(cmd)
			xy = []
			for line in program_call:
					xy=np.array(line.strip().split()[-2:]).astype(float)
			if len(xy) > 0:
				return xy

	else:
		cmd 	=("more %s | awk '{print $1,$2}' > %s" %(CAT, CAT.replace(CAT.split('.')[-1], 'reg')))
		os.system(cmd)
		cmd 	= ("sky2xy %s @%s | grep -v off | awk '{print $5, $6}'" %(FITS, CAT.replace(CAT.split('.')[-1], 'reg')))
		cat 	= os.popen(cmd)


		xy	= []
		for line in cat:
			xy.append(list(map(float, line.replace('\n', '').split())))

		return np.array(xy)

def xy2sky (FITSFILE,X,Y):

	'''
	Coordinate transformation: xy -> sky
	'''

	program_call = os.popen('xy2sky %s %s %s' %(FITSFILE, X, Y))
	sky = []

	for line in program_call:
		sky.append(line.strip().split()[:2])

	return sky
