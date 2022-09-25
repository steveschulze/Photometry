import pdb
from 	astropy import coordinates as coord
from 	astropy import wcs
from 	astropy.io import ascii, fits
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
	

def sky2xy_helper(FITS, RA, DEC):

        """
        Converts physical coordinates to WCS coordinates for
        calculations with wcstools' xy2sky.
        @param image_path: FITS image file name with path.
        @type image_path: str
        @param ra: RA coordinate of object.
        @type ra: string
        @param dec: DEC coordinate of object.
        @type dec: string
        @return: tuple
        """

        if FITS:
            hdulist = fits.open(FITS)#[0]
        else:
            print("FITS image has not been provided by the user!")
            raise SystemExit

        #fo = FitsOps(IMAGE)

        #header = [hdulist[ii].header for ii in range(len(hdulist))]
        #from astropy import table
        #header = table.vstack(header)

        #header = hdu.header

        flag_save = True

        for hdu in hdulist:
            _w = wcs.WCS(hdu.header)
            if _w.array_shape != None:
                if flag_save:
                    w = wcs.WCS(hdu.header)
                    naxis2, naxis1 = w.array_shape
                    flag_save = False

        if ":" in (str(RA) or str(DEC)):
            c = coord.SkyCoord('{0} {1}'.format(
                        RA, DEC), unit=(u.hour, u.deg),
                                 frame='icrs')
        else:
            c = coord.SkyCoord('{0} {1}'.format(
                        RA, DEC), unit=(u.deg, u.deg),
                                 frame='icrs')

        # target's X and Y coordinates

        t_x, t_y = w.wcs_world2pix(c.ra.degree, c.dec.degree, 1)

        if np.isnan(t_x):
            t_x, t_y = w.wcs_world2pix(c.dec.degree, c.ra.degree, 1)

        if naxis1 < t_x or naxis2 < t_y or t_x < 0 or t_y < 0:
            #print("Provided coordinates are out of frame!")
            return(False)
        else:
            return(float(t_x), float(t_y))

def sky2xy(FITS, RA=None, DEC=None, CAT=None):

    # For individual object
    if CAT == None:
        return sky2xy_helper(FITS, RA, DEC)

    # For a set of objects
    else:
        cat = ascii.read(CAT)
        cat_keys = cat.keys()

        xy = []
        for ii in range(len(cat)):
            _sky2xy = sky2xy_helper(FITS, cat[ii][cat_keys[0]], cat[ii][cat_keys[1]])

            if _sky2xy != False:
                xy.append(_sky2xy)

    return np.array(xy)

def xy2sky(FITS, X, Y, sep=":"):

        """
        Converts physical coordinates to WCS coordinates for STDOUT.
        @param file_name: FITS image file name with path.
        @type file_name: str
        @param x: x coordinate of object.
        @type x: float
        @param y: y coordinate of object.
        @type y: float
        @param sep: delimiter for HMSDMS format.
        @type sep: float
        @return: str
        """

        try:
            header = fits.getheader(FITS)
            w = wcs.WCS(header)
            astcoords_deg = w.wcs_pix2world([[X, Y]], 0)
            c = coord.SkyCoord(astcoords_deg * u.deg,
                                             frame='fk5')

            alpha = c.to_string(style='hmsdms', sep=sep, precision=3)[0]
            delta = c.to_string(style='hmsdms', sep=sep, precision=2)[0]

            return alpha.split(" ")[0], delta.split(" ")[1]
        
        except Exception as e:
            pass