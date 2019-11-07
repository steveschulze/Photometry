#!/usr/bin/env python

import 	argparse
from	astropy.io import fits, votable
import 	numpy as np
import	phot_routines
import	os
import	sys

parser				= argparse.ArgumentParser(description='Align WCS system of two images')

parser.add_argument('--ref-image',
					type = str, 		
					help = 'Reference image', 		
					required = True)

parser.add_argument('--new-image',
					type = str, 		
					help = 'New image', 			
					required = True)

parser.add_argument('--keep-temp',
					action = 'store_true',
					help = 'Keep temporary files',
					default  = False)

args				= parser.parse_args()

# Setting up environment for sextractor

phot_routines.setup_sextractor()

# Generate catalogue of the reference image

default_sex_param		= 'default.param'
output					= open(default_sex_param, 'w')

output.write('X_WORLD\nY_WORLD\nERRA_WORLD\nERRB_WORLD\nERRTHETA_WORLD\nMAG_AUTO\nMAGERR_AUTO')
output.close()

os.system('sex -c default.sex \
			-DETECT_MINAREA 2 -DETECT_THRESH 5  -ANALYSIS_THRESH 5 \
			-DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.000005 \
			-SATUR_LEVEL 64000 -CATALOG_NAME %s.ldac -CATALOG_TYPE FITS_LDAC %s' \
			%('ref_'+args.ref_image.replace('.fits', ''), args.ref_image) )

# Generate catalogue of the image that needs to aligned

output					= open(default_sex_param, 'w')

output.write('XWIN_IMAGE\nYWIN_IMAGE\nERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\nFLUX_AUTO\nFLUXERR_AUTO\nFLUX_RADIUS\nFLAGS\nFLAGS_WEIGHT')
output.close()

os.system('sex -c default.sex \
			-DETECT_MINAREA 2 -DETECT_THRESH 3 -ANALYSIS_THRESH 3 \
			-DEBLEND_NTHRESH 64 -DEBLEND_MINCONT 0.000005 \
			-SATUR_LEVEL 64000 -CATALOG_NAME %s.ldac -CATALOG_TYPE FITS_LDAC %s' \
			%('ima_'+args.new_image.replace('.fits', ''), args.new_image) )

os.system('ls %s.ldac > scamp_input.list' %('ima_'+args.new_image.replace('.fits', '')))

# Align astrometry

os.system('scamp -dd > default.scamp')

os.system('scamp -c default.scamp -ASTREF_CATALOG FILE -ASTREFCAT_NAME %s \
			-ASTREFMAG_KEY MAG_AUTO -MATCH Y -WRITE_XML Y -XML_NAME scamp.xml\
			-CHECKPLOT_DEV NULL \
			-DISTORT_DEGREES 1 -SOLVE_ASTROM Y -SOLVE_PHOTOM N -POSITION_MAXERR 0.15\
			-SN_THRESHOLDS 10.0,40.0 @scamp_input.list' %('ref_'+args.ref_image.replace('.fits', '.ldac')))

# Bug fix of in scamp.xml

scamp_log_input			= open('scamp.xml', 'r')
scamp_log_input_lines	= scamp_log_input.readlines()
scamp_log_input.close()

scamp_log_output		= open('scamp.xml', 'w')

for line in scamp_log_input_lines:
	if 'datatype="*"' in line:
		line=string.replace(line, 'datatype="*"', 'datatype="char"')
	
	scamp_log_output.write(line)

scamp_log_output.close()

# Create log file

scamp_log				= votable.parse("scamp.xml")
scamp_log_data			= scamp_log.get_first_table()
scamp_log_file			= scamp_log_data.array['Catalog_Name']
scamp_log_contrast		= scamp_log_data.array['XY_Contrast']
scamp_log_dx			= scamp_log_data.array['DX']
scamp_log_dy			= scamp_log_data.array['DY']

# Storing diagnotisc information

if not os.path.isdir('results_scamp'): os.system('mkdir results_scamp')

if os.path.exists('results_scamp/scamp.log'):
	output				= open('results_scamp/scamp.log', 'a')
else:
	output				= open('results_scamp/scamp.log', 'w')
	output.write('File\tXY_Contrast\tDX(arcsec)\tDY(arcsec)\n')

for i in range(len(scamp_log_file)):
	output.write('{file}\t{contrast:.2f}\t{dx:.3f}\t{dy:.3f}\n'.format(
			file	= str(scamp_log_file[i]).replace('ima_', '').replace('.ldac', '.fits'),
			contrast= scamp_log_contrast[i],
			dx		= float(scamp_log_dx[i]) * 3600,
			dy		= float(scamp_log_dy[i]) * 3600
			)
	)

output.close()

# Add keywords to header

for i in range(len(scamp_log_file)):

	file				= scamp_log_file[i].decode("utf-8") 

	hdu					= fits.open(file.replace('ima_', '').replace('.ldac', '.fits'))
	hdu_header			= hdu[0].header
	hdu_data			= hdu[0].data

	scamp_results_file	= open(file.replace('.ldac', '.head'), 'r')
	scamp_results_lines	= scamp_results_file.readlines()

	for line in scamp_results_lines[3:-2]:

		keyword			= np.array(line.split('='))[0].replace(" ", "")
		value 			= np.array(line.split('='))[1].split('/')[0].replace(" ", "")

		if value[0].isdigit() or value[1].isdigit():
		 	hdu_header[keyword] = float(value)
		else:
			hdu_header[keyword] = str(value).replace(" ", "").replace("'", "")

	fits.writeto(file.replace('ima_', '').replace('.ldac', '_astro.fits'), hdu_data, hdu_header, overwrite=True)

# Keep temporary files?

if not args.keep_temp:
	os.system('rm scamp_input.list')
	os.system('rm default.*')
	os.system('rm *ldac* *head*')

for i in range(len(scamp_log_file)):
	os.system('mv scamp.xml results_scamp/{file}_scamp.xml'.format(file=scamp_log_file[i].decode("utf-8").replace('.fits', '').replace('.ldac', '').replace('ima_', '')))