#!/usr/bin/env python

import	argparse
import	glob
import	numpy as np
import	os

parser = argparse.ArgumentParser(description='Programme to retrieve PanSTARRS images.')
parser.add_argument('--object',		type=str, help='Object name)', required=True)
parser.add_argument("--ra",	    	type=str, help='RA(J2000) of the object (HMS or DD)', required=True)
parser.add_argument('--dec', 		type=str, help='Dec(J2000) of the object (HMS or DD) (if negative do  " -12:20:20.2")', required=True)
parser.add_argument('--filter',		type=str, help='Which bands? (e.g. grizy)', default="grizy", required=False)
parser.add_argument('--width',		type=str, help='Image size in arcmin', default="20", required=False)
parser.add_argument('--outdir',		type=str, help='Output directory (default: ./)', default="./", required=False)

args 		= parser.parse_args()

random_number	= np.random.randint(1,1000)
temp_dir	= 'outdir_'+str(random_number)

os.system('mkdir %s' %(args.outdir + '/' + temp_dir))

cmd='panstamps --width=%s --filters=%s --downloadFolder=%s stack %s %s' %(args.width, args.filter, temp_dir, args.ra, args.dec)
print (cmd)
os.system(cmd)

files		= glob.glob(temp_dir+'/*')

[os.system("mv %s %s" %(x, args.outdir + args.object + '_PanSTARRS_' + x.split('_')[2] + '.fits')) for x in files]

os.system("rm -r %s" %temp_dir)