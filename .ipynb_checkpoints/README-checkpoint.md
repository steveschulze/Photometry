# Photometry tools

Tools to retrieve photometry catalogues, improve astrometry and perform aperture photometry.

## Prerequisites

Python libraries

```
python 2.7, astropy, astroquery, numpy, photutils, scipy, sewpy, and sip2pv.py (included)
```

```
astrometry.net, sextractor and wcs-tools
```

## Installing

Store the directory in your favourite folder and add the folder to your PATH and PYTHON_PATH environment.

```
export PATH=$PATH:/new/path
export PYTHON_PATH=$PYTHON_PATH:/new/path
```

Give field_calibration.py, improve_astro.py and photometry.py permission to be executed.

```
chmod u+x PROGRAMME.py
```

## Description of the tools

### improve_astrometry.py

#### How to call it?

```
usage: improve_astro.py [-h] --fits FITS --ra RA --dec DEC [--radius RADIUS]
                        [--downsample DOWNSAMPLE] [--tweak-order TWEAK_ORDER]

Improving WCS calibration with astrometry.net and converting SIP terms to
Sextractor format.

optional arguments:
  -h, --help            show this help message and exit
  --fits FITS           File name (Required)
  --ra RA               RA(J2000) of the Object (HMS and DD system allowed).
                        Note: if declination is negative, write "
                        -12:20:20.2"). (Required)
  --dec DEC             Dec (J2000) of the Object (HMS and DD system allowed).
                        (Required)
  --radius RADIUS       only search in indexes within 'radius' of the field
                        center (unit: deg; default: 0.125 deg)
  --downsample DOWNSAMPLE
                        downsample the image by factor <int> before running
                        source extraction (default: 2)
  --tweak-order TWEAK_ORDER
                        Polynomial order of SIP WCS corrections (default: 2)
```

#### How does it work?

It uses astrometry.net and sextractor to establish the field calibration of the image. Afterwards it uses sip2pv.py to transform the distortion keywords into the format sextractor understands. The new image has the suffix '_wcs'.


### field_calibration.py

#### How to call it?

```
usage: field_calibration.py [-h] --ra RA --dec DEC [--radius RADIUS]
                            [--outdir OUTDIR] [--type TYPE]

Retrieve photometric catalogues. 2MASS, PS1 and SDSS are the input sources.
Bessel catalogue is generated through colour equations. PS1 and SDSS cats are
in the AB system, whereas Bessel and 2MASS cats are in the Vega system.

optional arguments:
  -h, --help       show this help message and exit
  --ra RA          RA(J2000) of the Object (HMS and DD system allowed). Note:
                   if declination is negative, write " -12:20:20.2"). (Required)
                   keyword.
  --dec DEC        Dec (J2000) of the Object (HMS and DD system allowed).
                   (Required)
  --radius RADIUS  only search in indexes within 'radius' of the field center
                   (unit: arcmin; default: 20)
  --outdir OUTDIR  Output directory (default: results/)
  --type TYPE      Generate catalogues for optical/nir/all (default: all)
```

#### How does it work?

The photometric catalogues are build from the 2MASS point source catalogues and SDSS source catalogues. The PS1 point source catalogue was build following XXX. PS1 photometry was converted to the SDSS filters using the colour equations in XXX. Bessel photometry was derived following the [Lupton (2004)](http://classic.sdss.org/dr4/algorithms/sdssUBVRITransform.html) colour equations.

The files will be called ```PS1_FILTER.ascii```, ```PS1_SDSS_FILTER.ascii```, ```PS1_BESSEL_FILTER.ascii```, ```SDSS_FILTER.ascii```, ```SDSS_BESSEL_FILTER.ascii``` and ```2MASS_FILTER.ascii```. The columns in each file are: col1 = ra, col2 = dec, col3 = mag, and col4 = sigma_mag.

### photometry.py

#### How to call it?

```
usage: photometry_v6.py [-h] --ra RA --dec DEC --fits FITS
                        [--host-offset HOST_OFFSET] [--ref-cat REF_CAT]
                        [--ref-filter REF_FILTER] [--ref-file REF_FILE]
                        [--ref-image REF_IMAGE] [--ref-radius REF_RADIUS]
                        [--ana-thresh ANA_THRESH] [--ana-radius ANA_RADIUS]
                        [--ap-diam AP_DIAM [AP_DIAM ...]]
                        [--ap-diam-ul AP_DIAM_UL [AP_DIAM_UL ...]]
                        [--det-thresh DET_THRESH] [--corr-noise CORR_NOISE]
                        [--mag-cut MAG_CUT] [--mag-stdfaint MAG_STDFAINT]
                        [--mag-stdbright MAG_STDBRIGHT] [--tol TOL] [--auto]
                        [--bw] [--outdir OUTDIR]

Programme for aperture and forced photometry.

optional arguments:
  -h, --help            show this help message and exit
  --ra RA               RA(J2000) of the Object (HMS and DD system allowed).
                        Note: if declination is negative, write "
                        -12:20:20.2"). (Required)
  --dec DEC             Dec (J2000) of the Object (HMS and DD system allowed).
                        (Required)
  --fits FITS           Name of FITS file. (Required)
  --host-offset HOST_OFFSET
                        Host offset (unit: arcsec, default: 0.5)
  --ref-cat REF_CAT     Which reference catalogue should be used (SDSS, 2MASS,
                        PS1)?
  --ref-filter REF_FILTER
                        Filter of the refernce catalogue?
  --ref-file REF_FILE   Name of reference catalog file (if given overwrites '
                        --ref-cat' option)
  --ref-image REF_IMAGE
                        Name of reference image to run sextractor in dual
                        image mode?
  --ref-radius REF_RADIUS
                        Search radius in the reference catalogue query?
                        (default: 10)
  --ana-thresh ANA_THRESH
                        Analysis threshold. (default: 1 sigma)
  --ana-radius ANA_RADIUS
                        Analysis radius (unit: arcsec). (default: 90)
  --ap-diam AP_DIAM [AP_DIAM ...]
                        Aperture radii in FWHM (one or more). (default:
                        [1.,2.,3.,4.])
  --ap-diam-ul AP_DIAM_UL [AP_DIAM_UL ...]
                        Aperture radii in FWHM (one or more) (default: 2)
  --det-thresh DET_THRESH
                        Detection threshold. (default: 1 sigma)
  --corr-noise CORR_NOISE
                        Bootstrap photometric errors? (critical for correlated
                        noise, very time consuming). (default: False)
  --mag-cut MAG_CUT     Magnitude cut in ZP determination. Default: 12 mag
  --mag-stdfaint MAG_STDFAINT
                        Lower magnitude cut for secondary standards. (default:
                        0)
  --mag-stdbright MAG_STDBRIGHT
                        Upper magnitude cut for secondary standards. Default:
                        0
  --tol TOL             Tolerance of the cross-matching in arcsec. (default: 1)
                        arcsec
  --auto                Automatic mode? (default: False)
  --outdir OUTDIR       Output path. Default: 'results/'
```

#### How does it work?

The photometric catalogues are build from the 2MASS and SDSS point source catalogues. The PS1 point source catalogue was build following XXX. PS1 photometry was converted to the SDSS filters using the colour equations in XXX. Bessel photometry was derived following the [Lupton (2004)](http://classic.sdss.org/dr4/algorithms/sdssUBVRITransform.html) colour equations.

## Authors

* **Steve Schulze**

## How to reference?

The developement of these tools started in the paper [Schulze et al. (2018)](http://adsabs.harvard.edu/abs/2018MNRAS.473.1258S). If possible, please add a reference to this paper in your article.

Bibtex code

```
@ARTICLE{Schulze2018a,
   author = {{Schulze}, S. and {Kr{\"u}hler}, T. and {Leloudas}, G. and {Gorosabel}, J. and 
	{Mehner}, A. and {Buchner}, J. and {Kim}, S. and {Ibar}, E. and 
	{Amor{\'{\i}}n}, R. and {Herrero-Illana}, R. and {Anderson}, J.~P. and 
	{Bauer}, F.~E. and {Christensen}, L. and {de Pasquale}, M. and 
	{de Ugarte Postigo}, A. and {Gallazzi}, A. and {Hjorth}, J. and 
	{Morrell}, N. and {Malesani}, D. and {Sparre}, M. and {Stalder}, B. and 
	{Stark}, A.~A. and {Th{\"o}ne}, C.~C. and {Wheeler}, J.~C.},
    title = "{Cosmic evolution and metal aversion in superluminous supernova host galaxies}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1612.05978},
 keywords = {galaxies: evolution, galaxies: high-redshift, galaxies: luminosity function, mass function, galaxies: starburst, galaxies: star formation},
     year = 2018,
    month = jan,
   volume = 473,
    pages = {1258-1285},
      doi = {10.1093/mnras/stx2352},
   adsurl = {http://adsabs.harvard.edu/abs/2018MNRAS.473.1258S},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

## License

2-clause BSD