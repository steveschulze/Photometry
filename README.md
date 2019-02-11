# Photometry tools

Tools to align images, improve astrometry, perform aperture photometry and retrieve photometry catalogues.

## Prerequisites

Python: 3.6

Additional python libraries:

If your python installation is based on [astroconda](https://astroconda.readthedocs.io/en/latest/), you need to install these additional libraries:

```
astroquery, photutils, pysynphot, sympy, sewpy, and sip_tpv (included)
```

Futhermore, you need to install install these software packages

```
astrometry.net, Scamp, SExtractor and wcs-tools
```

## Installing

Store the directory in your favourite folder and add the folder to your PATH and PYTHON_PATH environment.

```
export PATH=$PATH:/new/path
export PYTHON_PATH=$PYTHON_PATH:/new/path
```

Give field_calibration.py, improve_astro.py, photometry.py, photometry_hst.py permission to be executed.

```
chmod u+x PROGRAMME.py
```

-To install Sewpy follow the instructions at [https://github.com/megalut/sewpy](https://github.com/megalut/sewpy).

-The first line in files field_calibration.py, improve_astro.py, photometry.py, photometry_hst.py is ```#!/usr/bin/env pythonw```. If you use linux, change ```pythonw``` to ```python```. With ```pythonw``` you can close plot windows with the keystroke &#8984;+w.

## Description of the tools

### align_images.py

#### How to call it?

Requires: ```Scamp```, ```SExtractor```

```
usage: align_images.py [-h] --ref-image REF_IMAGE --new-image NEW_IMAGE
                       [--keep-temp]

Align the WCS systems of two images

optional arguments:
  -h, --help            show this help message and exit
  --ref-image REF_IMAGE Reference image
  --new-image NEW_IMAGE New image
  --keep-temp           Keep temporary files
```

#### Example

```
python align_images.fits --ref-image SN2015bn_SDSS_r.fits --new-image
```

#### How does it work?

1) Identify objects in the reference image
2) Identify objects in the new image
3) Align the astrometry. Information about the quality of the alignment is saved in 'results_scamp/scamp.log'. The scamp log file is stored in 'results_scamp/{NEW_IMAGE}_scamp.xml'
4) Remove temporary files

#### Important note

The programm can fail if the images are very different from each other, e.g., a shallow transient image and an HST image.
If your browser does not display the scamp diagnostic file, reach these [notes](https://www.astromatic.net/2009/10/05/understanding-astromatic-metadata-files).

### improve_astrometry.py

#### How to call it?

Requires: ```sympy```, ```SExtractor```, ```sip_tpv```

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

#### Example

```
python improve_astro.py --fits SN2015bn_SDSS_r.fits --ra 173.423125 --dec 0.725972
```

#### How does it work?

It uses astrometry.net and sextractor to establish the field calibration of the image. Afterwards it uses sip2pv.py to transform the distortion keywords into the format sextractor understands. The new image has the suffix '_wcs'.

#### Important note

If you built ```astrometry.net``` with python 2.7, running ```astrometry.net``` in python 3 will possibly fail. There are two work arounds. 1) ```improve_astrometry.py``` is compatible with ```python 2.7```. 2) To run ```improve_astrometry.py``` in python 3, activate your python-3 environment and re-install ```astrometry.net```.

```astrometry.net``` requires an image with a size of at least 5x5 arcmin.

### field_calibration.py

#### How to call it?

Requires: ```astroquery```, ```sewpy```

(Also works in python 2.)

```
usage: field_calibration.py [-h] --ra RA --dec DEC [--radius RADIUS]
                            [--outdir OUTDIR] [--type TYPE]

Retrieve photometric catalogues. 2MASS, PS1 and SDSS are the input sources.
Bessel catalogues are generated through colour equations. PS1 and SDSS cats
are in the AB system, whereas Bessel and 2MASS cats are in the Vega system.

optional arguments:
  -h, --help       show this help message and exit
  --ra RA          RA(J2000) of the Object (HMS and DD system allowed). Note:
                   if declination is negative, write " -12:20:20.2"). (Required)
                   keyword.
  --dec DEC        Dec (J2000) of the Object (HMS and DD system allowed).
                   (Required)
  --radius RADIUS  only search in indexes within 'radius' of the field center
                   (unit: arcmin; default: 10)
  --outdir OUTDIR  Output directory (default: results/)
  --type TYPE      Generate catalogues for optical/nir/all (default: all)
```

### Example

```
python field_calibration.py --ra 173.423125 --dec 0.725972
```

#### How does it work?

The photometric catalogues are build from the 2MASS point source catalogues and SDSS/DR12 source catalogues. The PS1 point source catalogue was build following XXX. PS1 photometry was converted to the SDSS filters using the colour equations in XXX. Bessel photometry was derived following the [Lupton (2004)](http://classic.sdss.org/dr4/algorithms/sdssUBVRITransform.html) colour equations.

The files will be called ```PS1_PS1_FILTER.ascii```, ```PS1_SDSS_FILTER.ascii```, ```PS1_BESSEL_FILTER.ascii```, ```SDSS_FILTER.ascii```, ```SDSS_BESSEL_FILTER.ascii``` and ```2MASS_FILTER.ascii```. The columns in each file are: col1 = ra, col2 = dec, col3 = mag, and col4 = sigma_mag.

### photometry.py

#### How to call it?

Requires: ```astroquery```, ```phot_utils```, ```sewpy```, ```SExtractor```, ```wcs-tools```

```
usage: photometry.py [-h] --ra RA --dec DEC --fits FITS
                     [--host-offset HOST_OFFSET] [--ref-cat REF_CAT]
                     [--ref-filter REF_FILTER] [--ref-file REF_FILE]
                     [--ref-image REF_IMAGE] [--ref-radius REF_RADIUS]
                     [--ana-thresh ANA_THRESH]
                     [--ap-diam AP_DIAM [AP_DIAM ...]]
                     [--ap-diam-ul AP_DIAM_UL [AP_DIAM_UL ...]]
                     [--det-thresh DET_THRESH] [--gain GAIN]
                     [--back-size BACK_SIZE]
                     [--back-filtersize BACK_FILTERSIZE]
                     [--deblend-nthresh DEBLEND_NTHRESH]
                     [--deblend-mincont DEBLEND_MINCONT] [--mag-cut MAG_CUT]
                     [--mag-stdfaint MAG_STDFAINT]
                     [--mag-stdbright MAG_STDBRIGHT] [--maxstars MAXSTARS]
                     [--auto] [--bw] [--keeptemp] [--loglevel LOGLEVEL]
                     [--outdir OUTDIR] [--sex-loglevel SEX_LOGLEVEL]
                     [--tol TOL]

Programme for aperture photometry.

optional arguments:
  -h, --help            show this help message and exit
  --ra RA               RA(J2000) of the Object (HMS and DD system allowed).
                        Note: if declination is negative, write "
                        -12:20:20.2"). (required)
  --dec DEC             Dec (J2000) of the Object (HMS and DD system allowed)
                        (required)
  --fits FITS           Name of FITS file (required)
  --host-offset HOST_OFFSET
                        Host offset (default: 10 arcsec)
  --ref-cat REF_CAT     Which reference catalogue should be used (SDSS, 2MASS,
                        PS1)?
  --ref-filter REF_FILTER
                        Filter of the reference catalogue?
  --ref-file REF_FILE   Name of reference catalog file (if given overwrites '
                        --ref-cat' option)
  --ref-image REF_IMAGE
                        Name of reference image to run sextractor in dual
                        image mode?
  --ref-radius REF_RADIUS
                        Search radius in the reference catalogue query?
                        (default: 10, unit: arcmin)
  --ana-thresh ANA_THRESH
                        Analysis threshold (default: 1 sigma)
  --ap-diam AP_DIAM [AP_DIAM ...]
                        Aperture radii in FWHM (one or more) (default:
                        [1.,2.,3.,4.])
  --ap-diam-ul AP_DIAM_UL [AP_DIAM_UL ...]
                        Aperture radii in FWHM (one or more) (default: 2)
  --det-thresh DET_THRESH
                        Detection threshold (default: 1; unit: sigma)
  --gain GAIN           Gain keyword in the header (default: Gain; unit:
                        e-/ADU)
  --back-size BACK_SIZE
                        Background mesh: <size> (default: 64)
  --back-filtersize BACK_FILTERSIZE
                        Background filter: <size> (default: 3)
  --deblend-nthresh DEBLEND_NTHRESH
                        Number of deblending sub-thresholds (default: 64)
  --deblend-mincont DEBLEND_MINCONT
                        Minimum contrast parameter for deblending (default:
                        0.00001)
  --mag-cut MAG_CUT     Remove all objects brighter than a given magnitude
                        (default: 12)
  --mag-stdfaint MAG_STDFAINT
                        Lower magnitude cut for secondary standards (default:
                        0)
  --mag-stdbright MAG_STDBRIGHT
                        Upper magnitude cut for secondary standards (default:
                        0)
  --maxstars MAXSTARS   Maximum number of stars used for building the local
                        sequence (default: 200)
  --auto                Automatic mode? (default: False)
  --bw                  Text output in color (default: False)
  --keeptemp            Keep temporary files
  --loglevel LOGLEVEL   Logger level (default: INFO, possible values: DEBUG,
                        INFO, WARNING, ERROR, CRITICAL)
  --outdir OUTDIR       Output path efault: 'results/'
  --sex-loglevel SEX_LOGLEVEL
                        Sextractor logger level (default: WARNING, possible
                        values: DEBUG, INFO, WARNING, ERROR, CRITICAL)
  --tol TOL             Tolerance of the cross-matching in arcsec (default: 1)
```

#### Example

Step 1:

Go to the directory that contains the fits file

Step 2: You have no local star catalogue to constrain the flux scale of the image

```photometry.py --fits SN2015bn_SDSS_r.fits --ra 173.423125 --dec 0.725972 --ref-cat SDSS --ref-filter r```

or you would like to use a local star catalogue

```photometry.py --fits SN2015bn_SDSS_r.fits --ra 173.423125 --dec 0.725972 --ref-file results/SDSS_r.cat```

#### How does it work?

1) Select stars for the local sequence (either from the user catalogue or it downloads a catalogue from the VizieR database)
2) Measure the zeropoints for the different apertures
3) Measure the brightness of the science object and all other sources in the images. If no credible source was detected close to the specified coordinates, the programme will perform forced photometry at the specified coordinates.
4) Removing temporary files

You can find more information at [SN2015bn_SDSS_r_summary.pdf](Example/SN2015bn_SDSS_r_summary.pdf).

The programme generates several text files

| File  | Explanation |
| ----- | ------- |
| ./SN2015bn_SDSS_r.log                           | General log file|
| ./results/SN2015bn_SDSS_r_SDSS_r.cat            | SDSS point source catalogue| 
| ./results/SN2015bn_SDSS_r_all.log               | SExtractor log| 
| ./results/SN2015bn_SDSS_r_all_abs_cal.phot      | SExtractor output catalogue<br> all magnitudes are calibrated with the zeropoints measured from the image| 
| ./results/SN2015bn_SDSS_r_phot.log              | Summary of the photometry on the science source |
| ./results/SN2015bn_SDSS_r_zp.log                | Zeropoint summary |

and diagnostic plots

| File  | Explanation |
| ----- | ------- |
|./results/SN2015bn_SDSS_r_fwhm.pdf               | Distribution of the FWHMs of the stellar PSFs |
|./results/SN2015bn_SDSS_r_poststamp.pdf          | Poststage stamp of the source of interest |
|./results/SN2015bn_SDSS_r_std.pdf                | Local sequence: instrumental magnitude vs. tabulated magnitudes |
|./results/SN2015bn_SDSS_r_zp.pdf                 | Measured zeropoints for the different apertures <br> Datapoints: used in black; not used in grey <br> Lines: median solid, 1-sigma confidence interval dashed sigma-clipped region dotted |

#### Limitations

-The BW is currently not fully implemented.

-If you want to define different circular apertures, make sure that you always define 4.

-The dual-image mode is currently not activated.

### photometry_hst.py

#### How to call it?

Requires: ```phot_utils```, ```sewpy```, ```SExtractor```, ```pysynphot```

```
usage: photometry_hst.py [-h] --ra RA --dec DEC --fits FITS
                         [--ana-thresh ANA_THRESH] [--det-thresh DET_THRESH]
                         [--ap-diam AP_DIAM [AP_DIAM ...]]
                         [--ap-inner-annulus AP_INNER_ANNULUS]
                         [--ap-outer-annulus AP_OUTER_ANNULUS] [--auto] [--bw]
                         [--centroid] [--keeptemp] [--loglevel LOGLEVEL]
                         [--outdir OUTDIR] [--sex-loglevel SEX_LOGLEVEL]
                         [--tol TOL]

Programme for aperture photometry of HST images.

optional arguments:
  -h, --help            show this help message and exit
  --ra RA               RA(J2000) of the Object (HMS and DD system allowed).
                        Note: if declination is negative, write "
                        -12:20:20.2"). Required keyword.
  --dec DEC             Dec (J2000) of the Object (HMS and DD system allowed).
                        Required keyword.
  --fits FITS           Name of FITS file. Required keyword.
  --ana-thresh ANA_THRESH
                        Analysis threshold. Default: 3 sigma
  --det-thresh DET_THRESH
                        Detection threshold. Default: 3 sigma
  --ap-diam AP_DIAM [AP_DIAM ...]
                        Aperture radii in arcsec (one or more value allowed).
                        Default: [0.25, 0.5, 0.75, 1.00]
  --ap-inner-annulus AP_INNER_ANNULUS
                        Diameter of the inner annulus of the background.
                        Default: 1.5 x ap-diam
  --ap-outer-annulus AP_OUTER_ANNULUS
                        Diameter of the outer annulus of the background.
                        Default: 3 x ap-diam
  --auto                Automatic mode? Default: False
  --bw                  Screen output in B/W? Default: False
  --centroid            Centre on the most nearby object
  --keeptemp            Keep temporary files
  --loglevel LOGLEVEL   Logger level (default: INFO, possible values: DEBUG,
                        INFO, WARNING, ERROR, CRITICAL)
  --outdir OUTDIR       Output path. Default: 'results/'
  --sex-loglevel SEX_LOGLEVEL
                        Sextractor logger level (default: WARNING, possible
                        values: DEBUG, INFO, WARNING, ERROR, CRITICAL)
  --tol TOL             Tolerance of the cross-matching in arcsec. Default: 1
                        arcsec
```

#### Example

Step 1:

Go to the directory that contains the fits file

Step 2: 

```pythonw photometry_hst.py --ra 173.42306 --dec 0.72589793 --fits SN2015bn_F625W_drc.fits ```

#### How does it work?

1) Compute the ZP from the header keywords and applies aperture corrections for each aperture (using pysynphot)
2) Generate a general source catalogue to identify the object that is closests to your science object
3) Measure the brightness for the specified apertures
4) Perform curve of growth analysis
5) Removing temporary files

You can find more information at [SN2015bn_F625W_drc.pdf](Example/SN2015bn_F625W_drc.pdf).


| File  | Explanation |
| ----- | -------     |
| ./SN2015bn_F625W_drc.log                           | General log file |
| ./results/SN2015bn_F625W_drc.mag                   | Photometry log of the science object |
| ./results/SN2015bn_F625W_drc_centroid.log          | SExtractor output log |
| ./results/SN2015bn_F625W_drc_centroid.phot         | SExtractor output catalogue | 
| ./results/SN2015bn_F625W_drc_cog_stat.ascii        | Text file of the curve-of-growth analysis |
| ./results/SN2015bn_F625W_drc_phot.log              | Summary of the photometry on the science source |
| ./results/SN2015bn_F625W_drc_zp.log                | Zeropoint summary |

and diagnostic plots

| File  | Explanation |
| ----- | -------     |
| ./results/SN2015bn_F625W_drc_cog.pdf               | Curve of growth|
| ./results/SN2015bn_F625W_drc.pdf                   | Poststage stamp of the source of interest|

## Can you speed up the execution speed?

That's very easy. You can run each tool in parallel. For example

```
ls "*fits | xargs -n 1 -P 16 -I {} photometry.py --fits {} --ra 173.423125 --dec 0.725972 --ref-cat SDSS --ref-filter r
```

would run 16 parallal sessions of photometry.py.

## Authors

* **Steve Schulze**

## How to reference?

The developement of these tools started in the paper [Schulze et al. (2018)](http://adsabs.harvard.edu/abs/2018MNRAS.473.1258S). If possible, please add a reference to this paper in your article. 

Also cite the relevant papers for [astrometry.net](https://github.com/dstndstn/astrometry.net), SExtractor and
[sip_tpv](https://github.com/stargaser/sip_tpv).

Bibtex codes

```
@ARTICLE{Berin1996A,
   author = {{Bertin}, E. and {Arnouts}, S.},
    title = "{SExtractor: Software for source extraction.}",
  journal = {\aaps},
 keywords = {METHODS: DATA ANALYSIS, TECHNIQUES: IMAGE PROCESSING, GALAXIES: PHOTOMETRY},
     year = 1996,
    month = jun,
   volume = 117,
    pages = {393},
      doi = {10.1051/aas:1996164},
   adsurl = {http://adsabs.harvard.edu/abs/1996A%26AS..117..393B},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

```
@ARTICLE{Lang2010a,
   author = {{Lang}, D. and {Hogg}, D.~W. and {Mierle}, K. and {Blanton}, M. and 
	{Roweis}, S.},
    title = "{Astrometry.net: Blind Astrometric Calibration of Arbitrary Astronomical Images}",
  journal = {\aj},
archivePrefix = "arXiv",
   eprint = {0910.2233},
 primaryClass = "astro-ph.IM",
 keywords = {astrometry, catalogs, instrumentation: miscellaneous, methods: data analysis, methods: statistical, techniques: image processing},
     year = 2010,
    month = may,
   volume = 139,
    pages = {1782},
      doi = {10.1088/0004-6256/139/5/1782},
   adsurl = {http://adsabs.harvard.edu/abs/2010AJ....139.1782L},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

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
    pages = {1258},
      doi = {10.1093/mnras/stx2352},
   adsurl = {http://adsabs.harvard.edu/abs/2018MNRAS.473.1258S},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

```
@INPROCEEDINGS{Shupe2012a,
   author = {{Shupe}, D.~L. and {Laher}, R.~R. and {Storrie-Lombardi}, L. and 
	{Surace}, J. and {Grillmair}, C. and {Levitan}, D. and {Sesar}, B.
	},
    title = "{More flexibility in representing geometric distortion in astronomical images}",
booktitle = {Software and Cyberinfrastructure for Astronomy II},
     year = 2012,
   series = {\procspie},
   volume = 8451,
    month = sep,
      eid = {84511M},
    pages = {84511M},
      doi = {10.1117/12.925460},
   adsurl = {http://adsabs.harvard.edu/abs/2012SPIE.8451E..1MS},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## License

2-clause BSD