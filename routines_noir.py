from      astropy import stats, table, time
from      dl import authClient as ac, queryClient as qc
from      dl.helpers.utils import convert
from      getpass import getpass

def query_noir_datalab_stars(RA, DEC, RADIUS = 16):

    """Queries the LS DR9 at the NOIR Datalab and retrieves stars within a given radius

    Input:
        RA:     Right ascension (unit: deg)
        DEC:    Declination (unit: deg)
        RADIUS: Radius (unit: arcmin)

    Returns:
        Table:  Photometry catalogue (columsn: ra, dec, mag_[g,r,z], snr_[g,r,z], allmask_[g,r,z])
    """

    query = '''
    SELECT

    ra, dec, mag_g, mag_r, mag_z, snr_g, snr_r, snr_z, allmask_g, allmask_r, allmask_z

    FROM ls_dr9.tractor AS tractor

    WHERE 't' = q3c_join({ra}, {dec}, tractor.ra, tractor.dec, {radius:.1f}/60.0)
    AND   tractor.type = 'PSF'

    '''.format(ra=RA, dec=DEC, radius=RADIUS)
       
    query   = query.replace('\n', ' ')
    query   = qc.query(sql=query, timeout=30)
    
    result = convert(query, 'pandas')
    result = table.Table.from_pandas(result)
    
    return result

def query_noir_datalab_extended(RA, DEC, RADIUS = 16):

    """Queries the LS DR9 at the NOIR Datalab and retrieves stars within a given radius

    Input:
        RA:     Right ascension (unit: deg)
        DEC:    Declination (unit: deg)
        RADIUS: Radius (unit: arcmin)

    Returns:
        Table:  Photometry catalogue (columsn: ra, dec, mag_[g,r,z], snr_[g,r,z], allmask_[g,r,z])
    """

    query = '''
    SELECT

    ra, dec, mag_g, mag_r, mag_z, snr_g, snr_r, snr_z, allmask_g, allmask_r, allmask_z

    FROM ls_dr9.tractor AS tractor

    WHERE 't' = q3c_join({ra}, {dec}, tractor.ra, tractor.dec, {radius:.1f}/60.0)
    AND   NOT tractor.type = 'PSF'

    '''.format(ra=RA, dec=DEC, radius=RADIUS)

    query   = query.replace('\n', ' ')
    query   = qc.query(sql=query, timeout=30)
    
    result = convert(query, 'pandas')
    result = table.Table.from_pandas(result)
    
    return result

