from wds_astropy_table import parse_wdsweb
from wds_astropy_table import parse_sixth_orbital

from astropy.coordinates import SkyCoord
from astropy.time import Time

from astropy.table import Table

from astroquery.simbad import Simbad

import astropy.units as u
import numpy as np

import re
import sys

table_wds = parse_wdsweb('wdsweb_summ2.txt')
#table_orbital = parse_sixth_orbital('orb6orbits.txt')

# a list of (discoverer, components) pairs
targets = [
    ("STF 518", "BC"),
    ("AG 342", "AB"),
    ("STF 1645", "AB"),
    ("DUN 5", "AB")
]

def normalize(discoverer, components):
    return discoverer.replace(" ", "")+components

def lookup_wds(discoverer, components):
    normalized=normalize(discoverer, components)

    return table_wds[(table_wds['discoverer_normalized']==normalized) & (table_wds['components']==components)][0]

def extract_coords(wds_record):
    # "wds coordinate string"
    wcs = wds_record['j2000']

    # parse into format recognized by skycoord
    formatted_coord_string = wcs[0:2]+'h'+wcs[2:4]+'m'+wcs[4:9]+'s' + ' ' + wcs[9:12]+'d'+wcs[12:14]+'m'+wcs[14:18]+'s'
    
    # get the position angle and offset of secondary
    pa = wds_record['last_pa']*u.deg
    sep = wds_record['last_sep']*u.arcsec

    # create skycoord objects
    c_pri = SkyCoord(formatted_coord_string, unit=(u.hourangle, u.deg), obstime=Time('J2000.0'))
    c_sec = c_pri.directional_offset_by(pa,sep)

    return (c_pri, c_sec)

# extract the Gaia DR3 id from the Simbad ids string
def get_gaia_dr3_id(ids):
    regex = 'GAIA DR3 (\\d+)'
    result = re.search(regex, ids, re.IGNORECASE)

    return result.group(1) 


# lookup the wds records
wds_records = [lookup_wds(disc, comp) for (disc, comp) in targets]

#extract the coordinates
coords = [extract_coords(wds_record) for wds_record in wds_records]

# do a Simbad cone search for the primaries and secondaries



search_radius = 3*u.arcsec

Simbad.add_votable_fields('ids')
results_pri = Simbad.query_region(c_pri_j2000, radius=search_radius)
results_sec = Simbad.query_region(c_sec_j2000, radius=search_radius)

gaia_id_pri=get_gaia_dr3_id(results_pri['ids'][0])
gaia_id_sec=get_gaia_dr3_id(results_sec['ids'][0])
print(f"gaia_id_pri={gaia_id_pri}")
print(f"gaia_id_sec={gaia_id_sec}")

display(results_pri)
display(results_sec)
#results_t = Table()
#results_t = vstack([results_t, results])

print(c_pri_j2000.to_string())


