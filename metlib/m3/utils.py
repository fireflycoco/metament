from metlib.m3.m3_file import *
from mpl_toolkits.basemap import Basemap

def get_m3_map(m3_file_obj, border_thick=0.5, dot=False):
    """Args: m3_file_obj: a metlib.m3.m3_file.M3_File object 
        border_thick: extra grids around the data area for the map.
            0.5 is proper for pcolor type plots,
            0.0 is proper for contour type plots
        dot: dot grid (True) or cross grid (False)
    """
    # calculate corners
    imin = 0. - border_thick
    imax = m3_file_obj.NCOLS - 1. + border_thick
    jmin = 0. - border_thick
    jmax = m3_file_obj.NROWS - 1. + border_thick

    ll_lon, ll_lat = m3_file_obj.ij_to_lonlat(imin, jmin, dot=dot)
    ur_lon, ur_lat = m3_file_obj.ij_to_lonlat(imax, jmax, dot=dot)
    # calculate proper resolution
    if m3_file_obj.XCELL >= 27000.0:
        proper_reso = 'i'
    elif m3_file_obj.XCELL <= 3000.0:
        proper_reso = 'f'
    else:
        proper_reso = 'h'

    # setup a map
    m = Basemap(projection='lcc', 
            lat_0 = m3_file_obj.YCENT,
            lon_0 = m3_file_obj.XCENT,
            lat_1 = m3_file_obj.P_ALP,
            lat_2 = m3_file_obj.P_BET,
            llcrnrlat = ll_lat,
            llcrnrlon = ll_lon,
            urcrnrlat = ur_lat,
            urcrnrlon = ur_lon,
            resolution = proper_reso,
            rsphere=6370000.0
            )
    # return the map
    return m
