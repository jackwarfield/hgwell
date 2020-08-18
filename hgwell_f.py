from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

def open_images (filename):
    """
    Input
        filename:           filename for fits image
    Output
        im_list, wcs_list:  lists holding the SCI extensions and astropy.wcs objects from fits file
    """
    hdu_list = fits.open (filename)
    im_list, wcs_list = [], []
    for ext in hdu_list:
        if "SCI" in ext.name:
            im_list += [ext]
            wcs_list += [WCS(ext.header, hdu_list)] 
    return im_list, wcs_list

def find_stars (im, val_thr, pix_thr):
    """
    Input
        im:         fits image extension
        val_thr:    the minimum brightness to flag a pixel as a potential star
        pix_thr:    the size of the window used to sweep over flagged pixels and combine them as singular stars
    Output
        x,y,v:      lists containing the x pixel position, y pixel position, and mean brightness across combined
                    pixels for possible stars
    """
    data = im.data
    y_arr, x_arr, v_arr = [], [], []
    print ("Finding pixels brighter than %.2f." %val_thr)
    for i in range(len(data)):
        for j in range(len(data[0])):
            k = data[i][j]
            if k > val_thr:
                y_arr += [i]
                x_arr += [j]
                v_arr += [k]

    del_list = []
    print ("Locating weighted center of star with window of %.1f pixels." %pix_thr)
    for i in range(len(x_arr)):
        if i not in del_list:
            x1, y1, v1 = x_arr[i], y_arr[i], v_arr[i]
            xhold, yhold, vhold = [x1], [y1], [v1]
            for j in range(i+1, len(x_arr)):
                x2, y2, v2 = x_arr[j], y_arr[j], v_arr[j]
                if np.sqrt((xhold[-1]-x2)**2 + (yhold[-1]-y2)**2) < pix_thr:
                    xhold += [x2]
                    yhold += [y2]
                    vhold += [v2] # stronger weight to pull to the center of star
                    del_list += [j]
            x_arr[i] = np.average (xhold, weights=vhold)
            y_arr[i] = np.average (yhold, weights=vhold)
            v_arr[i] = np.average (vhold)

    x, y, v = [], [], []
    print ("Creating final table of object pixel positions and average values.")
    for i in range(len(x_arr)):
        if i not in del_list:
            x += [x_arr[i]]
            y += [y_arr[i]]
            v += [v_arr[i]]

    return x, y, v

def match_gaia (df, wcs, gaia_window):
    """
    Input
        df:             Pandas DataFrame with 'x', 'y', and 'g' columns for positions of possible stars and Gaia designations
        wcs:            astropy.wcs object for unit transformation from fits image pixel coordinates
        gaia_window:    the size of the window (in arcsec) for the Gaia query 
    Output
        df.g:           Pandas series of Gaia designations
    """
    from astroquery.gaia import Gaia
    window = gaia_window * u.arcsec
    for i in range(len(df)):
        x,y = df.loc[i,'x'],df.loc[i,'y'] 
        ra,dec = wcs.pixel_to_world_values (x,y)
        coord = SkyCoord (ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
        r = Gaia.query_object_async (coordinate=coord, width=window, height=window)
        try:
            df.loc[i,'g'] = r['designation'][0].decode('utf-8')
        except:
            pass
    return df.g

def testplot (im, df, out):
    """
    Input
        im:     fits image extension
        df:     Pandas DataFrame with 'x', 'y', and 'g' columns for positions of possible stars and Gaia designations
        out:    output filename
    Output
        void
    """
    import matplotlib.pyplot as plt
    height, width = len(im.data), len(im.data[0])
    fig = plt.figure (figsize=(20, 20*height/width))
    ax = fig.add_subplot (1,1,1)
    ax.imshow (im.data, cmap='gray', vmin=0, origin='lower')
    ax.scatter (df.x, df.y, s=20, color='none', edgecolor='red')
    ax.scatter (df[df.g.notna()].x, df[df.g.notna()].y, s=20, color='none', edgecolor='green')
    for x,y,z in zip(df[df.g.notna()].x, df[df.g.notna()].y, df[df.g.notna()].g):
        ax.text (x, y-60, s=z.replace ('Gaia DR2 ','Gaia DR2\n'), horizontalalignment='center', verticalalignment='center', color='green')
    fig.savefig (out, dpi=300, bbox_inches='tight')
    return
