from hgwell_f import *
import pandas as pd

import sys
# arguments
# filename val_threshold pix_threshold gaia_window(arcsec)
filename = sys.argv[1]
try:
    val_thr = float(sys.argv[2])
except:
    val_thr = 5.5e4
try:
    pix_thr = float(sys.argv[3])
except:
    pix_thr = 20
try:
    gaia_window = float(sys.argv[4])
except:
    gaia_window = 1.0

# open all SCI extensions into list, hold WCS info
im_list, wcs_list, date = open_images (filename)

x, y, v = [], [], []
print ("\n")
for i,im in enumerate(im_list):
    print ("====================")
    print ("Finding objects in " + im.name + " " + str(i+1) + " . . .")
    star_out = find_stars (im, val_thr, pix_thr)
    x += [star_out[0]]
    y += [star_out[1]]
    v += [star_out[2]]

print ("====================")
print ("Querying the Gaia archive via astroquery using list of %g objects . . ." %(len(x)))

df_list = []
for i,j,k in zip(x,y,v):
    df_list += [pd.DataFrame(data = {'x':i, 'y':j, 'v':k})]

# currently is very slow. will change to get all of the designation ids, and then
# do the epoch_prop in a single query, which I think will speed everything up
for df,wcs in zip(df_list, wcs_list):
    df = match_gaia (df, wcs, gaia_window, date)

print ("====================")
print ("Outputting tables and plots of Gaia matches!")

for i in range(len(im_list)):
    testplot (im_list[i], df_list[i], 'SCI'+str(i+1)+'.png')
    df_list[i].to_csv ('SCI'+str(i+1)+'.csv', index=False)
