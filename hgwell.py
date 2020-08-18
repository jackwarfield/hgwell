from hgwell_f import *
import pandas as pd

import sys
# arguments
# filename val_threshold pix_threshold gaia_window
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
    gaia_window = 7.2

# open all SCI extensions into list, hold WCS info
im_list, wcs_list = open_images (filename)

x, y, v = [], [], []
for i,im in enumerate(im_list):
    print ("\n====================")
    print ("Finding objects in " + im.name + " " + str(i+1) + ".")
    star_out = find_stars (im, val_thr, pix_thr)
    x += [star_out[0]]
    y += [star_out[1]]
    v += [star_out[2]]
print ("====================")


df_list = []
for i,j,k in zip(x,y,v):
    df_list += [pd.DataFrame(data = {'x':i, 'y':j, 'v':k, 'g':np.zeros(len(i))*np.nan})]

for df,wcs in zip(df_list, wcs_list):
    df['g'] = match_gaia (df, wcs, gaia_window)

for i in range(len(im_list)):
    testplot (im_list[i], df_list[i], 'SCI'+str(i)+'.png')
    df_list[i].to_csv ('SCI'+str(i)+'.csv', index=False)
