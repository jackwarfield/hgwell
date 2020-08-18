# hgwell
## <i>"Combine HST and Gaia, and do it Well!"</i>

This script takes in a `fits` image, attempts to locate stars, and looks for <i>Gaia</i> crossmatches at those locations.

`python3 hgwell.py <image filename> <star brightness threshold> <pixel window size> <Gaia window size>`

There is still much that needs to be improved!
Things immediately being worked on:
\[ ] better star-finding using PSF fit
