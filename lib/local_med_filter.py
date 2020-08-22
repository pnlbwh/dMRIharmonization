# ===============================================================================
# dMRIharmonization (2018) pipeline is written by-
#
# TASHRIF BILLAH
# Brigham and Women's Hospital/Harvard Medical School
# tbillah@bwh.harvard.edu, tashrifbillah@gmail.com
#
# ===============================================================================
# See details at https://github.com/pnlbwh/dMRIharmonization
# Submit issues at https://github.com/pnlbwh/dMRIharmonization/issues
# View LICENSE at https://github.com/pnlbwh/dMRIharmonization/blob/master/LICENSE
# ===============================================================================

import numpy as np


def local_med_filter(img, outLier, w=2):
    # img is the zero padded image
    # outLier is the mask over img>np.percentile(img,95)

    img_med = img.copy()

    indx, indy, indz = np.where(outLier > 0)

    for (ix, iy, iz) in zip(indx, indy, indz):
        neighborhood = img[ix - w:ix + w, iy - w:iy + w, iz - w:iz + w]
        img_med[ix, iy, iz] = np.median(neighborhood)

    return img_med
