import numpy as np

def local_med_filter(img, outLier, w=3):
    # img is the zero padded image
    # outLier is the mask over img>np.percentile(img,95)

    img_med= img.copy()

    indx, indy, indz= np.where(outLier>0)

    for (ix, iy, iz) in zip(indx, indy, indz):
        neighborhood= img[ix-w:ix+w, iy-w:iy+w, iz-w:iz+w]
        img_med[ix,iy,iz]= np.median(neighborhood)

    return img_med