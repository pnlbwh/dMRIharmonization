#!/usr/bin/env python

import numpy as np
import pyfftw

pyfftw.config.NUM_THREADS = -1
pyfftw.config.PLANNER_EFFORT= 'FFTW_ESTIMATE'

PI= np.float64(np.pi)
fftw_complex= 'complex128'
eps= 2.204e-16


def unring_1D(data, n, numlines, nsh, minW, maxW):

    inn = pyfftw.empty_aligned(n, dtype= fftw_complex)
    out = inn.copy()

    p = pyfftw.builders.fft(inn, n)
    p_inv = pyfftw.builders.ifft(inn, n)

    sh = pyfftw.empty_aligned(n *(2*nsh+1), dtype= fftw_complex)
    sh2 = sh.copy()
    
    nfac = 1./n
    
    shifts = np.zeros(2*nsh+1, dtype= 'int')
    shifts[0] = 0
    for j in range(nsh):
        shifts[j+1] = j+1
        shifts[1+nsh+j] = -(j+1)
    
    TV1arr= np.zeros(2*nsh+1, dtype= 'float64')
    TV2arr= TV1arr.copy()
    
    for k in range(numlines):
        p(data[k*n:(k+1)*n], sh[ :n])

        maxn = int(n/2-1 if n % 2 else (n-1)/2)
        
        for j in range(1, 2*nsh+1):

            phi = PI/n * shifts[j]/nsh
            u = [np.cos(phi), np.sin(phi)]
            e = [1, 0]

            sh[j*n] = sh[0]

            if not n%2:
                sh[j*n + n//2] = 0

            
            for l in range(maxn):

                tmp= np.float64(e[0])
                e[0] = u[0]*e[0] - u[1]*e[1]
                e[1] = u[1]*tmp + u[0]*e[1]

                L = l+1
                sh[j*n +L] = complex(e[0]*sh[L].real - e[1]*sh[L].imag, e[0]*sh[L].imag + e[1]*sh[L].real)

                L = n-1-l
                sh[j*n +L] = complex(e[0]*sh[L].real + e[1]*sh[L].imag, e[0]*sh[L].imag - e[1]*sh[L].real)

        
                
        for j in range(2*nsh+1):
            p_inv(sh[j*n:(j+1)*n], sh2[j*n:(j+1)*n])

            for j in range(2*nsh+1):

                TV1arr[j] = 0
                TV2arr[j] = 0
                l= 0

                for t in range(minW, maxW):

                    TV1arr[j] += abs(sh2[j*n + (l-t+n)%n ].real - sh2[j*n + (l-(t+1)+n)%n ].real)
                    TV1arr[j] += abs(sh2[j*n + (l-t+n)%n ].imag - sh2[j*n + (l-(t+1)+n)%n ].imag)
                    TV2arr[j] += abs(sh2[j*n + (l+t+n)%n ].real - sh2[j*n + (l+(t+1)+n)%n ].real)
                    TV2arr[j] += abs(sh2[j*n + (l+t+n)%n ].imag - sh2[j*n + (l+(t+1)+n)%n ].imag)

        
        
        for l in range(n):

            minTV = np.float64(999999999999)
            minidx= 0

            for j in range(2*nsh+1):
                
                if TV1arr[j] < minTV:
                    minTV = TV1arr[j]
                    minidx = j

                if TV2arr[j] < minTV:
                    minTV = TV2arr[j]
                    minidx = j
                
                TV1arr[j] += abs(sh2[j*n + (l-minW+1+n)%n ].real - sh2[j*n + (l-(minW)+n)%n ].real)
                TV1arr[j] -= abs(sh2[j*n + (l-maxW+n)%n ].real - sh2[j*n + (l-(maxW+1)+n)%n ].real)
                TV2arr[j] += abs(sh2[j*n + (l+maxW+1+n)%n ].real - sh2[j*n + (l+(maxW+2)+n)%n ].real)
                TV2arr[j] -= abs(sh2[j*n + (l+minW+n)%n ].real - sh2[j*n + (l+(minW+1)+n)%n ].real)
                
                TV1arr[j] += abs(sh2[j*n + (l-minW+1+n)%n ].imag - sh2[j*n + (l-(minW)+n)%n ].imag)
                TV1arr[j] -= abs(sh2[j*n + (l-maxW+n)%n ].imag - sh2[j*n + (l-(maxW+1)+n)%n ].imag)
                TV2arr[j] += abs(sh2[j*n + (l+maxW+1+n)%n ].imag - sh2[j*n + (l+(maxW+2)+n)%n ].imag)
                TV2arr[j] -= abs(sh2[j*n + (l+minW+n)%n ].imag - sh2[j*n + (l+(minW+1)+n)%n ].imag)
            

            a0r = np.float64(sh2[minidx*n + (l-1+n)%n ].real)
            a1r = np.float64(sh2[minidx*n + l].real)
            a2r = np.float64(sh2[minidx*n + (l+1+n)%n ].real)
            a0i = np.float64(sh2[minidx*n + (l-1+n)%n ].imag)
            a1i = np.float64(sh2[minidx*n + l].imag)
            a2i = np.float64(sh2[minidx*n + (l+1+n)%n ].imag)
            s   = np.float64(shifts[minidx])/nsh/2
            
            if s>0:
                data[k*n + l] =  complex((a1r*(1-s) + a0r*s)*nfac, (a1i*(1-s) + a0i*s)*nfac)

            else:
                s = -s
                data[k*n + l] =  complex((a1r*(1-s) + a2r*s)*nfac, (a1i*(1-s) + a2i*s)*nfac)


    return data



def unring_2d(data1, tmp2, dim_sz, nsh, minW, maxW):

    tmp1 = pyfftw.empty_aligned(dim_sz[0]*dim_sz[1], dtype= fftw_complex)
    data2 = tmp1.copy()

    p = pyfftw.builders.fft2(data1.reshape([dim_sz[1],dim_sz[0]]))
    p_inv = pyfftw.builders.ifft2(data1.reshape([dim_sz[1],dim_sz[0]]))

    p_tr = pyfftw.builders.fft2(data2.reshape([dim_sz[0],dim_sz[1]]))
    pinv_tr = pyfftw.builders.ifft2(data2.reshape([dim_sz[0],dim_sz[1]]))

    # p = pyfftw.FFTW(data1.reshape([dim_sz[1],dim_sz[0]]), data1.reshape([dim_sz[1],dim_sz[0]]),
    #                              direction= 'FFTW_FORWARD', flags= ('FFTW_ESTIMATE,'))
    # p_inv = pyfftw.FFTW(data1.reshape([dim_sz[1],dim_sz[0]]), tmp1.reshape([dim_sz[1],dim_sz[0]]),
    #                              direction= 'FFTW_BACKWARD', flags= ('FFTW_ESTIMATE,'))

    nfac = 1/np.float64(dim_sz[0]*dim_sz[1])

    for k in range(dim_sz[1]):
       for j in range(dim_sz[0]):
            data2[j*dim_sz[1]+k] = data1[k*dim_sz[0]+j]

    p(data1.reshape([dim_sz[1],dim_sz[0]]), tmp1.reshape([dim_sz[1],dim_sz[0]]))
    p_tr(data2.reshape([dim_sz[0],dim_sz[1]]), tmp2.reshape([dim_sz[0],dim_sz[1]]))

    for k in range(dim_sz[1]):
        ck = (1+np.cos(2*PI*(k/dim_sz[1])))*0.5 +eps
        for j in range(dim_sz[0]):
            cj = (1+np.cos(2*PI*j/dim_sz[0]))*0.5 +eps
            tmp1[k*dim_sz[0]+j] = nfac*(tmp1[k*dim_sz[0]+j] * ck) / (ck+cj)
            tmp2[j*dim_sz[1]+k] = nfac*(tmp2[j*dim_sz[1]+k] * cj) / (ck+cj)


    p_inv(tmp1.reshape([dim_sz[1],dim_sz[0]]), data1.reshape([dim_sz[1],dim_sz[0]]))
    pinv_tr(tmp2.reshape([dim_sz[0],dim_sz[1]]), data2.reshape([dim_sz[0],dim_sz[1]]))

    data1= unring_1D(data1,dim_sz[0],dim_sz[1],nsh,minW,maxW)
    data2= unring_1D(data2,dim_sz[1],dim_sz[0],nsh,minW,maxW)

    p(data1.reshape([dim_sz[1],dim_sz[0]]), tmp1.reshape([dim_sz[1],dim_sz[0]]))
    p_tr(data2.reshape([dim_sz[0],dim_sz[1]]), tmp2.reshape([dim_sz[0],dim_sz[1]]))

    for k in range(dim_sz[1]):
        for j in range(dim_sz[0]):
            tmp1[k*dim_sz[0]+j] = nfac*(tmp1[k*dim_sz[0]+j]  + tmp2[j*dim_sz[1]+k])


    p_inv(tmp1.reshape([dim_sz[1],dim_sz[0]]),tmp2.reshape([dim_sz[1],dim_sz[0]]))

    return tmp2


def pyFunction(data, params):

    nsh  = int(params[2])
    minW = int(params[0])
    maxW = int(params[1])

    dim_sz= np.array(data.shape)
    if len(dim_sz)==2:
        data= data[ :, :, np.newaxis, np.newaxis]
    elif len(dim_sz)==3:
        data= data[ :, :, :,np.newaxis]
    dim_sz = np.array(data.shape)

    output_vol= np.zeros(data.shape)

    data_complex = pyfftw.empty_aligned(dim_sz[0]*dim_sz[1], dtype=fftw_complex)
    res_complex = data_complex.copy()
    for t in range(dim_sz[3]):
        for z in range(dim_sz[2]):
            print(f'Unringing slice {z} of gradient {t} ...')
            for x in range(dim_sz[0]):
                for y in range(dim_sz[1]):
                    data_complex[dim_sz[0] * y + x] = complex(data[y, x, z, t],0)

            res_complex= unring_2d(data_complex, res_complex, dim_sz, nsh, minW, maxW)
            for x in range(dim_sz[0]):
                for y in range(dim_sz[1]):
                    output_vol[y, x, z, t] = res_complex[dim_sz[0] * y + x].real

    return output_vol.squeeze()

if __name__=='__main__':

    # data = np.random.randint(0, 10, (10, 10, 5))
    # params = [1, 3, 20]
    # pyFunction(data, params)

    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        from dipy.io.image import load_nifti, save_nifti
        import nibabel as nib

    lowResImgPath= '/home/tb571/Downloads/Harmonization-Python/test_data/test_a/connectom/connectom_a_dwi_bse.nii.gz'
    higResImgPath= lowResImgPath.split('.')[0] + '_py_unringed' + '.nii.gz'
    data, affine= load_nifti(lowResImgPath)
    params=[1, 3, 20]
    unring_data= pyFunction(data, params)
    save_nifti(higResImgPath, unring_data, affine)

    from subprocess import check_call
    higResImgPath= lowResImgPath.split('.')[0] + '_cc_unringed' + '.nii.gz'
    check_call(['unring.a64', lowResImgPath, higResImgPath])

    img= nib.load(higResImgPath)
    print(img.header)
    print((img.get_data()-unring_data).sum())






