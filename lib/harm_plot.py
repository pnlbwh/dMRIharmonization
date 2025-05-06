#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from conversion import read_imgs_masks
from os.path import isfile, basename

XERR=0.1
ELINEWIDTH=3
CAPSIZE=5
CAPTHICK=3
FMT='cD'

def harm_plot(ydata, labels, outPrefix):
    '''
    :param ydata: list of [y1, y2, y3, ...] where each yi is a list
    :param labels: list of strings
    :param outPrefix:
    :return:
    '''

    labels= list(labels)

    num_series= len(ydata)

    iter_obj= [i for i in range(num_series)]

    # errorbar plot
    plt.figure(1)
    plt.grid(True)
    for i in iter_obj:
        num_sub= len(ydata[i])
        x= list(i*np.ones((num_sub,)))
        y= ydata[i]
        plt.plot(x, y, 'r*')

        plt.errorbar([i], np.mean(y), xerr=XERR, yerr=np.std(y),
                     ecolor='k', capsize=CAPSIZE, capthick=CAPTHICK, elinewidth=ELINEWIDTH, fmt=FMT)

    plt.xticks(iter_obj, labels)
    plt.title('Comparison of meanFA before and after harmonization')
    plt.ylabel('meanFA over IIT_mean_FA_skeleton')
    plt.savefig(outPrefix+'_ebarplot.png')
    # plt.show()


    # box plot
    # plt.figure(2)
    # plt.grid(True)
    # for i in iter_obj:
    #     num_sub= len(ydata[i])
    #     x = list(i * np.ones((num_sub,)))
    #     y = ydata[i]
    #     plt.plot(x, y, 'r*')
    #
    # plt.boxplot(ydata, labels=labels, positions=iter_obj,
    #             boxprops=dict(linewidth=4),
    #             medianprops=dict(linewidth=4),
    #             whiskerprops=dict(linewidth=2))
    #
    #
    # plt.title('Comparison of boxplot before and after harmonization')
    # plt.ylabel('meanFA over IIT_mean_FA_skeleton')
    # plt.savefig(outPrefix+'_boxplot.png')
    # plt.show()
    # return (outPrefix+'_ebarplot.png', outPrefix+'_boxplot.png')


    return outPrefix+'_ebarplot.png'


def generate_csv(imgs, site_means, outPrefix):

    try:
        imgs, _= read_imgs_masks(imgs)
    except:
        pass

    stat = {'subject': [basename(f) for f in imgs], 'meanFA': site_means}
    df = pd.DataFrame(stat)
    statFile = outPrefix + '_stat.csv'

    df.to_csv(statFile, index=False)


if __name__=='__main__':

    sub=['hi','hello','go','come']
    ref_mean= list(0.05*np.random.randn(20)+0.5) # [0.46, 0.49, 0.44, 0.40]
    target_mean_before= list(0.05*np.random.randn(20)+0.45) # [0.42, 0.58, 0.43, 0.66]
    target_mean_after= list(0.05*np.random.randn(20)+0.5) # [0.5 , 0.45, 0.40, 0.55]
    labels=['Reference','Target_before','Target_after']

    harm_plot([ref_mean, target_mean_before, target_mean_after], labels, '/tmp/abc')
    # harm_plot([ref_mean], ['Reference'], '/tmp/abc')

