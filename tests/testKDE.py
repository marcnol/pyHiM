#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 11:52:27 2020

@author: marcnol
"""

import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
import numpy as np
from scipy.stats import norm
from sklearn.neighbors import KernelDensity


from sklearn import mixture
import matplotlib.mlab
from pylab import *
from scipy.optimize import leastsq, curve_fit


def make_data(N, f=0.3, shift=1, baseline=2, rseed=1):
    rand = np.random.RandomState(rseed)
    x = rand.randn(N)

    x[int(f * N) :] += shift
    x += baseline
    return x


def make_data1(N, baseline=2, errorAmp=1, rseed=1):
    rand = np.random.RandomState(rseed)
    x = rand.randn(int(N))
    error = errorAmp * (np.random.rand(int(N)) - 0.5)
    x = x + error
    x += baseline
    x = np.abs(x)
    return x


def kdeFit(x, x_d, bandwidth=1, kernel="gaussian"):

    kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
    kde.fit(x[:, None])

    logprob = kde.score_samples(x_d[:, None])

    return logprob, kde


# def double_gaussian( x, params ):
#     (c1, mu1, sigma1, c2, mu2, sigma2) = params
#     res =   c1 * np.exp( - (x - mu1)**2.0 / (2.0 * sigma1**2.0) ) \
#           + c2 * np.exp( - (x - mu2)**2.0 / (2.0 * sigma2**2.0) )
#     return res

# def double_gaussian_fit( params ):
#     fit = double_gaussian( x, params )
#     return (fit - y_proc)

# def _1gaussian(x, amp1,cen1,sigma1):
#     return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))

# def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
#     return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
#             amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))

f = 0.05
N = [100, 10000]
shift = 3
baselines = [0.1, 3]
bandwidth = 0.2
threshold = 2
errorAmp = 1

x1 = np.concatenate(
    (
        make_data1(f * N[0], errorAmp=errorAmp, baseline=baselines[0]),
        make_data1((1 - f) * N[0], errorAmp=errorAmp, baseline=baselines[1]),
    )
)
x2 = np.concatenate(
    (
        make_data1(f * N[1], errorAmp=errorAmp, baseline=baselines[0]),
        make_data1((1 - f) * N[1], errorAmp=errorAmp, baseline=baselines[1]),
    )
)
x = [x1, x2]


bins = np.linspace(0, 10, 30)
x_d = np.linspace(0, 10, 100)

fig, ax = plt.subplots(
    1, 2, figsize=(12, 4), sharex=True, sharey=True, subplot_kw={"xlim": (0, 10), "ylim": (-0.02, 0.6)}
)
fig.subplots_adjust(wspace=0.05)
offsets = [0.0, 0.0]
for i, offset, xi in zip(range(len(offsets)), offsets, x):
    ax[i].hist(xi, bins=bins + offset, normed=True, alpha=0.5, color="r")
    ax[i].plot(xi, np.full_like(xi, -0.01), "|k", markeredgewidth=0.5)

    logprob, kde = kdeFit(xi, x_d, bandwidth=bandwidth)
    ax[i].fill_between(x_d, np.exp(logprob) * 1.3, alpha=0.3)
    contactProbability = np.nonzero(xi < threshold)[0].shape[0] / xi.shape[0]
    ax[i].text(7, 0.5, "p = {}\nmean = {:.2f}".format(contactProbability, np.mean(xi)))
    print("subplot: {} | contactP = {}".format(i, contactProbability))


#%%

from sklearn.model_selection import GridSearchCV, LeaveOneOut

# from sklearn.cross_validation import LeaveOneOut

bandwidths = 10 ** np.linspace(-1, 1, 100)
grid = GridSearchCV(KernelDensity(kernel="gaussian"), {"bandwidth": bandwidths}, cv=LeaveOneOut())
grid.fit(x2[:, None])
grid.best_params_
# # Least squares fit. Starting values found by inspection.
# popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, x_d, y_array_2gauss, p0=[amp1, cen1, sigma1, amp2, cen2, sigma2])

# perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
# pars_1 = popt_2gauss[0:3]
# pars_2 = popt_2gauss[3:6]
# gauss_peak_1 = _1gaussian(x_array, *pars_1)
# gauss_peak_2 = _1gaussian(x_array, *pars_2)
