#!/usr/bin/env python3 

import matplotlib.pyplot as plt 
import pandas as pd 
import numpy as np 
from scipy import interpolate
from scipy.stats import gaussian_kde
import os, sys 
pwd = os.path.abspath(os.path.dirname(__file__))
tpath = os.path.abspath(os.path.join(pwd, "../tables"))

sumtable = pd.read_csv("events_26M_summary.csv")

fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.1, 0.1, 0.84, 0.84])

def get_weight(wmass, idx):
    dfs = sumtable[(sumtable['wmass'] == int(wmass)) & (sumtable['idx'] == int(idx))]
    dfs = dict(dfs.iloc[0])
    return dfs['xsect/events']

    
for mw in ['80362', '80367', '80372', '80377', '80382', '80387', '80392']:
    # df = pd.read_csv(os.path.join(tpath, "80362_1M.csv"))
    dfs = []

    for ii in range(26):
        idx = "{:003}".format(ii)
        print("Processing ",mw, idx)
        weight = get_weight(mw, idx)

        df = pd.read_csv(os.path.join(tpath, "{}/{}_26M_{}.csv".format(mw, mw, idx)))
        mrcb = np.array(df['mRC_B'].sort_values(ascending=True))
        weigt = np.ones(mrcb.shape[0])
        weigt = weigt * weight
        ddf = pd.DataFrame({
            "mRCmax":   mrcb, 
            "weight":   weigt 
        })
        print(ddf.shape)
        dfs.append(ddf)

    dff = pd.concat(dfs)
    dff = dff.sort_values(by="mRCmax", ascending=True)

    evs = np.histogram(dff['mRCmax'], range=[0, 120], bins=600, weights=dff['weight'])[0]
    xx = np.linspace(0, 120, 601)
    df = pd.DataFrame({
        "mRCmax":   evs, 
        "xminus":   xx[:-1],
        "xplus":    xx[1:]
    })
    df.to_csv("hist_mRCmax_{}_bin_02GeV.csv".format(mw), index=False)


    # wig = np.array(dff['weight'])
    # cdf = np.ones(dff.shape[0])
    # # cdf = np.insert(np.cumsum(wig), 0, 0)
    # cdf = np.cumsum(wig)
    # dff['CDF'] = cdf 
    # print(cdf.shape)

    # print(dff)

    # cutdf = dff[::400000]
    # print(cutdf.shape)
    # mrCB = np.concatenate(mrcbs)
    # cdf = np.concatenate(weights)
    # print("Finish concat")

    # dff = pd.DataFrame({"mRCmax":mrCB, "weight":cdf})
    # print(dff.shape())
    # # mrcb = np.array(df['mRC_B'].sort_values(ascending=True))

    # print(mrcb, mrcb.shape)
    # mrcb = np.insert(mrcb, 0, 0)

    # cdf = np.ones(mrcb.shape[0])
    # cdf = cdf * 1 
    # cdf = np.insert(np.cumsum(cdf), 0, 0)
    # cdf = cdf[:-1]
    # # print(cdf.shape, mrcb.shape, cdf[-1], mrcb[-1])
    # unique_mrcb, unique_indices = np.unique(mrcb, return_index=True)
    # unique_cdf = cdf[unique_indices]
    # cdf_fn = interpolate.interp1d(unique_mrcb, unique_cdf, kind='linear')
    # xx = np.linspace(0, 120, 1201)
    # yy = cdf_fn(xx)
    # pdf = np.zeros_like(xx)
    # pdf[1:] = yy[1:] - yy[:-1]
    # cdf_finefn = interpolate.interp1d(xx, yy, kind='cubic')
    # xx = np.linspace(0, 120, 241)
    # zz = cdf_finefn(xx)
    # pdf = np.zeros_like(xx)
    # pdf[1:] = zz[1:] - zz[:-1]
    # ax.plot(xx, pdf, '-')

    # pdf = np.zeros_like(np.array(cutdf["CDF"]))
    # print(pdf.shape)
    # pdf[1:] = (np.array(cutdf["CDF"])[1:] - np.array(cutdf["CDF"])[:-1]) / (np.array(cutdf['mRCmax'])[1:] - np.array(cutdf['mRCmax'])[:-1] )

    # print(pdf)






    # kde = gaussian_kde(cutdf['mRCmax'], bw_method="silverman", weights=pdf)
    # yy = kde(xx)
    # yy = yy/max(yy)
    # ax.plot(xx, yy, '.')

    # print(xx, yy)
    # ax.plot(cutdf['mRCmax'], pdf, '.-')
    # ax.plot(xx, zz, '-')
    # ax.hist(df['mRC_B'], bins=120, facecolor='white', edgecolor="b", linewidth=1.0)
# plt.show()