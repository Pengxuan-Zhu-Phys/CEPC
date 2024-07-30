#!/usr/bin/env python3 

import os
import pandas as pd 
import numpy as np

def run_task(name, path): 
    cmd = f"rivet -a CEPC_mRC:DFNAME={name} -o {name}.yoda {path} --pwd"
    print(cmd)
    os.system(cmd)

def read_lumin(path):
    with open(path) as f1: 
        cont = f1.readlines()
        nmc, lumi = "", ""
        for lin in cont: 
            if "| Simulation: requested number of events = " in lin: 
                nmc = lin.replace("| Simulation: requested number of events = ", "")
                # print(nmc)
            if "|             corr. to luminosity [fb-1] = " in lin: 
                lumi = lin.replace("|             corr. to luminosity [fb-1] = ", "")
                # print(lumi)
        if nmc != "" and lumi != "":
            nmc = int(nmc)
            lumi = float(lumi)
            # print("Total cross section is -> ", nmc/lumi, nmc)
    return nmc, lumi

    

if __name__=="__main__":
    wmss = ['80362', '80367', '80372', '80377', '80382', '80387', '80392']
    ids = ["01", "02", "03", "04", "05", "06", "07"]
    df = []
    for ii in range(7):
        # p01m = f"/Volumes/Buding_T5/Wmass/events_01M/{ids[ii]}/ww_l0e1e2_{ids[ii]}.hepmc" 
        # name = f"{wmss[ii]}_1M"
        # run_task(name=name, path=p01m)
        # idx = "{:003}".format(0)
        # name = f"{wmss[ii]}_26M_{idx}"
        # run_task(name=name, path=p26m)
        for jj in range(26): 
            idx = "{:003}".format(jj)
            # p26m = f"/Volumes/Buding_T5/Wmass/events_26M/{ids[ii]}/ww_l0e1e2_{ids[ii]}_{idx}.0.hepmc" 
            p26log = f"/Volumes/Buding_T5/Wmass/events_26M/{ids[ii]}/whizard_{idx}.log" 
            nmc, lumi = read_lumin(p26log)
            print(ids[ii], idx, "Total cross section is -> ", nmc/lumi, nmc, lumi)
            cross = np.float64(1. / lumi )
            df.append({
                "wmass":    wmss[ii], 
                "idx":  idx, 
                "xsect/events": cross
            })
    df = pd.DataFrame(df)
    print(df)
    df.to_csv("events_26M_summary.csv", index=False)




