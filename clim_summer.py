# Calculate summer (JJA) climatology over 1980-2018
# TR, ET, ED, TR/ET, Beta, surface RO, subsurface RO, RO, SM, root uptake Q, GPP
import os
import os.path
import subprocess
import numpy as np
import netCDF4 as nc

INDIR = "/scratch/06956/mguo/NoahMP_PHS/CN/experiment/clim_monthly/"
OUTDIR = "/scratch/06956/mguo/NoahMP_PHS/CN/experiment/clim_summer/"
TEMPLATE  = "/work2/06956/mguo/script/climatology/clim_monthly_TEM.cdl"

MASKFILE = '/work2/06956/mguo/script/mask/MASK_FOREST.nc'
mask = nc.Dataset(MASKFILE).variables['MASK_FOREST'][:]

lat = nc.Dataset("/scratch/06956/mguo/NoahMP_PHS/CN/output/process/clim_39/map_2m/clim.nc").variables['lat'][:]
lon = nc.Dataset("/scratch/06956/mguo/NoahMP_PHS/CN/output/process/clim_39/map_2m/clim.nc").variables['lon'][:]

os.system("mkdir " + OUTDIR)

phsnm = ["map","oak"]
soilm = [2,10]
month = ['06','07','08']

for i in range(len(phsnm)):
    for n in range(len(soilm)):
        phs = phsnm[i] + "_" + str(soilm[n]) + "m/"
        os.system("mkdir " + OUTDIR + phs)
        # TR ED ET Beta
        summer_TR = np.ma.zeros([1, 400, 700])
        summer_ED = np.ma.zeros([1, 400, 700])
        summer_ET = np.ma.zeros([1, 400, 700])
        summer_TR_ET = np.ma.zeros([1, 400, 700])
        summer_VB = np.ma.zeros([1, 400, 700])
        # RO
        summer_SFC_RO = np.ma.zeros([1, 400, 700])
        summer_UGD_RO = np.ma.zeros([1, 400, 700])
        summer_RO = np.ma.zeros([1, 400, 700])
        # GPP
        summer_GPP = np.ma.zeros([1, 400, 700])
        # SM Q
        summer_SM = np.ma.zeros([1, 400, 4, 700])
        summer_Q = np.ma.zeros([1, 400, 4, 700])
        for m in range(len(month)):
            input_file = INDIR + phs + month[m] + "/clim.nc"
            month_TR = nc.Dataset(input_file).variables['clim_TR'][:]
            month_ED = nc.Dataset(input_file).variables['clim_ED'][:]
            month_ET = nc.Dataset(input_file).variables['clim_ET'][:]
            month_TR_ET = nc.Dataset(input_file).variables['clim_TR_ET'][:]
            month_VB = nc.Dataset(input_file).variables['clim_VBTRAN'][:]
            month_SFC_RO = nc.Dataset(input_file).variables['clim_SFC_RO'][:]
            month_UGD_RO = nc.Dataset(input_file).variables['clim_UGD_RO'][:]
            month_RO = nc.Dataset(input_file).variables['clim_RO'][:]
            month_GPP = nc.Dataset(input_file).variables['clim_GPP'][:]
            month_SM = nc.Dataset(input_file).variables['clim_SM'][:]
            month_Q = nc.Dataset(input_file).variables['clim_Q'][:]

            summer_TR[0] += month_TR[0,:]
            summer_ED[0] += month_ED[0,:]
            summer_ET[0] += month_ET[0,:]
            summer_TR_ET[0] += month_TR_ET[0,:]
            summer_VB[0] += month_VB[0,:]
            summer_SFC_RO[0] += month_SFC_RO[0,:]
            summer_UGD_RO[0] += month_UGD_RO[0,:]
            summer_RO[0] += month_RO[0,:]
            summer_GPP[0] += month_GPP[0,:]

            for la in range(400):
                summer_SM[0, la, :] += month_SM[0, la, :]
                summer_Q[0, la, :] += month_Q[0, la, :]

        clim_TR = summer_TR/len(month)
        clim_TR.mask = mask 

        clim_ED = summer_ED/len(month)
        clim_ED.mask = mask 

        clim_ET = summer_ET/len(month)
        clim_ET.mask = mask 

        clim_TR_ET = summer_TR_ET/len(month)
        clim_TR_ET.mask = mask 

        clim_VBTRAN = summer_VB/len(month)
        clim_VBTRAN.mask = mask 

        clim_SFC_RO = summer_SFC_RO/len(month)
        clim_SFC_RO.mask = mask 

        clim_UGD_RO = summer_UGD_RO/len(month)
        clim_UGD_RO.mask = mask 

        clim_RO = summer_RO/len(month)
        clim_RO.mask = mask 

        clim_GPP = summer_GPP/len(month)
        clim_GPP.mask = mask 

        clim_SM = summer_SM/len(month)
        for j in range(4):
            tmp = clim_SM[0,:,j,:]
            tmp.mask = mask
            clim_SM[0,:,j,:] = tmp 

        clim_Q = summer_Q/len(month)
        for j in range(4):
            tmp = clim_Q[0,:,j,:]
            tmp.mask = mask
            clim_Q[0,:,j,:] = tmp 

        output = OUTDIR + phs + "clim.nc"
        os.system("rm " + output)
        subprocess.run(['ncgen', '-3', '-o', output, TEMPLATE])
        with nc.Dataset(output, 'a') as ncf:
            ncf.variables['lat'][:] = lat[:]
            ncf.variables['lon'][:] = lon[:]
            ncf.variables['clim_TR'][:] = clim_TR[:]
            ncf.variables['clim_ED'][:] = clim_ED[:]
            ncf.variables['clim_ET'][:] = clim_ET[:]
            ncf.variables['clim_TR_ET'][:] = clim_TR_ET[:]
            ncf.variables['clim_VBTRAN'][:] = clim_VBTRAN[:]
            ncf.variables['clim_SFC_RO'][:] = clim_SFC_RO[:]
            ncf.variables['clim_UGD_RO'][:] = clim_UGD_RO[:]
            ncf.variables['clim_RO'][:] = clim_RO[:]
            ncf.variables['clim_GPP'][:] = clim_GPP[:]
            ncf.variables['clim_SM'][:] = clim_SM[:]
            ncf.variables['clim_Q'][:] = clim_Q[:]