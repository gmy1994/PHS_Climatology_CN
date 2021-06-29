# Calculate the ralative difference (%) of monthly climatology for Oak (10m) vs. Maple
import os
import os.path
import subprocess
import numpy as np
import netCDF4 as nc

INDIR = "/scratch/06956/mguo/NoahMP_PHS/CN/experiment/clim_monthly/"
TEMPLATE = "./clim_monthly_TEM.cdl"
OUTDIR = "/scratch/06956/mguo/NoahMP_PHS/CN/experiment/clim_monthly/dif/"

os.system("mkdir " + OUTDIR)

lat = nc.Dataset("/scratch/06956/mguo/NoahMP_PHS/CN/output/process/clim_39/map_2m/clim.nc").variables['lat'][:]
lon = nc.Dataset("/scratch/06956/mguo/NoahMP_PHS/CN/output/process/clim_39/map_2m/clim.nc").variables['lon'][:]

month = ['01','02','03','04','05','06','07','08','09','10','11','12']

for m in range(len(month)):
    os.system("mkdir " + OUTDIR + month[m] + "/")

    oak = INDIR + "oak_10m/" + month[m] + "/" + "clim.nc"
    map = INDIR + "map_2m/" + month[m] + "/" + "clim.nc"

    vars = ["clim_ET","clim_TR","clim_ED","clim_TR_ET","clim_VBTRAN","clim_SFC_RO","clim_UGD_RO","clim_RO","clim_SM","clim_Q","clim_GPP"]

    output = OUTDIR + month[m] + "/" + "oak_vs_map.nc"
    os.system("rm " + output)
    subprocess.run(['ncgen', '-3', '-o', output, TEMPLATE])
    with nc.Dataset(output, 'a') as ncf:
        ncf.variables['lat'][:] = lat[:]
        ncf.variables['lon'][:] = lon[:]
        for var in vars:
            ncf.variables[var][:] = 100 * (nc.Dataset(oak).variables[var][:] - nc.Dataset(map).variables[var][:]) / nc.Dataset(oak).variables[var][:]

