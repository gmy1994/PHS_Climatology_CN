# Calculate the monthly climatology over 1980-2018
# TR, ET, TR/ET, Beta, surface RO, subsurface RO, RO, SM
import os
import os.path
import subprocess
import numpy as np
import netCDF4 as nc

INDIR = "/scratch/06956/mguo/NoahMP_PHS/CN/output/process/monthly/"
RO_DIR = "/scratch/06956/mguo/NoahMP_PHS/CN/output/process/rnoff/monthly/"
OUTDIR = "/scratch/06956/mguo/NoahMP_PHS/CN/experiment/clim_monthly/"
TEMPLATE  = "/work2/06956/mguo/script/climatology/clim_monthly_TEM.cdl"

MASKFILE = '/work2/06956/mguo/script/mask/MASK_FOREST.nc'
mask = nc.Dataset(MASKFILE).variables['MASK_FOREST'][:]

lat = nc.Dataset("/scratch/06956/mguo/NoahMP_PHS/CN/output/process/clim_39/map_2m/clim.nc").variables['lat'][:]
lon = nc.Dataset("/scratch/06956/mguo/NoahMP_PHS/CN/output/process/clim_39/map_2m/clim.nc").variables['lon'][:]

os.system("mkdir " + OUTDIR)

phsnm = ["map","oak"]
soilm = [2,10]
month = ['01','02','03','04','05','06','07','08','09','10','11','12']

for i in range(len(phsnm)):
    for n in range(len(soilm)):
        for m in range(len(month)):
            phs = phsnm[i] + "_" + str(soilm[n]) + "m/"
            os.system("mkdir " + OUTDIR + phs)
            os.system("mkdir " + OUTDIR + phs + month[m] + "/" )
        
            # TR ET Beta
            TR = np.ma.zeros([1, 400, 700])
            ET = np.ma.zeros([1, 400, 700])
            TR_ET = np.ma.zeros([1, 400, 700])
            VB = np.ma.zeros([1, 400, 700])

            # RO
            SFC_RO = np.ma.zeros([1, 400, 700])
            UGD_RO = np.ma.zeros([1, 400, 700])
            RO = np.ma.zeros([1, 400, 700])

            # SM
            SM = np.ma.zeros([1, 400, 4, 700])

            for y in range(1980,2018):
                input_file = INDIR + phs +"mon-" + str(y) + month[m] + ".nc"
                ETRAN = nc.Dataset(input_file).variables['ETRAN'][:]
                ECAN = nc.Dataset(input_file).variables['ECAN'][:]
                EDIR = nc.Dataset(input_file).variables['EDIR'][:]
                VBTRAN = nc.Dataset(input_file).variables['VBTRAN'][:]
                SOIL_M = nc.Dataset(input_file).variables['SOIL_M'][:]

                TR[0] += ETRAN[0,:]*86400
                ET[0] += (ETRAN[0,:] + ECAN[0,:] + EDIR[0,:])*86400
                TR_ET[0] += np.ma.true_divide(ETRAN[0,:], ETRAN[0,:] + ECAN[0,:] + EDIR[0,:])
                VB[0] += VBTRAN[0,:]
                for la in range(400):
                    SM[0, la, :] += SOIL_M[0, la, :]

                input_file = RO_DIR + phs +"mon-" + str(y) + month[m] + ".nc"
                UGDRNOFF = nc.Dataset(input_file).variables['UGDRNOFF'][:]
                SFCRNOFF = nc.Dataset(input_file).variables['SFCRNOFF'][:]

                SFC_RO[0]  += SFCRNOFF[0,:]*86400
                UGD_RO[0]  += UGDRNOFF[0,:]*86400
                RO[0]  += (SFCRNOFF[0,:] + UGDRNOFF[0,:])*86400

            clim_TR = TR/39
            clim_TR.mask = mask 

            clim_ET = ET/39
            clim_ET.mask = mask 

            clim_TR_ET = TR_ET/39
            clim_TR_ET.mask = mask 

            clim_VBTRAN = VB/39
            clim_VBTRAN.mask = mask 

            clim_SFC_RO = SFC_RO/39
            clim_SFC_RO.mask = mask 

            clim_UGD_RO = UGD_RO/39
            clim_UGD_RO.mask = mask 

            clim_RO = RO/39
            clim_RO.mask = mask 

            clim_SM = SM/39
            for j in range(4):
                tmp = clim_SM[0,:,j,:]
                tmp.mask = mask
                clim_SM[0,:,j,:] = tmp 

            output = OUTDIR + phs + month[m] + "/" + "clim.nc"
            os.system("rm " + output)
            subprocess.run(['ncgen', '-3', '-o', output, TEMPLATE])
            with nc.Dataset(output, 'a') as ncf:
                ncf.variables['lat'][:] = lat[:]
                ncf.variables['lon'][:] = lon[:]
                ncf.variables['clim_TR'][:] = clim_TR[:]
                ncf.variables['clim_ET'][:] = clim_ET[:]
                ncf.variables['clim_TR_ET'][:] = clim_TR_ET[:]
                ncf.variables['clim_VBTRAN'][:] = clim_VBTRAN[:]
                ncf.variables['clim_SFC_RO'][:] = clim_SFC_RO[:]
                ncf.variables['clim_UGD_RO'][:] = clim_UGD_RO[:]
                ncf.variables['clim_RO'][:] = clim_RO[:]
                ncf.variables['clim_SM'][:] = clim_SM[:]