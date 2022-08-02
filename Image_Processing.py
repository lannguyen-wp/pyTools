#------------------------------------------------------------------------------
#-- / Description -------------------------------------------------------------


#-- Script  : image_processing.py
#-- Version : 1.0
#-- Date    : 2020-01-15
#-- Created : @Lan.Nguyen (lnguyen1985@yahoo.com)

#-- List of tools:
#--     1. Projection: reproject a raster to another projection
#--     2. Subset:
#--     3. Extents_Obtain:
#--     4. Extents_Match:

#-- Description / -------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#-- / Variable Definition -----------------------------------------------------

#-- i_dir  : input directory
#-- o_dir  : output directory
#-- i_ras  : full directory of input raster
#-- o_ras  : full directory of output raster
#-- r_ras  : full directory of reference raster
#-- roi    : region-of-interest (raster | [coordinates, resolution])
#-- dtype  : data type of input/output raster
#-- vrange : [min_value, max_value] => valid value range of a raster   
#-- min_dp : min valid_data/total_data proportion
#-- flist  : list of filenames
#-- ds     : gdal dataset
#-- ar     : numpy array
#-- pj     : projection
#-- gt     : geotransform
#-- prows  : processing_rows - number of rows to be processed each time
#-- band   : band - key to be used in 'List_Dir' or part of filename
#-- window : [row, col] - number of rows and columns to extract data
#-- srs    : source features (raster/extent) to be processed
#-- xy     : an array contains x and y coordinates
#-- xycol  : column indexes of x(col) and y(row) from 'xy'
#-- sep    : separate (True) or combine (False) bands in dataset
#-- fid    : feature's id

#-- Variable Definition / -----------------------------------------------------
#------------------------------------------------------------------------------

import datetime, glob, math, os, sys
import numpy as np

try:
    import gdal, gdalconst, osr
except:
    from osgeo import gdal, gdalconst, osr

from scipy.optimize import leastsq
from pathlib import Path

usr_dir = str(Path.home())                      #- user directory
dbx_dir = usr_dir + '\\Dropbox'                 #- dropbox directory
lib_dir = dbx_dir + '\\Toolbox\\External_Libs'  #- library directory
log_dir = dbx_dir + '\\Toolbox\\Logs'           #- log directory
tst_dir = dbx_dir + '\\Toolbox\\Testing'        #- testing directory
tmp_dir = 'D:\\Temp'

sys.path.append(dbx_dir + '\\Toolbox')
import mytools as tls

#-- Project an image using a reference raster ---------------------------------
def Projection(i_ras, o_ras, r_ras, dtype=gdal.GDT_Int16):
    i_ds = gdal.Open(i_ras)
    i_pj = i_ds.GetProjection()
    i_pj_osr = osr.SpatialReference(i_pj)
    i_gt = i_ds.GetGeoTransform()
    i_ds.GetRasterBand(1).SetNoDataValue(-9999)
    r_ds = gdal.Open(r_ras)
    r_pj = r_ds.GetProjection()
    r_pj_osr = osr.SpatialReference(r_pj)
    r_gt = r_ds.GetGeoTransform()
    x_size = i_ds.RasterXSize
    y_size = i_ds.RasterYSize
    tx = osr.CoordinateTransformation (i_pj_osr, r_pj_osr)
    if int(gdal.VersionInfo()[0])>=3:
        (ulx, uly, ulz) = tx.TransformPoint(i_gt[3], i_gt[0])
        (lrx, lry, lrz) = tx.TransformPoint(i_gt[3] + i_gt[5]*y_size, 
                                            i_gt[0] + i_gt[1]*x_size)
    else:
        (ulx, uly, ulz) = tx.TransformPoint(i_gt[0], i_gt[3])
        (lrx, lry, lrz) = tx.TransformPoint(i_gt[0] + i_gt[1]*x_size,
                                            i_gt[3] + i_gt[5]*y_size)
    nc = int((lrx - ulx)/r_gt[1])
    nr = int((uly - lry)/(-1*r_gt[5]))
    xmin = r_gt[0] + round((ulx - r_gt[0])/r_gt[1], 0)*r_gt[1]
    ymax = r_gt[3] + round((uly - r_gt[3])/r_gt[5], 0)*r_gt[5]
    o_gt = (xmin, r_gt[1], 0, ymax, 0, r_gt[5])
    o_ds = gdal.GetDriverByName('GTiff').Create(o_ras, nc, nr, 1, dtype)
    o_ds.SetGeoTransform(o_gt)
    o_ds.SetProjection(r_pj)
    o_ds.GetRasterBand(1).SetNoDataValue(-9999)
    gdal.ReprojectImage(i_ds, o_ds, i_pj, r_pj, gdalconst.GRA_Bilinear)
#------------------------------------------------------------------------------


#-- Subset an image using ROI -------------------------------------------------
def Subset(i_ras, roi, data=False, vrange=[0,0], min_dp=0.1):
    vr = vrange; md = min_dp; output = [0]
    if (len(roi)==4) & (data==False):
        i_xmin = roi[0]; i_ymin = roi[1]
        i_xmax = roi[2]; i_ymax = roi[3]
    else:
        r_ds = gdal.Open(roi[0])
        r_gt = r_ds.GetGeoTransform()
        r_nr = r_ds.RasterYSize
        r_nc = r_ds.RasterXSize
        r_xmin, r_xres, r_xskew, r_ymax, r_yskew, r_yres = r_gt
        r_xmax = r_xmin + (r_nc * r_xres)
        r_ymin = r_ymax + (r_nr * r_yres)
    i_ds = gdal.Open(i_ras)
    i_gt = i_ds.GetGeoTransform()
    i_pj = i_ds.GetProjection()
    i_nr = i_ds.RasterYSize
    i_nc = i_ds.RasterXSize
    i_xmin, i_xres, i_xskew, i_ymax, i_yskew, i_yres = i_gt
    i_xmax = i_xmin + (i_nc * i_xres)
    i_ymin = i_ymax + (i_nr * i_yres)
    if ((i_xmin<r_xmax)&(i_xmax>r_xmin)&(i_ymin<r_ymax)&(i_ymax>r_ymin)):
        lx = int((max(r_xmin, i_xmin) - i_xmin) / r_xres)
        ly = int((i_ymax - min(r_ymax, i_ymax)) / r_xres)
        ux = int((min(r_xmax, i_xmax) - i_xmin) / r_xres)
        uy = int((i_ymax - max(r_ymin, i_ymin)) / r_xres)
        if data==False:
            output = [lx, ly, ux, uy]
        else:
            if round(i_xres, 4)==round(r_xres, 4):
                i_ar = i_ds.ReadAsArray()
            elif i_xres==(r_xres*int(i_xres/r_xres)):
                out_rows = i_ds.GetRasterBand(1).YSize * int(i_yres/r_yres)
                out_cols = i_ds.GetRasterBand(1).XSize * int(i_xres/r_xres)
                i_ar = i_ds.ReadAsArray(buf_xsize=out_cols, buf_ysize=out_rows)
                i_ar = i_ar[ly:uy, lx:ux]
                i_gt = [max(r_xmin, i_xmin), r_xres, i_xskew, \
                        min(r_ymax, i_ymax), i_yskew, r_yres]
            if md>0:
                unq, cnt = np.unique(i_ar, return_counts=True)
                pixs_valid = np.sum(cnt[(unq>=vr[0]) & (unq<=vr[1])])
                pixs_percent = pixs_valid*100/np.sum(cnt)
                i_ar[(i_ar<vr[0]) | (i_ar>vr[1])] = -9999
                if pixs_percent>=md:
                    output = [i_ar, i_gt, i_pj]
            else:
                output = [i_ar, i_gt, i_pj]
    return output
#------------------------------------------------------------------------------


#-- Obtain extents of a list of images ---------------------------------------
def Extents_Obtain(flist):
    if len(flist)==1:
        ds = gdal.Open(flist[0])
        nr = ds.RasterYSize
        nc = ds.RasterXSize
        xmin, xres, xskew, ymax, yskew, yres = ds.GetGeoTransform()
        xmax = xmin + (nc * xres)
        ymin = ymax + (nr * yres)
        extents = [xmin, ymin, xmax, ymax, xres]
    else:
        extents = np.zeros((len(flist), 5))
        count = 0
        for file in flist:
            ds = gdal.Open(file)
            nr = ds.RasterYSize
            nc = ds.RasterXSize
            xmin, xres, xskew, ymax, yskew, yres = ds.GetGeoTransform()
            xmax = xmin + (nc * xres)
            ymin = ymax + (nr * yres)
            extents[count, 0] = xmin
            extents[count, 1] = ymin
            extents[count, 2] = xmax
            extents[count, 3] = ymax
            extents[count, 4] = xres
            count += 1
    return extents
#------------------------------------------------------------------------------


#-- Match extent of ROI to a list of features -----------------------------
def Extents_Match(srs, roi, vrange=[0,0], min_dp=0):
    matching = []
    if len(roi)==1:
        s_xmin, s_ymin, s_xmax, s_ymax, s_xres = Extents_Obtain(roi)
    if len(roi)==5:
        s_xmin = roi[0]
        s_xmax = roi[1]
        s_ymin = roi[2]
        s_ymax = roi[3]
    for source in srs:
        if len(source)==1:
            t_xmin, t_ymin, t_xmax, t_ymax, t_xres = Extents_Obtain(source)
        if len(source)==5:
            t_xmin = source[0]
            t_xmax = source[1]
            t_ymin = source[2]
            t_ymax = source[3]
        if ((t_xmin<s_xmax)&(t_xmax>s_xmin)&(t_ymin<s_ymax)&(t_ymax>s_ymin)):
            if ((vrange==[0,0]) & (min_dp==0)):
                matching.append(1)
            else:
                op = Subset(i_ras=source, roi=roi, data=False, \
                            vrange=vrange, min_dp=min_dp)
                if len(op)==3:
                    matching.append(1)
                else:
                    matching.append(0)               
        else:
            matching.append(0)
    return matching

def Neighbors_C(mask, window):
    rs = window[0]
    cs = window[1]
    nr = mask.shape[0]
    nc = mask.shape[1]
    neighbor = np.zeros((nr-2, nc-2))
    for i in range(0, rs):
        for j in range(0, cs):
            neighbor = neighbor + mask[i:(nr-rs+i+1), j:(nc-cs+j+1)]
    return neighbor

def Neighbors_F(mask, value, window=[3,3]):
    neighbor = Neighbors_C(mask, window)
    if value>0:
        f_ar = ((neighbor>0) & (neighbor<value)) * 1
    if value<0:
        f_ar = (neighbor>=(-1 * value)) * 1
    f_ar = np.pad(f_ar, ((1,1),(1,1)), mode='constant', constant_values=1)
    nr = f_ar.shape[0]
    nc = f_ar.shape[1]
    f_ar[0,:] = mask[0,:]
    f_ar[:,0] = mask[:,0]
    f_ar[nr-1,:] = mask[nr-1,:]
    f_ar[:,nc-1] = mask[:,nc-1]
    f_ar[mask==0] = 0
    return f_ar

def Edge_Removal (array, edge_remove=[0,0]):
    from scipy import ndimage
    ar = array; er = edge_remove
    mask = (array==er[1])*1
    edge = Neighbors_F(mask, value=9)
    edge = Neighbors_F(edge, value=-3)
    buff = ndimage.morphology.binary_dilation(edge, iterations=er[0])
    ar[(buff==1) | (mask==1)] = -9999
    return ar
#------------------------------------------------------------------------------


#-- Create virtual dataset from an image stack --------------------------------
def Virtual_DS(i_dir, pattern='*.tif', \
               roi='', data_check=[0,0,0], \
               nodata=-9999, sep=True):
    tifs = glob.glob(i_dir + '\\' + pattern)     
    if roi=='':
        vt = tmp_dir + '\\virtualDS.vrt'
        ds = gdal.BuildVRT(vt, tifs, srcNodata=nodata, separate=sep)
    else:
        if len(roi)==5:
            xmin = roi[0]; ymin = roi[2]
            xmax = roi[1]; ymax = roi[3]
            xres = roi[4]
        else:
            xmin, xmax, ymin, ymax, xres = Extents_Obtain(roi)
        roi = (xmin, ymin, xmax, ymax)
        matching = Extents_Match(roi, tifs, data_check)
        in_list  = [tifs[i] for i in range(0, len(matching)) if matching[i]==1]    
        vt = tmp_dir + '\\virtualDS.vrt'
        ds = gdal.BuildVRT(vt, in_list,         \
                           srcNodata=nodata,    \
                           separate=sep,        \
                           outputBounds=roi)
        tifs = in_list
    tifs = np.array([os.path.basename(tif) for tif in tifs])
    return tifs, ds
#------------------------------------------------------------------------------
    

#-- Compute statistics for an image stack -------------------------------------
def st_Stats_Sliced(ds, vrange, stats, prows, band, o_dir):
    import warnings
    warnings.filterwarnings("ignore")
    lb = vrange[0]; ub = vrange[1]  
    xmin, xres, xskew, ymax, yskew, yres = ds.GetGeoTransform()
    pj = ds.GetProjection()
    nc = ds.RasterXSize
    nr = ds.RasterYSize
    blocks = list(range(0, nr, prows))
    blocks.append(nr); nb = len(blocks); zf = int(math.log10(nb))+1
    log_file = log_dir + '\st_Stats_Sliced_' + band + '.txt'
    if os.path.isfile(log_file)==False:
        f = open(log_file, 'w')
        f.close()
    for i in range(1, nb):
        now = datetime.datetime.now()
        now = now.strftime("%Y-%m-%d %H:%M")
        print('S--- ' + str(i).zfill(zf) + '/' + str(nb-1) + ' ---: ' + now)
        b_ar = ds.ReadAsArray(0, blocks[i-1], nc, blocks[i]-blocks[i-1])
        b_ar = b_ar.astype(float)
        b_xmin = xmin; b_ymax = ymax + blocks[i-1]*yres
        b_gt = [b_xmin, xres, xskew, b_ymax, yskew, yres] 
        b_ar[(b_ar<lb) | (b_ar>ub)] = np.nan
        #-calculate statistics
        f = open(log_file, 'r')
        fnames = f.readlines()
        f.close()
        fn = band + '_Min_' + str(i).zfill(zf) + '.tif'
        if (('Min' in stats) & ((fn + '\n') not in fnames)):
            print('Min'); tls.Log(log_dir, 'st_Stats_Sliced_' + band, fn)
            ar = np.nanmin(b_ar, axis=0); ar[np.isnan(ar)] = -9999
            fn = o_dir + '\\' +  fn
            Write_Raster(fn, ar, b_gt, pj)
        f = open(log_file, 'r')
        fnames = f.readlines()
        f.close()
        fn = band + '_Max_' + str(i).zfill(zf) + '.tif'
        if (('Max' in stats) & ((fn + '\n') not in fnames)):
            print('Max'); tls.Log(log_dir, 'st_Stats_Sliced_' + band, fn)
            ar = np.nanmax(b_ar, axis=0); ar[np.isnan(ar)] = -9999
            fn = o_dir + '\\' +  fn
            Write_Raster(fn, ar, b_gt, pj)
        f = open(log_file, 'r')
        fnames = f.readlines()
        f.close()
        fn = band + '_Avg_' + str(i).zfill(zf) + '.tif'
        if (('Avg' in stats) & ((fn + '\n') not in fnames)):
            print('Avg'); tls.Log(log_dir, 'st_Stats_Sliced_' + band, fn)
            ar = np.nanmean(b_ar, axis=0); ar[np.isnan(ar)] = -9999
            fn = o_dir + '\\' +  fn
            Write_Raster(fn, ar, b_gt, pj)
        f = open(log_file, 'r')
        fnames = f.readlines()
        f.close()
        fn = band + '_Std_' + str(i).zfill(zf) + '.tif'
        if (('Std' in stats) & ((fn + '\n') not in fnames)):
            print('Std'); tls.Log(log_dir, 'st_Stats_Sliced_' + band, fn)
            ar = np.nanstd(b_ar, axis=0); ar[np.isnan(ar)] = -9999
            fn = o_dir + '\\' +  fn
            Write_Raster(fn, ar, b_gt, pj)
        for stat in stats:
            if 'P' in stat:
                f = open(log_file, 'r')
                fnames = f.readlines()
                f.close()
                percentile = stat[1:3]
                fn = band + '_' + stat + '_' + str(i).zfill(zf) + '.tif'
                if ((fn + '\n') not in fnames):
                    print(stat); tls.Log(log_dir, 'st_Stats_Sliced_' + band, fn)
                    ar = np.nanpercentile(b_ar, int(percentile), axis=0)
                    ar[np.isnan(ar)] = -9999 
                    fn = o_dir + '\\' +  fn
                    Write_Raster(fn, ar, b_gt, pj)
        #-end statistical calculation
        now = datetime.datetime.now()
        now = now.strftime("%Y-%m-%d %H:%M")
        print('E--- ' + str(i).zfill(zf) + '/' + str(nb-1) + ' ---: ' + now)

def st_Stats_Combined(i_dir, o_dir, band, stat, o_file):
    fn = o_dir + '\\' + o_file
    if os.path.isfile(fn)==False:    
        tifs = glob.glob(i_dir + '\\*' + band + '_' + stat + '*.tif')
        vt = "/vsimem/stacked.vrt"
        ds = gdal.BuildVRT(vt, tifs, separate=False)
        ar = ds.ReadAsArray()
        gt = ds.GetGeoTransform()
        pj = ds.GetProjection() 
        Write_Raster(fn, ar, gt, pj)
#------------------------------------------------------------------------------
        
             
#-- Extract data from timeseries stack ----------------------------------------      
def st_Extract_Data(xy, xycol, window, fid, i_dir, bname, o_dir):
    c_dim = window[1]; c = xycol[0]               #- x direction
    r_dim = window[0]; r = xycol[1]               #- y direction
    tifs, ds = Virtual_DS(i_dir, pattern=bname)
    xmin, xres, xskew, ymax, yskew, yres = ds.GetGeoTransform()   
    rows = ((xy[:,r] - ymax)/yres) + 1; rows = rows.astype(int)
    cols = ((xy[:,c] - xmin)/xres) + 1; cols = cols.astype(int)   
    n = xy.shape[0]; zf = int(math.log10(n))+1
    if n==1:
        fn = o_dir + '\\' + bname + '_fnames_' + str(fid) + '.txt'
    else:
        fn = o_dir + '\\' + bname + '_fnames.txt'
    np.savetxt(fn, tifs, fmt='%s') 
    for j in range(0, n):           
        row = int(rows[j])
        col = int(cols[j])
        ar = ds.ReadAsArray(col-c_dim, row-r_dim, 2*c_dim + 1, 2*r_dim + 1)
        ar = ar.reshape((ar.shape[0], -1))
        if n==1:
            print(fid)
            fn = o_dir + '\\' + bname + '_' + str(fid) + '.txt'
        else:
            print(str(j).zfill(zf) + ' -- ' + str(n))
            fn = o_dir + '\\' + bname + '_' + str(j+1) + '.txt'
        np.savetxt(fn, ar, fmt='%i')
#------------------------------------------------------------------------------
  

#- Reclassify raster ----------------------------------------------------------
def Reclassify(iRas, fVal, tVal, nVal, oRas, nodata = 255):
    ds = gdal.Open(iRas)
    gt = ds.GetGeoTransform()
    pj = ds.GetProjection()
    ar = ds.ReadAsArray()
    n_ar = np.zeros((ar.shape[0], ar.shape[1])) + nodata
    for i in range(0, len(nVal)):
        fval = fVal.iloc[i]
        tval = tVal.iloc[i]
        nval = nVal.iloc[i]
        print(fval, tval, nval)
        n_ar[(ar>=fval) & (ar<=tval)] = nval
    Write_Raster(oRas, n_ar, gt, pj, dtype=gdal.GDT_Int16)
#------------------------------------------------------------------------------

      
#-- Write an 3D array to a Tif ------------------------------------------------
def Write_Raster(fn, ar, gt, pj, dtype=gdal.GDT_Int16):
    if len(ar.shape)==2:
        n = 1
        nr = ar.shape[0]
        nc = ar.shape[1]
    if len(ar.shape)==3:
        n = ar.shape[0]
        nr = ar.shape[1]
        nc = ar.shape[2]
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(fn, nc, nr, n, dtype, ['COMPRESS=LZW'])
    dataset.SetGeoTransform(gt)  
    dataset.SetProjection(pj)
    if n==1:
        dataset.GetRasterBand(1).WriteArray(ar)
    else:
        for i in range(0, n):
            dataset.GetRasterBand(i+1).WriteArray(ar[i])
    dataset = None
#------------------------------------------------------------------------------

           
###--- CURVE FITTING ----------------------------------------------------------
#------------------------------------------------------------------------------

###-Linear harmonic (Fourier series): a0 + SUM[ai*sin(n*t*i) + bi*cos(n*t*i)]
def Harmonic_1(params, time_ls, yData):
    p = params; times = 2*np.pi*time_ls/365
    yTh = np.zeros(len(times)) + p[0]
    nterms = int((len(p)-1)/2)
    for i in range(nterms):
        ai = p[i*2+1]; bi = p[i*2+2]
        yTh = yTh + ai*np.cos(times*(i+1)) + bi*np.sin(times*(i+1))
    errors = [(y-yt) for y, yt in zip(yData, yTh)]
    return errors

def Opt_H1(params, time_ls, yData):
    opt_params,_ = leastsq(Harmonic_1, params, args=(time_ls, yData))
    return opt_params

###-Non-linear harmonic model: p0 + p1*cos(n*t + p2 + p3*sin(n*t + p4))   
def Harmonic_2(params, time_ls, yData):
    p = params; times = 2*np.pi*time_ls/365
    yTh = p[0] + p[1]*np.cos(times + p[2] + p[3]*np.sin(times + p[4]))
    errors = [(y-yt) for y, yt in zip(yData, yTh)]
    return errors

def Opt_H2(params, time_ls, yData):
    opt_params,_ = leastsq(Harmonic_2, params, args=(time_ls, yData))
    return opt_params

###--- CURVE FITTING ----------------------------------------------------------
#------------------------------------------------------------------------------
    