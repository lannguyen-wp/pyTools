#------------------------------------------------------------------------------
#-- / Description -------------------------------------------------------------


#-- Script  : mytools.py
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
#-- valid  : data validation [vrange, min_data_proportion]
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

import datetime, os
import numpy as np

from pathlib import Path

usr_dir = str(Path.home())
dbx_dir = usr_dir + '\\Dropbox'
lib_dir = dbx_dir + '\\Toolbox\\External_Libs'
log_dir = dbx_dir + '\\Toolbox\\Logs'
test_dir = dbx_dir + '\\Toolbox\\Testing'


#-- Add information to a log file ---------------------------------------------
def Log(log_dir, log_file, log_info, overwrite=False):
    if log_file[-3:]=='txt':
        file = log_dir + '\\' + log_file
    else:
        file = log_dir + '\\' + log_file + '.txt'
    if overwrite==True:
        f = open(file, 'w')
    if overwrite==False:
        f = open(file, 'a')
    f.write(str(log_info) + '\n')
    f.close()


#-- Get DOYs from a list of file names ---------------------------------
def DOY_from_Fname(flist, s_pos, e_pos, date_fmt, path_mode='File'):
    if path_mode=='File':
        flist = flist
    if path_mode=='Full':
        flist = [os.path.basename(f) for f in flist]
    dates = [f[s_pos:e_pos+1] for f in flist] 
    dates = [datetime.datetime.strptime(date, date_fmt) for date in dates]
    doy_ls = [date.timetuple().tm_yday for date in dates]
    doy_ls = np.array(doy_ls)
    return doy_ls


###--- List files in a directory with particular patterns ---------------------  
def List_Path(i_dir, patterns=['All'], o_file='N/A', o_dir=log_dir, mode='File'):
    i_ls = os.listdir(i_dir)
    o_ls = []
    if patterns[0]=='All':
        o_ls = i_ls
    else:
        for pattern in patterns:
            o_ls = o_ls + [item for item in i_ls if pattern in item]
    if o_file!='N/A':
        o_file = o_dir + '\\' + o_file
        f = open(o_file, 'w')
        for item in o_ls:
            f.write(item + '\n')  
        f.close()
    if mode=='Full':
        o_ls = [i_dir + '\\' + item for item in o_ls]
    return o_ls


def _remove_temp(log_dir, text):
    ls = []
    f = open(log_dir + '\\Temp.txt','r')
    for line in f:
        if text not in line:
            ls.append(text)
    f.close()
    f = open(log_dir + '\\Temp.txt','w')
    for line in ls:
        f.write(line)
        f.close()
    