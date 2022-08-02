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

import datetime, dateutil.parser, fiona, geojson, json, math, os, sys, requests, shutil
import numpy  as np
import pandas as pd

try:
    import gdal, gdalconst, osr
except:
    from osgeo import gdal, gdalconst, osr

from sentinelsat.sentinel import SentinelAPI, read_geojson, geojson_to_wkt
from shapely.geometry import mapping
from shapely.geometry import shape
from six.moves.urllib.parse import urlencode
from six import string_types
from pathlib import Path
from tqdm import tqdm

usr_dir = str(Path.home())
dbx_dir = usr_dir + '\\Dropbox'
log_dir = dbx_dir + '\\ToolBox\\Logs\\Download_Sentinel'
lib_dir = dbx_dir + '\\Toolbox\\External_Libs'

gpt_dir = 'C:\\ESA_SNAP\\bin'
s2c_dir = 'C:\\SEN2COR'
tmp_dir = 'D:\\Temp'

sys.path.append(dbx_dir + '\\Toolbox')
import mytools as tls
import image_processing as ips

user, password = 'lnguyen1985', 'Meobeo5191.'
api = SentinelAPI(user, password, 'https://scihub.copernicus.eu/dhus')

###--- Search and download from Scihub ----------------------------------------
def Scihub_Sentinel_Query(sensor, sdate, edate, geojson, fprint=False):
    footprint = geojson_to_wkt(read_geojson(geojson))
    if sensor=='S1':
        query_response = api.query(footprint,
                             platformname = 'Sentinel-1',
                             sensoroperationalmode = 'IW',
                             producttype = 'GRD',
                             date = (sdate, edate))
    if sensor=='S2':
        query_response = api.query(footprint,
                             platformname = 'Sentinel-2',
                             instrumentshortname = 'MSI',
                             processinglevel = 'Level-1C',
                             date = (sdate, edate))
    query_df = api.to_dataframe(query_response)
    query_df.to_csv(log_dir + '\\Scihub_' + sensor + '_Scenes.csv')
    ##-- Create footprint shapefile
    if fprint==True:
        roi_name = os.path.basename(geojson)
        roi_name = roi_name.split('.')[0]
        out_dir = log_dir + '\\Scihub_' + roi_name + '_' + sensor + '.shp'
        Scihub_Footfrint(query_response, out_dir)
    
def Scihub_Footfrint(query_response, out_dir): 
    import shapely
    shapely.speedups.disable()
    json_collection = api.to_geojson(query_response)
    json_features = json_collection['features']
    ##-- Write to shapefile using fiona
    wgs84 = {'proj': 'longlat', 'ellps': 'WGS84', 'datum': 'WGS84', 'no_defs': True}
    myschema = {'geometry': 'Polygon',
               'properties': {'id': 'int', 'scene': 'str'}}
    with fiona.open(out_dir, 'w',
                    crs = wgs84,
                    driver = 'ESRI Shapefile',
                    schema = myschema) as op:
        for i in range(0, len(json_features)):
            #- Get information from geojson's feature collection
            json_item = json_features[i]
            json_item_geom = json_item['geometry']
            json_item_prop = json_item['properties']
            json_item_id   = json_item_prop['identifier']
            #- Convert to geojson.geometry to shapely.Polygon
            json_polygon = json.dumps(json_item_geom)
            g1 = geojson.loads(json_polygon)
            g2 = shape(g1)              
            #- Write polygon to shapefile
            op.write({'geometry': mapping(g2),
                      'properties': {'id': i, 'scene':json_item_id}})

def Scihub_Download(sensor, ras_dir, idx='all'):
    Scihub_Download_Verification(sensor, ras_dir)
    f1_file = log_dir + '\\Download_' + sensor + '.txt'
    if os.path.isfile(f1_file)==False:
        f1 = open(f1_file, "w") 
        f1.close() 
    f1 = open(f1_file, 'r')
    f1_ls = f1.readlines()
    f1.close()
    f2_file = log_dir + '\\Offline.txt'
    if os.path.isfile(f2_file)==False:
        f2 = open(f2_file, "w") 
        f2.close() 
    f2 = open(f2_file, 'r')
    f2_ls = f2.readlines()
    f2.close()
    flist = f1_ls + f2_ls
    flist = [f[:-1] for f in flist]
    now = datetime.datetime.now()
    now = now.strftime("%Y-%m-%d %H:%M")
    products = pd.read_csv(log_dir + '\\Scihub_' + sensor + '_Scenes.csv')
    index = []
    for i in range(0, len(products)):
        fne = products['title'][i]
        if fne in flist:
            index.append(1)
        else:
            index.append(0)
    index = np.array(index)
    sub_prods = products.iloc[index==0,:]
    print('START-------------------------------')
    print(now)
    if idx=="all":
        for i in range(0, len(sub_prods)):
            print(str(i+1) + '/' + str(sub_prods.shape[0]))         
            fne = sub_prods.iloc[i, 1]
            fne = fne.split('.')[0]
            fid = sub_prods.iloc[i, 0]
            product_info = api.get_product_odata(fid)
            if product_info['Online']:
                print('Product {} is online.'.format(fne))
                api.download(fid, ras_dir, checksum=True)
                tls.Log(log_dir, 'Download_' + sensor, fne)
                now = datetime.datetime.now()
                now = now.strftime("%Y-%m-%d %H:%M")
                print('---end: ' + now)
            else:
                print('Product {} is offline.'.format(fne))
                tls.Log(log_dir, 'Offline', fne)            
    else:
        fid = sub_prods.iloc[idx, 0]
        api.download(fid, ras_dir, checksum=True)
    now = datetime.datetime.now()
    now = now.strftime("%Y-%m-%d %H:%M")
    print('---------------------------------END')

def Scihub_Download_Verification(sensor, ras_dir):
    S_scenes = pd.read_csv(log_dir + '\\Scihub_' + sensor + '_Scenes.csv')
    S_inhand = tls.List_Path(ras_dir)
    S_inhand = [item.split('.')[0] for item in S_inhand]
    count_S = 0; index_S = []
    for i in range(0, S_scenes.shape[0]):
        scene = S_scenes['title'][i]
        if scene not in S_inhand:
            count_S += 1
            index_S.append(0)
        else:
            index_S.append(1)
    index_S = np.array(index_S)
    S_o = S_scenes.iloc[index_S==1,:]
    S = np.array(S_o['title'])
    S_ofile = log_dir + '\\Download_' + sensor + 'S.txt'
    if os.path.exists(S_ofile):
        os.remove(S_ofile)
    np.savetxt(log_dir + '\\Download_' + sensor + '.txt', S, fmt='%s')


###--- Search and download from Creodias --------------------------------------
API_URL = ('http://finder.creodias.eu/resto/api/collections/{collection}'
           '/search.json?maxRecords=1000')
ONLINE_STATUS_CODES = '34|37|0'

def Creodias_Query(collection, start_date=None, end_date=None,
          geometry=None, status=ONLINE_STATUS_CODES, **kwargs):
    query_url = _build_query(
        API_URL.format(collection=collection),
        start_date,
        end_date,
        geometry,
        status,
        **kwargs)
    query_response = {}
    while query_url:
        response = requests.get(query_url)
        response.raise_for_status()
        data = response.json()
        for feature in data['features']:
            query_response[feature['id']] = feature
        query_url = _get_next_page(data['properties']['links'])
    return query_response

def _build_query(base_url, sdate=None, edate=None, geometry=None, status=None, **kwargs):
    query_params = {}
    if sdate is not None:
        sdate = _parse_date(sdate)
        query_params['startDate'] = sdate.isoformat()
    if edate is not None:
        edate = _parse_date(edate)
        edate = _add_time(edate)
        query_params['completionDate'] = edate.isoformat()
    if geometry is not None:
        query_params['geometry'] = _parse_geometry(geometry)
    if status is not None:
        query_params['status'] = status
    for attr, value in sorted(kwargs.items()):
        value = _parse_argvalue(value)
        query_params[attr] = value
    url = base_url
    if query_params:
        url += f'&{urlencode(query_params)}'
    return url

def _get_next_page(links):
    for link in links:
        if link['rel'] == 'next':
            return link['href']
    return False

def _parse_date(date):
    if isinstance(date, datetime.datetime):
        return date
    elif isinstance(date, datetime.date):
        return datetime.datetime.combine(date, datetime.time())
    try:
        return dateutil.parser.parse(date)
    except ValueError:
        raise ValueError('Date {date} is not in a valid format. Use Datetime object or iso string')

def _add_time(date):
    if date.hour == 0 and date.minute == 0 and date.second == 0:
        date = date + datetime.timedelta(hours=23, minutes=59, seconds=59)
        return date
    return date

def _tastes_like_wkt_polygon(geometry):
    try:
        return geometry.replace(", ", ",").replace(" ", "", 1).replace(" ", "+")
    except Exception:
        raise ValueError('Geometry must be in well-known text format')

def _parse_geometry(geom):
    try:
        # If geom has a __geo_interface__
        return shape(geom).wkt
    except AttributeError:
        if _tastes_like_wkt_polygon(geom):
            return geom
        raise ValueError('geometry must be a WKT polygon str or have a __geo_interface__')

def _parse_argvalue(value):
    if isinstance(value, string_types):
        value = value.strip()
        if not any(
            value.startswith(s[0]) and value.endswith(s[1])
            for s in ["[]", "{}", "//", "()"]
        ):
            value.replace(" ", "+")
        return value
    elif isinstance(value, (list, tuple)):
        # Handle value ranges
        if len(value) == 2:
            value = "[{},{}]".format(*value)
            return value
        else:
            raise ValueError(
                "Invalid number of elements in list. Expected 2, received "
                "{}".format(len(value))
            )
    else:
        raise ValueError(
            "Additional arguments can be either string or tuple/list of 2 values"
        )

DOWNLOAD_URL = 'https://zipper.creodias.eu/download'
TOKEN_URL = 'https://auth.creodias.eu/auth/realms/DIAS/protocol/openid-connect/token'

def _get_token(username, password):
    token_data = {
        'client_id': 'CLOUDFERRO_PUBLIC',
        'username': username,
        'password': password,
        'grant_type': 'password'
    }
    response = requests.post(TOKEN_URL, data=token_data).json()
    try:
        return response['access_token']
    except KeyError:
        raise RuntimeError(f'Unable to get token. Response was {response}')

def Creodias_Download(uid='', url='', outfile='D:\Temp\Test.zip', show_progress=True):
    username = 'hoanglan.nguyen@ucalgary.ca',
    password = 'Meobeo5191.'
    if url == '':
        token = _get_token(username, password)
        url = f'{DOWNLOAD_URL}/{uid}?token={token}'
    outfile_temp = str(outfile) + '.incomplete'
    try:
        downloaded_bytes = 0
        with requests.get(url, stream=True, timeout=100) as req:
            with tqdm(unit='B', unit_scale=True, disable=not show_progress) as progress:
                chunk_size = 2 ** 20  # download in 1 MB chunks
                with open(outfile_temp, 'wb') as fout:
                    for chunk in req.iter_content(chunk_size=chunk_size):
                        if chunk:  # filter out keep-alive new chunks
                            fout.write(chunk)
                            progress.update(len(chunk))
                            downloaded_bytes += len(chunk)
        shutil.move(outfile_temp, outfile)
    finally:
        try:
            Path(outfile_temp).unlink()
        except OSError:
            pass   
###----------------------------------------------------------------------------


###-- Preprocessing S1 and S2 imagery -----------------------------------------
def Preprocessing_S1_IW(i_dir, o_dir):
    flist = os.listdir(i_dir)
    nf = len(flist)                             #- number-of-file in 'flist'
    zf = int(math.log10(nf))                    #- zero-fill
    xml_dir = tmp_dir + '\\S1_IW_GRDH_PP.xml'
    shutil.copy(lib_dir + '\\S1_IW_GRDH_PP.xml', xml_dir)
    count = 1
    for file in flist:
        print(str(count).zfill(zf) + '//' + str(nf) + ': ' + file)
        log_file = log_dir + '\\Preprocessing_S1_IW.txt'
        if os.path.isfile(log_file)==False:
            file = open(log_file, 'w')
            file.close()
        f = open(log_file, 'r')
        flist = f.readlines()
        f.close()
        if (file + '\n') in flist:
            print(file + ' was preprocessed !!!')
            count += 1
            continue
        else:
            tls.Log(log_dir, 'Preprocessing_S1_IW', file)
            os.chdir(gpt_dir)
            i_file = i_dir + '\\' + file
            o_file = o_dir + '\\' + file.split('.')[0] + '.dim'
            args = ' -Pifile=' + i_file + ' -Pofile=' + o_file
            cmd = 'gpt.exe ' + xml_dir + args
            os.system(cmd)
            path = o_dir + '\\' + file.split('.')[0] + '.data'
            images = ['Sigma0_VV', 'Sigma0_VH', 'localIncidenceAngle']
            for image in images:
                fn = path + '\\' + image
                ds = gdal.Open(fn + '.img')
                ar = ds.ReadAsArray()
                gt = ds.GetGeoTransform()
                pj = ds.GetProjection()
                if (image=='Sigma0_VV') | (image=='Sigma0_VH'):
                    ar = np.ma.log10(ar)*1000
                    ar = ar.filled(-9999)
                    ar = ar.astype(int)
                if image=='localIncidenceAngle':
                    ar = (ar * 100).astype(int)
                ips.Write_Raster(fn + '.tif', ar, gt, pj)
                ds = None
                os.remove(path + '\\' + image + '.img')
                os.remove(path + '\\' + image + '.hdr')
            shutil.rmtree(path + '\\vector_data')
            count += 1
            
def Preprocessing_S2_SR(i_dir, o_dir):
    flist = os.listdir(i_dir)
    flist = [f for f in flist if (f[-5:]=='.SAFE')]
    nf = len(flist)
    zf = int(math.log10(nf))
    count = 1
    for file in flist:
        print(str(count).zfill(zf) + '//' + str(nf) + ': ' + file)
        log_file = log_dir + '\\Preprocessing_S2_SR.txt'
        if os.path.isfile(log_file)==False:
            f = open(log_file, 'w')
            f.close()
        f = open(log_file, 'r')
        fls = f.readlines()
        f.close()
        fls = [f[:-1] for f in fls]
        if (file in fls):
            print('This tile was processed !!!')
            count += 1
            continue
        else:
            print('This tile is processing !!!')
            tls.Log(log_dir, 'Preprocessing_S2_SR', file)
            os.chdir(s2c_dir)
            i_file = i_dir + '\\' + file
            cmd = 'L2A_Process.bat --resolution 10 --output_dir '
            cmd = cmd + o_dir + ' ' + i_file
            os.system(cmd)          
            count += 1

def Preprocessing_2(sensor, r_ras, i_dir, o_dir, pattern):
    S2_bands = ['B02_10m', 'B03_10m', 'B04_10m', 'B08_10m', \
            'B05_20m', 'B06_20m', 'B07_20m', 'B8A_20m', \
            'B11_20m', 'B12_20m']
    S2_range = [[0, 12000]]*len(S2_bands) 
    
    S1_bands = ['Sigma0_VV', 'Sigma0_VH', 'localIncidenceAngle']
    S1_range = [[-3000, 1000], [-3000, 1000], [1000, 7000]]
    count = 1
    flist = os.listdir(i_dir)
    flist = [item for item in flist if item.split('.')[1]==pattern]
    for f in flist:
        print(str(count).zfill(3) + '|' + str(len(flist)) + ': ' + f)
        log_file = log_dir + '\\' + sensor + '_Preprocessed_2.txt'
        if os.path.isfile(log_file)==False:
            file = open(log_file, 'w')
            file.close()
        file = open(log_file, 'r')
        fnames = file.readlines()
        file.close()
        if (f + '\n') in fnames:
            print(f + ' was preprocessed !!!')
            count = count + 1
            continue
        else:
            tls.Log(log_dir, sensor + '_Preprocessed_2', f)
            now = datetime.datetime.now()
            now = now.strftime("%Y-%m-%d %H:%M")
            print('S---: ' + now)
            strs = f.split('_')
            if sensor=='S1':
                f1 = strs[0] + '_' + strs[4] + '_' + strs[8].split('.')[0]
                for band in S1_bands:
                    print(band)
                    i_ras = i_dir + '\\' + f + '\\' + band + '.tif'
                    o_ras = o_dir + '\\' + f1 + '_' + band + '.tif'
                    t_ras = tmp_dir  + '\\Temp_' + f1 + '.tif'
                    ips.Projection(i_ras, t_ras, r_ras)
                    lb = S1_range[S1_bands.index(band)][0]
                    ub = S1_range[S1_bands.index(band)][1]
                    op = ips.Subset(t_ras, roi=[r_ras], data=True, value_range=[lb,ub])
                    if len(op)==3:
                        ar, gt, pj = op
                        ips.Write_Raster(o_ras, ar, gt, pj)
                    os.remove(t_ras)
            if sensor=='S2':
                l1_dir = i_dir + '\\' + f + '\\GRANULE'
                l2_dir = l1_dir + '\\' + os.listdir(l1_dir)[0] + '\\IMG_DATA'
                i_ras = strs[5] + '_' + strs[2] + '_SCL_20m.jp2'
                i_ras = l2_dir + '\\R20m\\' + i_ras
                op = ips.SubsetI(i_ras, roi=[r_ras], data=True, value_range=[4,7])
                if len(op)==3:
                    s_ar, s_gt, s_pj = op
                    s_ar[(s_ar<=3) | (s_ar>=8)] = 0
                    s_ar[s_ar>0] = 1
                    for band in S2_bands:
                        print(band)
                        bres = band.split('_')[1]
                        ip = strs[5] + '_' + strs[2] + '_' + band + '.jp2'
                        i_ras = l2_dir + '\\R' + bres + '\\' + ip
                        o_ras = o_dir + '\\' + strs[5] + '_' + strs[2] + '_' + band.split('_')[0] + '.tif'
                        lb = S2_range[S2_bands.index(band)][0]
                        ub = S2_range[S2_bands.index(band)][1]
                        ar, gt, pj = ips.Subset(i_ras, roi=[r_ras], data=True, value_range=[lb,ub])
                        ar[s_ar==0] = 0
                        ips.Write_Raster(o_ras, ar, gt, pj)
            count = count + 1
            now = datetime.datetime.now()
            now = now.strftime("%Y-%m-%d %H:%M")
            print('---E: ' + now)
#------------------------------------------------------------------------------