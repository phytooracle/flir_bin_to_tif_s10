#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Michele Cosi, Holly Ellingson, Jeffrey Demieville
Date   : 2020-07-09
Purpose: Convert FLIR .bin files to .tif (Season 10)
"""

import argparse
import os
import sys
import logging
import json
import numpy as np
import glob
from terrautils.spatial import scanalyzer_to_utm, scanalyzer_to_latlon, geojson_to_tuples
from terrautils.formats import create_geotiff
import matplotlib.pyplot as plt
from osgeo import gdal, osr
import math
from numpy.matlib import repmat


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Season 10 flir2tif',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bin',
                        metavar='str',
                        help='Bin file to be converted to TIF')

    parser.add_argument('-m',
                        '--metadata',
                        help='Raw metadata file',
                        metavar='metadata',
                        type=str,
                        required=True)
                        #default='cleanmetadata_out')
    parser.add_argument('-z',
                        '--zoffset',
                        help='Z-axis offset',
                        metavar='z-offset',
                        type=int,
                        default=0.76)

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory where .tif files will be saved',
                        metavar='str',
                        type=str,
                        default='flir2tif_out')

    args = parser.parse_args()

    if '/' not in args.outdir:
        args.outdir = args.outdir + '/'

    return args


# --------------------------------------------------
def get_boundingbox(metadata, z_offset):

    with open(metadata) as f:
        meta = json.load(f)['lemnatec_measurement_metadata']

    loc_gantry_x = float(meta['sensor_fixed_metadata']['location in camera box x [m]'])
    loc_gantry_y = float(meta['sensor_fixed_metadata']['location in camera box y [m]'])
    loc_gantry_z = float(meta['sensor_fixed_metadata']['location in camera box z [m]'])

    gantry_x = float(meta['gantry_system_variable_metadata']['position x [m]']) + loc_gantry_x
    gantry_y = float(meta['gantry_system_variable_metadata']['position y [m]']) + loc_gantry_y
    gantry_z = float(meta['gantry_system_variable_metadata']['position z [m]']) + z_offset + loc_gantry_z#offset in m

    fov_x, fov_y = float(meta['sensor_fixed_metadata']['field of view x [m]']), float(meta['sensor_fixed_metadata']['field of view y [m]'])
    
    img_height, img_width = 640, 480
    
    B = gantry_z
    A_x = np.arctan((0.5*float(fov_x))/2)
    A_y = np.arctan((0.5*float(fov_y))/2)
    L_x = 2*B*np.tan(A_x)
    L_y = 2*B*np.tan(A_y)

    x_n = gantry_x + (L_x/2)
    x_s = gantry_x - (L_x/2)
    y_w = gantry_y + (L_y/2)
    y_e = gantry_y - (L_y/2)

    bbox_nw_latlon = scanalyzer_to_latlon(x_n, y_w)
    bbox_se_latlon = scanalyzer_to_latlon(x_s, y_e)

    # TERRA-REF
    lon_shift = 0.000020308287

    # Drone
    lat_shift = 0.000018292 #0.000015258894
    b_box =  ( bbox_se_latlon[0] - lat_shift,
                bbox_nw_latlon[0] - lat_shift,
                bbox_nw_latlon[1] + lon_shift,
                bbox_se_latlon[1] + lon_shift)

    return b_box, img_height, img_width


# --------------------------------------------------
def flirRawToTemperature(rawData, meta):

    R = 16923.6
    B = 1434.6
    F = 1
    J0 = 4124
    J1 = 70.0704

    X = 1.9
    a1 = 0.006569
    b1 = -0.002276
    a2 = 0.01262
    b2 = -0.00667

    H2O_K1 = 1.56
    H2O_K2 = 0.0694
    H2O_K3 = -0.000278
    H2O_K4 = 0.000000685

    K0 = 273.15
    H = 0.1

    if 'shutter temperature [K]' in meta['sensor_variable_metadata']:
        T = int(float(meta['sensor_variable_metadata']['shutter temperature [K]']) - 273.15) 
    else: 
        T = int(float(meta['sensor_variable_metadata']['shutter_temperature_K']) - 273.15)
    
    D = 2.5
    E = 0.98

    im = rawData

    # AmbTemp = T
    # AtmTemp = T
    AmbTemp = T + K0
    AtmTemp = T + K0

    H2OInGperM2 = H*math.exp(H2O_K1 + H2O_K2*T + H2O_K3*math.pow(T, 2) + H2O_K4*math.pow(T, 3))
    #H2OInGperM2 = H*math.exp(H2O_K1 + H2O_K2*T + H2O_K3*(T**2) + H2O_K4*(T**3))
    a1b1sqH2O = (a1+b1*math.sqrt(H2OInGperM2))
    a2b2sqH2O = (a2+b2*math.sqrt(H2OInGperM2))
    exp1 = math.exp(-math.sqrt(D/2)*a1b1sqH2O)
    exp2 = math.exp(-math.sqrt(D/2)*a2b2sqH2O)

    tao = X*exp1 + (1-X)*exp2

    obj_rad = im*E*tao

    theo_atm_rad = (R*J1/(math.exp(B/AtmTemp)-F)) + J0
    atm_rad = repmat((1-tao)*theo_atm_rad, 640, 480)

    theo_amb_refl_rad = (R*J1/(math.exp(B/AmbTemp)-F)) + J0
    amb_refl_rad = repmat((1-E)*tao*theo_amb_refl_rad, 640, 480)

    corr_pxl_val = obj_rad + atm_rad + amb_refl_rad

    pxl_temp = B/np.log(R/(corr_pxl_val-J0)*J1+F)

    return pxl_temp


# --------------------------------------------------
def main():
    """Create TIF here"""

    args = get_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    bin_file = args.bin
    if bin_file is not None:
        with open(args.metadata, 'r') as mdf:

            full_md = json.load(mdf)['lemnatec_measurement_metadata']#['content']
            extractor_info = None

            if full_md:
                if bin_file is not None:
                    out_file = os.path.join(args.outdir, bin_file.split('/')[-1].replace(".bin", ".tif"))

                    #gps_bounds_bin = geojson_to_tuples(
                        #full_md['spatial_metadata']['flirIrCamera']['bounding_box'])
                    gps_bounds_bin, img_height, img_width = get_boundingbox(args.metadata, args.zoffset)
                    
                    raw_data = np.fromfile(bin_file, np.dtype('<u2')).reshape(
                        [480, 640]).astype('float')
                    raw_data = np.rot90(raw_data, 3)

                    tc = flirRawToTemperature(raw_data, full_md)
                    create_geotiff(tc, gps_bounds_bin, out_file, None,
                                True, extractor_info, None, compress=True)

                    print(f'Done. See output in {args.outdir}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
