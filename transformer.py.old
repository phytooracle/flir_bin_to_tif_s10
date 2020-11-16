#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Michele Cosi, Holly Ellingson, Jeffrey Demieville
Note   : Parts of this code was initially developed by the AgPipeline and TERRA-REF teams.
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
from terrautils.spatial import scanalyzer_to_utm, geojson_to_tuples
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
                        help='Cleaned metadata file',
                        metavar='metadata',
                        type=str,
                        required=True)
                        #default='cleanmetadata_out')

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
def flirRawToTemperature(rawData, meta):
    
    # Normalizing the tempearature measurements according to the Season 10 calibration
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
    
    # Obtain shutter temperature from metadata
    K0 = 273.15
    H = 0.1
    T = int(float(meta['sensor_variable_metadata']['shutter_temperature_K']) - 273.15)
    # T = 22
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

    pxl_temp = B/np.log(R/(corr_pxl_val-J0)*J1+F)#.astype(int)
    #pxl_temp = pxl_temp.astype(int)

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

            full_md = json.load(mdf)['content']
            extractor_info = None

            if full_md:
                if bin_file is not None:
                    out_file = os.path.join(args.outdir, bin_file.split('/')[-1].replace(".bin", ".tif"))

                    gps_bounds_bin = geojson_to_tuples(
                        full_md['spatial_metadata']['flirIrCamera']['bounding_box'])
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
