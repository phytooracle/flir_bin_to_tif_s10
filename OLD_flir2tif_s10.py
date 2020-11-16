#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Michele Cosi, Holly Ellingson, Jeffrey Demieville
Date   : 2020-07-09
Purpose: Convert FLIR .bin files to .tif
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


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Season 10 flir2tif',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('dir',
                        metavar='str',
                        help='Directory containing bin and metadata files')

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory where .tif files will be saved',
                        metavar='str',
                        type=str,
                        default='flir2tif_out')

    args = parser.parse_args()

    if '/' not in args.dir:
        args.dir = args.dir + '/'

    if '/' not in args.outdir:
        args.outdir = args.outdir + '/'

    return args



# --------------------------------------------------
def get_calibrate_param(metadata):
    calibparameter = calibParam()
    try:
        if 'terraref_cleaned_metadata' in metadata:
            fixedmd = metadata['sensor_fixed_metadata']
            # added this to pull in shutter temperature
            variablemd = metadata['sensor_variable_metadata']
            if fixedmd['is_calibrated'] == 'True':
                return calibparameter
            else:
                calibparameter.calibrated = False
                calibparameter.calibrationR = float(fixedmd['calibration_R'])
                calibparameter.calibrationB = float(fixedmd['calibration_B'])
                calibparameter.calibrationF = float(fixedmd['calibration_F'])
                calibparameter.calibrationJ1 = float(fixedmd['calibration_J1'])
                calibparameter.calibrationJ0 = float(fixedmd['calibration_J0'])
                calibparameter.calibrationa1 = float(
                    fixedmd['calibration_alpha1'])
                calibparameter.calibrationa2 = float(
                    fixedmd['calibration_alpha2'])
                calibparameter.calibrationX = float(fixedmd['calibration_X'])
                calibparameter.calibrationb1 = float(
                    fixedmd['calibration_beta1'])
                calibparameter.calibrationb2 = float(
                    fixedmd['calibration_beta2'])
                # To-Do: verify variable names to match what is present in cleaned metadata
                calibparameter.shuttertemperature = float(
                    variablemd['shutter_temperature_[K]'])
                return calibparameter
    except KeyError as err:
        return calibparameter


# --------------------------------------------------
def flirRawToTemperature(rawData, calibP):

    R = calibP.calibrationR
    B = calibP.calibrationB
    F = calibP.calibrationF
    J0 = calibP.calibrationJ0
    J1 = calibP.calibrationJ1

    X = calibP.calibrationX
    a1 = calibP.calibrationa1
    b1 = calibP.calibrationb1
    a2 = calibP.calibrationa2
    b2 = calibP.calibrationb2

    H2O_K1 = 1.56
    H2O_K2 = 0.0694
    H2O_K3 = -0.000278
    H2O_K4 = 0.000000685

    H = 0.1
    shutter_temp = calibP['content']['sensor_variable_metadata']['shutter_temperature_K']
    T = float(shutter_temp) 
    D = 2.5
    E = 0.98

    K0 = 273.15

    im = rawData

    AmbTemp = T
    AtmTemp = T

    H2OInGperM2 = H*math.exp(H2O_K1 + H2O_K2*T + H2O_K3*math.pow(T, 2) + H2O_K4*math.pow(T, 3))
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
    pxl_temp = pxl_temp.astype(int)
   
    return pxl_temp


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()
    files_path = args.dir

    bin_files = glob.glob(f'{files_path}*/*.bin')

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    for bin_file in bin_files:
        meta_file_cleaned = bin_file.replace(
            '_ir.bin', '_metadata_cleaned.json')

        with open(meta_file_cleaned, 'r') as mdf:
            full_md = json.load(mdf)
            extractor_info = None

            if full_md:
                if bin_file is not None:
                    out_file = os.path.join(args.outdir, bin_file.split('/')[-1].replace(".bin", ".tif"))
                    gps_bounds_bin = geojson_to_tuples(
                        full_md['content']['spatial_metadata']['flirIrCamera']['bounding_box'])
                    raw_data = np.fromfile(bin_file, np.dtype('<u2')).reshape(
                        [480, 640]).astype('float')
                    raw_data = np.rot90(raw_data, 3)
                    tc = flirRawToTemperature(raw_data, full_md)
                    print(tc)

                    create_geotiff(tc, gps_bounds_bin, out_file, None,
                                False, extractor_info, None, compress=True)


# --------------------------------------------------
if __name__ == '__main__':
    main()
