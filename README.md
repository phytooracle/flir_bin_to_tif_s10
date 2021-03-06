# Flir bin2tif for s10 (Lettuce)

Re-calibrated transformer: converts bin to tif for gantry files (s10). Updated for geo-correction.

## Inputs

Raw metadata, raw flir compressed (.bin) files

## Outputs

Decompressed flir geotiffs from `.bin` to `.tif`.

## Arguments and Flags
- **Required Arguments:** 
    - **Directory of flir files to decompress:** 'bin'
    - **raw metadata files:** '-m', '--metadata' 

- **Optional Arguments**
    - **Output directory:** '-o', '--outdir', default='flir2tif_out/'

## Adapting the Script
- Metadata has to match name of input flir bin such as `example_0.bin`, `example_0_metadata.bin`    
                                        
## Executing example (using singularity)
`singularity run -B $(pwd):/mnt --pwd /mnt/ docker://phytooracle/flir_bin_to_tif_s10 -m <metadata.json> <bin_dir>`
