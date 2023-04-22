import glob
import pandas as pd
from snakemake.utils import Paramspace



dicoms = (
    pd.read_csv(config["dicoms"], 
                sep="\t",
                dtype={"fullpath": str, 
                      "sub": str, 
                      "ses": str,
                      "target": str, 
                      "scan": str} 
                )
    .set_index("fullpath",drop=False)
 )

paramspace = Paramspace(dicoms,filename_params='*')

# def get_converted_nii():
#     converted_nii = expand(
#        f"{base}/{sub}/{ses}/dwi/{scan}",
#         zip,
#         base=config["BIDSdatasetpath"],

#     )
#     return converted_nii

# def get_nii_target(wildcards):

    