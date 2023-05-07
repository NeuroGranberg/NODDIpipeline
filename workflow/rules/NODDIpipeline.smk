import glob
import os
import pandas as pd
from pathlib import Path
import snakemake.io as snake

# temporary dictionarys for running isolated debugging
config = {'dicoms': 'config/dicomstest.tsv', 
          'BIDSdatasetpath': '/home/oscar/testdicom/testdicom'}

dicoms = (
    pd.read_csv(config["dicoms"], 
                sep="\t",
                dtype={"fullpath": str, 
                      "sub": str, 
                      "ses": str,
                      "target": str, 
                      "scan": str} 
                )
 )

dicoms["acq"]=dicoms["scan"].str.replace(r"[()_\\n]","",regex=True)
suffixes={"dwi":"dwi.nii.gz","fmap":"PA_epi.nii.gz"}
dicoms["suffix"]=dicoms["target"].map(suffixes.get)
dicoms = dicoms.replace('\n','',regex=True)
dicompath = dicoms.loc[:,('sub','ses','scan','target','acq','suffix')]
lookupdir = pd.Series(dicompath.scan.values,index=dicompath.acq).to_dict()

#lookup helpers for rule inputs
def get_dicomtarget(wildcards):
    path=f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/" + lookupdir[wildcards.acq]
    iter=glob.iglob(path+"/*.dcm")
    return Path(next(iter,None)).name

#workaround for single wildcard object, also no multiple matches required in dicompath.tsv
#write decorator or return a function with single wildcard, ie currying
def get_acq_from_sub_ses(wildcards,acqmask,suffixstr,bidsdir):
    pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    mask=pd.Series(pattern,dtype='string').str.contains(acqmask)
    return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/{bidsdir}/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_{suffixstr}"

def get_acq_from_sub_ses_b1000(wildcards):
     return get_acq_from_sub_ses(wildcards,acqmask="b1000",suffixstr="dwi.nii.gz",bidsdir="dwi")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"b1000")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_dwi.nii.gz"

def get_acq_from_sub_ses_b1000_bvec(wildcards):
    return get_acq_from_sub_ses(wildcards,acqmask="b1000",suffixstr="dwi.voxel_space.bvecs",bidsdir="dwi")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"b1000")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_dwi.voxel_space.bvecs"

def get_acq_from_sub_ses_b1000_bvals(wildcards):
    return get_acq_from_sub_ses(wildcards,acqmask="b1000",suffixstr="dwi.bvals",bidsdir="dwi")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"b1000")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_dwi.bvals"

def get_acq_from_sub_ses_b3000(wildcards):
    return get_acq_from_sub_ses(wildcards,acqmask="b3000",suffixstr="dwi.nii.gz",bidsdir="dwi")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"b3000")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_dwi.nii.gz"

def get_acq_from_sub_ses_b3000_bvec(wildcards):
    return get_acq_from_sub_ses(wildcards,acqmask="b3000",suffixstr="dwi.voxel_space.bvecs",bidsdir="dwi")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"b3000")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_dwi.voxel_space.bvecs"

def get_acq_from_sub_ses_b3000_bvals(wildcards):
    return get_acq_from_sub_ses(wildcards,acqmask="b3000",suffixstr="dwi.bvals",bidsdir="dwi")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"b3000")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/dwi/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_dwi.bvals"

def get_acq_from_sub_ses_PA(wildcards):
    return get_acq_from_sub_ses(wildcards,acqmask="PA",suffixstr="PA_epi.nii.gz",bidsdir="fmap")
    #pattern=dicompath.loc[(dicompath["sub"] == wildcards.sub) & (dicompath["ses"] == wildcards.ses), "acq"]
    #mask=pd.Series(pattern,dtype='string').str.contains(r"PA")
    #return f"{config['BIDSdatasetpath']}/{wildcards.sub}/{wildcards.ses}/fmap/{wildcards.sub}_{wildcards.ses}_acq-{pattern[mask].values[0]}_PA_epi.nii.gz"


testw = snake.Wildcards()
print(testw)

#build targets
def get_converted_nii(): 
    opattern = config['BIDSdatasetpath']+"/{sub}/{ses}/{target}/{sub}_{ses}_acq-{acq}_{suffix}"
    output_patterns =[
    opattern.format(**dict(i for i in row.items())) 
    for _, row in dicompath.iterrows()]
    return output_patterns

def get_NODDI_protocols():
    opattern = f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/protocol.prtcl"
    output_patterns =[
    opattern.format(**dict(i for i in row.items())) 
    for _, row in dicompath.iterrows()]
    return output_patterns

def get_preprocessing():
    opattern = f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/bvals"
    output_patterns =[
    opattern.format(**dict(i for i in row.items())) 
    for _, row in dicompath.iterrows()]
    return output_patterns


rule all:
    input:
        get_NODDI_protocols()

#bvec and bvals are side effects
# rule nii_conversion:
#     input:
#          lambda wc: f"{config['BIDSdatasetpath']}/{{sub}}/{{ses}}/dwi/" + lookupdir[wc.acq]
#     output:
#          f"{config['BIDSdatasetpath']}/{{sub}}/{{ses}}/{{target}}/{{sub}}_{{ses}}_acq-{{acq}}_{{suffix}}"
#     wildcard_constraints:
#         acq="[^_]+"
#     params:
#         dicomtarget=get_dicomtarget
#     shell:
#         "cd {input} && mri_convert -ot nii {params.dicomtarget} {output}"


#dcm2niix -z y -f {output} {input}

#checkpoint needed to revaluate DAG to capture bvec, bvals, ie run rule all sequentially calling 
#get_converted_nii() followed by get_NODDI_protocols()

rule image_for_topup:
    input:
        b1000=get_acq_from_sub_ses_b1000,
        PA=get_acq_from_sub_ses_PA
    output: f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/Diff_32_0_AP_PA.nii.gz"
    shadow:"shallow"
    shell:
        r"""
        fslroi {input.b1000} Diff_32_0_AP.nii.gz 0 2 
        fslroi {input.b1000} Diff_32_0_APfix.nii.gz 0 1 
        fslroi {input.PA} Diff_0_G3_PA0.nii.gz 0 1 
        fslroi {input.PA} Diff_0_G3_PA1.nii.gz 1 1 
        flirt -in Diff_0_G3_PA0.nii.gz -ref Diff_32_0_APfix.nii.gz -out Diff_0_G3_PAfixed0.nii.gz 
        flirt -in Diff_0_G3_PA1.nii.gz -ref Diff_32_0_APfix.nii.gz -out Diff_0_G3_PAfixed1.nii.gz 
        fslmerge -t {output} Diff_32_0_AP.nii.gz Diff_0_G3_PAfixed0.nii.gz Diff_0_G3_PAfixed1.nii.gz 
        rm Diff_32_0_AP.nii.gz Diff_32_0_APfix.nii.gz Diff_0_G3_PA0.nii.gz Diff_0_G3_PA1.nii.gz Diff_0_G3_PAfixed0.nii.gz Diff_0_G3_PAfixed1.nii.gz 
        """

#side effects
rule topup:
    input:
        b0appa= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/Diff_32_0_AP_PA.nii.gz",
        acqparams= "config/acqparams.txt"
    output:
        image= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/my_unwarped_images.nii.gz",
        topupim= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/my_topup_results_fieldcoef.nii.gz",
        topuptxt= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/my_topup_results_movpar.txt"
    shadow:"shallow"
    shell:
        r"""
        topup --imain={input.b0appa} --datain={input.acqparams} --config=b02b0.cnf --out=my_topup_results --iout=my_unwarped_images
        mv my_unwarped_images.nii.gz {output.image}
        mv my_topup_results_fieldcoef.nii.gz {output.topupim}
        mv my_topup_results_movpar.txt {output.topuptxt}
        """

#side effects
rule brainmask:
    input:
        im= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/my_unwarped_images.nii.gz"
    output:
        mask= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/brain_mask.nii.gz"
    params:
        bet_thres= "0.3"
    shadow:"shallow"
    shell:
        r"""
        fslmaths {input.im} -Tmean MeanB0.nii.gz
        bet MeanB0.nii.gz brain -m -R -f {params.bet_thres}
        rm MeanB0.nii.gz brain.nii.gz
        mv brain_mask.nii.gz {output.mask}
        """
# fails bash strict mode ie set -euo pipefail
rule merge_shells:
    input:
        b1000=  get_acq_from_sub_ses_b1000,
        b1000bvecs= get_acq_from_sub_ses_b1000_bvec,
        b1000bvals= get_acq_from_sub_ses_b1000_bvals,
        b3000= get_acq_from_sub_ses_b3000,
        b3000bvecs= get_acq_from_sub_ses_b3000_bvec,
        b3000bvals= get_acq_from_sub_ses_b3000_bvals
    output:
        dwi= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.nii.gz",
        bvecs= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.bvecs",
        bvals= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.bvals",
        index= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/index.txt"
    shadow:"shallow"
    shell:
        r"""    
        fslmerge -t {output.dwi} {input.b1000} {input.b3000}
        cat {input.b1000bvecs} {input.b3000bvecs} > {output.bvecs}
        cat {input.b1000bvals} {input.b3000bvals} > {output.bvals}
        indx=""
        lines=`wc -l < {output.bvals}`
        for ((i=1;i<=$lines;i+=1)); do indx="$indx 1";done
        echo $indx > {output.index}      
        """  
        #fix wildcard string substitution
    # run:
    #     shell("fslmerge -t {output.dwi} {input.b1000} {input.b3000}")
    #     shell("cat {input.b1000bvecs} {input.b3000bvecs} > {output.bvecs}")
    #     shell("cat {input.b1000bvals} {input.b3000bvals} > {output.bvals}")
    #     with open({output.bvals}, 'r') as f:
    #         lines = sum(1 for line in f)
    #     indx = " ".join("1" for i in range(lines))
    #     with open({output.index}, 'w') as f:
    #         f.write(indx)

#side effects
rule eddy:
    input:
        imain= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.nii.gz",
        mask= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/brain_mask.nii.gz",
        acqp= "config/acqparams.txt",
        index= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/index.txt",
        bvecs= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.bvecs",
        bvals= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.bvals",
        topupim= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/my_topup_results_fieldcoef.nii.gz",
        topuptxt= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/my_topup_results_movpar.txt",
    output: 
        image= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/eddy_corrected_data.nii.gz",
        bvecs= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/eddy_rotated_bvecs"
    shadow:"shallow"
    shell:
        r"""
        cp config/acqparams.txt {config[BIDSdatasetpath]}/derivatives/NODDI_MDT/{wildcards.sub}/{wildcards.ses}/preproc/acqparams.txt
        cd {config[BIDSdatasetpath]}/derivatives/NODDI_MDT/{wildcards.sub}/{wildcards.ses}/preproc/
        eddy_cuda10.2 --imain={input.imain} --mask={input.mask} --acqp=acqparams.txt --index={input.index} --bvecs={input.bvecs} --bvals={input.bvals} --topup=my_topup_results --out=eddy_corrected_data1
        mv eddy_corrected_data1.nii.gz {output.image}
        mv eddy_corrected_data1.eddy_rotated_bvecs {output.bvecs}
        """

AWK_CMD =r"""{printf "%i ",$0}"""
rule bvals_for_mdt:
     input:
         bvals= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/DWI.bvals",
         script= "workflow/scripts/transpose.sh" 
     output:
         bvals= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/bvals"     
     shadow:"shallow"
     shell:
         r"""
         awk {AWK_CMD:q} {input.bvals}> bvals1
         {input.script} bvals1 > bvals2
         {input.script} bvals2 > {output.bvals}
         rm bvals1 bvals2
         """
# #move run section to script to avoid import mdt
rule mdt_fit:
     input:
         bvecs= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/eddy_rotated_bvecs",
         bvals= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/bvals",
         dwi= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/eddy_corrected_data.nii.gz",
         mask= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/preproc/brain_mask.nii.gz"
     output:
         protocol= f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/protocol.prtcl",
         outputdir= directory(f"{config['BIDSdatasetpath']}/derivatives/NODDI_MDT/{{sub}}/{{ses}}/NODDI")

     run:
         import mdt
         protocol = mdt.create_protocol(
         bvecs= {input.bvecs}, 
         bvals= {input.bvals},
         out_file= {output.protocol})

         input_data = mdt.load_input_data(
             {input.dwi},
             {output.protocol},
             {input.mask},
             noise_std=None,
             gradient_deviations=None,
             extra_protocol={})

         mdt.fit_model(
             'NODDI',
             input_data,
             {output.outputdir},
             recalculate=True,
             double_precision=False,
             cl_device_ind=[0],
             use_cascaded_inits=True,
             method='Powell',
             optimizer_options={'patience': 2})