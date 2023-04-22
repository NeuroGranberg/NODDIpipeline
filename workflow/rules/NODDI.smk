configfile:

rule all:
    input:
        bvecs: "{indir}/eddy_rotated_bvecs"
        bvals: "{indir}/bvals"
        dwi: "{indir}/eddy_corrected_data.nii.gz"
        mask: "{indir}/brain_mask.nii.gz"


rule image_for_topup:
    input:
        b1000: "{indir}/test_ep2d_diff_2mm_hcp_32_b1000_AP.nii"
    output: temp("{indir}/Diff_32_0_AP_PA.nii.gz")
    shadow:"shallow"
    log: "{indir}/snakelog"
    shell:
    r"""
    fslroi {input.b1000} Diff_32_0_AP.nii.gz 0 2 &>> {log}
    fslroi {input.b1000} Diff_32_0_APfix.nii.gz 0 1 &>> {log}
    fslroi {input.b1000} Diff_0_G3_PA0.nii.gz 0 1 &>> {log}
    fslroi {input.b1000} Diff_0_G3_PA1.nii.gz 1 1 &>> {log}
    flirt -in Diff_0_G3_PA0.nii.gz -ref Diff_32_0_APfix.nii.gz -out Diff_0_G3_PAfixed0.nii.gz &>> {log}
    flirt -in Diff_0_G3_PA1.nii.gz -ref Diff_32_0_APfix.nii.gz -out Diff_0_G3_PAfixed1.nii.gz &>> {log}
    fslmerge -t {output} Diff_32_0_AP.nii.gz Diff_0_G3_PAfixed0.nii.gz Diff_0_G3_PAfixed1.nii.gz &>> {log}
    rm Diff_32_0_AP.nii.gz Diff_32_0_APfix.nii.gz Diff_0_G3_PA0.nii.gz Diff_0_G3_PA1.nii.gz Diff_0_G3_PAfixed0.nii.gz Diff_0_G3_PAfixed1.nii.gz &>> {log}
    """

#side effects
rule topup:
    input:
        b0appa: "{indir}/Diff_32_0_AP_PA.nii.gz"
        acqparams: "{indir}/acqparams.txt"
    output:
        image: "{indir}/my_unwarped_images.nii.gz"
        topupim: "{indir}/my_topup_results_fieldcoeff.nii.gz"
        topuptxt: "{indir}/my_topup_results_movpar.txt"
    shadow:"shallow"
    log: "{indir}/snakelog"
    shell:
    r"""
    topup --imain={input.b0appa} --datain={input.acqparams} --config=b02b0.cnf --out=my_topup_results --iout=my_unwarped_images
    mv my_unwarped_images.nii.gz {output.image}
    mv my_topup_results_fieldcoeff.nii.gz {output.topupim}
    mv my_topup_results_movpar.txt {output.topuptxt}
    """

#side effects
rule brainmask:
    input:
        im: "{indir}/my_unwarped_images.nii.gz"
    output:
        mask: "{indir}/brain_mask.nii.gz"
    params:
        bet_thres:"0.3"
    shadow:"shallow"
    log: "{indir}/snakelog"
    shell:
    r"""
    fslmaths {input.im} -Tmean MeanB0.nii.gz
    bet MeanB0.nii.gz brain -m -R -f {params.bet_thres}
    rm MeanB0.nii.gz brain.nii.gz
    mv brain_mask.nii.gz {output.mask}
    """

rule merge_shells:
    input:
        b1000:"{indir}/test_ep2d_diff_2mm_hcp_32_b1000_AP.nii"
        b1000bvecs: "{indir}/test_ep2d_diff_2mm_hcp_32_b1000_AP.voxel_space.bvecs"
        b1000bvals: "{indir}/test_ep2d_diff_2mm_hcp_32_b1000_AP.bvals"
        b3000:"{indir}/test_ep2d_diff_2mm_hcp_73_b3000_AP.nii"
        b3000bvecs: "{indir}/test_ep2d_diff_2mm_hcp_73_b3000_AP.voxel_space.bvecs"
        b3000bvals: "{indir}/test_ep2d_diff_2mm_hcp_73_b3000_AP.bvals"
    output:
        dwi: temp("{indir}/DWI.nii.gz")
        bvecs: temp("{indir}/DWI.bvecs")
        bvals: temp("{indir}/DWI.bvals")
        index: temp("{indir}/index.txt")
    log: "{indir}/snakelog"
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
#side effects
rule eddy:
    input:
        imain: "{indir}/DWI.nii.gz"
        mask: "{indir}/brain_mask.nii.gz"
        acqp: "{indir}/acqparams.txt"
        index: "{indir}/index.txt"
        bvecs: "{indir}/DWI.bvecs"
        bvals: "{indir}/DWI.bvals"
        topupim: "{indir}/my_topup_results_fieldcoeff.nii.gz"
        topuptxt: "{indir}/my_topup_results_movpar.txt"
    output: 
        image: "{indir}/eddy_corrected_data.nii.gz"
        bvecs: "{indir}/eddy_rotated_bvecs"
    shadow:"shallow"
    log: "{indir}/snakelog"
    shell:
    r"""
    eddy_cuda10.2 --imain={input.imain} --mask={input.mask} --acqp={input.acqp} --index={input.index} --bvecs={input.bvecs} --bvals={input.bvals} --topup=my_topup_results --out=eddy_corrected_data
    mv eddy_corrected_data.nii.gz {output.image}
    mv eddy_corrected_data.eddy_rotated_bvecs {output.bvecs}
    """

AWK_CMD =r"""{printf "%i ",$0}"""
rule bvals_for_mdt:
    input:
        bvals: "{indir}/DWI.bvals"
        script: "$scripts/transpose.sh" 
    output:
        bvals: "{indir}/bvals"     
    log: "{indir}/snakelog"
    shadow:"shallow"
    shell:
    r"""
    awk {AWK_CMD:q} {input.bvals}> bvals1
    {input.script} bvals1 > bvals2
    {input.script} bvals2 > {output.bvals}
    rm bvals1 bvals2
    """

# rule mdt_fit:
#     input:
#         bvecs: "{indir}/eddy_rotated_bvecs"
#         bvals: "{indir}/bvals"
#         dwi: "{indir}/eddy_corrected_data.nii.gz"
#         mask: "{indir}/brain_mask.nii.gz"
#     # output:
    #     protocol: temp("/protocol.prtcl")
    #     outputdir: directory("/outputdir")
    # log: 'run.log'
    # run:
    #     import mdt
    #     protocol = mdt.create_protocol(
    #     bvecs= {input.bvecs}, 
    #     bvals= {input.bvals},
    #     out_file= {output.protocol})

    #     input_data = mdt.load_input_data(
    #         {input.dwi},
    #         {output.protocol},
    #         {input.mask},
    #         noise_std=None,
    #         gradient_deviations=None,
    #         extra_protocol={})

    #     mdt.fit_model(
    #         'NODDI',
    #         input_data,
    #         {output.outputdir},
    #         recalculate=True,
    #         double_precision=False,
    #         cl_device_ind=[0],
    #         use_cascaded_inits=True,
    #         method='Powell',
    #         optimizer_options={'patience': 2})