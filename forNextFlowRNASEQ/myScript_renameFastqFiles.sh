
# run this file from /data2/hanlab/richardvan/tubby_pipeline

# making a variable for our main pipeline path
TUBBY_PIPELINE_RAWDATA=/data2/han_lab/richardvan/tubby_pipeline/data/raw_fastq
TUBBY_PIPELINE_RENAMED=/data2/han_lab/richardvan/tubby_pipeline/data/renamed_fastq

set -e  #exits the script when an error occurs for safety

### use associative array to handle this
# (source) https://stackoverflow.com/questions/28725333/looping-over-pairs-of-values-in-bash
# (tutorial code)
        # animals=(dog cat mouse)
        # declare -A size=(
        #   [dog]=big
        #   [cat]=medium
        #   [mouse]=small
        # )
        # declare -A sound=(
        #   [dog]=barks
        #   [cat]=purrs
        #   [mouse]=cheeps
        # )
        # for animal in "${animals[@]}"; do
        #   echo "$animal is ${size[$animal]} and it ${sound[$animal]}"
        # done

samples=(sample01 sample02 sample03 sample04 sample05 sample06 \
         sample07 sample08 sample09 sample10 sample11 sample12 \
         sample13 sample14 sample15 sample16 sample17 sample18 \
         sample19 sample20 sample21 sample22 sample23 sample24 \
        )
declare -A raw_fastq=(
    [sample01]="1Tub24-R1_S4"
    [sample02]="2Tub24-R2_S6"
    [sample03]="3Tub24-R3_S3"
    [sample04]="13-38Nora_S12"
    [sample05]="14-39Nora_S11"
    [sample06]="15-40Nora_S15"
    [sample07]="18-46Nora_S3"
    [sample08]="8WT24-R2_S2"
    [sample09]="9WT24-R3_S1"
    [sample10]="10WT3m-R1_S16"
    [sample11]="11WT3m-R2_S14"
    [sample12]="12WT3m-R3_S8"
    [sample13]="13Tub24-H1_S17"
    [sample14]="14Nora_S15"
    [sample15]="15Nora_S9"
    [sample16]="17-42Nora_S10"
    [sample17]="17Nora_S17"
    [sample18]="18Tub3m-H5_S2"
    [sample19]="19WT24-H4_S11"
    [sample20]="16-41Nora_S13"
    [sample21]="21Nora_S4"
    [sample22]="22WT3m-H1_S10"
    [sample23]="23Nora_S2"
    [sample24]="24Nora_S13"
    )
declare -A renamed_fastq=(
    [sample01]="tb_24d_ret1"
    [sample02]="tb_24d_ret2"
    [sample03]="tb_24d_ret3"
    [sample04]="tb_3mo_ret1"
    [sample05]="tb_3mo_ret2"
    [sample06]="tb_3mo_ret3"
    [sample07]="wt_24d_ret1"
    [sample08]="wt_24d_ret2"
    [sample09]="wt_24d_ret3"
    [sample10]="wt_3mo_ret1"
    [sample11]="wt_3mo_ret2"
    [sample12]="wt_3mo_ret3"
    [sample13]="tb_24d_hyp1"
    [sample14]="tb_24d_hyp2"
    [sample15]="tb_24d_hyp3"
    [sample16]="tb_3mo_hyp1"
    [sample17]="tb_3mo_hyp2"
    [sample18]="tb_3mo_hyp3"
    [sample19]="wt_24d_hyp1"
    [sample20]="wt_24d_hyp2"
    [sample21]="wt_24d_hyp3"
    [sample22]="wt_3mo_hyp1"
    [sample23]="wt_3mo_hyp2"
    [sample24]="wt_3mo_hyp3"
    )

#### copy files from data/raw_fastq to data/renamed_fastq
for sample in ${samples[@]}
do
    # echo $sample raw is ${raw_fastq[$sample]} 
    # echo $sample renamed is ${renamed_fastq[$sample]} 
    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L001_R1_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R1.fastq.gz 

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L002_R1_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L2_R1.fastq.gz 

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L003_R1_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L3_R1.fastq.gz 

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L004_R1_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L4_R1.fastq.gz  

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L001_R2_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L1_R2.fastq.gz 

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L002_R2_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L2_R2.fastq.gz 

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L003_R2_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L3_R2.fastq.gz 

    cp -av ${TUBBY_PIPELINE_RAWDATA}/${raw_fastq[$sample]}_L004_R2_001.fastq.gz \
           ${TUBBY_PIPELINE_RENAMED}/${renamed_fastq[$sample]}_L4_R2.fastq.gz 
done


# for filename in ${TUBBY_PIPELINE_RAWDATA}/*_L001_R1_001.fastq.gz
# do
#     base=$(basename $filename _L1_R1.fastq.gz)
#     echo $base
#     echo $filename
# done



