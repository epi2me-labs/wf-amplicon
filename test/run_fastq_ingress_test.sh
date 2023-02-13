#!/bin/bash
set -exo pipefail

get-test_data-from-aws () {
    # get glibc for aws-cli
    apk del libc6-compat
    apk --no-cache add binutils curl
    # Started to get breakages with glibc (2.35-r0), so fix to previous good version
    #- "GLIBC_VER=$(curl -s https://api.github.com/repos/sgerrand/alpine-pkg-glibc/releases/latest | grep tag_name | cut -d : -f 2,3 | tr -d \\\",' ' )"
    GLIBC_VER="2.34-r0"
    curl -sL https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub -o /etc/apk/keys/sgerrand.rsa.pub
    curl -sLO https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${GLIBC_VER}/glibc-${GLIBC_VER}.apk
    curl -sLO https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${GLIBC_VER}/glibc-bin-${GLIBC_VER}.apk
    apk add --no-cache glibc-${GLIBC_VER}.apk glibc-bin-${GLIBC_VER}.apk

    # get aws-cli
    curl -s "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
    unzip -q awscliv2.zip

    # get test data
    aws/dist/aws s3 cp --recursive --quiet \
        "$S3_TEST_DATA" \
        test_data_from_S3
}

fastq=$1
wf_output_dir=$2
sample_sheet=$3

# `fastq` and `wf_output_dir` are required
if ! [[ $# -eq 2 || $# -eq 3 ]]; then
    echo "Provide 2 or 3 arguments!" >&2
    exit 1
fi

# get test data from s3 if required
if [[ $fastq =~ ^s3:// ]]; then
    get-test_data-from-aws
    fastq="$PWD/test_data_from_S3/${fastq#*test_data/}"
    [[ -n $sample_sheet ]] &&
        sample_sheet="$PWD/test_data_from_S3/${sample_sheet#*test_data/}"
fi

# add CWD if paths are relative
[[ ( $fastq != /* ) ]] && fastq="$PWD/$fastq"
[[ ( $wf_output_dir != /* ) ]] && wf_output_dir="$PWD/$wf_output_dir"
[[ ( -n $sample_sheet ) && ( $sample_sheet != /* ) ]] &&
    sample_sheet="$PWD/$sample_sheet"

# add flags to parameters
fastq="--fastq $fastq"
wf_output_dir="--wf-output-dir $wf_output_dir"
[[ -n $sample_sheet ]] && sample_sheet="--sample_sheet $sample_sheet"

# get container hash from config
img_hash=$(grep 'container_sha.\?=' nextflow.config | grep -oE 'sha[0-9,a-f,A-F]+')

# run test
docker run -v "$PWD":"$PWD" \
    ontresearch/wf-template:"$img_hash" \
    python "$PWD"/test/test_fastq_ingress.py $fastq $wf_output_dir $sample_sheet
