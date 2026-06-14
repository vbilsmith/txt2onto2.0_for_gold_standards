#!/bin/bash
#SBATCH --job-name=download_job
#SBATCH --nodes=1
#SBATCH --time=08:00:00


# Download a file using wget
wget -O /work/pi_hrando_smith_edu/gold_standards/BioSentVec_PubMed_MIMICIII-bigram_d700.bin "https://ftp.ncbi.nlm.nih.gov/pub/lu/Suppl/BioSentVec/BioSentVec_PubMed_MIMICIII-bigram_d700.bin"

echo "Download complete!"
