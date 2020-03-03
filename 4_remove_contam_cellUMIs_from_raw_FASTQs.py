#!/usr/bin/env python3.6
import sys,time,re,os
import numpy  as np
import subprocess
import python_functions as pf
import gzip
from multiprocessing import Pool

def remove_contam_reads_from_vdj_fastq(fh_tup):
  sample,fh,fastq_dir = fh_tup
  master_results_dir  = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/All_cleaned_vdj_fastqs/"
  sample_results_dir  = os.path.join(master_results_dir,sample) + '/'
  removed_reads_dir   = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/removed_reads/"
  decontam_script     = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/scripts_decontam2020/fastq_cleaner.pl"
  
  if os.path.isdir(master_results_dir) == False:
    os.mkdir(master_results_dir)
  if os.path.isdir(sample_results_dir) == False:
    os.mkdir(sample_results_dir)
  if os.path.isdir(removed_reads_dir) == False:
    os.mkdir(removed_reads_dir)
  
  print(">>> Processing:\t",sample)
  
  # b. remove contam read IDs from sample FASTQS
  for root,dirs,filenames in os.walk(fastq_dir):
    for file in filenames:
      # Check if file already exists
      output_file_gz = os.path.join(sample_results_dir,file)
      if os.path.isfile(output_file_gz):
        continue
      
      fastq_fh = os.path.join(root,file)

      ## INPUT of perl script: fastq_fh,fh,sample_results_dir
      print("\tDecontaminating =>\t", file)
      subprocess.call([decontam_script,fastq_fh,fh,sample_results_dir,removed_reads_dir])
      

  
  

if __name__ == '__main__':
  # INPUT
  THREADS           = 17
  contam_reads_dir  = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/reads_to_remove/"
  new_cellranger_output = "/data/ryan/10x/2020_vdj_fastqs/2020_vdj_results/BCR/"
  old_cellranger_output = "/data/rdandekar/Reruns_cellranger/All_VDJ_cellranger_results_Dec2018/"
  
  # Contaminated FASTQs
  sc_vdj_fastq_dir1     = "/data/ryan/10x/dec_2018_igseq/IG1/"
  sc_vdj_fastq_dir2     = "/data/ryan/10x/dec_2018_igseq/IG2/"
  sc_vdj_fastq_dir3     = "/data/ryan/10x/dec_2018_igseq/IG2_Cat/"
  # New
  sc_vdj_fastq_dir4 = "/data/ryan/10x/2020_vdj_fastqs/RSMW02/"
  
  # output of samples that are already cleaan
  samples_already_cleaned = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_new_samples_only/round1/samples_already_clean_pre_decontam.txt"

  # 1. Create dictionary of all samples needed for the analysis
  target_samples_hash = pf.AutoVivification()
  target_samples_bcr = os.listdir(bcr_cellranger_output)
  target_samples_tcr = os.listdir(tcr_cellranger_output)
  
  for sample in target_samples_bcr:
    target_samples_hash[sample] += 1
  
  for sample in target_samples_tcr:
    target_samples_hash[sample] += 1
  
  # 2. Create dictionary of sample -> cellranger results directory
  sample_fastq_dir_hash = pf.AutoVivification()
  
  dir1_samples = os.listdir(sc_vdj_fastq_dir1)

  for sample in dir1_samples:
    if sample in target_samples_hash.keys():
      sample_dir = sc_vdj_fastq_dir1 + sample + "/"
      sample_fastq_dir_hash[sample] = sample_dir
  
  
  # 3. Create list of files with contaminated reads and their corresponding FASTQ directories
  sample_file_list = os.listdir(contam_reads_dir)
  fh_list = [contam_reads_dir + filename for filename in sample_file_list]
  fh_list_tuple = []
  for fh in fh_list:
    filename  = fh.split('/')[-1]
    sample    = filename.split('.')[0]
    fastq_dir = sample_fastq_dir_hash[sample]
    fh_list_tuple.append((sample,fh,fastq_dir))
    ### DEBUG: only run one sample
    #break
    ###
  
  
  ### SANITY CHECK ###
  # make sure you have accounted for every sample correctly
  clean_samples = []
  final_sample_list = [s for s,r,f in fh_list_tuple]
  for sample in target_samples_hash:
    if sample not in final_sample_list:
      clean_samples.append(sample)
  
  OUT = open(samples_already_cleaned,'w')
  for sample in clean_samples:
    OUT.write("%s\n" % (sample))
    
  OUT.close()
      
  ###

  
  # 4. remove contaminated reads from raw FASTQs
  
  # Single Thread
  #for fh_tup in fh_list_tuple:
  #  remove_contam_reads_from_vdj_fastq(fh_tup)
  
  # Multi-threading
  with Pool(THREADS) as p:
    p.map(remove_contam_reads_from_vdj_fastq,fh_list_tuple)
  
  
  print("\n>>> DONE! <<<\n")
    


