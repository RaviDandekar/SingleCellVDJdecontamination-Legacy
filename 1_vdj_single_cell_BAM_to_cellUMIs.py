#!/usr/bin/env python3.6
import sys,time,re,os
import numpy  as np
import pandas as pd
import pysam
import python_functions as pf
from multiprocessing import Pool


def extract_cellumi(x):
  sample = x.split('/')[-3]
  print(sample)
  
  # make results file
  results_dir = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/CellUMI_files/"
  sample_cellumi_fh = results_dir + sample + '.csv'
  
  print("CREATING:\t",sample_cellumi_fh)
  
  OUT = open(sample_cellumi_fh,'w')
  OUT.write("read,cell_barcode,umi_seq\n")
  
  samfile = pysam.AlignmentFile(x, "rb")
  for read in samfile.fetch():
    read_id = read.qname
    tag_array = read.tags
    cell_barcode = ""
    umi          = ""
  
    for k,v in tag_array:
      if k == "CB":
        cell_barcode = v
      elif k == "UB":
        umi = v
      else:
        continue
    entry = "%s,%s,%s\n" % (read_id,cell_barcode,umi)
    OUT.write(entry)
  
  OUT.close()
  

if __name__ == '__main__':
  THREADS = 10
  
  # Create list of BAM files
  bam_list = []
  sample_hash = pf.AutoVivification()
  vdj_dir_bcr_new     = "/data/ryan/10x/2020_vdj_fastqs/2020_vdj_results/BCR/"
  vdj_dir_bcr_old     = "/data/rdandekar/Reruns_cellranger/All_VDJ_cellranger_results_Dec2018/"
  results_dir = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/CellUMI_files/"
  if os.path.isdir(results_dir) == False:
    os.mkdir(results_dir)
    
  dir_samples_bcr_new = os.listdir(vdj_dir_bcr_new)
  for sample in dir_samples_bcr_new:
    bam_file = vdj_dir_bcr_new + sample + "/outs/all_contig.bam"
    sample_hash[sample] = bam_file
    bam_list.append(bam_file)
  
  dir_samples_bcr_old = os.listdir(vdj_dir_bcr_old)
  for sample in dir_samples_bcr_old:
    bam_file = vdj_dir_bcr_old + sample + "/outs/all_contig.bam"
    sample_hash[sample] = bam_file
    bam_list.append(bam_file)
  
  
  ## Extract cell barcode UMI info from BAM files
  with Pool(THREADS) as p:
    p.map(extract_cellumi,bam_list)

