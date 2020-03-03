#!/usr/bin/env python3.6
import sys,time,re,os
import numpy  as np
import pandas as pd
import python_functions as pf

if __name__ == '__main__':
  cell_hash = pf.AutoVivification()
  
  vdj_cellumi_dir         = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/CellUMI_files/"
  results_contam_cells_fh = "/data/rdandekar/rprojects/sc_analysis_2019/contamination_analysis/vdj_contamination_analysis/decontam_Jan2020/results_all_BCR/round1/FLAGGED_contam_cellumi.VDJ.txt"
  
  for root,dirs,filenames in os.walk(vdj_cellumi_dir):
	  for file in filenames:
	    sample = file.split('.')[0]
	    fh     = os.path.join(root,file)
	    print(">> PROCESSING:\t",sample, "\t", fh)
	    start = time.time()
	    
	    sample_df = pd.read_csv(fh)
	    sample_df['cellumi'] = sample_df['cell_barcode'].map(str) + '_' + sample_df['umi_seq'].map(str)
	    cellumi_array = np.array(sample_df['cellumi'])
	    cellumi_array = np.unique(cellumi_array,return_counts=True)
	    for i in range(len(cellumi_array[0])):
	      cellumi = cellumi_array[0][i]
	      count   = cellumi_array[1][i]
	      cell_hash[cellumi][sample] += count
	    end = time.time()
	    print("Time (min) =  ", (end - start)/60)
      
  print("\n>> Aggregate all contamination and write to output file")
  OUT = open(results_contam_cells_fh,'w')
  for cellumi in cell_hash:
    count = 0
    sample_str = ""
    for sample in cell_hash[cellumi]:
      count += 1
      read_count = cell_hash[cellumi][sample]
      sample_str += sample +  '=' + str(read_count) + ':'
    sample_str = sample_str[:-1]
    entry = "%s\t%s\n" % (cellumi,sample_str)
    if count > 1:
      OUT.write(entry)
    
  OUT.close()
    
      
	    
	    


