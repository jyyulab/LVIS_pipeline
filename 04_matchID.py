import pyfastx
import simplesam
import os

os.chdir('/research/projects/yu3grp/IO_JY/yu3grp/LVXSCID/patients_scATACseq/multiome_P1');
bam_file='./03_chimeric/P1_scMulti_ATAC_S1_pe.mated.filter.bam';
out_sam_file='./04_match_CB/P1_scMulti_ATAC_S1_pe.mated.filter_wCB.sam'
cellID_file='./04_match_CB/P1_scMulti_ATAC_S1_pe.mated.filter_R2.fastq';


#fa = pyfastx.Fastx('./LVX_SCID_P1_S1_L001_pe.mated.filter2.bam_readbarcode')
fa = pyfastx.Fastx(cellID_file)

barcodes = {}
for name,seq,qual,comment in fa:
	barcodes[name] = seq

barcode_tag = 'CB'

with simplesam.Reader(open(bam_file)) as in_bam:
  with simplesam.Writer(open(out_sam_file, 'w'), in_bam.header) as out_sam:
    for read in in_bam:
      #read[umi_tag] = barcodes[read.qname][0]
      read[barcode_tag] = barcodes[read.qname];
      out_sam.write(read)

