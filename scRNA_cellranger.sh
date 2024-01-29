cd /media/ysl/
cellranger count --id=run_count_N35 \
   --fastqs=/media/ysl/scRNA_rawdata/N35 \
   --sample=N35 \
   --transcriptome=/media/ysl/yard/run_cellranger_count/refdata-cellranger-mm10-1.2.0

cellranger count --id=run_count_W35 \
   --fastqs=/media/ysl/scRNA_rawdata/W35 \
   --sample=N35 \
   --transcriptome=/media/ysl/yard/run_cellranger_count/refdata-cellranger-mm10-1.2.0
