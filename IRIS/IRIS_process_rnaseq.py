import argparse,os,sys,time,logging,datetime

def STAR_alignment(readsFilesRNA, starGenomeDir, dbLength, outDir):
	unsort=True
	gz=readsFilesRNA.split(',')[0].endswith('gz')
	rfs=readsFilesRNA.split(',')
	file_num=len(rfs)
	if file_num%2!=0:
		exit('need pairs')
	readsFiles_split=','.join(rfs[0:file_num/2])+' '+','.join(rfs[file_num/2:file_num])
	if unsort:
		if gz:
			cmd1='STAR --genomeDir '+starGenomeDir+' --twopassMode Basic --readFilesIn '+readsFiles_split+' --runThreadN 6 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang '+str(dbLength)+' --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMheaderHD @HD VN:1.4 --alignEndsType EndToEnd --readFilesCommand zcat --outFileNamePrefix '+outDir+'/'
		else:
			cmd1='STAR --genomeDir '+starGenomeDir+' --twopassMode Basic --readFilesIn '+readsFiles_split+' --runThreadN 6 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang '+str(dbLength)+' --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMheaderHD @HD VN:1.4 --alignEndsType EndToEnd --outFileNamePrefix '+outDir+'/'   
		logging.debug('[RNA-seq] Running STAR alignment: '+cmd1+'\n')
	else:   
		if gz:
			cmd1='STAR --genomeDir '+starGenomeDir+' --twopassMode Basic --readFilesIn '+readsFiles_split+' --runThreadN 6 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang '+str(dbLength)+' --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --alignEndsType EndToEnd --readFilesCommand zcat --outFileNamePrefix '+outDir+'/'
		else:
			cmd1='STAR --genomeDir '+starGenomeDir+' --twopassMode Basic --readFilesIn '+readsFiles_split+' --runThreadN 6 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang '+str(dbLength)+' --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --alignEndsType EndToEnd --outFileNamePrefix '+outDir+'/' 
		logging.debug('[RNA-seq] Running STAR alignment: '+cmd1+'\n')
	os.system(cmd1)
	
def cufflinks(gtf, sampleID):
  if os.path.exists(sampleID+'/Aligned.sortedByCoord.out.bam'):
      cmd2='cufflinks --multi-read-correct -p 16 --GTF '+gtf+' -o '+sampleID+'/cufflinks '+sampleID+'/Aligned.sortedByCoord.out.bam'
  else:
      cmd_sort='samtools sort '+sampleID+'/Aligned.out.bam -o '+sampleID+'/Aligned.sortedByCoord.out.bam'
      os.system(cmd_sort)
      logging.debug('[RNA-seq] Running command sort: '+cmd_sort+'\n')

      cmd2='cufflinks --multi-read-correct -p 16 --GTF '+gtf+' -o '+sampleID+'/cufflinks '+sampleID+'/Aligned.sortedByCoord.out.bam'
	
  logging.debug('[RNA-seq] Running command 2: '+cmd2+'\n')
  os.system(cmd2)

def samtools_sort(sampleID):
	cmd_sort='samtools sort '+sampleID+'/Aligned.out.bam -o '+sampleID+'/Aligned.sortedByCoord.out.bam'
	logging.debug('[RNA-seq] Running command sort: '+cmd_sort+'\n')
	os.system(cmd_sort)

def kallisto(readsFilesRNA, sampleID ,kallisto_index):
	index =  kallisto_index# '/u/nobackup/yxing/NOBACKUP/ywang/reference_genome/kallisto_index/Homo_sapiens.GRCh37.88.cdna.all.fa_kallisto.idx'
	cmd_kallisto = 'kallisto quant -i '+index+' -b 100 -t 6 -o '+sampleID+'/kallisto '+' '.join(readsFilesRNA.split(','))
	logging.debug('[RNA-seq] Running kallisto: '+cmd_kallisto)
	os.system(cmd_kallisto)
	cmd_kallisto_mk='echo '+str(datetime.datetime.now())+' >'+sampleID+'/kallisto/finished' 
	logging.debug('[RNA-seq] Running kallisto: '+cmd_kallisto_mk)
	os.system(cmd_kallisto_mk)

def main(args):
	sampleID=args.sampleID_outdir.rstrip('/')
	dbLength=args.db_length
	all_the_way=True
	if args.mapping or args.quant or args.sort:
		all_the_way=False

	os.system('mkdir -p '+sampleID)
	logging.basicConfig(level=logging.DEBUG,
						format='%(asctime)s %(message)s',
						filename=sampleID+'/IRIS_process_rnaseq.log'+ str(datetime.datetime.now())+'.txt' ,
						filemode='w')

	### STAR
	if all_the_way or args.mapping:
		logging.debug('[RNA-seq] # Start STAR 2pass alignment.')
		print '[RNA-seq] # Start STAR 2pass alignment.'

		if os.path.exists(sampleID+'/Aligned.sortedByCoord.out.bam')==False and os.path.exists(sampleID+'/Aligned.out.bam')==False:
			STAR_alignment(args.readsFilesRNA, args.starGenomeDir, dbLength, sampleID)
			if os.path.exists(sampleID+'/Aligned.sortedByCoord.out.bam')==False and os.path.exists(sampleID+'/Aligned.out.bam')==False:
				sys.exit('[RNA-seq] # An Error Occured. STAR Incomplete. Exit!')
		else:
			logging.debug('[RNA-seq] # Skipped STAR 2pass alignment.')
			print '[RNA-seq] # Skipped STAR 2pass alignment.'

	### samtools sortedByCoord
	if all_the_way or args.quant or args.sort:
		logging.debug('[RNA-seq] # Start samtools sortedByCoord.')
		print '[RNA-seq] # Start samtools sortedByCoord.'

		if os.path.exists(sampleID+'/Aligned.sortedByCoord.out.bam')==True:
			logging.debug('[RNA-seq] # Skipped samtools sortedByCoord. BAM is already sortedByCoord')
		else:
			if os.path.exists(sampleID+'/Aligned.out.bam'):
				samtools_sort(sampleID)
			else:
				sys.exit('[RNA-seq] # An Error Occured. Unsorted bam does not exist.')  

	### Quantification
	if all_the_way or args.quant:
		logging.debug('[RNA-seq] # Start Gene Expression Quantification.')
		print '[RNA-seq] # Start Gene Expression Quantification.'
		# if os.path.exists(sampleID+'/kallisto/abundance.h5')==True:
		# 	logging.debug('[RNA-seq] # Skipped GE Quantification.')
		# 	print '[RNA-seq] # Skipped GE Quantification.'
		# else:
		# 	kallisto(args.readsFilesRNA,sampleID, kallisto_index)
		# 	if os.path.exists(sampleID+'/kallisto/finished')==False:
		# 		logging.debug('[RNA-seq] # An Error Occured. GE Quantification Incomplete. Exit!')
		# 		sys.exit('[RNA-seq] # An Error Occured. GE Quantification Incomplete. Exit!')

		if os.path.exists(sampleID+'/cufflinks/genes.fpkm_tracking')==False:
			cufflinks(args.gtf,sampleID)
			if os.path.exists(sampleID+'/cufflinks/genes.fpkm_tracking')==False:
				sys.exit('[RNA-seq] # An Error Occured. cufflinks Incomplete. Exit!')
		else:
			logging.debug('[RNA-seq] # Skipped GE Quantification.')
			print '[RNA-seq] # Skipped GE Quantification.'


	logging.debug('[RNA-seq] # Completed.')
	print '[RNA-seq] # Completed.'

if __name__ == '__main__':
	main()