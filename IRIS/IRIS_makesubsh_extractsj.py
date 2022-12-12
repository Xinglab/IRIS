import sys, csv, glob, os, argparse

def parseMappingLog(log_fin):
	read_length=''
	for l in open(log_fin):
		if l.find('Average input read length')!=-1:
			ls=l.strip().split('|')
			read_length=ls[1].strip().strip(' ')
	if read_length=='':
		exit('can not find log file for read length')
	return str(int(round(int(read_length)/2)))

def main(args):
	#extractSJ_path=args.extractSJ_path
	BAM_prefix=args.BAM_prefix
	bam_folder_list=args.bam_folder_list
	rl=args.rmats_used_read_length
	parserl=False
	if rl=='':
		parserl=True
		print('checking read legnth')
	else:
		print('use user specified read length')
	gtf=args.gtf
	task_name= args.task_name
	genome_fasta= args.genome_fasta
	task_dir = args.task_dir
	if not task_dir:
		task_dir = os.getcwd()

	list_name='cmdlist.extract_sjc.'+task_name
	list_name = os.path.join(task_dir, list_name)

	n=0
	fout=open(list_name,'w')
	for bam_folder in open(bam_folder_list):
		n+=1
		bam_folder=bam_folder.strip()
		if parserl:
			rl=parseMappingLog(bam_folder+'/Log.final.out')
		fout.write('IRIS extract_sjc -i '+bam_folder+'/'+BAM_prefix+'.bam -g '+gtf+' -r '+rl+' -f '+genome_fasta+' -o '+bam_folder+'/SJcount.txt \n')
	fout.close()

	sh_file_name = 'subsh.extract_sjc.{}.sh'.format(task_name)
	sh_file_name = os.path.join(task_dir, sh_file_name)
	fout_qsub=open(sh_file_name,'w')

	cmd='sbatch --array=1-{} {}'.format(str(n), sh_file_name)

	fout_qsub.write('#!/bin/bash\n#SBATCH --job-name=extract_sjc\n#SBATCH --mem=5G\n#SBATCH -t 15:00:00\n')
	fout_qsub.write('export s=`sed -n ${SLURM_ARRAY_TASK_ID}p '+list_name+'`\necho $s\n$s')
	fout_qsub.close()
	print(cmd)
if __name__ == '__main__':
	main()
