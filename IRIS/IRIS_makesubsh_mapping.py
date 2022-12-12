import sys,glob,os
from . import config

def write_task_script(out_prefix, fin, label_string, starGenomeDir, gtf, task_dir):
	fq_dir=fin
	fqs=glob.glob(fq_dir+'/*')
	r1=[]
	r2=[]
	for fq in fqs:
		if fq.find('1'+label_string+'f')!=-1:
			r1.append(os.path.abspath(fq))
		elif fq.find('2'+label_string+'f')!=-1:
			r2.append(os.path.abspath(fq))
	if len(r1)!=len(r2) or len(r1)==0 or len(r2)==0:
		print 'file name can not be recognize'
		return '','' 

	out_dir=out_prefix.rstrip('/')+'.aln'
	sample_name=out_dir.split('/')[-1].split('.')[0]
	task_script_base1 = 'STARmap.{}.sh'.format(sample_name)
	task_script1 = os.path.join(task_dir, task_script_base1)
	task_script_base2 = 'Cuffquant.{}.sh'.format(sample_name)
	task_script2 = os.path.join(task_dir, task_script_base2)
	fout1=open(task_script1,'w')
	fout2=open(task_script2,'w')
	fq_path=','.join(sorted(r1)+sorted(r2))

	cmd1='IRIS process_rnaseq --starGenomeDir '+starGenomeDir+' --gtf '+gtf+' --mapping --sort -p '+out_dir+' '+fq_path
	fout1.write('#!/bin/bash\n'+cmd1+'\n')
	cmd2='IRIS process_rnaseq --starGenomeDir '+starGenomeDir+' --gtf '+gtf+' --quant -p '+out_dir+' '+fq_path
	fout2.write('#!/bin/bash\n'+cmd2+'\n')
	return task_script1,task_script2

def main(args):
	starGenomeDir=args.starGenomeDir
	gtf=args.gtf
	fastq_folder_dir=args.fastq_folder_dir.rstrip('/')
	fastq_folder_list=glob.glob(fastq_folder_dir+'/*')
	out_dir=args.outdir.rstrip('/')
	task_name=args.data_name
	label_string=args.label_string
	os.system('mkdir -p '+out_dir)
	task_dir=args.task_dir
	if not os.path.exists(task_dir):
		os.makedirs(task_dir)

	for folder in fastq_folder_list:
		print folder
		fn1,fn2=write_task_script(out_dir+'/'+folder.split('/')[-1], folder, label_string, starGenomeDir, gtf, task_dir)
		if fn1=='':
			continue


if __name__ == '__main__':
	main()
