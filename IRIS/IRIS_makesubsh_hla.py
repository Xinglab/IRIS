import sys,glob,os
from . import config


def write_task_script(out_prefix, fin, label_string, task_dir):
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
		print '[Error] File name can not be recognized'
		return

	out_dir=out_prefix.rstrip('/')
	os.system('mkdir -p '+out_dir)
	sample_name=out_dir.split('/')[-1].split('.')[0]	
	task_script_base = 'seq2hla.{}.sh'.format(sample_name)
	task_script = os.path.join(task_dir, task_script_base)
	fout=open(task_script,'w')
	suffix=r1[0].split(label_string)[1]
	path='/'.join(r1[0].split('/')[:-1])
	fq_path1=' '.join(sorted(r1))
	cmd_fq1='cat '+fq_path1+' > '+path+'/fq1.'+suffix
	fout.write('#!/bin/bash\n'+cmd_fq1+'\n')

	fq_path2=' '.join(sorted(r2))
	cmd_fq2='cat '+fq_path2+' > '+path+'/fq2.'+suffix
	fout.write(cmd_fq2+'\n')

	cmd1='seq2HLA -1 '+path+'/fq1.'+suffix+' -2 '+path+'/fq2.'+suffix+' -r '+out_dir+'/'+sample_name+' > '+out_dir+'/seq2hla.log 2>&1'
	fout.write(cmd1+'\n')
	fout.write('rm '+path+'/fq1.'+suffix+'\nrm '+path+'/fq2.'+suffix+'\n')
	fout.close()


def main(args):
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
		write_task_script(out_dir+'/'+folder.split('/')[-1], folder, label_string, task_dir)


if __name__ == '__main__':
	main()
