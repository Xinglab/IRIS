import sys, csv, glob, os
from . import config
def main(args):

	bam_dir= args.bam_dir
	file_list=glob.glob(bam_dir.rstrip('/')+'/*/*.bam')
	length=str(args.read_length)
	gtf=args.gtf
	task_name=args.data_name
	rMATS_path=args.rMATS_path#/u/home/s/shiehshi/rMATS-2017-3-15/rmats.py
	list_name='cmdlist.rMATS_prep_'+task_name

	i=0
	fout=open(list_name,'w')
	for fin in file_list:
		fin_name='/'.join(fin.split('/')[:-1]).rstrip('/')
		i+=1
		fout_local=open(fin_name+'/bam_list.txt','w')
		fout_local.write(fin)
		fout_local.close()
		fout.write('python '+args.rMATS_path+' --b1 '+fin_name+'/bam_list.txt --od '+fin_name+' --tmp '+fin_name+'.tmp --anchorLength 1 --readLength '+length+' --gtf '+gtf+' -t paired --task prep --nthread 8 --statoff\n')
	fout.close()

	fout_qsub=open('qsub.rMATSturboPrep.'+task_name+'.sh','w')
	cmd='qsub -t 1-'+str(i)+':1 qsub.rMATSturboPrep.'+task_name+'.sh'
	fout_qsub.write('#!/bin/bash\n#$ -N rmats_prep\n#$ -S /bin/bash\n#$ -R y\n#$ -l '+config.QSUB_RMATS_PREP_CONFIG+'\n#$ -V\n#$ -cwd\n#$ -j y\n#$ -m bea\n')
	fout_qsub.write('export s=`sed -n ${SGE_TASK_ID}p '+list_name+'`\necho $s\n$s')
	fout_qsub.close()
	print cmd
	
if __name__ == '__main__':
	main()
