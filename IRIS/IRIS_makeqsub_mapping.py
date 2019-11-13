import sys,glob,os
from . import config
#python make.shellsubmiter.py Out_dir Fastq_dir
# two inputs: 1)output folder/sh name/prefix 2)dir of all RNAseq files (if more than 2, recognize and sort by R1 R2 for star input)
def makeSubmit(out_prefix, fin, label_string, starGenomeDir, gtf):
	fq_dir=fin
	fqs=glob.glob(fq_dir+'/*')
	#abs_path=os.path.abspath(sys.argv[2])
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
        fout1=open('submit.STARmap.'+sample_name+'.sh','w')
        fout2=open('submit.Cuffquant.'+sample_name+'.sh','w')
	fq_path=','.join(sorted(r1)+sorted(r2))
	cmd1='IRIS process_rnaseq --starGenomeDir '+starGenomeDir+' --gtf '+gtf+' --mapping --sort -p '+out_dir+' '+fq_path
	fout1.write(cmd1+'\n')
	cmd2='IRIS process_rnaseq --starGenomeDir '+starGenomeDir+' --gtf '+gtf+' --quant -p '+out_dir+' '+fq_path
	fout2.write(cmd2+'\n')
	return 'submit.STARmap.'+sample_name+'.sh','submit.Cuffquant.'+sample_name+'.sh'

def main(args):
	starGenomeDir=args.starGenomeDir
	gtf=args.gtf
	fastq_folder_dir=args.fastq_folder_dir
	fastq_folder_list=glob.glob(fastq_folder_dir+'')
	out_dir=args.out_dir.rstrip('/')
	task_name=args.data_name
	label_string=args.label_string
	os.system('mkdir -p '+out_dir)

	fout1=open('cmdlist.STAR.'+task_name,'w')
	fout2=open('cmdlist.Cufflinks.'+task_name,'w')
	i=0
	for folder in fastq_folder_list:
		print folder
		fn1,fn2=makeSubmit(out_dir+'/'+folder.split('/')[-1], folder, label_string, starGenomeDir)
		if fn1=='':
			continue
		i+=1
		fout1.write(fn1+'\n')
		fout2.write(fn2+'\n')
	fout1.close()
	fout2.close()

	fout_qsub1=open('qsub.STARmapping.'+task_name+'.sh','w')
	cmd='qsub -t 1-'+str(i)+':1 qsub.STARmapping.'+task_name+'.sh'
	fout_qsub1.write('#!/bin/bash\n#$ -N STARmapping\n#$ -S /bin/bash\n#$ -R y\n#$ -l '+config.QSUB_ALIGNMENT_CONFIG+'\n#$ -V\n#$ -cwd\n#$ -j y\n#$ -m bea\n')
	fout_qsub1.write('export s=`sed -n ${SGE_TASK_ID}p '+'cmdlist.STAR.'+task_name+'`\necho $s\nbash $s')
	fout_qsub1.close()
	print cmd

	fout_qsub2=open('qsub.Cufflinks.'+task_name+'.sh','w')
	cmd='qsub -t 1-'+str(i)+':1 qsub.Cufflinks.'+task_name+'.sh'
	fout_qsub2.write('#!/bin/bash\n#$ -N Cufflinks\n#$ -S /bin/bash\n#$ -R y\n#$ -l '+config.QSUB_EXPRESSION_CONFIG+'\n#$ -V\n#$ -cwd\n#$ -j y\n#$ -m bea\n')
	fout_qsub2.write('export s=`sed -n ${SGE_TASK_ID}p '+'cmdlist.Cufflinks.'+task_name+'`\necho $s\nnbash $s')
	fout_qsub2.close()
	print cmd

if __name__ == '__main__':
	main()