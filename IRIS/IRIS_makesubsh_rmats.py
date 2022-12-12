import sys, csv, glob, os
from . import config

def writeShell(rMATS_path, fin_name, folder_name, bam_dir, read_length_argument, gtf, novelSS, task_name, task_dir):
	fout_local=open(folder_name+'/bam_list.txt','w')
	fout_local.write(fin_name)
	fout_local.close()

	sample_name=folder_name.split('/')[-1].split('.')[0]
	task_script_base = 'rMATS_prep.{}.sh'.format(sample_name)
	task_script = os.path.join(task_dir, task_script_base)
	fout=open(task_script,'w')
	fout.write('#!/bin/bash\n')
	novelSS_str=''
	if novelSS:
		novelSS_str='--novelSS '
	# TODO the '|| true' at the end of this command ignores a
	# failure return code from the python command.
	# rMATS produces the desired output file despite the error return.
	# A future version of rMATS may fix this behavior.
	fout.write('python {} --b1 {}/bam_list.txt --od {} --tmp {}/{}.RL{}/{}.tmp --anchorLength 1 --readLength {} --gtf {} -t paired --task prep --nthread 8 --statoff {}|| true\n'.format(rMATS_path, folder_name, folder_name, bam_dir, task_name, read_length_argument, sample_name, read_length_argument, gtf, novelSS_str))
	fout.close()

def organizeReadLength(rMATS_path, file_list_mapping, gtf, novelSS, bam_prefix, task_name, task_dir):
	rl_dict={}
	folder_names={}
	for fin_name in file_list_mapping:
		for l in open(fin_name):
			if l.find('Average input read length |')!=-1:
				map_rl=int(round(float(l.split('Average input read length |')[-1].strip())/2,0))
				rl_dict['/'.join(fin_name.split('/')[:-1])]=map_rl
				folder_names[map_rl]=''
				break
	bam_dir='/'.join(file_list_mapping[0].split('/')[:-2])
	for folder_name in folder_names:
		os.system('mkdir -p '+bam_dir+'/'+task_name+'.RL'+str(folder_name))
	for folder_name in rl_dict:
		writeShell(rMATS_path, folder_name+'/'+bam_prefix+'.bam', folder_name, bam_dir, str(rl_dict[folder_name]), gtf, novelSS, task_name, task_dir)

def main(args):
	gtf=args.gtf
	task_name=args.data_name
	rMATS_path=args.rMATS_path
	bam_dir= args.bam_dir.rstrip('/')
	bam_prefix=args.bam_prefix
	novelSS=args.novelSS
	task_dir=args.task_dir
	if not os.path.exists(task_dir):
		os.makedirs(task_dir)
	if args.read_length:
		read_length= int(args.read_length)
	
	print 'preparing rMATS-turbo prep directories'
	if args.read_length==False:
		mapping_log_file_list=glob.glob(bam_dir+'/*/Log.final.out')
		organizeReadLength(rMATS_path, mapping_log_file_list, gtf, novelSS, bam_prefix, task_name, task_dir) #relocated based on the read length
	else:
		mapping_bam_list=glob.glob(bam_dir+'/*/'+bam_prefix+'.bam')
		os.system('mkdir -p '+bam_dir+'/'+task_name+'.RL'+str(read_length))
		for fin_name in mapping_bam_list:
			folder_name= '/'.join(fin_name.split('/')[:-1])
			writeShell(rMATS_path, fin_name, folder_name, bam_dir, str(read_length), gtf, novelSS, task_name, task_dir)


if __name__ == '__main__':
	main()
