import sys, os

def loadFinlist(fin_list_input):
	fin_list=[]
	for l in open(fin_list_input):
		fin_list.append(l.strip())
	return fin_list
def readSJfile_STAR(fin, SJ_dict):
	for l in open(fin):
		ls=l.strip().split('\t')
		sj=':'.join(ls[0:3])
		SJ_dict[sj]=ls[5] #add strand info
	return SJ_dict

def readSJfile(fin, SJ_dict):
	for l in open(fin):
		ls=l.strip().split('\t')
		sj=ls[0]
		SJ_dict[sj]='' #add strand info
	return SJ_dict

def index_SJMatrix(fn, outdir, delim):
	out_fp = outdir+'/'+fn.split('/')[-1]+'.idx'
	line_formatter =  "{id}\t{offset}\n"
	offset = 0	
	with open(outdir+'/'+fn, 'r') as fin:
		with open(out_fp, 'w') as fout:
			offset += len(fin.readline()) 
			for line in fin:
				ele = line.strip().split(delim)
				eid = ele[0]
				fout.write( line_formatter.format(id=eid, offset=offset) )
				offset += len(line)
	return

def main(args):
	
	fname_pos=args.sample_name_field#2
	fin_list_input=args.file_list_input
	fin_list=loadFinlist(fin_list_input)
	data_name=args.data_name
	db_dir=args.iris_db_path.rstrip('/')
	os.system('mkdir -p '+db_dir+'/'+data_name+' '+db_dir+'/'+data_name+'/sjc_matrix')


	if os.path.exists(db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.txt'):
		print '[INFO] Output '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.txt'+' exists. Only perform indexing.'
		index_SJMatrix('SJ_count.'+data_name+'.txt', db_dir+'/'+data_name+'/sjc_matrix/', '\t')
		exit('[INFO] Index finished.')
		
	fname_dict={}
	for fn in fin_list:
		name=fn.split('/')[-fname_pos].split('.aln')[0]
		if name in fname_dict:
			print name
			exit('dup name'+fn+' '+name)
		fname_dict[name]=''
	fname_list=fname_dict.keys()
	print '[INFO] Done checking file names.'

	SJ_dict={}
	for i,fin in enumerate(fin_list):
		if i%100==0:
			print i
		SJ_dict=readSJfile(fin, SJ_dict)

	SJ_list=sorted(SJ_dict.keys())

	fout_SJ=open(db_dir+'/'+data_name+'/sjc_matrix/SJ_coordiate.'+data_name+'.txt','w')
	for s in SJ_list:
		fout_SJ.write(s+'\t'+SJ_dict[s]+'\n')
	fout_SJ.close()

	print '[INFO] Done summarizing SJ coordinates.' 

	batch=100000
	for b in range(0, len(SJ_list), batch):
		print b
		batch_SJ_list=SJ_list[b:min(b+batch,len(SJ_list))]
		batch_SJ_dict=dict.fromkeys(batch_SJ_list,0) # can't use this to store values. will only store the last input
		batch_SJ_count={}
		for fin in fin_list:
			fname=fin.split('/')[-fname_pos].split('.aln')[0]
			for l in open(fin):
				ls=l.strip().split('\t')
				# sj=':'.join(ls[:3])# STAR
				# count=ls[6]# STAR
				sj=ls[0]
				count=ls[1]
				if sj in batch_SJ_dict:
					if sj not in batch_SJ_count:
						batch_SJ_count[sj]={}
					batch_SJ_count[sj][fname]=count
					continue
		fout_intermediate=open(db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.batch_'+str(b)+'.txt','w')
		for k in sorted(batch_SJ_count.keys()):
			sj_line=[k]
			for sample in fname_list:
				if sample in batch_SJ_count[k]:
					sj_line.append(batch_SJ_count[k][sample])
				else:
					sj_line.append('0') ##It's okay to change to 0 later
			fout_intermediate.write('\t'.join(sj_line)+'\n')
		fout_intermediate.close()
	fout_head=open(db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.header.txt','w')
	fout_head.write('\t'.join(['SJ']+fname_list)+'\n')
	fout_head.close()
	cmd='cat '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.batch_*.txt > '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.txt_tmp'
	cmd_merge='cat '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.header.txt '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.txt_tmp > '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.txt'
	cmd_rm='rm '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.batch_*.txt'
	cmd_rm2='rm '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.header.txt'
	cmd_rm3='rm '+db_dir+'/'+data_name+'/sjc_matrix/SJ_count.'+data_name+'.txt_tmp'
	print cmd
	os.system(cmd)
	os.system(cmd_merge)
	print cmd_rm
	os.system(cmd_rm)
	os.system(cmd_rm2)
	os.system(cmd_rm3)

	print '[INFO] Matrix finished.'
	index_SJMatrix('SJ_count.'+data_name+'.txt', db_dir+'/'+data_name+'/sjc_matrix/', '\t')
	print '[INFO] Index finished.'

if __name__ == '__main__':
	main()
