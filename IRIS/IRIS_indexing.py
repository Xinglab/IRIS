import numpy as np
import os,sys,argparse
from scipy import stats
import statsmodels.stats.weightstats as smw

def index_PsiMatrix(fn,outdir,delim):
	out_fp = outdir+'/'+fn.split('/')[-1]+'.idx'
	line_formatter =  "{id}\t{offset}\n"
	offset = 0
	with open(fn, 'r') as fin:
		with open(out_fp, 'w') as fout:
			offset += len(fin.readline()) 
			for line in fin:
				ele = line.strip().split(delim)
				eid = ':'.join([ele[0].split('_')[0].split('.')[0]]+ele[1:8])
				fout.write( line_formatter.format(id=eid, offset=offset) )
				offset += len(line)
	return

def main(args):
	fin=args.splicing_matrix
	data_name=args.data_name
	db_dir=args.db_dir.rstrip('/')

	#prepare files/folders in IRIS db directory
	os.system('mkdir -p '+db_dir+'/'+data_name+' '+db_dir+'/'+data_name+'/splicing_matrix')
	new_dir_fin=db_dir+'/'+data_name+'/splicing_matrix/splicing_matrix.SE.cov10.'+data_name+'.txt'
	os.system('mv '+fin+' '+new_dir_fin)
	#create index in IRIS db directory
	index_PsiMatrix(new_dir_fin,db_dir+'/'+data_name+'/splicing_matrix','\t')
	print '[INFO] Finished. Created matrix: '+new_dir_fin

if __name__ == '__main__':
	main()
