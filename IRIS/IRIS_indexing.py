import numpy as np
import os,sys,argparse
from scipy import stats
import statsmodels.stats.weightstats as smw

def index_PsiMatrix(fn, outdir, delim, splicing_event_type, out_fp):
	line_formatter =  "{id}\t{offset}\n"
	offset = 0
	col_index=10
	if splicing_event_type=='SE':#handle SE and other types of AS events
		col_index=8
	with open(fn, 'r') as fin:
		with open(out_fp, 'w') as fout:
			offset += len(fin.readline()) 
			for line in fin:
				ele = line.strip().split(delim)
				eid = ':'.join([ele[0].split('_')[0].split('.')[0]]+ele[1:col_index])
				fout.write( line_formatter.format(id=eid, offset=offset) )
				offset += len(line)
	return

def main(args):
	fin=args.splicing_matrix
	splicing_event_type=args.splicing_event_type
	cov_cutoff=args.cov_cutoff
	data_name=args.data_name
	outdir=args.outdir.rstrip('/')
	out_fp = outdir+'/'+fin.split('/')[-1]+'.idx'
	#create index in the current directory
	index_PsiMatrix(fin,outdir,'\t',splicing_event_type, out_fp)
	print '[INFO] Finished. Created matrix: '+out_fp

if __name__ == '__main__':
	main()
