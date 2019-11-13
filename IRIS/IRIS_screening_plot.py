import numpy as np
import os,sys,glob,argparse
from scipy import stats
import statsmodels.stats.weightstats as smw
import matplotlib
matplotlib.use('agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from . import config
#If max and min is less than 5% filtered

#OUTPUT abs PSI for query group
def read_PsiMatrix(fn, delim,srr_list=None):
	with open(fn, 'r') as f:
		ele = f.readline().strip().split(delim)
		header = {ele[i].split('.')[0]:i for i in range(len(ele)) }
		if srr_list is None:
			srr_list = [ x.split('.')[0] for x in ele[8::] ]
		for line in f:
			ele = line.strip().split(delim)
			eid = ':'.join([ ele[header['chr']], ele[header['strand']], 
				ele[header['upstreamEE']], ele[header['exonStart_0base']], ele[header['exonEnd']], ele[header['downstreamES']] ])
			psi_list = np.empty(len(srr_list))
			psi_list[:] = np.nan
			for i in range(len(srr_list)):
				srr = srr_list[i]
				psi_list[i] = ele[ header[srr] ]
			yield ( eid, srr_list, psi_list )


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

def read_PsiMatrix_index(fn,outdir):
	index = {}
	if not os.path.isfile(outdir+'/'+fn.split('/')[-1]+'.idx'):
		raise Warning("indexing %s.."%fn)
		index_PsiMatrix(fn,outdir)
	with open(outdir+'/'+fn.split('/')[-1]+'.idx', 'r') as fidx:
		for line in fidx:
			ele = line.strip().split()
			index[ele[0]] = int(ele[1])
	return index


def fetch_PsiMatrix(eid, fn,outdir,delim,index=None):
	if index is None:
		index = read_PsiMatrix_index(fn,outdir)
	if eid not in index:
		raise Exception('Exon "%s" not in index'%eid)
	with open(fn, 'r') as f:
		ele = f.readline().strip().split(delim)
		header = np.asarray([ x.split('.')[0] if x.startswith('SRR') else x for x in ele ])
		f.seek(index[eid], 0)
		data = np.asarray(f.readline().strip().split(delim))
	return (header, data)

def fileLength(fin):
	i=sum(1 for l in open(fin))
	return i

def indiviualPlot(psi_df, event_id, ref_list):
	plt.figure(figsize=(15,9))
	sns_plot=sns.violinplot(data=psi_df, inner="box",cut=0)
	#sns_plot.set_yticks(np.arange(0,1,0.2))	
	sns_plot.set(ylim=(0, 1))
	sns.despine(offset=10, trim=True)
	sns_plot.set_ylabel('Percent-Spliced-In',fontweight='bold',fontsize=14)
	sns_plot.set_xticklabels(ref_list,rotation=20,ha='right',fontsize=12)
	sns_plot.set_title(event_id,fontsize=15,fontweight='bold')
	sns_plot.figure.savefig(event_id+".pdf")

def loadParametersRow(filter_para, ref_list):
	if len(filter_para.split(' '))==6:
		filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_list =filter_para.split(' ')[1:]
		filter_cutoff_pval=float(filter_cutoff_pval)
		filter_cutoff_dpsi=float(filter_cutoff_dpsi)
		filter_cutoff_foc=float(filter_cutoff_foc)
		filter_group_cutoff=int(filter_group_cutoff)
		filter_list=filter_list.split(',')		
		ref_list+=filter_list
	else:
		filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_list =['','','','',[]]
	return filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_list, ref_list


## screening
def main(args):
	index={}
	fin_list={}
	para_fin=args.parameter_fin
	fin_plot_query=args.event_list #to the txt file
	step=int(args.step)
	has_header=args.header

	out_prefix,db_dir,filter1_para,filter2_para,filter3_para,test_mode,use_ratio,blacklist_path,mappability_path,ref_genome=[l.strip() for l in open(para_fin)]
	ref_list=[out_prefix]
	filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc, filter1_group_cutoff, filter1, ref_list =loadParametersRow(filter1_para, ref_list)
	filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter2_group_cutoff, filter2, ref_list =loadParametersRow(filter2_para, ref_list)
	filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, filter3_group_cutoff, filter3, ref_list =loadParametersRow(filter3_para, ref_list)
	if filter1==[] and filter2==[] and filter3==[] and test_mode!='summary':
		exit("no filtering required in para file. exit!")
	group_test=False if test_mode!='group' else True
	individual_test=False if test_mode!='personalized' else True
	summary_file=False if test_mode!='summary' else True
	if [group_test,individual_test,summary_file]==[False,False,False]:
		exit('Need to choose one mode.exit!')
	single_plot=False
	group_plot=True
	if group_test==individual_test:
		exit('can only choose one mode')

	db_dir=db_dir
	ref_list=list([out_prefix]+filter1+filter2+filter3)

	for group_name in ref_list:##Load IRIS reference panels
		fin_list[group_name]=db_dir+'/'+group_name+'/splicing_matrix/splicing_matrix.SE.cov10.'+group_name+'.txt'

	for group in fin_list.keys():
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_PsiMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))
	skipped=0
	filered=0

	has={}
	for group in ref_list:
		if group!=out_prefix:
			has[group]=True
	
	if group_plot:
		from matplotlib.backends.backend_pdf import PdfPages

	fin_plot_query_len=fileLength(fin_plot_query)
	start_point=0
	if has_header:
		start_point=1
	for r_start in xrange(start_point,fin_plot_query_len,step):
		if group_plot:
			pdf = PdfPages(fin_plot_query+'.'+str(r_start)+'.pdf')
			fig = plt.figure(figsize=(8.8,12))
			fig.suptitle('', fontsize=14,fontweight='bold')
			sns.set(style="white", color_codes=True)
		for j,k in enumerate(open(fin_plot_query)):
			ks=k.strip().split('\t')
			event_id=ks[0]
			psi={}
			if j < r_start:
				continue
			if j > r_start+step-1:
				break
			psi[out_prefix]=map(float,fetch_PsiMatrix(event_id,fin_list[out_prefix],'.','\t',index[out_prefix])[1][8:])
			
			for group in ref_list:
				try:
					psi[group]=map(float,fetch_PsiMatrix(event_id,fin_list[group],'.','\t',index[group])[1][8:])
				except:
					has[group]=False
					psi[group]=[]
			cat_psi=[]
			for i in psi:
				cat_psi+=psi[i]
			if abs(max(cat_psi)-min(cat_psi))<0.05:#if change less than 5% skipped.
				filered+=1
				continue

			psi_df = pd.DataFrame.from_dict(psi, orient='index').transpose()[ref_list]
			if single_plot:
				indiviualPlot(psi_df, event_id, ref_list)
			if group_plot:
				print step, j
				ax_i = plt.subplot2grid((step,11), (j-r_start*1,0), colspan=10, rowspan=1)
				##SHOULD SET TO area. use width just for prelim
				sns_plot=sns.violinplot(data=psi_df,ax=ax_i,inner="box",cut=0,scale='width',linewidth=1.5)
				sns_plot.set_yticks(np.arange(0,1.1,0.5))     
				sns_plot.set(xticklabels=[])
				sns.despine(offset=0, trim=False)
				sns_plot.text(15.7, 0.4, event_id.split(':')[1], horizontalalignment='left', size='large', color='black',fontweight='bold')
		if group_plot:
			pdf.savefig(fig)
			pdf.close()


if __name__ == '__main__':
	main()

