import numpy as np
import sys
import os, glob, pyBigWig, argparse
from scipy import stats
import statsmodels.stats.weightstats as smw
from . import config
import warnings
warnings.filterwarnings("ignore")

#python retreive_SJ_info.py test_sj.para SE event_list_test.txt sj_info
def read_SJMatrix_index(fn,outdir):
	index = {}
	for line in open(outdir+'/'+fn.split('/')[-1]+'.idx', 'r'):
		ele = line.strip().split()
		index[ele[0]] = int(ele[1])
	return index

def fetch_SJMatrix(eid, fn, delim, index, head_only):
	with open(fn, 'r') as f:
		if head_only:
			ele = f.readline().strip().split(delim)
			retrieved_text = np.asarray([ x.split('.aln')[0] for x in ele ])
		else:
			f.seek(index[eid], 0)
			retrieved_text = np.asarray(f.readline().strip().split(delim))
	return retrieved_text


def loadParametersRow(filter_para, panel_list):
	filter_cutoffs=''
	if filter_para.strip()!='':
		filter_cutoffs = map(float,filter_para.strip().split(' ')[0].split(','))
		filter_panel_list = filter_para.strip().split(' ')[1].split(',')	
		panel_list+=filter_panel_list
	else:
		filter_panel_list =[]
	return filter_cutoffs, filter_panel_list, panel_list



def readEventRow(row, header_line):
	if header_line=='' or header_line==False:
		rs=row.strip().split('\t')
		return rs
	else:
		rs=row.strip().split('\t')
		return dict(zip(header_line, rs))

def convertAS2SJevent(line_dict, splicing_event_type):#only for this script. Diff
	as_event=line_dict['as_event']
	event_s=as_event.split(':')
	if splicing_event_type=='SE':	
		event_row_list=[event_s[2]+':'+str(int(event_s[6])+1)+':'+event_s[4], event_s[2]+':'+str(int(event_s[5])+1)+':'+event_s[7], event_s[2]+':'+str(int(event_s[6])+1)+':'+event_s[7]]
	
	elif splicing_event_type=='A5SS':# Only use one junction for inc. Need to improve by updating db later
		event_row_list=[line_dict['chr']+':'+str(int(event_s[5])+1)+':'+event_s[8], line_dict['chr']+':'+str(int(event_s[7])+1)+':'+event_s[8]]	
		
	elif splicing_event_type=='A3SS': # Only use one junction for inc. Need to improve by updating db later
		event_row_list=[line_dict['chr']+':'+str(int(event_s[9])+1)+':'+event_s[4],line_dict['chr']+':'+str(int(event_s[9])+1)+':'+event_s[6]]	

	else:
		exit('[Error] Splicing event type not supported. Exiting.')
	return event_row_list, as_event

def summarizeSJ2ASevent(event_list_fin, splicing_event_type, sig_junction, outdir, out_prefix):
	fout_summary_fname=outdir+'/SJ.'+out_prefix+'.'+splicing_event_type+'.summary_by_sig_event.txt'
	fout_summary=open(fout_summary_fname,'w')
	for event_idx, event_row in enumerate(open(event_list_fin)):
		if event_idx==0:
			header_list=readEventRow(event_row,'')
			continue
		line_dict=readEventRow(event_row, header_list)
		event_row_list, as_event=convertAS2SJevent(line_dict, splicing_event_type)
		as_event_result=[]
		as_event_result_list=[]
		for k in event_row_list:
			if k not in sig_junction:
				as_event_result.append(False)
			else:
				as_event_result.append(True)
				as_event_result_list.append(k)
		if as_event_result[0]==as_event_result[1]==as_event_result[2]==True:
			fout_summary.write(as_event+'\tAll junctions\t'+';'.join(as_event_result_list)+'\n')
		elif as_event_result[0]==as_event_result[1]==as_event_result[2]==False:
			continue
		else:
			if as_event_result[0]==as_event_result[1]!=as_event_result[2]:
				fout_summary.write(as_event+'\tOnly alternative junctions\t'+';'.join(as_event_result_list)+'\n')
			else:
				fout_summary.write(as_event+'\tOther combination\t'+';'.join(as_event_result_list)+'\n')
	fout_summary.close()
	return fout_summary_fname




def main(args):
	###Loading Parameters####
	para_fin=args.parameter_file
	splicing_event_type=args.splicing_event_type
	if splicing_event_type!='SE':
		exit('[Error] Invalid AS event type.')
	event_list_fin=args.screening_result_event_list
	outdir=args.outdir.rstrip('/')
	
	os.system('mkdir -p '+outdir)
	fetching_sj_col=1
	out_prefix,db_dir,filter1_para,filter2_para,filter3_para=[l.strip() for l in open(para_fin)][:5]
	db_dir=db_dir.rstrip('/')
	if os.path.isdir(db_dir+'_sjc'): #automatically use db_sjc if in the same dir. Otherwise, use the user input db_dir
		db_dir=db_dir+'_sjc'
	panel_list=[out_prefix]

	filter1_cutoffs, filter1_panel_list, panel_list = loadParametersRow(filter1_para, panel_list)
	filter2_cutoffs, filter2_panel_list, panel_list = loadParametersRow(filter2_para, panel_list)
	filter3_cutoffs, filter3_panel_list, panel_list = loadParametersRow(filter3_para, panel_list)
	tumor_dict=dict.fromkeys(filter2_panel_list,'')
	tumor_dict[out_prefix]=''
	pvalue_cutoff_normal=''; pvalue_cutoff_tumor=''
	filter1_group_cutoff=''; filter2_group_cutoff=''; filter3_group_cutoff='';
	if filter1_cutoffs!='':
		pvalue_cutoff_normal,filter1_group_cutoff=filter1_cutoffs[3:]
	if filter2_cutoffs!='':
		pvalue_cutoff_tumor,filter2_group_cutoff=filter2_cutoffs[3:]
	if filter3_cutoffs!='':
		pvalue_cutoff_normal,filter3_group_cutoff=filter3_cutoffs[3:]

	inc_read_cov_cutoff=int(args.inc_read_cov_cutoff)#2
	event_read_cov_cutoff=int(args.event_read_cov_cutoff)#10

	##Load IRIS reference panels to 'fin_list', 'index'
	index={}
	fin_list={}
	for group_name in panel_list:
		fin_list[group_name]=db_dir+'/'+group_name+'/sjc_matrix/SJ_count.'+group_name+'.txt'
	for group in fin_list.keys():
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_SJMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))

	print('[INFO] Done loading index '+' '.join(panel_list))
	tot=config.file_len(event_list_fin)-1
	fout_ijc=open(outdir+'/'+event_list_fin.split('/')[-1]+'.ijc_info.txt','w')
	header_line=[]
	sample_size={}

	header_line=['ijc_ratio', 'mean_ijc_by_group', 'percent_sample_imbalanced', 'sample_imbalanced_by_group']
	header_list=[]
	print('[INFO] Retrieving inclusion junction info')
	for event_idx, event_row in enumerate(open(event_list_fin)):
		config.update_progress(event_idx/(0.0+tot))
		if event_idx==0:
			header_list=readEventRow(event_row,'')
			fout_ijc.write(event_row.strip()+'\t'+'\t'.join(header_line)+'\n')
			continue
		line_dict=readEventRow(event_row, header_list)
		event_row_list, as_event=convertAS2SJevent(line_dict, splicing_event_type)

		forms={}
		if splicing_event_type=='SE':
			inc1, inc2, skp= event_row_list
			forms={'inc1':inc1,'inc2':inc2,'skp':skp}
		else:
			exit('[Error] Invalid AS event type.')
	#Initiate psi matrix by each row to 'sj' 

		sj={}
		for form in forms:
			k=forms[form]
			sj[form]={}
			for group in panel_list:
				random_key=index[group].keys()[0]
				sample_names=map(str,fetch_SJMatrix(random_key,fin_list[group],'\t',index[group],True)[fetching_sj_col:])
				sample_size[group]=len(sample_names)
				if k in index[group]:
					sj[form][group]=map(int,fetch_SJMatrix(k,fin_list[group],'\t',index[group], False)[fetching_sj_col:])
				else:
					sj[form][group]=[0]*sample_size[group]

		imbalance={}
		total_sample={}
		inc_ratio={}
		i1_list={}
		i2_list={}
		for group in panel_list:
			i1=map(int,sj['inc1'][group])
			i2=map(int,sj['inc2'][group])
			s=map(int,sj['skp'][group])
			i1_list[group]=[]
			i2_list[group]=[]
			imbalance[group]=0
			total_sample[group]=0
			for pidx in range(0, len(sj['inc1'][group])):
				if i1[pidx]+i2[pidx]<inc_read_cov_cutoff or i1[pidx]+i2[pidx]+s[pidx]<event_read_cov_cutoff:
					continue
				total_sample[group]+=1
				ratio=i1[pidx]/(i2[pidx]+0.01)
				if ratio>2 or ratio<0.5:
					imbalance[group]+=1
				i1_list[group].append(i1[pidx])
				i2_list[group].append(i2[pidx])

		i_by_group='|'.join([str(round(np.mean(i1_list[group]),1))+','+str(round(np.mean(i2_list[group]),1)) for group in panel_list])
		i1_sum=sum(sum(i1_list[g]) for g in panel_list)
		i2_sum=sum(sum(i2_list[g]) for g in panel_list)
		ratio=round(i1_sum/(i2_sum+0.01),3)

		imbalance_by_group='|'.join([str(imbalance[group])+','+str(total_sample[group]) for group in panel_list])
		imb_count=sum(imbalance[group] for group in panel_list)
		tot_count=sum(total_sample[group] for group in panel_list)
		imb_perc='-'
		if tot_count!=0:
			imb_perc=round(imb_count/(tot_count+0.0),3)
		fout_ijc.write(event_row.rstrip()+'\t'+'\t'.join([as_event,str(ratio), i_by_group, str(imb_perc),imbalance_by_group])+'\n')

	fout_ijc.close()

		
if __name__ == '__main__':
	main()
