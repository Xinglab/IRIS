import numpy as np
import sys
import os, glob, pyBigWig, argparse
from scipy import stats
import statsmodels.stats.weightstats as smw
from . import config
import warnings
warnings.filterwarnings("ignore")

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

def convert2SJevent(line_dict, splicing_event_type):
	if splicing_event_type=='SE':
		event_row_list=[line_dict['chr']+':'+str(int(line_dict['upstreamEE'])+1)+':'+line_dict['exonStart'],line_dict['chr']+':'+str(int(line_dict['exonEnd'])+1)+':'+line_dict['downstreamES'], line_dict['chr']+':'+str(int(line_dict['upstreamEE'])+1)+':'+line_dict['downstreamES']]
	elif splicing_event_type=='A5SS':# Only use one junction for inc. Need to improve by updating db later
		event_row_list=[line_dict['chr']+':'+str(int(line_dict['longExonEnd'])+1)+':'+line_dict['flankingES'],line_dict['chr']+':'+str(int(line_dict['shortEE'])+1)+':'+line_dict['flankingES']]	
	elif splicing_event_type=='A3SS': # Only use one junction for inc. Need to improve by updating db later
		event_row_list=[line_dict['chr']+':'+str(int(line_dict['flankingEE'])+1)+':'+line_dict['longExonStart'],line_dict['chr']+':'+str(int(line_dict['flankingEE'])+1)+':'+line_dict['shortES']]	
	else:
		exit('splicine event type not supported. Exiting.')
	return event_row_list

def convert2SJASevent(line_dict, splicing_event_type):
	if splicing_event_type=='SE':
		event_row_list=[line_dict['chr']+':'+str(int(line_dict['upstreamEE'])+1)+':'+line_dict['exonStart'],line_dict['chr']+':'+str(int(line_dict['exonEnd'])+1)+':'+line_dict['downstreamES'], line_dict['chr']+':'+str(int(line_dict['upstreamEE'])+1)+':'+line_dict['downstreamES']]	
		as_event=line_dict['AC'].strip('"').split('.')[0]+':'+line_dict['GeneName'].strip('"')+':'+line_dict['chr']+':'+line_dict['strand']+':'+line_dict['exonStart']+':'+line_dict['exonEnd']+':'+line_dict['upstreamEE']+':'+line_dict['downstreamES']
	elif splicing_event_type=='A5SS':# Only use one junction for inc. Need to improve by updating db later
		event_row_list=[line_dict['chr']+':'+str(int(line_dict['longExonEnd'])+1)+':'+line_dict['flankingES'],line_dict['chr']+':'+str(int(line_dict['shortEE'])+1)+':'+line_dict['flankingES']]	
		as_event=line_dict['AC'].strip('"').split('.')[0]+':'+line_dict['GeneName'].strip('"')+':'+line_dict['chr']+':'+line_dict['strand']+':'+line_dict['longExonStart']+':'+line_dict['longExonEnd']+':'+line_dict['shortES']+':'+line_dict['shortEE']+':'+line_dict['flankingES']+':'+line_dict['flankingEE']
	elif splicing_event_type=='A3SS': # Only use one junction for inc. Need to improve by updating db later
		event_row_list=[line_dict['chr']+':'+str(int(line_dict['flankingEE'])+1)+':'+line_dict['longExonStart'],line_dict['chr']+':'+str(int(line_dict['flankingEE'])+1)+':'+line_dict['shortES']]	
		as_event=line_dict['AC'].strip('"').split('.')[0]+':'+line_dict['GeneName'].strip('"')+':'+line_dict['chr']+':'+line_dict['strand']+':'+line_dict['longExonStart']+':'+line_dict['longExonEnd']+':'+line_dict['shortES']+':'+line_dict['shortEE']+':'+line_dict['flankingES']+':'+line_dict['flankingEE']
	else:
		exit('splicine event type not supported. Exiting.')
	return event_row_list, as_event

def loadSigJunction(fin):
	sig_junction={}
	for i,l in enumerate(open(fin)):
		if i==0:
			continue
		sig_junction[l.strip().split('\t')[0]]=''
	return sig_junction

def summarizeSJ2ASevent(event_list_fin, splicing_event_type, sig_junction, outdir, out_prefix):
	fout_summary_fname=outdir+'/SJ.'+out_prefix+'.'+splicing_event_type+'.summary_by_sig_event.txt'
	fout_summary=open(fout_summary_fname,'w')
	for event_idx, event_row in enumerate(open(event_list_fin)):
		if event_idx==0:
			header_list=readEventRow(event_row,'')
			continue
		line_dict=readEventRow(event_row, header_list)
		event_row_list, as_event=convert2SJASevent(line_dict, splicing_event_type)
		as_event_result=[]
		as_event_result_list=[]
		for k in event_row_list:
			if k not in sig_junction:
				as_event_result.append(False)
			else:
				as_event_result.append(True)
				as_event_result_list.append(k)

		if splicing_event_type == 'SE':
			all_true = as_event_result[0]==as_event_result[1]==as_event_result[2]==True
			all_false = as_event_result[0]==as_event_result[1]==as_event_result[2]==False
			only_alt = as_event_result[0]==as_event_result[1]!=as_event_result[2]
		else:
			all_true = as_event_result[0]==as_event_result[1]==True
			all_false = as_event_result[0]==as_event_result[1]==False
			only_alt = False
		if all_true:
			fout_summary.write(as_event+'\tAll junctions\t'+';'.join(as_event_result_list)+'\n')
		elif all_false:
			continue
		else:
			if only_alt:
				fout_summary.write(as_event+'\tOnly alternative junctions\t'+';'.join(as_event_result_list)+'\n')
			else:
				fout_summary.write(as_event+'\tOther combination\t'+';'.join(as_event_result_list)+'\n')
	fout_summary.close()
	return fout_summary_fname

def main(args):
	###Loading Parameters####
	para_fin=args.parameter_file
	splicing_event_type=args.splicing_event_type
	event_list_fin=args.event_list_file
	use_existing_test_result=args.use_existing_test_result

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

	tumor_read_cov_cutoff=int(args.tumor_read_cov_cutoff)#5
	normal_read_cov_cutoff=int(args.normal_read_cov_cutoff)#2
	# if filter1_panel_list==[] and filter2_panel_list==[] and filter3_panel_list==[] and test_mode[0]!='summary':
	# 	exit("[Error] No filtering required in parameteres file. exit!")

	##Load IRIS reference panels to 'fin_list', 'index'
	index={}
	fin_list={}
	for group_name in panel_list:
		fin_list[group_name]=db_dir+'/'+group_name+'/sjc_matrix/SJ_count.'+group_name+'.txt'
	for group in fin_list.keys():
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_SJMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))


	tot=config.file_len(event_list_fin)-1
	if tot==0:
		exit('[Ended] no test performed because no testable events. Check input or filtering parameteres.') #Modified 2021
	print '[INFO] IRIS screening - started. Total input events:', tot+1
	if use_existing_test_result==False:
		fout_sj_count=open(outdir+'/SJ.'+out_prefix+'.'+splicing_event_type+'.test_all.txt','w')
		header_line=[]
		sample_size={}
		for group in panel_list:
			random_key=index[group].keys()[0]
			sample_names=map(str,fetch_SJMatrix(random_key,fin_list[group],'\t',index[group],True)[fetching_sj_col:])
			sample_size[group]=len(sample_names)
			if group==out_prefix:
				header_line+=[out_prefix+'_carrier_number', out_prefix+'_fraction']
				continue
			header_line+=[group+'_carrier_number', group+'_fraction', group+'_pvalue']
		fout_sj_count.write('Junction\t'+'\t'.join(header_line)+'\n')
	
	header_line=[]
	fout_sj_sig_fname=outdir+'/SJ.'+out_prefix+'.'+splicing_event_type+'.test_sig.txt'
	fout_sj_sig=open(fout_sj_sig_fname,'w')
	for group in panel_list:
		if group==out_prefix:
			header_line+=[out_prefix+'_carrier_number', out_prefix+'_fraction']
			continue
		header_line+=[group+'_carrier_number', group+'_fraction', group+'_pvalue']
	fout_sj_sig.write('Junction\t'+'\t'.join(header_line)+'\n')

	if use_existing_test_result==False:
		header_list=[]
		junction_dict={}
		for event_idx, event_row in enumerate(open(event_list_fin)):
			if event_idx==0:
				header_list=readEventRow(event_row,'')
				continue
			line_dict=readEventRow(event_row, header_list)
			event_row_list=convert2SJevent(line_dict, splicing_event_type)
			for k in event_row_list:
				if k not in junction_dict:
					junction_dict[k]=''
				else:
					continue
				config.update_progress(event_idx/(0.0+tot))
				
				#Initiate psi matrix by each row to 'sj' 
				sj={}
				for group in panel_list:
					if k in index[group]:
						sj[group]=map(int,fetch_SJMatrix(k,fin_list[group],'\t',index[group], False)[fetching_sj_col:])
					else:
						sj[group]=[0]*sample_size[group]
				write_sj_list=[]
				significant_normal_match=0
				significant_normal=0
				significant_tumor=0
				prevalence={}
				prevalence_test=''
				for group in panel_list:
					if group in tumor_dict:
						prevalence[group]=sum(v>=tumor_read_cov_cutoff for v in sj[group])
					else:
						prevalence[group]=sum(v>=normal_read_cov_cutoff for v in sj[group])

					if group==out_prefix:
						write_sj_list=[prevalence[group],prevalence[group]/(0.0+sample_size[group])]
						continue
					else:
						if group in filter2_panel_list: # filter2 always require filer1!!! if no filter1, all references should be defined in filter3
							prevalence_test=stats.fisher_exact([[prevalence[group],sample_size[group]-prevalence[group]],[prevalence[filter1_panel_list[0]], sample_size[filter1_panel_list[0]]-prevalence[filter1_panel_list[0]]]], alternative='greater')
						else:
							prevalence_test=stats.fisher_exact([[prevalence[out_prefix],sample_size[out_prefix]-prevalence[out_prefix]],[prevalence[group], sample_size[group]-prevalence[group]]], alternative='greater')
						write_sj_list+=[prevalence[group], prevalence[group]/(0.0+sample_size[group]), prevalence_test[1]]
						#determine difference of a junction
						if group in filter1_panel_list:
							if prevalence_test[1]<=pvalue_cutoff_normal:
								significant_normal_match+=1
						elif group in filter2_panel_list:
							if prevalence_test[1]<=pvalue_cutoff_tumor:
								significant_tumor+=1
						else:
							if prevalence_test[1]<=pvalue_cutoff_normal:
								significant_normal+=1 
				if (significant_normal_match>=filter1_group_cutoff or filter1_group_cutoff=='') and (significant_tumor>=filter2_group_cutoff or filter2_group_cutoff=='') and (significant_normal>=filter3_group_cutoff or filter3_group_cutoff==''):
					fout_sj_sig.write(k+'\t'+'\t'.join(map(str,write_sj_list))+'\n')
				fout_sj_count.write(k+'\t'+'\t'.join(map(str,write_sj_list))+'\n')	
		fout_sj_count.close()
		fout_sj_sig.close()

	else:
		print 'Use existing testing result.'
		fout_sj_count_name=outdir+'/SJ.'+out_prefix+'.'+splicing_event_type+'.test_all.txt'
		for i, l in enumerate(open(fout_sj_count_name)):
			if i==0:
				header=l.strip().split('\t')
				group_list=map(lambda x: x.split('_carrier_number')[0], header[3::3])
			else:
				ls=l.strip().split('\t')
				prevalence_value= map(int,ls[3::3])# don't do map because '-'
				percent_value = map(float,ls[4::3])
				p_value= map(float,ls[5::3])
				significant_normal_match=0
				significant_normal=0
				significant_tumor=0
				#determine difference of a junction
				for j,group in enumerate(group_list):
					if group in filter1_panel_list:
						if p_value[j]<=pvalue_cutoff_normal:
							significant_normal_match+=1
					elif group in filter2_panel_list:
						if p_value[j]<=pvalue_cutoff_tumor:
							significant_tumor+=1
					else:
						if p_value[j]<=pvalue_cutoff_normal:
							significant_normal+=1
				if (significant_normal_match>=filter1_group_cutoff or filter1_group_cutoff=='') and (significant_tumor>=filter2_group_cutoff or filter2_group_cutoff=='') and (significant_normal>=filter3_group_cutoff or filter3_group_cutoff==''):
					fout_sj_sig.write(l.strip()+'\n')
		fout_sj_sig.close()

	sig_junction=loadSigJunction(fout_sj_sig_fname)
	fout_summary_fname=summarizeSJ2ASevent(event_list_fin, splicing_event_type, sig_junction, outdir, out_prefix)


if __name__ == '__main__':
	main()
