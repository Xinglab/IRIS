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


def get_column_sums(input_path):
    sums = list()
    with open(input_path, 'rt') as handle:
        for line_i, line in enumerate(handle):
            columns = line.strip().split('\t')
            if line_i == 0:
                headers = columns
                num_columns_to_sum = len(headers) - 1
                sums.extend([0] * num_columns_to_sum)
                continue

            for column_i, col in enumerate(columns):
                if column_i == 0:
                    continue

                sum_i = column_i - 1
                col_val = float(col)
                sums[sum_i] += col_val

    return sums


def factors_from_sums(sums):
    million = 1000000
    factors = list()
    for col_sum in sums:
        if col_sum == 0:
            factor = None
        else:
            factor = million / col_sum

        factors.append(factor)

    return factors


def mult_cols_by_factors(input_path, output_path, factors):
    with open(input_path, 'rt') as in_handle:
        with open(output_path, 'wt') as out_handle:
            for line_i, line in enumerate(in_handle):
                columns = line.strip().split('\t')
                if line_i == 0:
                    out_handle.write(line)
                    continue

                new_cols = list()
                for col_i, col in enumerate(columns):
                    if col_i == 0:
                        new_cols.append(col)
                        continue

                    factor = factors[col_i - 1]
                    if factor is None:
                        new_cols.append('NA')
                        continue

                    col_val = float(col)
                    with_factor = factor * col_val
                    str_with_factor = str(with_factor)
                    if str_with_factor == '0.0':
                        str_with_factor = '0'

                    new_cols.append(str_with_factor)

                out_handle.write('{}\n'.format('\t'.join(new_cols)))


def make_cpm_factor(SJ_count_input_path, CPM_factor_output_path):
    #output_path = name + '_SJC_factors.txt'
    column_sums = get_column_sums(SJ_count_input_path)
    factors = factors_from_sums(column_sums)
    fout = open(CPM_factor_output_path, 'w')
    with open(SJ_count_input_path, 'rt') as f:
        header = f.readline().strip().split('\t')[1:]
        fout.write('\t'.join(header) + '\n')
        fout.write('\t'.join([str(x) for x in factors]) + '\n')
    fout.close()


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


def statTest(g1,g2, direction, non_parametric):
	if direction != 'equivalence':
		if non_parametric:
			pvalue=stats.mannwhitneyu(g1,g2, alternative=direction)[1]
		else:
			pvalue=smw.ttest_ind(g1,g2, alternative=direction, usevar='unequal')[1]
			#pvalue=smw.ttest_ind(g1,g2, alternative=direction)[1]
	else:
		threshold_tost = 0.05
		pvalue=smw.ttost_ind(g1,g2,-threshold_tost,threshold_tost,usevar='unequal')[0] #equivalence test
	return pvalue

def statTest_minSampleCount(g1,g2, direction, non_parametric):#Only enabled when filters out by min_sample_count. With enough sample, using the default setting assume equal var for both groups is ok.
	if direction != 'equivalence':
		if non_parametric:
			pvalue=stats.mannwhitneyu(g1,g2, alternative=direction)[1]
		else:
			pvalue=smw.ttest_ind(g1,g2, alternative=direction)[1]
			#pvalue=smw.ttest_ind(g1,g2, alternative=direction)[1]
	else:
		threshold_tost = 0.05
		pvalue=smw.ttost_ind(g1,g2,-threshold_tost,threshold_tost)[0] #equivalence test
	return pvalue

def groupTest(g1,g2, non_parametric=True, direction='greater', min_sample_count=False):
	if non_parametric==False:
		direction='larger'
	g1=np.array(g1)
	g2=np.array(g2)
	g1=g1[~np.isnan(g1)]
	g2=g2[~np.isnan(g2)]
	delta_cpm=np.nanmean(g1)-np.nanmean(g2)
	g2_cpm_mean=np.nanmean(g2)	
	if g2_cpm_mean==0:
		g2_cpm_mean=0.00001#pseudo count for fold change calculation
	foc_cpm=round(np.nanmean(g1)/g2_cpm_mean,2)
	if min_sample_count:
		pvalue = statTest_minSampleCount(g1, g2, direction, non_parametric)
	else:	
		pvalue = statTest(g1, g2, direction, non_parametric)
	return [pvalue, delta_cpm, g2_cpm_mean]

def loadSigJunction(fin):
	sig_junction={}
	for i,l in enumerate(open(fin)):
		if i==0:
			continue
		sig_junction[l.strip().split('\t')[0]]=''
	return sig_junction

def summarizeSJ2ASevent(event_list_fin, splicing_event_type, sig_junction, outdir, out_prefix):
	fout_summary_fname=outdir+'/CPM.'+out_prefix+'.'+splicing_event_type+'.summary_by_sig_event.txt'
	fout_summary=open(fout_summary_fname,'w')
	for event_idx, event_row in enumerate(open(event_list_fin)):
		if event_idx==0:
			header_list=readEventRow(event_row,'')
			continue
		line_dict=readEventRow(event_row, header_list)
		event_row_list, as_event=convert2SJASevent(line_dict, splicing_event_type)
		as_event_result=[]
		as_event_result_list=[]
		as_event_result_info=[]
		cpm_info=[]
		for k in event_row_list:
			if k not in sig_junction:
				as_event_result.append(False)
			else:
				as_event_result.append(True)
				as_event_result_list.append(k)
				as_event_result_info.append(sig_junction[k])
				cpm_info.append(float(sig_junction[k].split(':')[0]))
		if as_event_result[0]==as_event_result[1]==as_event_result[2]==True:
			fout_summary.write(as_event+'\tAll junctions\t'+';'.join(as_event_result_list)+'\t'+';'.join(as_event_result_info)+'\t'+str(max(cpm_info))+'\n')
		elif as_event_result[0]==as_event_result[1]==as_event_result[2]==False:
			continue
		else:
			if as_event_result[0]==as_event_result[1]!=as_event_result[2]:
				fout_summary.write(as_event+'\tOnly alternative junctions\t'+';'.join(as_event_result_list)+'\t'+';'.join(as_event_result_info)+'\t'+str(max(cpm_info))+'\n')
			else:
				fout_summary.write(as_event+'\tOther combination\t'+';'.join(as_event_result_list)+'\t'+';'.join(as_event_result_info)+'\t'+str(max(cpm_info))+'\n')
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
		pvalue_cutoff_normal,filter1_group_cutoff=[filter1_cutoffs[0],filter1_cutoffs[4]]
	if filter2_cutoffs!='':
		pvalue_cutoff_tumor,filter2_group_cutoff=[filter2_cutoffs[0],filter2_cutoffs[4]]
	if filter3_cutoffs!='':
		pvalue_cutoff_normal,filter3_group_cutoff=[filter3_cutoffs[0],filter3_cutoffs[4]]

	#tumor_read_cov_cutoff=int(args.tumor_read_cov_cutoff)#5
	#normal_read_cov_cutoff=int(args.normal_read_cov_cutoff)#2
	
	
	##Load IRIS reference panels to 'fin_list', 'index'
	index={}
	fin_list={}
	factor_list={}
	for group_name in panel_list:
		fin_list[group_name]=db_dir+'/'+group_name+'/sjc_matrix/SJ_count.'+group_name+'.txt'
	for group in fin_list.keys():
		#Check CPM factor file. If not exist, create one.
		if not os.path.isfile(fin_list[group]+'.factors.tsv'):
			print '[INFO] Factor file doesn\'t exist. Creating the library size factor file. '+fin_list[group]+'.factors.tsv'
			make_cpm_factor(fin_list[group], fin_list[group]+'.factors.tsv')
		else:
			print '[INFO] Using existing factor file. '+fin_list[group]+'.factors.tsv'
		#Load CPM factor to dict
		for fi,l in enumerate(open(fin_list[group]+'.factors.tsv')):
			if fi==1:		
				factor_list[group]=l.strip().split('\t')
				
		#Check index file
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_SJMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))


	tot=config.file_len(event_list_fin)-1
	if tot==0:
		exit('[Terminated] No test performed because no testable events. Check input or filtering parameteres.') #Modified 2021
	print '[INFO] IRIS screening - started. Total input events:', tot+1
	if use_existing_test_result==False:
		fout_cpm_count=open(outdir+'/CPM.'+out_prefix+'.'+splicing_event_type+'.test_all.txt','w')
		header_line=[]
		sample_size={}
		for group in panel_list:
			random_key=index[group].keys()[0]
			sample_names=map(str,fetch_SJMatrix(random_key,fin_list[group],'\t',index[group],True)[fetching_sj_col:])
			sample_size[group]=len(sample_names)
			if group==out_prefix:
				header_line+=[out_prefix+'_CPM']
				continue
			header_line+=[group+'_CPM', group+'_change', group+'_pvalue']
		fout_cpm_count.write('Junction\t'+'\t'.join(header_line)+'\n')
	
	header_line=[]
	fout_cpm_sig_fname=outdir+'/CPM.'+out_prefix+'.'+splicing_event_type+'.test_sig.txt'
	fout_cpm_sig=open(fout_cpm_sig_fname,'w')
	for group in panel_list:
		if group==out_prefix:
			header_line+=[out_prefix+'_CPM']
			continue
		header_line+=[group+'_CPM', group+'_change', group+'_pvalue']
	fout_cpm_sig.write('Junction\t'+'\t'.join(header_line)+'\n')
	
	sig_junction={}
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
				cpm={}
				cpm_test=''
				for group in panel_list:
					#CPM = SJ count * factor 
					#factor=10^6/sum SJ count
					cpm[group]=[round(float(v)*float(factor_list[group][i]),2) for i,v in enumerate(sj[group])]
					
					if group==out_prefix:
						cpm_mean_input=np.nanmean(cpm[group])
						write_sj_list=[cpm_mean_input]
						if set(cpm[group])==set([0]):#Handle all-zero cases. 1. mannwhiteney U will have error; 2. skip to speed up this one-sided test 
							write_sj_list=['NA','NA',1]*(len(panel_list)-1)
							break 
						continue
					else:
						if group in filter2_panel_list: # filter2 always require filer1!!! if no filter1, all references should be defined in filter3
							cpm_test=groupTest(cpm[group],cpm[filter1_panel_list[0]])
						else:
							cpm_test=groupTest(cpm[out_prefix],cpm[group])

						write_sj_list+=[cpm_test[2], cpm_test[1], cpm_test[0]]
						#determine difference of a junction
						if group in filter1_panel_list:
							if cpm_test[0]<=pvalue_cutoff_normal:
								significant_normal_match+=1
						elif group in filter2_panel_list:
							if cpm_test[0]<=pvalue_cutoff_tumor:
								significant_tumor+=1
						else:
							if cpm_test[0]<=pvalue_cutoff_normal:
								significant_normal+=1 
				if (significant_normal_match>=filter1_group_cutoff or filter1_group_cutoff=='') and (significant_tumor>=filter2_group_cutoff or filter2_group_cutoff=='') and (significant_normal>=filter3_group_cutoff or filter3_group_cutoff==''):
					fout_cpm_sig.write(k+'\t'+'\t'.join(map(str,write_sj_list))+'\n')
					sig_junction[k]='|'.join(map(str, [write_sj_list[0],significant_normal_match,significant_tumor,significant_normal]))
				fout_cpm_count.write(k+'\t'+'\t'.join(map(str,write_sj_list))+'\n')	
		fout_cpm_count.close()
		fout_cpm_sig.close()

	else:
		print 'Use existing testing result.'
		fout_cpm_count_name=outdir+'/CPM.'+out_prefix+'.'+splicing_event_type+'.test_all.txt'
		for i, l in enumerate(open(fout_cpm_count_name)):
			if i==0:
				header=l.strip().split('\t')
				group_list=map(lambda x: x.split('_CPM')[0], header[2::3])
			else:
				ls=l.strip().split('\t')
				if ls[1]=='NA':
					continue
				cpm_value= map(float,ls[2::3])# don't do map because '-'
				change_value = map(float,ls[3::3])
				p_value= map(float,ls[4::3])
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
					fout_cpm_sig.write(l.strip()+'\n')
					sig_junction[ls[0]]=':'.join(map(str, [round(float(ls[1]),2),significant_normal_match,significant_tumor,significant_normal]))
		fout_cpm_sig.close()

	#sig_junction=loadSigJunction(fout_cpm_sig_fname)
	fout_summary_fname=summarizeSJ2ASevent(event_list_fin, splicing_event_type, sig_junction, outdir, out_prefix)


if __name__ == '__main__':
	main()
