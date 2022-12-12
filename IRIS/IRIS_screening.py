import numpy as np
import os, glob, pyBigWig, argparse
from scipy import stats
import statsmodels.stats.weightstats as smw
from . import config
import warnings
warnings.filterwarnings("ignore")


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
				ele[header['upstreamEE']], ele[header['exonStart_0base']], ele[header['exonEnd']], ele[header['downstreamES']]])
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
	for line in open(outdir+'/'+fn.split('/')[-1]+'.idx', 'r'):
		ele = line.strip().split()
		index[ele[0]] = int(ele[1])
	return index

def fetch_PsiMatrix(eid, fn, outdir, delim, index=None):
	if index is None:
		index = read_PsiMatrix_index(fn,outdir)
	with open(fn, 'r') as f:
		ele = f.readline().strip().split(delim)
		header = np.asarray([ x.split('.')[0] if x.startswith('SRR') else x for x in ele ])
		f.seek(index[eid], 0)
		data = np.asarray(f.readline().strip().split(delim))
	return (header, data)

def openTestingFout(outdir, out_prefix, splicing_event_type, summary_file, panel_list, test_mode, fin_name=''):
	header=['as_event','meanPSI','Q1PSI','Q3PSI']
	if test_mode=='group':
		header_prefix=['_pVal','_deltaPSI','_tumorFC']
	if test_mode=='personalized':
		header_prefix=['_modifiedPctl','_deltaPSI','_tumorFC']
	fout_name=outdir+'/'+out_prefix+'.'+splicing_event_type+'.test.all_'+fin_name+'.txt'
	fout=open(fout_name,'w')
	if summary_file==False:
		header+=['\t'.join(map(lambda x:ref+x ,header_prefix)) for ref in panel_list if ref!=out_prefix]
	fout.write('\t'.join(header)+'\n')
	return fout, fout_name

def openScreeningFout(outdir, out_prefix, splicing_event_type, fout_name):
	fout=open(outdir+'/'+out_prefix+'.'+splicing_event_type+'.'+fout_name+'.txt','w')
	fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('as_event','meanPSI','Q1PSI','Q3PSI','deltaPSI','fc_of_tumor_isoform','tissue_matched_normal_panel','tumor_panel','normal_panel','tag','mappability','mappability_tag'))
	return fout

def writeSummaryFile(out_prefix, splicing_event_type, db_dir, index, fout, fetching_data_col): #no screening, preparing for MS search
	fin_name=db_dir+'/'+out_prefix+'/splicing_matrix/splicing_matrix.'+splicing_event_type+'.cov10.'+group_name+'.txt'
	index[out_prefix]=read_PsiMatrix_index(fin_name,'/'.join(fin_name.split('/')[:-1]))
	for j,k in enumerate(index[out_prefix]):
		psi_event=map(float,fetch_PsiMatrix(k,fin_name,'.','\t',index[out_prefix])[1][fetching_data_col:])
		query_mean=[np.nanmean(psi_event),np.nanpercentile(psi_event,25),np.nanpercentile(psi_event,75)]
		result=[k]+query_mean
		fout.write('\t'.join(map(str,result))+'\n')
	fout.close()

def loadMappability(bigwig_fin):
	bw_map=pyBigWig.open(bigwig_fin)
	d=45
	return bw_map,d

def getMappability(splicing_event,bw_map,d):
	arr=splicing_event.split(':')
	chrom,strand,start,end,up,down=arr[2],arr[3],int(arr[4]),int(arr[5]),int(arr[6]),int(arr[7])
	up_mean=bw_map.stats("%s"%chrom,up-d,up,type="mean")[0]
	down_mean=bw_map.stats("%s"%chrom,down,down+d,type="mean")[0]
	if abs(start-end)<2*d:
		target_mean=bw_map.stats("%s"%chrom,start,end,type="mean")[0]
	else:
		target_left=bw_map.stats("%s"%chrom,start,start+d,type="mean")[0]
		target_right=bw_map.stats("%s"%chrom,end-d,end,type="mean")[0]
		target_mean=(target_right+target_left)/2
	if strand=='-':  #switch the order
		li=[up_mean,down_mean]
		up_mean,down_mean=li[1],li[0]

	mappability=[str(up_mean),str(target_mean),str(down_mean)]
	return mappability

def getDirection(filter1_panel_list, psi, test, non_parametric, screening_type):
	psi_primary=[]
	deltaPSI_primary=[]
	for primary_group in filter1_panel_list:
		if test[primary_group]!=['-']*3:
			psi_primary+=psi[primary_group]
			deltaPSI_primary.append(test[primary_group][1])
	if psi_primary==[]: #in case tissue-matched normal doesn't have data
		return [],[]
	else:
		direction='greater' if non_parametric else 'larger'
		if np.median(deltaPSI_primary)<=0:
			direction='less' if non_parametric else 'smaller'
		return psi_primary, direction
		
def calcTumorFormFoc(delta_psi, mean_psi):
	if delta_psi>0:
		return mean_psi/(mean_psi - delta_psi+10**-8)
	elif delta_psi<0:
		return (1- mean_psi)/ (1- mean_psi+ delta_psi+10**-8)
	else:
		return 1

def one2N(p1, g1, test_type):	
	p1=np.array(p1)
	g1=np.array(g1)
	Q1=np.nanpercentile(g1,25)
	Q3=np.nanpercentile(g1,75)
	delta_psi=float(p1-np.nanmean(g1))
	tumor_foc=calcTumorFormFoc(delta_psi, float(p1))
	IQR=Q3-Q1
	low_bound=max(Q1-1.5*IQR,0)
	high_bound=min(Q3+1.5*IQR,1)
	points_add=np.append(g1,p1)
	empirical_percentile = stats.percentileofscore(points_add,p1)
	if empirical_percentile>50:
		empirical_percentile=100- empirical_percentile	
	outlier_val=float(max(low_bound-p1,p1-high_bound))
	if test_type=='sig' and outlier_val>=0:
		return [empirical_percentile,delta_psi,tumor_foc]
	elif test_type=='equ' and outlier_val<=0:
		return [50-empirical_percentile, delta_psi, tumor_foc]
	else:
		return [np.nan,delta_psi,tumor_foc]

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

def groupTest(g1,g2, non_parametric=False, direction='two-sided', min_sample_count=False):
	g1=np.array(g1)
	g2=np.array(g2)
	g1=g1[~np.isnan(g1)]
	g2=g2[~np.isnan(g2)]
	delta_psi=np.nanmean(g1)-np.nanmean(g2)
	tumor_foc=calcTumorFormFoc(delta_psi,np.nanmean(g1))
	if min_sample_count:
		pvalue = statTest_minSampleCount(g1, g2, direction, non_parametric)
	else:	
		pvalue = statTest(g1, g2, direction, non_parametric)
	return [pvalue, delta_psi, tumor_foc]

def performTest(set_matched_tumor, has, j, group, screening_type_list, psi, out_prefix, non_parametric, test, filter1_panel_list, psi_primary, direction, min_sample_count): #PSI value-based screen allow two-sided or one-sided tests. Different from SJ count or CPM-based screen, where only one-sided test is needed.
	screening_type = screening_type_list[j]
	redirect_output = False
	test_result = ['-']*3 #For missing in non-eesential tests/comparisons
	if screening_type == 'association':
		if has[group]:
			test_result = groupTest(psi[out_prefix],psi[group], non_parametric,"two-sided", min_sample_count)
		return test_result

	has_matched_tumor = False if psi_primary == [] else True #for clarity

	if screening_type == 'recurrence':
		if set_matched_tumor and has_matched_tumor:#set_matched_tumor is redundent here. kept for future implemtation of additional output type.
			if has[group]:
				test_result = groupTest(psi[group],psi_primary, non_parametric, direction, min_sample_count)
		else:
			if has[group]:
				test_result = groupTest(psi[out_prefix],psi[group], non_parametric, "equivalence", min_sample_count)#No or equivalent testing -John's use case		
		return test_result

	if screening_type == 'association_high':
		if set_matched_tumor and has_matched_tumor:
			if has[group]:
				test_result = groupTest(psi[out_prefix],psi[group], non_parametric, direction, min_sample_count)
		else:
			if has[group]:
				test_result = groupTest(psi[out_prefix],psi[group], non_parametric, "two-sided", min_sample_count) #Two-sided testing 		
		return test_result

def summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, pval, deltaPSI, foc, screening_type_list):
	association_passed,recurrence_passed,specificity_positive,specificity_negative, specificity_testable=[0,0,0,0,0]
	primary_result, primary_result_foc=[[],[]] #take care of multiple tiisue matched norm
	deltapsi_list_voted,foc_list_voted=[[],[]] #if no tissue-matched norm, use median
	for i,group_type in enumerate(screening_type_list):
		if pval[i]=='-': #This is important - skip all missing, which is not useful for summarizing but not affecting consistancy.
			continue
		if screening_type_list[i]=='association': 
			primary_result.append(float(deltaPSI[i])) 
			primary_result_foc.append(float(foc[i]))
			if float(pval[i])<=filter1_cutoff_pval and abs(float(deltaPSI[i]))>=filter1_cutoff_dpsi and float(foc[i])>=filter1_cutoff_foc:
					association_passed+=1
					continue
		if screening_type_list[i]=='recurrence':
			if float(pval[i])<=filter2_cutoff_pval and abs(float(deltaPSI[i]))>=filter2_cutoff_dpsi:
					recurrence_passed+=1
					continue
		if screening_type_list[i]=='association_high':# TODO: judge set/has, then run
			specificity_testable+=1
			if float(pval[i])<=filter3_cutoff_pval and float(foc[i])>=filter3_cutoff_foc:
				deltapsi_list_voted.append(float(deltaPSI[i]))
				foc_list_voted.append(float(foc[i]))
				if float(deltaPSI[i])>=filter3_cutoff_dpsi:
					specificity_positive+=1
					continue
				if float(deltaPSI[i])<=-filter3_cutoff_dpsi:
					specificity_negative+=1
					continue
	return association_passed,recurrence_passed,specificity_positive,specificity_negative,specificity_testable, np.median(primary_result), np.median(primary_result_foc), deltapsi_list_voted,foc_list_voted

def defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, set_matched_tumor, specificity_panel_len, association_passed, recurrence_passed, specificity_positive, specificity_negative, specificity_testable, primary_result, primary_result_foc, deltapsi_list_voted, foc_list_voted, use_ratio):
	tag=[]
 	if association_passed>=filter1_group_cutoff:#Improvement? current: 0>='' is false
 		tag.append('associated')
 	if recurrence_passed>=filter2_group_cutoff:
 		tag.append('recurrent') 
	if set_matched_tumor:#TODO-FUTURE: take care of set-yes has-no redirected events
		if primary_result>0:
			tissue_specificity=specificity_positive
			ratio=False if use_ratio==False else specificity_positive/(specificity_testable+10**-8)>=filter3_group_cutoff/(specificity_panel_len+0.0)
			if specificity_positive>=filter3_group_cutoff or ratio:
				tag.append('high_assoc')
		else:
			tissue_specificity=specificity_negative
			ratio=False if use_ratio==False else specificity_negative/(specificity_testable+10**-8)>=filter3_group_cutoff/(specificity_panel_len+0.0)
			if specificity_negative>=filter3_group_cutoff or ratio:
				tag.append('high_assoc')
	else: 
		tissue_specificity=max(specificity_positive,specificity_negative)
		ratio=False if use_ratio==False else tissue_specificity/(specificity_testable+10**-8)>=filter3_group_cutoff/(specificity_panel_len+0.0)
		if tissue_specificity>=filter3_group_cutoff or ratio:
			tag.append('high_assoc')
	primary_deltapsi=primary_result if (primary_result!=0 and set_matched_tumor) else np.median(deltapsi_list_voted)#TODO
	primary_foc=primary_result_foc if (primary_result!=0 and set_matched_tumor) else np.median(foc_list_voted)

	return primary_deltapsi, primary_foc, tissue_specificity, tag

def mappability_write(k, bw_map, calc_length):
	mappability_list=getMappability(k, bw_map, calc_length)
	mappability_tag='PASS'
	min_map_score=min(map(float, mappability_list))
	if min_map_score<0.8:
		mappability_tag='MID'
		if min_map_score<0.6:#from 0.7
			mappability_tag='LOW'
	return mappability_tag, mappability_list

def loadBlacklistEvents(fin):
	BlacklistEvents={}
	for l in open(fin):
		ls=l.strip().split('\t')
		des=ls[0].split(':')
		ensg=des[0].split('_')[0].split('.')[0]
		des=':'.join([ensg]+des[1:]).strip(':')
		BlacklistEvents[des]=''
	return BlacklistEvents

def translationCMD(ref_genome, gtf, outdir, out_prefix, splicing_event_type, all_orf, ignore_annotation, remove_early_stop, fout_name):
	#uSE name space 
	argument_line=''
	if all_orf:
		argument_line+=' --all-orf'
	if ignore_annotation:
		argument_line+=' --ignore-annotation'
	if remove_early_stop:
		argument_line+=' --remove-early-stop'
	cmd_translation='IRIS translate '+outdir+'/'+out_prefix+'.'+splicing_event_type+'.'+fout_name+'.txt '+' -o '+outdir+'/'+splicing_event_type+'.'+fout_name+' -g '+ref_genome+' -t '+splicing_event_type+argument_line+' --gtf '+gtf
	print '[INFO] Working on translating: '+fout_name
	os.system(cmd_translation)

def loadParametersRow(filter_para, panel_list):
	filter_cutoffs=''
	if filter_para.strip()!='':
		filter_cutoffs = map(float,filter_para.strip().split(' ')[0].split(','))
		filter_panel_list = filter_para.strip().split(' ')[1].split(',')	
		panel_list+=filter_panel_list
	else:
		filter_panel_list =[]
	return filter_cutoffs, filter_panel_list, panel_list


def main(args):
	###Loading Parameters####
	index={}
	fin_list={}
	para_fin=args.parameter_fin
	splicing_event_type=args.splicing_event_type
	fetching_data_col=8 if splicing_event_type == 'SE' else 10
	out_prefix,db_dir,filter1_para,filter2_para,filter3_para,test_mode,use_ratio,blacklist_path,mappability_path,ref_genome=[l.strip() for l in open(para_fin)]
	panel_list=[out_prefix]
	test_mode=test_mode.split(' ')
	use_ratio=True if use_ratio=='True' else False
	blacklist_events={}
	min_sample_count=args.min_sample_count
	if min_sample_count:
		min_sample_count=int(min_sample_count)
        if blacklist_path!='':
            #if blacklist_path=='BRAIN_BLACKLIST_PATH':
	    blacklist_events=loadBlacklistEvents(blacklist_path)
	bw_map,calc_length=loadMappability(mappability_path)

	all_orf=args.all_orf
	ignore_annotation=args.ignore_annotation
	remove_early_stop=args.remove_early_stop
	use_existing_test_result=args.use_existing_test_result
	
	filter1_cutoffs, filter1_panel_list, panel_list = loadParametersRow(filter1_para, panel_list)
	filter2_cutoffs, filter2_panel_list, panel_list = loadParametersRow(filter2_para, panel_list)
	filter3_cutoffs, filter3_panel_list, panel_list = loadParametersRow(filter3_para, panel_list)
	
	if filter1_cutoffs!='':
		filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc, filter1_group_cutoff=filter1_cutoffs[:3]+[filter1_cutoffs[4]]
	else:
		filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc, filter1_group_cutoff=['','','','']
	if filter2_cutoffs!='':
		filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter2_group_cutoff=filter2_cutoffs[:3]+[filter2_cutoffs[4]]
	else:
		filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter2_group_cutoff=['','','','']
	if filter3_cutoffs!='':
		filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, filter3_group_cutoff=filter3_cutoffs[:3]+[filter3_cutoffs[4]]
	else:
		filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, filter3_group_cutoff=['','','','']
	if filter1_panel_list==[] and filter2_panel_list==[] and filter3_panel_list==[] and test_mode[0]!='summary':
		exit("[Error] No filtering required in parameteres file. exit!")
	
	non_parametric=False
	if len(test_mode)>1:
		non_parametric=True if test_mode[1]=='nonparametric' else False

	group_test=False if test_mode[0]!='group' else True
	individual_test=False if test_mode[0]!='personalized' else True
	summary_file=False if test_mode[0]!='summary' else True
	if [group_test,individual_test,summary_file]==[False,False,False]:
		exit('[Error] Need to choose one mode.exit!')
	association_panel_len=len(filter1_panel_list)
	recurrence_panel_len=len(filter2_panel_list)
	specificity_panel_len=len(filter3_panel_list)
	panel_count=sum(1 for i in [association_panel_len, recurrence_panel_len, specificity_panel_len] if i!=0)
	screening_type_list=['association']*association_panel_len+['recurrence']*recurrence_panel_len+['association_high']*specificity_panel_len
	set_matched_tumor= True if screening_type_list[0] == 'association' else False
        
        if args.translating:
            gtf=args.gtf
            if os.path.exists(gtf)==False:
                exit('[Error] No gtf file provided for translation. exit!')

	###Create Folders/Output####
	outdir=args.outdir.rstrip('/')
	os.system('mkdir -p  '+outdir)
	db_dir=db_dir.rstrip('/')

	if use_existing_test_result==False:
		fout_direct, fout_direct_name= openTestingFout(outdir, out_prefix, splicing_event_type, summary_file, panel_list, test_mode[0], 'guided')
		fout_redirect, fout_redirect_name=openTestingFout(outdir, out_prefix, splicing_event_type, summary_file, panel_list, test_mode[0], 'voted')
		if summary_file:
			writeSummaryFile(out_prefix, splicing_event_type, db_dir, index, fout_direct, fetching_data_col)
			exit()
		fout_filtered=open(outdir+'/'+out_prefix+'.'+splicing_event_type+'.notest.txt','w')
	else:
		fout_direct_name=outdir+'/'+out_prefix+'.'+splicing_event_type+'.test.all_guided.txt'
		fout_redirect_name=outdir+'/'+out_prefix+'.'+splicing_event_type+'.test.all_voted.txt'

	fout_primary=openScreeningFout(outdir, out_prefix, splicing_event_type, 'tier1')
	fout_prioritized=openScreeningFout(outdir,out_prefix, splicing_event_type, 'tier2tier3')

	##Load IRIS reference panels
	for group_name in panel_list:
		fin_list[group_name]=db_dir+'/'+group_name+'/splicing_matrix/splicing_matrix.'+splicing_event_type+'.cov10.'+group_name+'.txt'
	for group in fin_list.keys():
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_PsiMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))


	## Load and perform test by row/event
	if use_existing_test_result==False:
		has={}
		tot=len(index[out_prefix])-1
		print '[INFO] IRIS screen - started. Total input events:', tot+1
		
		for event_idx,k in enumerate(index[out_prefix]):
			config.update_progress(event_idx/(0.0+tot))
			
			#Initiate  
			for group in panel_list:
				if group!=out_prefix:
					has[group]=True
			psi={}
			has_count=0
			for group in panel_list:
				if k in index[group]:
					psi[group]=map(float,fetch_PsiMatrix(k,fin_list[group],'.','\t',index[group])[1][fetching_data_col:])
					has_count+=1
				else:
					has[group]=False
			#Filtering
			cat_psi=[]
			for i in psi:
				cat_psi+=psi[i]
			if abs(max(cat_psi)-min(cat_psi))<0.05:#if change less than 5% skipped and no comparison available
				fout_filtered.write('[Low Range]{}\t{}\t{}\n'.format(k,str(abs(max(cat_psi)-min(cat_psi))),str(has_count))) 
				continue
			if k in blacklist_events:
				fout_filtered.write('[Blacklisted]{}\t{}\t{}\n'.format(k,'-',str(has_count))) 
				continue
			if  has_count<=1:
				fout_filtered.write('[Unique in Input]{}\t{}\t{}\n'.format(k,'-',str(has_count)))
				continue
			if min_sample_count:
				sample_count=np.count_nonzero(~np.isnan(psi[out_prefix]))
				if sample_count<min_sample_count:
					fout_filtered.write('[Low Sample Count in Input]{}\t{}\t{}\n'.format(k,str(sample_count),str(has_count)))
					continue
			#Testing and store/output results(both group and individual mode)	
			if group_test:
				test={}
				redirect=False
				psi_primary=''
				direction=''
				query_mean=map(lambda k:round(k,2),[np.nanmean(psi[out_prefix]),np.nanpercentile(psi[out_prefix],25),np.nanpercentile(psi[out_prefix],75)])
				for j,group in enumerate(panel_list[1:]):
					if screening_type_list[j]!='association' and psi_primary=='':#redirect or find one-sided test direction
						psi_primary, direction= getDirection(filter1_panel_list, psi, test, non_parametric, screening_type_list[j])
						redirect = True if psi_primary == [] else False #redirect if tumor-matched normal is missing
					test[group]=performTest(set_matched_tumor, has, j, group, screening_type_list, psi, out_prefix, non_parametric, test, filter1_panel_list, psi_primary, direction, min_sample_count)				
				result=[k]+query_mean+['\t'.join(map(str,test[t])) for t in panel_list if t!=out_prefix]
				
				if redirect:
					fout_redirect.write('\t'.join(map(str,result))+'\n') #to redicted; summarize direction and calculate FDR differently.
				else:
					fout_direct.write('\t'.join(map(str,result))+'\n')

			elif individual_test:###TODO
				test={}
				query_mean=[psi[out_prefix][0],'-','-']
				for j,group in enumerate(panel_list[1:]):
					test[group]=['-']*3
					if has[group]:
						test[group]=one2N(psi[out_prefix],psi[group],screening_type_list[j])

				result=[k]+query_mean+['\t'.join(map(str,test[t])) for t in panel_list if t!=out_prefix]
				fout_direct.write('\t'.join(map(str,result))+'\n')###TODO
			
			else:
				exit('[Error] no test performed')
		fout_filtered.close()
		fout_redirect.close()
		fout_direct.close()

	##Summarize results and priortize screening candidates
	testing_intermediate_file = fout_direct_name if set_matched_tumor else fout_redirect_name

	tot=config.file_len(testing_intermediate_file)-1
	if tot==0:
		exit('[Ended] no test performed because no testable events. Check input or filtering parameteres.') #Modified 2021
	print '[INFO] IRIS screen - summarizing. Total events from last step:', tot
	for event_idx,l in enumerate(open(testing_intermediate_file)):
		config.update_progress(event_idx/(0.0+tot))
		if event_idx==0:
			continue
		ls=l.strip().split('\t')
		pval= ls[4::3]# don't do map because '-'
		deltaPSI = ls[5::3]
		foc= ls[6::3]

		association_passed,recurrence_passed,specificity_positive,specificity_negative,specificity_testable,primary_result, primary_result_foc, deltapsi_list_voted,foc_list_voted=summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, pval, deltaPSI, foc, screening_type_list)

		primary_deltapsi, primary_foc, tissue_specificity, tag = defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, set_matched_tumor, specificity_panel_len, association_passed,recurrence_passed,specificity_positive,specificity_negative,specificity_testable,primary_result, primary_result_foc, deltapsi_list_voted,foc_list_voted,use_ratio)

		if tag!=[]:
			if tag[0]=='associated':
				mappability_tag, mappability_list=mappability_write(ls[0], bw_map, calc_length)
				fout_primary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ls[0],ls[1],ls[2],ls[3],primary_deltapsi, primary_foc, str(association_passed)+'/'+str(association_panel_len),str(recurrence_passed)+'/'+str(recurrence_panel_len),str(tissue_specificity)+'/'+str(specificity_panel_len),';'.join(tag), ';'.join(mappability_list),mappability_tag))

		if panel_count==len(tag):
			#if mappability_tag=='':#only if not in associated, then load mappability here
			mappability_tag, mappability_list=mappability_write(ls[0], bw_map, calc_length) #Modifed 2021
			fout_prioritized.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ls[0],ls[1],ls[2],ls[3],primary_deltapsi, primary_foc, str(association_passed)+'/'+str(association_panel_len),str(recurrence_passed)+'/'+str(recurrence_panel_len),str(tissue_specificity)+'/'+str(specificity_panel_len),';'.join(tag),';'.join(mappability_list),mappability_tag))

		mappability_tag, mappability_list=['',''] #clear
			
		# elif individual_test:
		# 	test={}
		# 	query_mean=[psi[out_prefix][0],'-','-']
		# 	for j,group in enumerate(panel_list[1:]):
		# 		test[group]=['-']*3
		# 		if has[group]:
		# 			test[group]=one2N(psi[out_prefix],psi[group],screening_type_list[j])

		# 	result=[k]+query_mean+['\t'.join(map(str,test[t])) for t in panel_list if t!=out_prefix]
		# 	fout_direct.write('\t'.join(map(str,result))+'\n')

		# 	pval=[test[t][0] for t in panel_list if t!=out_prefix]
		# 	deltaPSI=[test[t][1] for t in panel_list if t!=out_prefix]# select deltapsi col of screening result
		# 	foc=[test[t][2] for t in panel_list if t!=out_prefix]
			
		# 	differential,equal,positive,negative,testable,primary_result, primary_result_foc, deltapsi_list,foc_list=summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, association_panel_len, recurrence_panel_len, pval, deltaPSI, foc, screening_type_list)
			
		# 	primary_deltapsi, primary_foc, tissue_specificity, tag = defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, association_panel_len, specificity_panel_len, differential,equal,positive,negative,testable,primary_result, primary_result_foc, deltapsi_list,foc_list,use_ratio)
		# if tag!=[]:
		# 	if tag[0]=='associated':
		# 		mappability_tag, mappability_list=mappability_write(ls[0], bw_map, calc_length)
		# 		fout_primary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ls[0],ls[1],ls[2],ls[3],primary_deltapsi, primary_foc, str(differential)+'/'+str(association_panel_len),str(equal)+'/'+str(recurrence_panel_len),str(tissue_specificity)+'/'+str(specificity_panel_len),';'.join(tag), ';'.join(mappability_list),mappability_tag))

		# if panel_count==len(tag):
		# 	if mappability_tag=='':#only if not in associated, then load mappability here
		# 		mappability_tag, mappability_list=mappability_write(ls[0], bw_map, calc_length)
		# 	fout_prioritized.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ls[0],ls[1],ls[2],ls[3],primary_deltapsi, primary_foc, str(differential)+'/'+str(association_panel_len),str(equal)+'/'+str(recurrence_panel_len),str(tissue_specificity)+'/'+str(specificity_panel_len),';'.join(tag),';'.join(mappability_list),mappability_tag))

		# mappability_tag, mappability_list=['',''] #clear
	fout_primary.close()
	fout_prioritized.close()

	##Translation
	if args.translating:

		translationCMD(ref_genome, gtf, outdir, out_prefix, splicing_event_type, all_orf, ignore_annotation, remove_early_stop, 'tier1')
		translationCMD(ref_genome, gtf, outdir, out_prefix, splicing_event_type, all_orf, ignore_annotation, remove_early_stop, 'tier2tier3')

if __name__ == '__main__':
	main()
