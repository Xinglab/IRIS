import numpy as np
import sys
import os, glob, pyBigWig, argparse
from scipy import stats
import statsmodels.stats.weightstats as smw
from . import config
import warnings
warnings.filterwarnings("ignore")

def read_PsiMatrix_index(fn,outdir):
	index = {}
	for line in open(outdir+'/'+fn.split('/')[-1]+'.idx', 'r'):
		ele = line.strip().split()
		index[ele[0]] = int(ele[1])
	return index

def fetch_PsiMatrix(eid, fn, delim, index=None):
	with open(fn, 'r') as f:
		ele = f.readline().strip().split(delim)
		header = np.asarray([ x.split('.aln')[0] for x in ele ])
		f.seek(index[eid], 0)
		data = np.asarray(f.readline().strip().split(delim))
	return (header, data)

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
	if filter_para.strip()!='':
		para, filter_panel_list=filter_para.split(' ')
		filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_cutoff_pval_PT, filter_group_cutoff =para.split(',')
		filter_cutoff_pval=float(filter_cutoff_pval)
		filter_cutoff_dpsi=float(filter_cutoff_dpsi)
		filter_cutoff_foc=float(filter_cutoff_foc)
		filter_group_cutoff=int(filter_group_cutoff)
		filter_panel_list=filter_panel_list.split(',')		
		panel_list+=filter_panel_list
	else:
		filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_panel_list =['','','','',[]]
	return filter_cutoff_pval, filter_cutoff_dpsi, filter_cutoff_foc, filter_group_cutoff, filter_panel_list, panel_list

def openTestingFout(outdir, out_prefix, splicing_event_type, summary_file, panel_list, fin_name=''):
	header=['as_event','meanPSI','Q1PSI','Q3PSI']
	header_prefix=['_pVal','_deltaPSI','_tumorFC']
	fout_name=outdir+'/'+out_prefix+'.'+splicing_event_type+'.test_JCRA.all_'+fin_name+'.txt'
	fout=open(fout_name,'w')
	if summary_file==False:
		header+=['\t'.join(map(lambda x:ref+x ,header_prefix)) for ref in panel_list if ref!=out_prefix]
	fout.write('\t'.join(header)+'\n')
	return fout, fout_name

def openScreeningFout(outdir, out_prefix, splicing_event_type, fout_name):
	fout=open(outdir+'/'+out_prefix+'.'+splicing_event_type+'.'+fout_name+'.txt','w')
	fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('as_event','meanPSI','Q1PSI','Q3PSI','deltaPSI','fc_of_tumor_isoform','tissue_matched_normal_panel','tumor_panel','normal_panel','tag','mappability','mappability_tag','novel_ss_info'))
	return fout

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

def loadGTF(gtf):
	exon_start_dict={}
	exon_end_dict={}
	for l in open(gtf):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		if ls[2]=='exon':
			chrom=ls[0]
			if chrom.startswith('chr')==False:
				chrom='chr'+chrom
			exon_start_dict[ls[6]+':'+chrom+':'+ls[3]]=''
			exon_end_dict[ls[6]+':'+chrom+':'+ls[4]]=''
	return exon_start_dict, exon_end_dict

def selectJunction_forGTF(AS_coord, deltaPSI_c2n, cut_off, if_select_all, splicing_event_type): #For A3 and A5, inc1/2 are not considered as current sj db doesn't capture that(and it will not always generate novel sequence). May update later.
	if splicing_event_type == 'SE':
		# end_pos, start_pos; because rMATS is 0-based for start exon, it should +1 when compare to GTF
		skp = (AS_coord[2], str(int(AS_coord[3])+1),'skp') 
		inc1 = (AS_coord[2], str(int(AS_coord[0])+1),'inc1')
		inc2 = (AS_coord[1], str(int(AS_coord[3])+1),'inc2')
	elif splicing_event_type == 'A3SS':
		skp = (AS_coord[5], str(int(AS_coord[2])+1),'skp')
		inc1 = (AS_coord[5],str(int(AS_coord[0])+1),'inc1')
		inc2 = ''#(AS_coord[2],AS_coord[2],'inc2')
	elif splicing_event_type == 'A5SS':
		skp = (AS_coord[3], str(int(AS_coord[4])+1),'skp')
		inc1 = ''#(AS_coord[3],AS_coord[3],'inc1')
		inc2 = (AS_coord[1],str(int(AS_coord[4])+1),'inc2')
	elif splicing_event_type == 'RI':
		skp = ''#(AS_coord[3], AS_coord[4],'skp')#????
		inc1 = ''#(AS_coord[3],AS_coord[3],'inc1')
		inc2 = ''#(AS_coord[4],AS_coord[4],'inc2')
	if if_select_all:
		return [skp, inc1, inc2]
	if float(deltaPSI_c2n) < cut_off:  # tumor skipping
		return [skp]
	else:
		return [inc1, inc2]

def findNovelSS(event_name, deltaPSI_c2n, cut_off, splicing_event_type, exon_start_dict, exon_end_dict):# This is a conservative def of novelSS than rMATS4.1 (0.4% events less- complex cases )
	es=event_name.split(':')
	chrom=es[2]
	strand=es[3]
	AS_coord=es[4:]
	selected=selectJunction_forGTF(AS_coord,deltaPSI_c2n, cut_off, False, splicing_event_type)
	info=[]
	for j in selected:
		if j!='':
			check1=strand+':'+chrom+':'+j[0] not in exon_end_dict
			check2=strand+':'+chrom+':'+j[1] not in exon_start_dict
			if check1 or check2:
				info.append(j[2])##
	if info==[]:
		info.append('none')
	return ';'.join(info)

def readEventRow(row, header_line):
	if header_line=='' or header_line==False:
		rs=row.strip().split('\t')
		return rs
	else:
		rs=row.strip().split('\t')
		return dict(zip(header_line, rs))

def convert2SJevent(line_dict, splicing_event_type):#match to SJ db (different from gtf and rMATS - 0-based start; +1 for end exon)
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

def getDirection(panel_dict, psi, test, non_parametric):
	psi_primary=[]
	deltaPSI_primary=[]
	for primary_group in panel_dict: #check if matching normal
		if test[primary_group]!=['-']*3:
			psi_primary+=psi[primary_group]
			deltaPSI_primary.append(test[primary_group][1])
	if psi_primary==[]: #in case tissue-matched normal doesn't have data
		return [],['']
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

def performTest(group, matching_norm_dict, tumor_dict, normal_dict, psi, out_prefix, non_parametric, psi_primary, direction):
	redirect_output = False
	test_result = ['-']*3 #For missing in non-eesential tests/comparisons
	has_matched_tumor = False if psi_primary == [] else True #for clarity
	if group in matching_norm_dict:
		if psi[group]!=[]:
			test_result = groupTest(psi[out_prefix],psi[group], non_parametric,"two-sided")
		return test_result

	elif group in tumor_dict:
		if has_matched_tumor:#set_matched_tumor is redundent here. kept for future implemtation of additional output type.
			if psi[group]!=[]:
				test_result = groupTest(psi[group],psi_primary, non_parametric, direction)
		else:
			if psi[group]!=[]:
				test_result = groupTest(psi[out_prefix],psi[group], non_parametric, "equivalence")#No or equivalent testing 
		return test_result

	elif group in normal_dict:
		if has_matched_tumor:
			if psi[group]!=[]:
				test_result = groupTest(psi[out_prefix],psi[group], non_parametric, direction)
		else:
			if psi[group]!=[]:
				test_result = groupTest(psi[out_prefix],psi[group], non_parametric, "two-sided") #Two-sided testing 
		return test_result
	else:
		exit('error in group.')


def summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, pval, deltaPSI, foc, panel_list, matching_norm_dict, tumor_dict, normal_dict):
	association_passed,recurrence_passed,specificity_positive,specificity_negative, specificity_testable=[0,0,0,0,0]
	primary_result, primary_result_foc=[[],[]] #take care of multiple tiisue matched norm
	deltapsi_list_voted,foc_list_voted=[[],[]] #if no tissue-matched norm, use median
	for i,group in enumerate(panel_list[1:]):
		if pval[i]=='-': #This is important - skip all missing, which is not useful for summarizing but not affecting consistancy.
			continue
		if group in matching_norm_dict: 
			primary_result.append(float(deltaPSI[i])) 
			primary_result_foc.append(float(foc[i]))
			if float(pval[i])<=filter1_cutoff_pval and abs(float(deltaPSI[i]))>=filter1_cutoff_dpsi and float(foc[i])>=filter1_cutoff_foc:
					association_passed+=1
					continue
		elif group in tumor_dict:
			if float(pval[i])<=filter2_cutoff_pval and abs(float(deltaPSI[i]))>=filter2_cutoff_dpsi:
					recurrence_passed+=1
					continue
		elif group in normal_dict:# TODO: judge set/has, then run
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

def defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, matching_norm_dict, specificity_panel_len, association_passed, recurrence_passed, specificity_positive, specificity_negative, specificity_testable, primary_result, primary_result_foc, deltapsi_list_voted, foc_list_voted, use_ratio):
	tag=[]
	if filter1_group_cutoff!='':
 		if association_passed>=filter1_group_cutoff:#Improvement? current: 0>='' is false
 			tag.append('associated')
 	if filter2_group_cutoff!='':
 		if recurrence_passed>=filter2_group_cutoff:
 			tag.append('recurrent') 
	if matching_norm_dict!={}:#TODO-FUTURE: take care of set-yes has-no redirected events
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
	primary_deltapsi=primary_result if primary_result!=0 else np.median(deltapsi_list_voted)#TODO
	primary_foc=primary_result_foc if primary_result!=0 else np.median(foc_list_voted)
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

def translationCMD(ref_genome, gtf, outdir, out_prefix, splicing_event_type, all_orf, ignore_annotation, remove_early_stop, fout_name, find_novel_ss):
	#uSE name space 
	argument_line=''
	if all_orf:
		argument_line+=' --all-orf'
	if ignore_annotation:
		argument_line+=' --ignore-annotation'
	if remove_early_stop:
		argument_line+=' --remove-early-stop'
	if find_novel_ss:
		argument_line+=' --check-novel'
	cmd_translation='IRIS translate '+outdir+'/'+out_prefix+'.'+splicing_event_type+'.'+fout_name+'.txt '+' -o '+outdir+'/'+splicing_event_type+'.'+fout_name+' -g '+ref_genome+' -t '+splicing_event_type+argument_line+' --gtf '+gtf
	print '[INFO] Working on translating: '+fout_name
	os.system(cmd_translation)

def main(args):
	###Loading Parameters####
	para_fin=args.parameter_fin
	splicing_event_type=args.splicing_event_type
	event_list_fin=args.event_list_fin
	fetching_sj_col=1
	out_prefix,db_dir,filter1_para,filter2_para,filter3_para,test_mode,use_ratio,blacklist_path,mappability_path,ref_genome=[l.strip() for l in open(para_fin)]
	panel_list=[out_prefix]
	test_mode=test_mode.split(' ')
	use_ratio=True if use_ratio=='True' else False
	blacklist_events={}
        if blacklist_path!='':
	    blacklist_events=loadBlacklistEvents(blacklist_path)
	bw_map,calc_length=loadMappability(mappability_path)
	deltaPSI_cutoff=args.deltaPSI_cut_off
	find_novel_ss=True if args.report_known_and_novelss_tumor_junction==False else False
	if find_novel_ss:
		gtf=args.gtf
		exon_start_dict, exon_end_dict= loadGTF(gtf)

	all_orf=args.all_orf 
	ignore_annotation=args.ignore_annotation 
	remove_early_stop=args.remove_early_stop 
	use_existing_test_result=args.use_existing_test_result

	filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc, filter1_group_cutoff, filter1_panel_list, panel_list = loadParametersRow(filter1_para, panel_list)
	filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter2_group_cutoff, filter2_panel_list, panel_list = loadParametersRow(filter2_para, panel_list)
	filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, filter3_group_cutoff, filter3_panel_list, panel_list = loadParametersRow(filter3_para, panel_list)
	if filter1_panel_list==[] and filter2_panel_list==[] and filter3_panel_list==[]:
	# if filter1_panel_list==[] and filter2_panel_list==[] and filter3_panel_list==[] and test_mode[0]!='summary':
		exit("[Error] No filtering required in parameteres file. exit!")
	summary_file=False

	non_parametric=False
	if len(test_mode)>1:
		non_parametric=True if test_mode[1]=='nonparametric' else False
	
	##this block is 2020 code; different and improved from screening.py
	matching_norm_dict=dict.fromkeys(filter1_panel_list,'')
	tumor_dict=dict.fromkeys(filter2_panel_list,'')
	normal_dict=dict.fromkeys(filter3_panel_list,'')
	tumor_dict[out_prefix]=''
	association_panel_len=len(matching_norm_dict)
	recurrence_panel_len=len(tumor_dict)-1
	specificity_panel_len=len(normal_dict)

	panel_count=sum(1 for i in [association_panel_len, recurrence_panel_len, specificity_panel_len] if i!=0)

        if args.translating:
            gtf=args.gtf
            if os.path.exists(gtf)==False:
                exit('[Error] No gtf file provided for translation. exit!')
	###Create Folders/Output####
	outdir=args.outdir.rstrip('/')
	os.system('mkdir -p '+outdir)
	db_dir=db_dir.rstrip('/')

	if use_existing_test_result==False:
		fout_direct, fout_direct_name= openTestingFout(outdir, out_prefix, splicing_event_type, summary_file, panel_list, 'guided')
		fout_redirect, fout_redirect_name=openTestingFout(outdir, out_prefix, splicing_event_type, summary_file, panel_list, 'voted')
	else:
		fout_direct_name=outdir+'/'+out_prefix+'.'+splicing_event_type+'.test_JCRA.all_guided.txt'
		fout_redirect_name=outdir+'/'+out_prefix+'.'+splicing_event_type+'.test_JCRA.all_voted.txt'

	fout_primary=openScreeningFout(outdir, out_prefix, splicing_event_type, 'tier1_JCRA')
	fout_prioritized=openScreeningFout(outdir,out_prefix, splicing_event_type, 'tier2tier3_JCRA')
	
	##Load IRIS reference panels to 'fin_list', 'index'
	index={}
	fin_list={}
	for group_name in panel_list:
		fin_list[group_name]=db_dir+'/'+group_name.split('_')[0]+'/SJ_count.'+group_name+'.txt'
	for group in fin_list.keys():
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_SJMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))

	tot=config.file_len(event_list_fin)-1
	print '[INFO] IRIS screen_novelss - started. Total input events:', tot+1
	## Load and perform test by row/event
	sample_size={}
	for group in panel_list:
		random_key=index[group].keys()[0]
		sample_names=map(str,fetch_SJMatrix(random_key,fin_list[group],'\t',index[group],True)[fetching_sj_col:])
		sample_size[group]=len(sample_names)

	if use_existing_test_result==False:
		header_list=[]
		junction_dict={}
		for event_idx, event_row in enumerate(open(event_list_fin)):
			if event_idx==0:
				header_list=readEventRow(event_row,'')
				continue
			line_dict=readEventRow(event_row, header_list)
			event_row_list, as_event=convert2SJevent(line_dict, splicing_event_type)
			config.update_progress(event_idx/(0.0+tot))		
			
			sj={}
			for eid, k in enumerate(event_row_list):
				sj[eid]={}
				# if k not in junction_dict:
				# 	junction_dict[k]=''
				# else:
				# 	continue

				#Initiate psi matrix by each row to 'sj' 
				for group in panel_list:
					if k in index[group]:
						sj[eid][group]=map(int,fetch_SJMatrix(k,fin_list[group],'\t',index[group], False)[fetching_sj_col:])
					else:
						sj[eid][group]=[0]*sample_size[group]
			psi={}
			for group in panel_list:
				psi[group]=[]
				for sid in xrange(0,sample_size[group]):
					inc1=sj[0][group][sid]
					inc2=sj[1][group][sid]
					skp=sj[2][group][sid]
					if inc1+ inc2+ skp>=10:
						psi[group].append(((inc1+inc2)/2.0)/(((inc1+inc2)/2.0)+skp))

			test={}
			redirect=False
			psi_primary=''
			direction=''
			query_mean=map(lambda x:round(x,2),[np.nanmean(psi[out_prefix]),np.nanpercentile(psi[out_prefix],25),np.nanpercentile(psi[out_prefix],75)])
			for group in panel_list[1:]:
				if group not in matching_norm_dict and direction=='':#redirect or find one-sided test direction
					psi_primary, direction= getDirection(matching_norm_dict, psi, test, non_parametric)
					redirect = True if psi_primary == [] else False #redirect if tumor-matched normal is missing
				test[group]=performTest(group, matching_norm_dict, tumor_dict, normal_dict, psi, out_prefix, non_parametric, psi_primary, direction)
			result=[as_event]+query_mean+['\t'.join(map(str,test[t])) for t in panel_list if t!=out_prefix]
			
			if redirect:
				fout_redirect.write('\t'.join(map(str,result))+'\n') #to redicted; summarize direction and calculate FDR differently.
			else:
				fout_direct.write('\t'.join(map(str,result))+'\n')

		fout_redirect.close()
		fout_direct.close()

	testing_intermediate_file = fout_direct_name if matching_norm_dict!={} else fout_redirect_name
	
	tot=config.file_len(testing_intermediate_file)-1
	print '[INFO] IRIS screen_novelss - summarizing. Total events from last step:', tot
	for event_idx,l in enumerate(open(testing_intermediate_file)):
		config.update_progress(event_idx/(0.0+tot))
		if event_idx==0:
			continue
		ls=l.strip().split('\t')
		pval= ls[4::3]# don't do map because '-'
		deltaPSI = ls[5::3]
		foc= ls[6::3]

		association_passed,recurrence_passed,specificity_positive,specificity_negative,specificity_testable,primary_result, primary_result_foc, deltapsi_list_voted,foc_list_voted=summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, pval, deltaPSI, foc, panel_list, matching_norm_dict, tumor_dict, normal_dict)
		
		primary_deltapsi, primary_foc, tissue_specificity, tag = defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, matching_norm_dict, specificity_panel_len, association_passed, recurrence_passed, specificity_positive, specificity_negative, specificity_testable, primary_result, primary_result_foc, deltapsi_list_voted, foc_list_voted, use_ratio)
		novel_ss_status='disabled'
		if find_novel_ss:
			novel_ss_status=findNovelSS(ls[0], primary_deltapsi, deltaPSI_cutoff, splicing_event_type, exon_start_dict, exon_end_dict)
		if novel_ss_status=='none':
			continue
		if tag!=[]:
			if tag[0]=='associated':
				mappability_tag, mappability_list=mappability_write(ls[0], bw_map, calc_length)
				fout_primary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ls[0],ls[1],ls[2],ls[3],primary_deltapsi, primary_foc, str(association_passed)+'/'+str(association_panel_len),str(recurrence_passed)+'/'+str(recurrence_panel_len),str(tissue_specificity)+'/'+str(specificity_panel_len),';'.join(tag), ';'.join(mappability_list),mappability_tag, novel_ss_status))
		if panel_count==len(tag):#Modified 2021
			mappability_tag, mappability_list=mappability_write(ls[0], bw_map, calc_length) #Modified 2021
			fout_prioritized.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ls[0],ls[1],ls[2],ls[3],primary_deltapsi, primary_foc, str(association_passed)+'/'+str(association_panel_len),str(recurrence_passed)+'/'+str(recurrence_panel_len),str(tissue_specificity)+'/'+str(specificity_panel_len),';'.join(tag),';'.join(mappability_list),mappability_tag, novel_ss_status))
		mappability_tag, mappability_list=['',''] #clear

	fout_primary.close()
	fout_prioritized.close()

	##Translation
	if args.translating:
		translationCMD(ref_genome, gtf, outdir, out_prefix, splicing_event_type, all_orf, ignore_annotation, remove_early_stop, 'tier1_JCRA', find_novel_ss)
		translationCMD(ref_genome, gtf, outdir, out_prefix, splicing_event_type, all_orf, ignore_annotation, remove_early_stop, 'tier2tier3_JCRA', find_novel_ss)


if __name__ == '__main__':
	main()
