
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

def openTestingFout(outdir, out_prefix, summary_file, ref_list, test_mode):
	header=['as_event','meanPSI','Q1PSI','Q3PSI']
	if test_mode=='group':
		header_prefix=['_pVal','_deltaPSI','_tumorFC']
	if test_mode=='personalized':
		header_prefix=['_modifiedPctl','_deltaPSI','_tumorFC']
	fout=open(outdir+'/'+out_prefix+'.test.all.txt','w')
	if summary_file==False:
		header+=['\t'.join(map(lambda x:ref+x ,header_prefix)) for ref in ref_list if ref!=out_prefix]
	fout.write('\t'.join(header)+'\n')
	return fout

def openScreeningFout(outdir, out_prefix, fout_name):
	fout=open(outdir+'/'+out_prefix+'.'+fout_name+'.txt','w')
	fout.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('as_event','meanPSI','Q1PSI','Q3PSI','deltaPSI','fc_of_tumor_isoform','tissue_matched_normal_panel','tumor_panel','normal_panel','tag','mappability','mappability_tag'))
	return fout

def writeSummaryFile(out_prefix, db_dir, index, fout): #no screening, preparing for MS search
	fin_name=db_dir+'/'+out_prefix+'/splicing_matrix/splicing_matrix.SE.cov10.'+group_name+'.txt'
	index[out_prefix]=read_PsiMatrix_index(fin_name,'/'.join(fin_name.split('/')[:-1]))
	for j,k in enumerate(index[out_prefix]):
		psi_event=map(float,fetch_PsiMatrix(k,fin_name,'.','\t',index[out_prefix])[1][8:])
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

#def groupTest(g1,g2, test_type, threshold_tost=0.05):
def groupTest(g1,g2, test_type, direction='two-sided'):
	g1=np.array(g1)
	g2=np.array(g2)
	g1=g1[~np.isnan(g1)]
	g2=g2[~np.isnan(g2)]
	delta_psi=np.nanmean(g1)-np.nanmean(g2)
	tumor_foc=calcTumorFormFoc(delta_psi,np.nanmean(g1))
	if test_type=='sig':
		t1=stats.ttest_ind(g1,g2)[1]
	elif test_type=='equ':
	 	t1=smw.ttest_ind(g1,g2,alternative='two-sided')[1]
	# 	t1=smw.ttost_ind(g1,g2,-threshold_tost,threshold_tost,usevar='unequal')[0] #equalvalence test
	return [t1, delta_psi, tumor_foc]

def getDirection(filter1, psi, test):
	psi_primary=[]
	deltaPSI_primary=[]
	for primary_group in filter1:
		if test[primary_group]!=['-']*3:
			psi_primary+=psi[primary_group]
			deltaPSI_primary+=test[primary_group]
	if psi_primary==[]: #in case tissue-matched normal doesn't have data
		return [],[]
	else:
		direction='larger'
		if np.median(deltaPSI_primary)<=0:
			direction='smaller'
		return psi_primary, direction
	
def summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, primary, tumor_rec, pval, deltaPSI, foc, testing_type_index):
	differential,equal,positive,negative=[0,0,0,0]
	testable=0
	primary_result, primary_result_foc=[[],[]] #take care of multiple tiisue matched norm
	deltapsi_list,foc_list=[[],[]] #if no tissue-matched norm, use median
	for i,group_type in enumerate(testing_type_index):
		if pval[i]=='-':
			continue
		if i<primary:
			primary_result.append(float(deltaPSI[i]))
			primary_result_foc.append(float(foc[i]))
			if float(pval[i])<=filter1_cutoff_pval and abs(float(deltaPSI[i]))>=filter1_cutoff_dpsi and float(foc[i])>=filter1_cutoff_foc:
					differential+=1
					continue
		if testing_type_index[i]=='equ':
			if float(pval[i])<=filter2_cutoff_pval and abs(float(deltaPSI[i]))>=filter2_cutoff_dpsi:
					equal+=1
					continue
		if testing_type_index[i]=='sig' and i>=(primary+tumor_rec):
			testable+=1
			if float(pval[i])<=filter3_cutoff_pval and float(foc[i])>=filter3_cutoff_foc:
				deltapsi_list.append(float(deltaPSI[i]))
				foc_list.append(float(foc[i]))
				if float(deltaPSI[i])>=filter3_cutoff_dpsi:
					positive+=1
					continue
				if float(deltaPSI[i])<=-filter3_cutoff_dpsi:
					negative+=1
					continue
	return differential,equal,positive,negative,testable, np.median(primary_result), np.median(primary_result_foc), deltapsi_list,foc_list

def defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, primary, norm_tissue, differential, equal, positive, negative, testable, primary_result, primary_result_foc, deltapsi_list,foc_list, use_ratio):
	tag=[]
 	if differential>=filter1_group_cutoff:
 		tag.append('associated')
 	if equal>=filter2_group_cutoff:
 		tag.append('recurrent') 
	if primary==0:
		tissue_specificity=max(positive,negative)
		ratio=False if use_ratio==False else tissue_specificity/(testable+10**-8)>=filter3_group_cutoff/(norm_tissue+0.0)
		if tissue_specificity>=filter3_group_cutoff or ratio:
			tag.append('specific')
	else: 
		if primary_result>0:
			tissue_specificity=positive
			ratio=False if use_ratio==False else positive/(testable+10**-8)>=filter3_group_cutoff/(norm_tissue+0.0)
			if positive>=filter3_group_cutoff or ratio:
				tag.append('specific')
		else:
			tissue_specificity=negative
			ratio=False if use_ratio==False else negative/(testable+10**-8)>=filter3_group_cutoff/(norm_tissue+0.0)
			if negative>=filter3_group_cutoff or ratio:
				tag.append('specific')
	primary_deltapsi=primary_result if primary_result!=0 else np.median(deltapsi_list)
	primary_foc=primary_result_foc if primary_result!=0 else np.median(foc_list)
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

def translationCMD(ref_genome, outdir, out_prefix, fout_name):
	#uSE name space 
	#Namespace(outdir='test1', parameter_fin='/u/home/p/panyang/bigdata-nobackup/Glioma_test/GBM.prioritze.par', subcommand='screening', translating=False)
	cmd_translation='IRIS translation '+outdir+'/'+out_prefix+'.'+fout_name+'.txt '+' -o '+outdir+'/'+fout_name+' -g '+ref_genome
	print '[INFO] Working on translating: '+fout_name
	os.system(cmd_translation)

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


def main(args):
	###Loading Parameters####
	index={}
	fin_list={}
	para_fin=args.parameter_fin
	out_prefix,db_dir,filter1_para,filter2_para,filter3_para,test_mode,use_ratio,blacklist_path,mappability_path,ref_genome=[l.strip() for l in open(para_fin)]
	ref_list=[out_prefix]
	use_ratio=True if use_ratio=='True' else False
	if blacklist_path=='':
		blacklist_path=config.BRAIN_BLACKLIST_PATH
	blacklist_events=loadBlacklistEvents(blacklist_path)
	bw_map,calc_length=loadMappability(mappability_path)

	filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc, filter1_group_cutoff, filter1, ref_list =loadParametersRow(filter1_para, ref_list)
	filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter2_group_cutoff, filter2, ref_list =loadParametersRow(filter2_para, ref_list)
	filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, filter3_group_cutoff, filter3, ref_list =loadParametersRow(filter3_para, ref_list)
	if filter1==[] and filter2==[] and filter3==[] and test_mode!='summary':
		exit("[Error] No filtering required in parameteres file. exit!")
	group_test=False if test_mode!='group' else True
	individual_test=False if test_mode!='personalized' else True
	summary_file=False if test_mode!='summary' else True
	if [group_test,individual_test,summary_file]==[False,False,False]:
		exit('[Error] Need to choose one mode.exit!')
	primary=len(filter1)
	tumor_rec=len(filter2)
	norm_tissue=len(filter3)
	filter_count=sum(1 for i in [primary, tumor_rec, norm_tissue] if i!=0)
	testing_type_index=['sig']*primary+['equ']*tumor_rec+['sig']*norm_tissue
	
	db_dir=db_dir.rstrip('/')
	outdir=args.outdir.rstrip('/')
	os.system('mkdir -p  '+outdir)

	###Create Folders/Output####
	fout= openTestingFout(outdir, out_prefix, summary_file, ref_list, test_mode)
	if summary_file:
		writeSummaryFile(out_prefix, db_dir, index, fout)
		exit()
	
	fout_filtered=open(outdir+'/'+out_prefix+'.notest.txt','w')
	fout_primary=openScreeningFout(outdir, out_prefix, 'primary')
	fout_prioritized=openScreeningFout(outdir,out_prefix, 'prioritized')

	for group_name in ref_list:##Load IRIS reference panels
		fin_list[group_name]=db_dir+'/'+group_name+'/splicing_matrix/splicing_matrix.SE.cov10.'+group_name+'.txt'

	for group in fin_list.keys():
		if not os.path.isfile(fin_list[group]+'.idx'):
			exit('[Error] Need to index '+fin_list[group])
		index[group]=read_PsiMatrix_index(fin_list[group],'/'.join(fin_list[group].split('/')[:-1]))

	has={}
	#for j,k in enumerate(['ENSG00000083520:DIS3:chr13:-:73345041:73345126:73343050:73345218','ENSG00000110075:PPP6R3:chr11:+:68350510:68350597:68343511:68355265']):	
	tot=len(index[out_prefix])-1
	print '[INFO] IRIS screening started. Total input events:', tot+1
	for event_idx,k in enumerate(index[out_prefix]):
		config.update_progress(event_idx/(0.0+tot))
		
		for group in ref_list:#Initiate
			if group!=out_prefix:
				has[group]=True
		psi={}
		has_count=0
		for group in ref_list:
			if k in index[group]:
				psi[group]=map(float,fetch_PsiMatrix(k,fin_list[group],'.','\t',index[group])[1][8:])
				has_count+=1
			else:
				has[group]=False

		cat_psi=[]
		for i in psi:
			cat_psi+=psi[i]
		if abs(max(cat_psi)-min(cat_psi))<0.05:#if change less than 5% skipped and no comparison available
			fout_filtered.write('[LowVar]{}\t{}\t{}\n'.format(k,str(abs(max(cat_psi)-min(cat_psi))),str(has_count))) 
			continue
		if k in blacklist_events:
			fout_filtered.write('[Blacklisted]{}\t{}\t{}\n'.format(k,'-',str(has_count))) 
			continue
		if  has_count<=1:
			fout_filtered.write('[NoTest]{}\t{}\t{}\n'.format(k,'-',str(has_count)))
			continue
			
		if group_test:
			test={}
			query_mean=[np.nanmean(psi[out_prefix]),np.nanpercentile(psi[out_prefix],25),np.nanpercentile(psi[out_prefix],75)]
			for j,group in enumerate(ref_list[1:]):
				test[group]=['-']*3
				if has[group]:
					if testing_type_index[j]=='equ':		
						psi_primary, direction= getDirection(filter1, psi, test)
						if psi_primary==[]: #in case tissue-matched normal doesn't have data
							continue 
						test[group]=groupTest(psi[group],psi_primary, testing_type_index[j], direction)
					else:
						test[group]=groupTest(psi[out_prefix],psi[group],testing_type_index[j])

			result=[k]+query_mean+['\t'.join(map(str,test[t])) for t in ref_list if t!=out_prefix]
			fout.write('\t'.join(map(str,result))+'\n')

			## summarize test result and prioritze
			pval=[test[t][0] for t in ref_list if t!=out_prefix]
			deltaPSI=[test[t][1] for t in ref_list if t!=out_prefix]# select deltapsi col of screening result
			foc=[test[t][2] for t in ref_list if t!=out_prefix]

			differential,equal,positive,negative,testable,primary_result, primary_result_foc, deltapsi_list,foc_list=summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, primary, tumor_rec, pval, deltaPSI, foc, testing_type_index)
			
			primary_deltapsi, primary_foc, tissue_specificity, tag = defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, primary, norm_tissue, differential,equal,positive,negative,testable,primary_result, primary_result_foc, deltapsi_list,foc_list,use_ratio)

			if tag!=[]:
				if tag[0]=='associated':
					mappability_tag, mappability_list=mappability_write(k, bw_map, calc_length)
					fout_primary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(k,query_mean[0],query_mean[1],query_mean[2],primary_deltapsi, primary_foc, str(differential)+'/'+str(primary),str(equal)+'/'+str(tumor_rec),str(tissue_specificity)+'/'+str(norm_tissue),';'.join(tag), ';'.join(mappability_list),mappability_tag))

			if filter_count==len(tag):
				if mappability_tag=='':#only if not in associated, then load mappability here
					mappability_tag, mappability_list=mappability_write(k, bw_map, calc_length)
				fout_prioritized.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(k,query_mean[0],query_mean[1],query_mean[2],primary_deltapsi, primary_foc, str(differential)+'/'+str(primary),str(equal)+'/'+str(tumor_rec),str(tissue_specificity)+'/'+str(norm_tissue),';'.join(tag),';'.join(mappability_list),mappability_tag))

			mappability_tag, mappability_list=['',''] #clear
			
		elif individual_test:
			test={}
			query_mean=[psi[out_prefix][0],'-','-']
			for j,group in enumerate(ref_list[1:]):
				test[group]=['-']*3
				if has[group]:
					test[group]=one2N(psi[out_prefix],psi[group],testing_type_index[j])

			result=[k]+query_mean+['\t'.join(map(str,test[t])) for t in ref_list if t!=out_prefix]
			fout.write('\t'.join(map(str,result))+'\n')

			pval=[test[t][0] for t in ref_list if t!=out_prefix]
			deltaPSI=[test[t][1] for t in ref_list if t!=out_prefix]# select deltapsi col of screening result
			foc=[test[t][2] for t in ref_list if t!=out_prefix]
			
			differential,equal,positive,negative,testable,primary_result, primary_result_foc, deltapsi_list,foc_list=summarizeTestResult(filter1_cutoff_pval, filter1_cutoff_dpsi, filter1_cutoff_foc,filter2_cutoff_pval, filter2_cutoff_dpsi, filter2_cutoff_foc, filter3_cutoff_pval, filter3_cutoff_dpsi, filter3_cutoff_foc, primary, tumor_rec, pval, deltaPSI, foc, testing_type_index)
			
			primary_deltapsi, primary_foc, tissue_specificity, tag = defineTumorEvents(filter1_group_cutoff,filter2_group_cutoff,filter3_group_cutoff, primary, norm_tissue, differential,equal,positive,negative,testable,primary_result, primary_result_foc, deltapsi_list,foc_list,use_ratio)

			if tag!=[]:
				if tag[0]=='associated':
					mappability_tag, mappability_list=mappability_write(k, bw_map, calc_length)
					fout_primary.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(k,query_mean[0],query_mean[1],query_mean[2],primary_deltapsi, primary_foc, str(differential)+'/'+str(primary),str(equal)+'/'+str(tumor_rec),str(tissue_specificity)+'/'+str(norm_tissue),';'.join(tag), ';'.join(mappability_list),mappability_tag))

			if filter_count==len(tag):
				if mappability_tag=='':#only if not in associated, then load mappability here
					mappability_tag, mappability_list=mappability_write(k, bw_map, calc_length)
				fout_prioritized.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(k,query_mean[0],query_mean[1],query_mean[2],primary_deltapsi, primary_foc, str(differential)+'/'+str(primary),str(equal)+'/'+str(tumor_rec),str(tissue_specificity)+'/'+str(norm_tissue),';'.join(tag),';'.join(mappability_list),mappability_tag))

			mappability_tag, mappability_list=['',''] #clear
		else:
			exit('[Error] no test performed')
	fout.close()
	fout_filtered.close()
	fout_primary.close()
	fout_prioritized.close()

	if args.translating:
		translationCMD(ref_genome, outdir, out_prefix, 'primary')
		translationCMD(ref_genome, outdir, out_prefix, 'prioritized')

if __name__ == '__main__':
	main()