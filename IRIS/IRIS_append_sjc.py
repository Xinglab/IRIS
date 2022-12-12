import sys, glob, os
from . import config
#append annotation to all screening result and epitope result in a specified screening folder
def loadPTsummary(fin):
	PT_event={}
	for l in open(fin):
		ls=l.strip().split('\t')
		PT_event[ls[0]]=ls[1]+'|'+ls[2]
	return PT_event

def loadIJCsummary(fin):
	IJC_info={}
	for i,l in enumerate(open(fin)):
		if i==0:
			continue
		ls=l.strip().split('\t')
		IJC_info[ls[0]]=ls[-4]+'\t'+ls[-3]+'\t'+ls[-2]+'\t'+ls[-1]
	return IJC_info

def annotateRAbyPTEvent(fin, event_col, PT_event, fout_fname, IJC_info):
	fout=open(fout_fname,'w')
	for n,l in enumerate(open(fin)):
		if n==0:
			annotation='sjc_test_summary'
			if IJC_info!={}:
				annotation+='\t'+'\t'.join(['ijc_ratio', 'mean_ijc_by_group', 'percent_sample_imbalanced', 'sample_imbalanced_by_group'])
			fout.write(l.strip()+'\t'+annotation+'\n')
			continue
		ls=l.strip().split('\t')
		annotation='-'
		if PT_event.has_key(ls[event_col]):
			annotation=PT_event[ls[event_col]]
		if IJC_info!={}:
			if IJC_info.has_key(ls[event_col]):
				annotation+='\t'+IJC_info[ls[event_col]]
			else:
				annotation+='\t'+'\t'.join(['-']*4)
		fout.write(l.strip()+'\t'+annotation+'\n')
	fout.close()

def main(args):
	PTsummary=args.sjc_summary
	PT_event=loadPTsummary(PTsummary)
	print('[INFO] Number of event in SJ count screen summary:'+str(len(PT_event)))
	splicing_event_type=args.splicing_event_type
	screening_result_dir=args.outdir.rstrip('/')

	epitope_file_junction_file_list= glob.glob(screening_result_dir+'/'+splicing_event_type+'.*/epitope_summary.junction-based.txt')
	epitope_file_peptide_file_list= glob.glob(screening_result_dir+'/'+splicing_event_type+'.*/epitope_summary.peptide-based.txt')
	extracellular_as_file_list= glob.glob(screening_result_dir+'/*'+splicing_event_type+'.*.ExtraCellularAS.txt')
	screening_file_list=glob.glob(screening_result_dir+'/*.'+splicing_event_type+'.tier1.txt')+glob.glob(screening_result_dir+'/*.'+splicing_event_type+'.tier2tier3.txt')

	#optional, if ijc imbalance
	add_ijc_info=args.add_ijc_info
	use_existing_result=args.use_existing_result#use existing ijc result
	para_file=args.parameter_file
	event_list_file=args.screening_result_event_list
	inc_read_cov_cutoff=args.inc_read_cov_cutoff#2
	event_read_cov_cutoff=args.event_read_cov_cutoff#10
	IJC_info={}
	if add_ijc_info or use_existing_result:
		if para_file=='' or event_list_file=='':
			exit('[Error] Specify parameters and event list file for retrieving inclusion junction information.')
		
	if add_ijc_info and use_existing_result==False:
		cmd='IRIS annotate_ijc -p '+para_file+' --splicing-event-type '+splicing_event_type+' -e '+event_list_file+' -o '+screening_result_dir+' --inc-read-cov-cutoff '+str(inc_read_cov_cutoff)+' --event-read-cov-cutoff '+str(event_read_cov_cutoff)
		print('[INFO] Annotating inclusion junction info to '+event_list_file)
		print(cmd)
		os.system(cmd)
	if use_existing_result or add_ijc_info:
		file_path=screening_result_dir+'/'+event_list_file.split('/')[-1]+'.ijc_info.txt'
		if os.path.exists(file_path)==False:
			exit('[Error] Result IJC annotation file not found in path '+file_path)
		print('[INFO] Loading inclusion junction info from '+file_path)
		IJC_info=loadIJCsummary(file_path)

	for fin_fname in screening_file_list:
		print '[INFO] Integrating SJC test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 0, PT_event, fin_fname+'.integratedSJC.txt',IJC_info)

	for fin_fname in extracellular_as_file_list:
		print '[INFO] Integrating SJC test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 0, PT_event, fin_fname+'.integratedSJC.txt',IJC_info)

	for fin_fname in epitope_file_peptide_file_list:
		print '[INFO] Integrating SJC test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 1, PT_event, fin_fname+'.integratedSJC.txt',IJC_info)

	for fin_fname in epitope_file_junction_file_list:
		print '[INFO] Integrating SJC test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 0, PT_event, fin_fname+'.integratedSJC.txt',IJC_info)


if __name__ == '__main__':
	main()
