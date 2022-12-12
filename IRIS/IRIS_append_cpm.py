import sys, glob, os
from . import config
#append annotation to all screening result and epitope result in a specified screening folder
def loadPTsummary(fin):
	PT_event={}
	for l in open(fin):
		ls=l.strip().split('\t')
		PT_event[ls[0]]=ls[1]+'|'+ls[2]+'|'+ls[3]+'|'+ls[4]
	return PT_event


def annotateRAbyPTEvent(fin, event_col, CPM_event, fout_fname):
	fout=open(fout_fname,'w')
	for n,l in enumerate(open(fin)):
		if n==0:
			annotation='cpm_test_summary'
			fout.write(l.strip()+'\t'+annotation+'\n')
			continue
		ls=l.strip().split('\t')
		annotation='-'
		if CPM_event.has_key(ls[event_col]):
			annotation=CPM_event[ls[event_col]]
		fout.write(l.strip()+'\t'+annotation+'\n')
	fout.close()

def main(args):
	CPMsummary=args.cpm_summary
	CPM_event=loadPTsummary(CPMsummary)
	print('[INFO] Number of event in CPM screen summary:'+str(len(CPM_event)))
	splicing_event_type=args.splicing_event_type
	screening_result_dir=args.outdir.rstrip('/')

	epitope_file_junction_file_list= glob.glob(screening_result_dir+'/'+splicing_event_type+'.*/epitope_summary.junction-based.txt')
	epitope_file_peptide_file_list= glob.glob(screening_result_dir+'/'+splicing_event_type+'.*/epitope_summary.peptide-based.txt')
	extracellular_as_file_list= glob.glob(screening_result_dir+'/*'+splicing_event_type+'.*.ExtraCellularAS.txt')
	screening_file_list=glob.glob(screening_result_dir+'/*.'+splicing_event_type+'.tier1.txt')+glob.glob(screening_result_dir+'/*.'+splicing_event_type+'.tier2tier3.txt')
	

	for fin_fname in screening_file_list:
		print '[INFO] Integrating CPM test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 0, CPM_event, fin_fname+'.integratedCPM.txt')

	for fin_fname in extracellular_as_file_list:
		print '[INFO] Integrating CPM test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 0, CPM_event, fin_fname+'.integratedCPM.txt')

	for fin_fname in epitope_file_peptide_file_list:
		print '[INFO] Integrating CPM test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 1, CPM_event, fin_fname+'.integratedCPM.txt')

	for fin_fname in epitope_file_junction_file_list:
		print '[INFO] Integrating CPM test result to', fin_fname
		annotateRAbyPTEvent(fin_fname, 0, CPM_event, fin_fname+'.integratedCPM.txt')


if __name__ == '__main__':
	main()
