import sys, argparse, os ,datetime,logging, uuid
#import subprocess, csv, re, numpy
#import multiprocessing as mp
# from Bio.Blast import NCBIXML
from subprocess import Popen, PIPE

ID = str(uuid.uuid4()).split('-')[0]

def loadSampleMHC(f_s_in):
	hla_allele_list=f_s_in.split(',')
	return hla_allele_list

def parsePred(stdout):
	ref_dict={}
	ref_out=csv.DictReader(stdout,delimiter='\t')
	for r in ref_out:
		ic50=[k for k in r.keys() if k.find('ic50')!=-1]
		p = re.compile('\d+(\.\d+)?')
		ic50_value=[float(r[k]) for k in ic50 if p.match(r[k]) != None]
		ref_dict[r['allele']+'_'+r['seq_num']+'_'+r['start']+'_'+r['length']]=[numpy.median(ic50_value),r['peptide']]
	return ref_dict

def localIEDBCommand(iedb_path, hla_allele, epitope_len, outdir, JC_pep_fasta, enst_id, form, IEDB_model='IEDB_recommended'):
	# cmd = 'python '+iedb_path+'/predict_binding.py '+IEDB_model+' '+hla_allele+' '+epitope_len+' '+outdir+'/tmp/prot/'+JC_pep_fasta +' > '+outdir+'/tmp/pred/'+JC_pep_fasta.split('.prot')[0]+'.'+IEDB_model+'.'+hla_allele.replace('*','_').replace(':','_')+'.'+epitope_len+'.pred.txt'
	# print cmd
	# os.system(cmd)
	pred_out=outdir+'/tmp/pred/'+enst_id+'/'+form+'.'+JC_pep_fasta.split('.prot')[0]+'.'+IEDB_model+'.'+hla_allele.replace('*','_').replace(':','_')+'.'+epitope_len+'.pred.txt'
	fout=open(pred_out,'w')
	response = Popen(['python',iedb_path+'/predict_binding.py',IEDB_model,hla_allele,epitope_len,outdir+'/tmp/prot.compared/'+form+'/'+form+'.'+JC_pep_fasta], stdout=fout,shell = False)
	#response = Popen(['/u/home/p/panyang/local/bin/python',iedb_path+'/predict_binding.py',IEDB_model,hla_allele,epitope_len,outdir+'/tmp/prot.compared/'+form+'/'+form+'.'+JC_pep_fasta], stdout=fout,shell = False)
	fout.close()
	return response

def pep2antigen(JC_pep_fasta, enst_id, hla_allele_list,epitope_len_list, iedb_path, outdir):
	predicting=[]
	n=0
	#file_list=[]
	for hla_allele in hla_allele_list:
		for epitope_len in epitope_len_list:
			#file_list.append(outdir+'/tmp/'+JC_pep_fasta+' '+hla_allele.replace('*','_').replace(':','_')+'.'+epitope_len+'_iedb.txt')
			predicting.append(mp.Process(target=localIEDBCommand,args=(iedb_path, hla_allele, epitope_len, outdir, JC_pep_fasta, enst_id)))
			predicting[n].start()
			n+=1
	for i in range(0,n):
		predicting[i].join()
	#parsed_dict=[parsePred(open(tmp_file)) for tmp_file in file_list]
	#return parsed_dict
def pep2antigen_single(JC_pep_fasta, enst_id, form, hla_allele_list,epitope_len_list, iedb_path, outdir):
	predicting=[]
	for hla_allele in hla_allele_list:
		for epitope_len in epitope_len_list:
			response = localIEDBCommand(iedb_path, hla_allele, epitope_len, outdir, JC_pep_fasta, enst_id, form)
			predicting.append(response)

	for response in predicting:
		response.wait()


def main(args):
	#Define parameters
	fin=args.junction_pep_input
	outdir=args.outdir.rstrip('/')

	iedb_path=args.iedb_local
	if args.iedb_local!=False:
		iedb_path=args.iedb_local.rstrip('/')
	hla_allele_list=loadSampleMHC(args.hla_allele_list)
	if hla_allele_list==[]:
		sys.exit("# No HLA Alleles. Exit.")
	epitope_len_list=map(int,args.epitope_len_list.split(','))
	if min(epitope_len_list)<8:
		sys.exit("# The request epitope length is too small. Exit.")
	epitope_len_list=map(str,epitope_len_list)

	fs=fin.split('/')[-1].split('.')
	JC_pep_fasta='.'.join(fs[1:])
	form=fs[0]
	if os.stat(outdir+'/tmp/prot.compared/'+form+'/'+form+'.'+JC_pep_fasta).st_size != 0:
		enst_id=JC_pep_fasta.split('_')[0]
		os.system('mkdir -p '+outdir+'/tmp/pred/'+enst_id)
		pred_result_mut=pep2antigen_single(JC_pep_fasta, enst_id, form, hla_allele_list, epitope_len_list, iedb_path, outdir)
 

if __name__ == '__main__':
	main()
