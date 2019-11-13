import sys, argparse, os ,datetime,logging, uuid
import multiprocessing as mp
from . import config

ID = str(uuid.uuid4()).split('-')[0]

def loadFrame(fin):
	nomap=0
	mis=0
	exonFrameDict={}
	for l in open(fin):
		ls=l.strip().split('|')
		if len(ls)>5: #no hits
			if (float(ls[1])/3)*0.4>float(ls[11]): #identities too less
				mis+=1
				continue	
			exon_start_end=ls[0].split('-') ##TODO: name two variables, exonFrameDict key can be "'5end|'+exon_start_coord", the same for end. much clearner variable name in downstream.
			if exon_start_end[0] not in exonFrameDict:
				exonFrameDict[exon_start_end[0]]=[]
			exonFrameDict[exon_start_end[0]].append((exon_start_end[1],[ls[1]]+ls[6:9]))
			if exon_start_end[1] not in exonFrameDict:
				exonFrameDict[exon_start_end[1]]=[]
			exonFrameDict[exon_start_end[1]].append((exon_start_end[0],[ls[1]]+ls[6:9]))
		else:
			nomap+=1
	return exonFrameDict, nomap, mis 

def findJCRange(AS_chrom,JC_coord,exonFrameDict):
	JC_region_s=''# 5 prime (considered +/-)
	JC_region_e=''# 3 prime (considered +/-)
	JC_region_d=''
	evidence_level=[False,False]
	if JC_coord[0] in exonFrameDict:# check the left coord of the juction
		for exon_info in exonFrameDict[JC_coord[0]]:
			if exon_info[0].startswith(AS_chrom+':'):# check if it is the right pair
				exon_len,frame,align_s,align_e=exon_info[1]
				if int(frame)>0: # if the left coord is the start, then find the ideal starting poistion
					JC_region_d='+'
					JC_region_len=min(int(align_e)-int(align_s),30)
					JC_region_s=(int(JC_coord[0])-(int(exon_len)-int(align_e))-JC_region_len,JC_coord[0])
					evidence_level[0]=True
					break ### TODO: allow multiple starting
				else: # if left coord is the end coord, frame doesn't matter, keep proper length(30-33)
					JC_region_d='-'
					if JC_region_e=='':
						JC_region_e=(max(int(JC_coord[0])-33,int(JC_coord[0])-int(exon_len)), JC_coord[0])# max means min here
						evidence_level[1]=True
						break ### allowing loop will give more restrained junc region. but will be slower. disable for now.
					# else:
					# 	JC_region_e=(max(JC_region_e[0],int(JC_coord[0])-int(exon_len)),JC_coord[0]) 		 
			else:
				continue
	if AS_chrom+':'+JC_coord[1] in exonFrameDict:# check the right coord of the juction & check if it is the right pair
		for exon_info in exonFrameDict[AS_chrom+':'+JC_coord[1]]:
			exon_len,frame,align_s,align_e=exon_info[1]
			if int(frame)>0:# if right coord is the end coord, frame doesn't matter, keep proper length(30-33)
				if JC_region_d=='-':
					print 'frame conflict.', JC_coord;JC_region_e='';evidence_level[1]=False;break
				if JC_region_d=='':# when left coord has no hit
					JC_region_d='+'
				JC_region_e=(JC_coord[1],min(int(JC_coord[1])+33,int(exon_info[0])))
				evidence_level[1]=True
				break ### allowing loop will give more restrained junc region. but will be slower. disable for now.
			else:# right with be the start. then find the ideal starting poistion
				if JC_region_d=='+':
					print 'frame conflict.', JC_coord;JC_region_s='';evidence_level[0]=False;break
				if JC_region_d=='':# when left coord has no hit
					JC_region_d='-'
				if int(align_e)<=32:
					JC_region_s=(JC_coord[1],int(JC_coord[1])+int(align_e))
				elif int(align_e)>32:
					shift=int(align_e)%3
					JC_region_len=30+shift
					JC_region_s=(JC_coord[1],int(JC_coord[1])+JC_region_len)
				evidence_level[0]=True
				break ### TO-DO: allow multiple starting
	if JC_region_s!='' and JC_region_e=='':
		if JC_region_d=='+':
			JC_region_e=(JC_coord[1],int(JC_coord[1])+30)
		else:
			JC_region_e=(int(JC_coord[0])-30,JC_coord[0]) 
	return JC_region_s,JC_region_e, JC_region_d, evidence_level

def selectJC(AS_coord,deltaPSI_c2n,cut_off, select_all):
	if select_all:
		return [(AS_coord[2], AS_coord[3],'skp'),(AS_coord[2],AS_coord[0],'inc1'),(AS_coord[1],AS_coord[3],'inc2')]
	if float(deltaPSI_c2n)<cut_off:# tumor skipping
		return [(AS_coord[2], AS_coord[3],'skp')]
	else:
		return [(AS_coord[2],AS_coord[0],'inc1'),(AS_coord[1],AS_coord[3],'inc2')]

def translateNuc(seq, orf_option=3, rm_early_stop=False):
	codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

	prot_seqs = []
	for start_site in range(orf_option):
		cds = seq[start_site:]
		prot_seq = ''
		for n in range(0, len(cds), 3):
			current_triple = cds[n:n+3]

			if len(current_triple)!=3:
				break
			if current_triple.upper() in codontable:
				letter = codontable[current_triple.upper()]
				if letter=='-':
					prot_seq += letter
					break
				if current_triple == current_triple.upper():
					prot_seq += letter
				else:
					prot_seq += letter.lower()
			else:
				prot_seq += '*'
				print 'unknown codon',current_triple.upper() 
				break
		if prot_seq.find('_')!=-1:
			if rm_early_stop:
				continue #totally remove stop codon containing seq
			else:
				prot_seq = prot_seq.split('_')[0] #trim aa after stop codon
		prot_seqs.append((prot_seq,start_site))

	return prot_seqs

def JC2pep(JC_region_bed,ref_genome,outdir,info):#getfasta of the coord and translated based on the codon table #Another To-do: map gtf direction back to uniprot2gtf file.
	JC_region_fasta=JC_region_bed.split('.')[0]+'.fa'
	cmd='bedtools getfasta -fi '+ref_genome+' -bed '+outdir+'/tmp/junction/'+JC_region_bed+' -fo '+outdir+'/tmp/rna/'+JC_region_fasta+' -s -name'
	os.system(cmd)# DO NOT SPLIT!!!!
	JC_nuc_seq_ID=''
	JC_nuc_seq=''
	JC_pep_fasta=info+'.'+JC_region_fasta.split('.')[0]+'.prot.fa'
	fout_JC_pep=open(outdir+'/tmp/prot/'+JC_pep_fasta,'w')
	for row in open(outdir+'/tmp/rna/'+JC_region_fasta):
		if row.startswith('>'):
			if JC_nuc_seq_ID!='':
				JC_nuc_seq_ID+='|'+row.strip()[1:]
			else:
				JC_nuc_seq_ID=row.strip()[1:]
			continue
		if JC_nuc_seq!='':
			first_half_len=len(JC_nuc_seq)
			JC_nuc_seq+=row.strip().upper()
			JC_prot_seq=translateNuc(JC_nuc_seq,1)[0][0]
			if len(JC_prot_seq)<8 or len(JC_prot_seq)*3<=first_half_len:# skip seq finds PTC before junction site
				JC_nuc_seq_ID=''
				JC_nuc_seq=''
				continue
			junction_pos=first_half_len/3+1
			JC_prot_seq_up=JC_prot_seq[:junction_pos-1]+JC_prot_seq[junction_pos-1].lower()+JC_prot_seq[junction_pos:]
			if junction_pos>11:#long upstream junction peptides, often due to half-mapped but quality BLAST aligment(known exon but not/partially known in protein annotation)
				diff=junction_pos-11
				JC_nuc_seq_ID=JC_nuc_seq_ID+'|trim_'+str(diff)
				JC_prot_seq_up=JC_prot_seq_up[diff:]
			#print JC_nuc_seq, JC_prot_seq_up, JC_prot_seq, first_half_len, junction_pos
			fout_JC_pep.write('>{}\n{}\n'.format(JC_nuc_seq_ID,JC_prot_seq_up))
			JC_nuc_seq_ID=''
			JC_nuc_seq=''
			continue
		JC_nuc_seq=row.strip().upper()
	fout_JC_pep.close()
	return JC_pep_fasta

def loadSeq(fin):
	dic={}
	name=''
	form=''
	direction=''
	for l in open(fin):
		if l.startswith('>'):
			name=l.strip()[1:]
			form_field=name.split(':')[5]#May be different for different Bedtools#Feb19
			if form_field.startswith('inc1'):
				form='inc1'
			elif form_field.startswith('inc2'):
				form='inc2'
			elif form_field.startswith('skp'):
				form='skp'
			else:
				exit('error for the prot files.')
			direction=name.split(':')[3]
			continue
		dic[form]=[name,l.strip()]
		name=''
		form=''
	return dic,direction

def compSeq(seq1,seq2,seq3):
	true_junction1=0
	true_junction2=0
	short_seq2=False
	short_seq1=False
	short_seq3=False
	if seq2!='':
		for i,n in enumerate(list(seq1)):
			if len(seq2)>i:# else when seq2 is shorter
				if n.upper()==seq2[i].upper():#format to up case
					if i+1==len(seq1):
						short_seq1=True
						true_junction1=i+1
						break
					continue
				else:
					true_junction1=i
					break
			else:
				short_seq2=True
				true_junction1=i
				break
		if short_seq1==False: #have differences, then output
			seq1=seq1[:true_junction1]+seq1[true_junction1].lower()+seq1[true_junction1+1:]
			if short_seq2==False:
				seq2=seq2[:true_junction1]+seq2[true_junction1].lower()+seq2[true_junction1+1:]
			else: #seq2 has no difference
				seq2=seq2.upper()
		else: #have NO difference in seq1
			seq1=seq1.upper()
			if len(seq2)>len(seq1):#if seq2 is longer
				seq2=seq2[:true_junction1]+seq2[true_junction1].lower()+seq2[true_junction1+1:]
			else: #seq1 seq2 have no difference
				seq2=seq2.upper()
	if seq3!='': #it will skip if seq1 is the same of seq2
		rev_seq1=seq1[::-1]
		rev_seq3=seq3[::-1]
		pass_js2=False#if the difference occur after junction site, THIS CASE WILL ONLY in seq3.
		for j,n in enumerate((list(rev_seq1))):
			if len(seq3)>j:# check if seq3 is long enough
				if n.upper()==rev_seq3[j].upper():#format to up case
					if rev_seq3[j].islower():# if seq3 has a lower letter--junction site
						pass_js2=True
					if j+1==len(seq1):
						short_seq1=True
						true_junction2=i+1
						break
					continue
				else:
					true_junction2=j
					break
			else:
				short_seq3=True
				true_junction2=j
				break
		if short_seq1==False: #have differences, then output
			if short_seq3==False:
				if pass_js2:#if not passed seq3 js, no need to mark
					rev_seq3=rev_seq3[:true_junction2]+rev_seq3[true_junction2].lower()+rev_seq3[true_junction2+1:]
					seq3=rev_seq3[::-1]
				if seq1!=seq1.upper():#only label unique aa when seq1 passed seq2's comparison
					if true_junction2>len(seq1)-true_junction1-1:#if not passed seq1 seq2 js site on seq1 js, no need to mark
						rev_seq1=rev_seq1[:true_junction2]+rev_seq1[true_junction2].lower()+rev_seq1[true_junction2+1:]
						seq1=rev_seq1[::-1]
			else:#seq3 has no difference
				seq3=seq3.upper()
		else:#seq1 has no difference
			seq1==seq1.upper()
			if len(seq3)>len(seq1):#if seq3 is longer
				if pass_js2:
					rev_seq3=rev_seq3[:true_junction2]+rev_seq3[true_junction2].lower()+rev_seq3[true_junction2+1:]
				seq3=rev_seq3[::-1]
			else:# seq1 and seq3 has no difference
				seq3==seq3.upper()
	return seq1, seq2, seq3

def compPepFile(fin,outdir,select_form):
	comp_result_dic={}
	dic,direction=loadSeq(outdir+'/tmp/prot/'+fin)
	if len(dic)>1:
		if 'skp' in dic: #where the comparison is needed
			seq2=''
			seq3=''
			if dic.has_key('inc1'):
				seq2=dic['inc1'][1]
			if dic.has_key('inc2'):
				seq3=dic['inc2'][1]
			
			if direction=='+':
				comp_result=compSeq(dic['skp'][1],seq2,seq3)
				if 'skp' in dic:
					dic['skp'][1]=comp_result[0]
				if 'inc1' in dic:
					dic['inc1'][1]=comp_result[1]
				if 'inc2' in dic:
					dic['inc2'][1]=comp_result[2]

			else:
				comp_result=compSeq(dic['skp'][1],seq3,seq2)
				if 'skp' in dic:
					dic['skp'][1]=comp_result[0]
				if 'inc1' in dic:
					dic['inc1'][1]=comp_result[2]
				if 'inc2' in dic:
					dic['inc2'][1]=comp_result[1]
	if select_form!=2:
		fout_skp=open(outdir+'/tmp/prot.compared/skp/skp.'+fin.split('/')[-1],'w')
		if 'skp' in dic:
			fout_skp.write('>{}\n{}\n'.format(dic['skp'][0],dic['skp'][1]))
		fout_skp.close()
	if select_form!=1:
		fout_inc=open(outdir+'/tmp/prot.compared/inc/inc.'+fin.split('/')[-1],'w')
		if 'inc1' in dic:
			fout_inc.write('>{}\n{}\n'.format(dic['inc1'][0],dic['inc1'][1]))
		if 'inc2' in dic:
			fout_inc.write('>{}\n{}\n'.format(dic['inc2'][0],dic['inc2'][1]))
		fout_inc.close()

def AS2pep(AS_chrom, AS_coord,AS_direction, deltaPSI_c2n,cuf_off, select_all, exonFrameDict,ref_genome,outdir,info):
	JC_coord_list=selectJC(AS_coord,deltaPSI_c2n,cuf_off, True) #under select_all=True
	select_form=len(selectJC(AS_coord,deltaPSI_c2n,cuf_off, select_all))
	JC_region_bed='_'.join([AS_chrom,AS_direction]+AS_coord)+'.JCregion.bed'
	fout_JC_region=open(outdir+'/tmp/junction/'+JC_region_bed,'w') #The JC file created. with all forms included regardless of selecting setting. Required for comparison.
	for JC_coord in JC_coord_list:
		JC_region=findJCRange(AS_chrom,JC_coord[:2], exonFrameDict)
		if JC_region[3][0]:
			fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[0][0],JC_region[0][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':uniprotFrame:'+AS_direction+':form:'+JC_coord[2],'.',JC_region[2]))
			fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[1][0],JC_region[1][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':uniprotFrame:'+AS_direction+':form:'+JC_coord[2],'.',JC_region[2]))
		#else: #no known frame; implement 3 orf search latter
			#fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[0][0],JC_region[0][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':unknownFrame:'+AS_direction,'.',JC_region[2]))
			#fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[1][0],JC_region[1][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':unknownFrame:'+AS_direction,'.',JC_region[2]))
	fout_JC_region.close() 
	JC_pep_fasta=JC2pep(JC_region_bed,ref_genome,outdir,info)# Junction peptides file created. With all forms included regardless of form selection. 
	compPepFile(JC_pep_fasta, outdir, select_form) # Compared Junction peptides file created. Form selection STARTS HERE.

def main(args):
	orf_mapping_file=config.ORF_MAP_PATH
	ref_genome=args.ref_genome #ref_genome='~/nobackup-yxing/references.annotations/37.chr/ucsc.hg19.fasta'
	fin=args.as_input
	outdir=args.outdir.rstrip('/')
	select_all=args.no_tumor_form_selection
	deltaPSI_cut_off=float(args.deltaPSI_cut_off)
	deltaPSI_column=int(args.deltaPSI_column)-1

	os.system('mkdir -p '+outdir+'/tmp')
	os.system('mkdir -p '+outdir+'/tmp/junction '+outdir+'/tmp/rna '+outdir+'/tmp/prot '+outdir+'/tmp/pred')
	os.system('mkdir -p '+outdir+'/tmp/prot.compared '+outdir+'/tmp/prot.compared/skp '+outdir+'/tmp/prot.compared/inc')
	exonFrameDict,nomap, mis =loadFrame(orf_mapping_file)
	print '[INFO] Total exon-orf loaded',len(exonFrameDict), nomap, mis

	#Select AS version, region and frame, translate to peptides.
	events_processed=0
	tot=config.file_len(fin)-1
	for i,l in enumerate(open(fin)):
		if l.startswith('ENSG')==False:
			continue
		config.update_progress(i/(0.0+tot))
		# if i%500==0:
		# 	print i
		ls=l.strip('\n').split('\t')
		events_processed+=1
		des=ls[0].split(':')
		info='_'.join(des[:2]).strip('_').replace('/','+')
		AS2pep(des[2],des[4:8],des[3],ls[deltaPSI_column],deltaPSI_cut_off,False,exonFrameDict,ref_genome,outdir,info)
	print '[INFO] Total processed',events_processed

if __name__ == '__main__':
	main()
