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
		if len(ls)>5: #skip no hits
			if (float(ls[1])/3)*0.4>float(ls[11]): #skip identities too less
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

def loadGTFMicroexonInfo(gtf):#laod gtf and exon start and end and legnth, used to check microexon
	micro_exon={}	
	for l in open(gtf):
		if l.startswith('#'):
			continue
		ls=l.strip().split('\t')
		if ls[2]=='exon':
			exon_name=ls[0]+':'+ls[3]+'-'+ls[4]
			if exon_name in micro_exon:
				continue
			length=int(ls[4])-int(ls[3])
			if length<=30:
				micro_exon[exon_name]=int(ls[4])-int(ls[3])
	microexon_start={}#by junction
	microexon_end={}	
	for e in micro_exon:
		chrom=e.split(':')[0]
		es, ee=e.split(':')[1].split('-')
		es=chrom+':'+es
		ee=chrom+':'+ee
		if es not in microexon_start:
			microexon_start[es]=micro_exon[e]
		else:
			if micro_exon[e]<microexon_start[es]:
				microexon_start[es]=micro_exon[e]
		if ee not in microexon_end:
			microexon_end[ee]=micro_exon[e]
		else:
			if micro_exon[e]<microexon_end[ee]:
				microexon_end[ee]=micro_exon[e]
	return microexon_start,microexon_end
		


def findJCRange(AS_chrom, JC_coord, AS_direction, exonFrameDict, all_orf, microexon_start, microexon_end, ignore_annotation):
	##select the range of junction sequence based on 1) +/-30bp 2)known in annotation 3)ORF/frame in annotation 4)not exceed exon length
	##Step 1: assigned value based on the annotation
	JC_region_s=''# 5 prime (considered +/-)
	JC_region_e=''# 3 prime (considered +/-)
	JC_region_d=''
	evidence_level=[False,False]#This is used when generating bed file annotation column.
	if ignore_annotation==False:

		if JC_coord[0] in exonFrameDict:# check the left coord of the juction
			for exon_info in exonFrameDict[JC_coord[0]]:
				if exon_info[0].startswith(AS_chrom+':'):# check if it is the right pair
					exon_len,frame,align_s,align_e=exon_info[1]
					if int(frame)>0: # if the left coord is the start, then find the ideal starting poistion
						JC_region_d='+'
						JC_region_len=min(int(align_e)-int(align_s),30)
						JC_region_s=(int(JC_coord[0])-(int(exon_len)-int(align_e))-JC_region_len,JC_coord[0])
						evidence_level[0]=True
						break ### Limitation: allow multiple starting
					else: # if left coord is the end coord, frame doesn't matter, keep proper length(30-33)
						JC_region_d='-'
						if JC_region_e=='':
							JC_region_e=(max(int(JC_coord[0])-33,int(JC_coord[0])-int(exon_len)), JC_coord[0])# max means min here
							evidence_level[1]=True
							break ### allowing loop will give more restrained junc region. but will be slower. disable for now.
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
					break ### Limitation: allow multiple starting
		if JC_region_s!='' and JC_region_e=='':#check if JC region end is beyond exon end. handle microexons
			if JC_region_d=='+':
				range_bp=30#default	
				if AS_chrom+':'+str(int(JC_coord[1])+1) in microexon_start:#check if downstream exon start pos is a  microexon start
					range_bp=min(microexon_start[AS_chrom+':'+str(int(JC_coord[1])+1)],30)	
				JC_region_e=(JC_coord[1],int(JC_coord[1])+range_bp+1)#+/-1 from gtf
			else:
				range_bp=30#default				
				if AS_chrom+':'+JC_coord[0] in microexon_end:#check if upstream exon(downstream in translation) start pos is microexon end
					range_bp=min(microexon_end[AS_chrom+':'+JC_coord[0]],30)	
				JC_region_e=(int(JC_coord[0])-range_bp-1,JC_coord[0]) # For events partially in the ORF annotation: Extend the end based on the known start frame #Control Microexon!!
	##Step 2:  assign 30bp +/-. This is for 3 ORF search. Can be used when need/force 3 ORF search. TODO: check up/down exon length 
	if (all_orf and evidence_level[0]==False) or ignore_annotation: # For events not in the ORF annotation: use all 3 orf. Note: use 'evidence_level' instead of 'JC_region_s' to avoid mini exon/wrong annotation.
		JC_region_d = AS_direction
		if JC_region_d=='+':
			JC_region_s, JC_region_e= [(int(JC_coord[0])-30,JC_coord[0]),(JC_coord[1],int(JC_coord[1])+30)]
		else:
			JC_region_s, JC_region_e= [(JC_coord[1],int(JC_coord[1])+30),(int(JC_coord[0])-30,JC_coord[0])]
	return JC_region_s, JC_region_e, JC_region_d, evidence_level

#handle multiple AS types by taking 4 or 6 coordinate and select antigen-deriving junctions
def selectJunction(AS_coord, deltaPSI_c2n, cut_off, if_select_all, splicing_event_type): 
	if splicing_event_type == 'SE':
		skp = (AS_coord[2], AS_coord[3],'skp')
		inc1 = (AS_coord[2],AS_coord[0],'inc1')
		inc2 = (AS_coord[1],AS_coord[3],'inc2')
	elif splicing_event_type == 'A3SS':
		skp = (AS_coord[5], AS_coord[2],'skp')
		inc1 = (AS_coord[5],AS_coord[0],'inc1')
		inc2 = (AS_coord[2],AS_coord[2],'inc2')
	elif splicing_event_type == 'A5SS':
		skp = (AS_coord[3], AS_coord[4],'skp')
		inc1 = (AS_coord[3],AS_coord[3],'inc1')
		inc2 = (AS_coord[1],AS_coord[4],'inc2')
	elif splicing_event_type == 'RI':
		skp = (AS_coord[3], AS_coord[4],'skp')
		inc1 = (AS_coord[3],AS_coord[3],'inc1')
		inc2 = (AS_coord[4],AS_coord[4],'inc2')

	if if_select_all:
		return [skp, inc1, inc2]
	if float(deltaPSI_c2n) < cut_off:  # tumor skipping
		return [skp]
	else:
		return [inc1, inc2]

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

def JC2pep(JC_region_bed, ref_genome, outdir, info, remove_early_stop, all_orf, pep_dir_prefix):#getfasta of the coord and translated based on the codon table #Another To-do: map gtf direction back to uniprot2gtf file.
	orf_range=1 if all_orf==False else 3
	JC_region_fasta=JC_region_bed.split('.')[0]+'.fa'
	cmd='bedtools getfasta -fi '+ref_genome+' -bed '+outdir+'/tmp/junction/'+JC_region_bed+' -fo '+outdir+'/tmp/rna/'+JC_region_fasta+' -s -name'
	os.system(cmd)# DO NOT SPLIT!!!!
	JC_nuc_seq_ID=''
	JC_nuc_seq=''
	JC_pep_fasta_name=info+'.'+JC_region_fasta.split('.')[0]
	fout_JC_pep=open(outdir+'/tmp/'+pep_dir_prefix+'/'+JC_pep_fasta_name+'.prot.fa','w')
	if all_orf:
		fout_JC_pep_1=open(outdir+'/tmp/'+pep_dir_prefix+'/'+JC_pep_fasta_name+'_orf1.prot.fa','w')
		fout_JC_pep_2=open(outdir+'/tmp/'+pep_dir_prefix+'/'+JC_pep_fasta_name+'_orf2.prot.fa','w')
	
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
			JC_prot_seq_list=translateNuc(JC_nuc_seq, orf_range, remove_early_stop)
			for JC_prot_seq_pair in JC_prot_seq_list:
				JC_prot_seq,JC_prot_orf_start=JC_prot_seq_pair
				if len(JC_prot_seq)<8 or len(JC_prot_seq)*3<=first_half_len:# skip seq finds PTC before junction site
					continue
				junction_pos=first_half_len/3+1
				JC_prot_seq_up=JC_prot_seq[:junction_pos-1]+JC_prot_seq[junction_pos-1].lower()+JC_prot_seq[junction_pos:]
				if junction_pos>11:#long upstream junction peptides, often due to half-mapped but quality BLAST aligment(known exon but not/partially known in protein annotation)
					diff=junction_pos-11
					JC_nuc_seq_ID=JC_nuc_seq_ID+'|trim_'+str(diff)
					JC_prot_seq_up=JC_prot_seq_up[diff:]
				if JC_prot_orf_start==0:
					fout_JC_pep.write('>{}\n{}\n'.format(JC_nuc_seq_ID,JC_prot_seq_up))
				elif JC_prot_orf_start==1:
					fout_JC_pep_1.write('>{}\n{}\n'.format('FrameAdd1'.join(JC_nuc_seq_ID.split('Frame')),JC_prot_seq_up))
				elif JC_prot_orf_start==2:
					fout_JC_pep_2.write('>{}\n{}\n'.format('FrameAdd2'.join(JC_nuc_seq_ID.split('Frame')),JC_prot_seq_up))
			JC_nuc_seq_ID=''
			JC_nuc_seq=''
			continue
		JC_nuc_seq=row.strip().upper()
	fout_JC_pep.close()
	if all_orf:
		fout_JC_pep_1.close()
		fout_JC_pep_2.close()
	return JC_pep_fasta_name

def loadSeq(fin):
	dic={}
	name=''
	form=''
	direction=''
	for l in open(fin):
		if l.startswith('>'):
			name=l.strip()[1:]
			form_field=name.split(':')[5]#May different for different Bedtools. Feb19
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
						true_junction2=j+1
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

def compPepFile(fin,outdir,pep_dir_prefix,select_form, novel_info):
	comp_result_dic={}
	dic,direction=loadSeq(outdir+'/tmp/'+pep_dir_prefix+'/'+fin)
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
	if novel_info==[]:
		novel_info=dic
	if select_form!=2:
		fout_skp=open(outdir+'/tmp/'+pep_dir_prefix+'.compared/skp/skp.'+fin.split('/')[-1],'w')
		if 'skp' in dic and 'skp' in novel_info:
			fout_skp.write('>{}\n{}\n'.format(dic['skp'][0],dic['skp'][1]))
		fout_skp.close()
	if select_form!=1:
		fout_inc=open(outdir+'/tmp/'+pep_dir_prefix+'.compared/inc/inc.'+fin.split('/')[-1],'w')
		if 'inc1' in dic and 'inc1' in novel_info:
			fout_inc.write('>{}\n{}\n'.format(dic['inc1'][0],dic['inc1'][1]))
		if 'inc2' in dic and 'inc2' in novel_info:
			fout_inc.write('>{}\n{}\n'.format(dic['inc2'][0],dic['inc2'][1]))
		fout_inc.close()

def AS2pep(AS_chrom, AS_coord, AS_direction, deltaPSI_c2n, cuf_off, if_select_all, all_orf, pep_dir_prefix, microexon_start, microexon_end, ignore_annotation, remove_early_stop, splicing_event_type, exonFrameDict, ref_genome, outdir, info, novel_info):
	JC_coord_list = selectJunction(AS_coord, deltaPSI_c2n, cuf_off, True, splicing_event_type) #under if_select_all=True
	select_form = len(selectJunction(AS_coord, deltaPSI_c2n, cuf_off, if_select_all, splicing_event_type))
	JC_region_bed = '_'.join([AS_chrom,AS_direction]+AS_coord)+'.JCregion.bed'
	fout_JC_region = open(outdir+'/tmp/junction/'+JC_region_bed,'w') #The JC file created. with all forms included regardless of selecting setting. Required for comparison.
	for JC_coord in JC_coord_list:
		JC_region=findJCRange(AS_chrom, JC_coord[:2], AS_direction, exonFrameDict, all_orf, microexon_start, microexon_end, ignore_annotation)
		if JC_region[3][0]:#check if the starting end orf is known in annotation
			fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[0][0],JC_region[0][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':uniprotFrame:'+JC_region[2]+':form:'+JC_coord[2],'.',JC_region[2]))
			fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[1][0],JC_region[1][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':uniprotFrame:'+JC_region[2]+':form:'+JC_coord[2],'.',JC_region[2]))
		elif (all_orf or ignore_annotation): #3 orf frame;
			fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[0][0],JC_region[0][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':unknownFrame:'+AS_direction+':form:'+JC_coord[2],'.',JC_region[2]))
			fout_JC_region.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(AS_chrom,JC_region[1][0],JC_region[1][1],AS_chrom+':'+JC_coord[0]+'|'+JC_coord[1]+':unknownFrame:'+AS_direction+':form:'+JC_coord[2],'.',JC_region[2]))
	fout_JC_region.close() 
	JC_pep_fasta_name=JC2pep(JC_region_bed, ref_genome, outdir, info, remove_early_stop, all_orf, pep_dir_prefix)# Junction peptides file created. With all forms included regardless of form selection. 
	compPepFile(JC_pep_fasta_name+'.prot.fa', outdir, pep_dir_prefix, select_form, novel_info) # Compared Junction peptides file created. Form selection STARTS HERE.
	if all_orf:
		compPepFile(JC_pep_fasta_name+'_orf1.prot.fa', outdir, pep_dir_prefix, select_form, novel_info)
		compPepFile(JC_pep_fasta_name+'_orf2.prot.fa', outdir, pep_dir_prefix, select_form, novel_info)

def main(args):
	orf_mapping_file=config.ORF_MAP_PATH
	ref_genome=args.ref_genome
	fin=args.as_input
	splicing_event_type=args.splicing_event_type
	all_orf=args.all_reading_frames
	gtf=args.gtf
	microexon_start, microexon_end=loadGTFMicroexonInfo(gtf)
	ignore_annotation=args.ignore_annotation
	remove_early_stop=args.remove_early_stop
	outdir=args.outdir.rstrip('/')
	if_select_all=args.no_tumor_form_selection
	deltaPSI_cut_off=float(args.deltaPSI_cut_off)
	deltaPSI_column=int(args.deltaPSI_column)-1
	check_novel=args.check_novel
	pep_dir_prefix='prot'
	if all_orf:
		pep_dir_prefix='prot_allorf'
	os.system('mkdir -p '+outdir+'/tmp')
	os.system('mkdir -p '+outdir+'/tmp/junction '+outdir+'/tmp/rna '+outdir+'/tmp/'+pep_dir_prefix+' '+outdir+'/tmp/pred')
	os.system('mkdir -p '+outdir+'/tmp/'+pep_dir_prefix+'.compared '+outdir+'/tmp/'+pep_dir_prefix+'.compared/skp '+outdir+'/tmp/'+pep_dir_prefix+'.compared/inc')
	exonFrameDict,nomap, mis =loadFrame(orf_mapping_file)
	print '[INFO] Total exon-orf loaded',len(exonFrameDict), nomap, mis

	#Select AS version, region and frame, translate to peptides.
	events_processed=0
	tot=config.file_len(fin)-1
	header=[]
	for i,l in enumerate(open(fin)):
		if l.startswith('ENSG')==False:
			header=l.strip().split('\t')
			continue
		config.update_progress(i/(0.0+tot))
		ls=l.strip().split('\t')
		ld=dict(zip(header,ls))
		events_processed+=1
		des=ls[0].split(':')
		info='_'.join(des[:2]).strip('_').replace('/','+')
		novel_info=[]
		if check_novel:
			novel_info=ld['novel_ss_info'].split(';')
		AS2pep(des[2], des[4:], des[3], ls[deltaPSI_column], deltaPSI_cut_off, False, all_orf, pep_dir_prefix, microexon_start, microexon_end, ignore_annotation, remove_early_stop, splicing_event_type, exonFrameDict, ref_genome, outdir, info, novel_info)
	print '[INFO] Total processed',events_processed

if __name__ == '__main__':
	main()
