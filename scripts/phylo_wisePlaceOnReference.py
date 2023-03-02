#! /usr/bin/env python3
from Bio import SeqIO
from querier3.utils import *
from querier3.wise2 import *
from argparse import ArgumentParser
import os, os.path, sys, subprocess, shutil, json, math
from ete3 import *

if not 'WISECONFIGDIR' in os.environ:
    os.environ["WISECONFIGDIR"]="/home/gkoczyk/miniconda2/envs/phyloplace/share/wise2/wisecfg/"
argparser = ArgumentParser("Automated phylogenetic placement of short DNA sequences (amplicons/reads) on a reference protein phylogeny using WISE2 estimated protein fragments")
argparser.add_argument('query_fn', help='Input file (FASTA) with short sequences', type=os.path.abspath)
argparser.add_argument('out_dirn', help='Name of output directory  [DEFAULT: %(default)s]', type=str, default='output', nargs='?')
argparser.add_argument('ref_prot_tree_fn', help='Reference protein phylogeny [DEFAULT: %(default)s]',type=os.path.abspath, default="/data/reference.nwk", nargs='?')
argparser.add_argument('ref_prot_ali_fn', help='Reference protein sequence alignment  [DEFAULT: %(default)s]', type=os.path.abspath, default="/data/reference.align.fasta", nargs='?')
argparser.add_argument('ref_prot_raw_fn', help='Reference protein sequences unaligned  [DEFAULT: %(default)s]', type=os.path.abspath, default="/data/reference.unali.fasta", nargs='?')
argparser.add_argument('ref_hmm_fn', help='Reference HMMer v2 model corresponding to reference alignment (or untrimmed alignment)  [DEFAULT: %(default)s]',type=os.path.abspath, default="/data/reference.hmm", nargs='?')
argparser.add_argument('--annots_fn', help='Name of annotations TSV file for visualization  [DEFAULT: %(default)s]', type=os.path.abspath, default="/data/reference.annots.tsv")
argparser.add_argument('--clades_fn', help='Name of clade TSV file for classification  [DEFAULT: %(default)s]', type=os.path.abspath, default="/data/clades.tsv")
argparser.add_argument('--far_dmnd_fn', help='Name of diamond database for further classification of (mainly) outliers [DEFAULT: %(default)s]', type=os.path.abspath, default="/data/uniprot_sprot.dmnd")



#DIAMOND_EXEC='diamond'
MAX_DIAMOND_EVALUE='1e-10'
DIAMOND_CMD='diamond blastx -q {query_fn} -d {db_fn} -f 6 qseqid sseqid pident evalue bitscore qframe stitle -k 1 -e {evalue} -o {out_fn} --quiet' #>/dev/null 2>/dev/null '
DIAMOND_MAKEDB_CMD='diamond makedb --in {input_fn} -d {db_fn} --quiet >/dev/null 2>/dev/null '
ISOFORM_SEP="####"
def degap(s):
    return ''.join(l for l in s if l!='-')
def makeDiamondDb(ali_fn, out_db_fn):
    with tmpfile(give_name=True) as (tmp_fh, tmp_fn):
        with open(tmp_fn, 'w') as wfh:
            for seq in SeqIO.parse(ali_fn, 'fasta'):
                print(">" + seq.description, file=wfh)
                print(degap(str(seq.seq)), file=wfh)
        diamond_makedb_cmdstr = DIAMOND_MAKEDB_CMD.format(input_fn=tmp_fn, db_fn=out_db_fn)
        subprocess.run(diamond_makedb_cmdstr, shell=True)
        #print(diamond_makedb_cmdstr)
        #input('WAIT')

def runDiamond(query_fn, db_fn, out_fn, evalue):
    diamond_cmdstr = DIAMOND_CMD.format(query_fn=query_fn, db_fn=db_fn, evalue=evalue, out_fn=out_fn)
    #print(diamond_cmdstr)
    #sys.exit(0)
    return subprocess.run(diamond_cmdstr, shell=True)
@handlify(mode='r')
def readDiamondTsv(fh):
    for line in fh:
        fields = line.split('\t')
        yield StructObject(query=fields[0], hit=fields[1], percid=float(fields[2]), evalue=float(fields[3]), bitscore=float(fields[4]), qframe=int(fields[5]), stitle=fields[6].strip())

GENEWISEDB_CMD='genewisedb -init local -pfam {hmm_fn} -dnadb {query_fn} -quiet -block 200 -genes -aln {max_size} -kbyte {max_mem} -pthr_no {ncpu} -pthread -cut 10. > {out_fn}'
def runGenewiseDb(query_fn, hmm_fn, out_fn, max_size=10000000, max_mem=10000000, ncpu=4):
    genewisedb_cmdstr = GENEWISEDB_CMD.format(query_fn=query_fn, hmm_fn=hmm_fn, out_fn=out_fn, max_size=max_size, max_mem=max_mem, ncpu=ncpu)
    #print(diamond_cmdstr)
    #sys.exit(0)
    return subprocess.run(genewisedb_cmdstr, shell=True)

GENEWISEDBPROT_CMD='genewisedb -init local -prodb {protdb_fn} -dnadb {query_fn} -quiet -block 200 -genes -aln {max_size} -kbyte {max_mem} -cut 10. > {out_fn}'
def runGenewiseDbUsingProts(query_fn, protdb_fn, out_fn, max_size=10000000, max_mem=10000000, ncpu=4):
    genewisedb_cmdstr = GENEWISEDBPROT_CMD.format(query_fn=query_fn, protdb_fn=protdb_fn, out_fn=out_fn, max_size=max_size, max_mem=max_mem, ncpu=ncpu)
    #print(genewisedb_cmdstr)
    #sys.exit(0)
    return subprocess.run(genewisedb_cmdstr, shell=True)

def runGenewiseDbForSeqs(seqs, hmm_fn, out_fn, max_size=10000000, max_mem=10000000, ncpu=4):
    hmm_fn = os.path.abspath(hmm_fn)
    out_fn = os.path.abspath(out_fn)
    with tempDir(change_dir=True) as dirn:
        tmp_query_fn = 'tmp.fasta'
        with open(tmp_query_fn, 'w') as wfh:
            for seq_id in seqs:
                print(">%s\n%s" % (seq_id, seqs[seq_id]), file=wfh)
        r = runGenewiseDb(tmp_query_fn, hmm_fn, out_fn, max_size=max_size, max_mem=max_mem, ncpu=4)
    return r

def runGenewiseDbForSeqsUsingProts(seqs, protdb_fn, out_fn, diamondr=None, max_size=10000000, max_mem=10000000, ncpu=4):
    protdb_fn = os.path.abspath(protdb_fn)
    out_fn = os.path.abspath(out_fn)
    with tempDir(change_dir=True) as dirn:
        tmp_query_fn = 'tmp.fasta'
        with open(tmp_query_fn, 'w') as wfh:
            for seq_id in seqs:
                print(">%s\n%s" % (seq_id, seqs[seq_id]), file=wfh)
        if diamondr:
            tmpdb_fn = 'tmpdb.fasta'
            dbids = set()
            for seq_id in seqs:
                if seq_id in diamondr:
                    dbids.add(diamondr[seq_id].hit)
            with open(tmpdb_fn, 'w') as wfh:
                for prot in SeqIO.parse(protdb_fn, 'fasta'):
                    if prot.id in dbids:
                        print(">%s\n%s" % (prot.id, str(prot.seq)), file=wfh)
            #print(open(tmp_query_fn, 'r').read(), "\n", open(tmpdb_fn, 'r').read())
            r = runGenewiseDbUsingProts(tmp_query_fn, tmpdb_fn, out_fn, max_size=max_size, max_mem=max_mem, ncpu=4)
            #print("DONE")
        else:
            r = runGenewiseDbUsingProts(tmp_query_fn, protdb_fn, out_fn, max_size=max_size, max_mem=max_mem, ncpu=4)
    return r

#genewisedb_cmdstr = GENEWISEDB_CMD.format(query_fn=query_fn, hmm_fn=hmm_fn, out_fn=out_fn, max_size=max_size, max_mem=max_mem, ncpu=ncpu)
    #print(diamond_cmdstr)
    #sys.exit(0)
    #return subprocess.run(genewisedb_cmdstr, shell=True)


@handlify(mode='rt')
def parseJplace(fh):
    r = json.loads(fh.read())
    fnames = r.get('fields',[])
    placements = r.get('placements',[])
    for rec in placements:
        names = rec['n']
        ps = rec['p']
        v = StructObject(seqid=names[0], seqids=names, placements=[])
        for p in ps:
            v.placements.append(StructObject(**dict(zip(fnames, p))))
        yield v
    
#############################
# TODO: encapsulate proper MAFFT alignment
#MAFFT_CMD= 'mafft --thread {ncpu} --addfragments {input_fn} --keeplength {orig_ali_fn} > {output_fn}'
MAFFT_CMD= 'mafft --6merpair --thread {ncpu} --addfragments {input_fn} --keeplength {orig_ali_fn} > {output_fn} 2>/dev/null'
def runMafft(input_fn, orig_ali_fn, out_fn, seqids, ncpu=4):
    with tmpfile(give_name=True) as (tmp_fh, tmp_fn):
        mafft_cmdstr = MAFFT_CMD.format(input_fn=input_fn, orig_ali_fn=orig_ali_fn, output_fn=tmp_fn, ncpu=ncpu)
        #print(mafft_cmdstr)
        #raw("WAIT")
        subprocess.run(mafft_cmdstr, shell=True)
        with open(out_fn, 'w') as wfh:
            for seq in SeqIO.parse(tmp_fn, 'fasta'):
                #if seq.description.startswith('_R_'):
                #    seq.description = seq.description[3:]
                #    seq.id = seq.id[3:]
                if seq.id in seqids:
                    print(">" + seq.description, file=wfh)
                    print(str(seq.seq), file=wfh)
    #mafft_cmdstr = MAFFT_CMD.format(input_fn=input_fn, orig_ali_fn=orig_ali_fn, output_fn=out_fn, ncpu=ncpu)
    #print(mafft_cmdstr)
    #raw("WAIT")
    #subprocess.run(mafft_cmdstr, shell=True)
    
    
############################
EPANG_CMD='epa-ng -t {ref_tree_fn} -s {ref_ali_fn} -q {query_ali_fn} -m {model} --redo -w {out_dirn} >/dev/null'
def runEpaNg(query_ali_fn, ref_ali_fn, ref_tree_fn, out_dirn, model='LG'):
    #epa-ng -t rooted_protein_tree.nwk -s reference_set_protein.mafft -q x.fasta -m LG
    epang_cmdstr = EPANG_CMD.format(query_ali_fn=query_ali_fn, ref_ali_fn=ref_ali_fn, ref_tree_fn=ref_tree_fn, model=model, out_dirn=out_dirn)
    subprocess.run(epang_cmdstr, shell=True)
    #print(os.getcwd())

GAPPA_GRAFT_CMD='gappa examine graft --jplace-path {input_fn}  --allow-file-overwriting --out-dir {out_dirn} > /dev/null'
def runGappaGraft(input_fn, out_dirn):
    gappa_graft_cmdstr = GAPPA_GRAFT_CMD.format(input_fn=input_fn, out_dirn=out_dirn)
    subprocess.run(gappa_graft_cmdstr, shell=True)
    
def translateMsaFn(in_fn, out_fn):
    with open(out_fn, 'w') as wfh:
        for seq in SeqIO.parse(in_fn, 'fasta'):
            print(">%s" % (seq.description), file=wfh)
            trseq = ''
            for codon in [ seq[i:i+3] for i in range(0, len(seq), 3) ]:
                ctxt = str(codon.seq)
                #print(ctxt)
                if ctxt=='---':
                    trseq+='-'
                elif '-' in ctxt:
                    trseq+='X'
                else:
                    trseq+=str(codon.translate().seq)
                #trseq = str(seq.translate(gap='-').seq)
            print(trseq, file=wfh)
def translateExons(seq, exons):
    result = None
    for exon in exons:
        subseq = seq[exon.xstart-1:exon.xend]
        if exon.frame<0:
            subseq = subseq.reverse_complement()
        if result is None:
            result = subseq
        else:
            result+=subseq
    # Trim to multiple of 3
    result = result[:int(math.floor(len(result)/3)*3)]
    return str(result.translate().seq)

if __name__=='__main__':
    args = argparser.parse_args()
    args.out_dirn = os.path.join('/current', args.out_dirn) 
    if args.clades_fn:
        clade_data = {rec.clade_name:rec for rec in  parseTsvList(args.clades_fn) }
    else:
        clade_data = {}
    # If DIAMOND dbs do not exist - initialise
    dmnd_ref_prot_fn = args.ref_prot_ali_fn + '.dmnd'
    if not os.path.exists(dmnd_ref_prot_fn):
        makeDiamondDb(args.ref_prot_ali_fn, dmnd_ref_prot_fn)

    # Create the output directory if needed
    os.makedirs(args.out_dirn, exist_ok=True)
    # Final output file names
    summary_fn = os.path.join(args.out_dirn, 'summary.tsv')
    final_diamond_fn = os.path.join(args.out_dirn, 'diamond.tsv')
    final_far_diamond_fn = os.path.join(args.out_dirn, 'fardb_diamond.tsv')
    final_wise_fn = os.path.join(args.out_dirn,'wise2.txt')
    final_filt_query_fn = filt_query_fn = os.path.join(args.out_dirn, 'filtered.fasta')
    final_query_translated_fn = os.path.join(args.out_dirn, 'query_translated.fasta')
    final_query_prot_ali_fn = os.path.join(args.out_dirn, 'query_aligned_protein.fasta')
    # This one is by default - if we want another we have to rewrite runEpaNg to move files around
    final_epa_jplace_fn = os.path.join(args.out_dirn, 'epa_result.jplace')    
    final_epa_newick_fn = os.path.join(args.out_dirn, 'epa_result.newick')
    # Everything runs in a temporary dir
    with tempDir(change_dir=True) as tmp_dirn:
        # Quickly curate input to avoid quotes which mess up JavaScript
        tmp_query_fn = 'query.fasta'
        all_counts = {}
        with open(tmp_query_fn, 'w') as wfh:
            for seq in SeqIO.parse(args.query_fn, 'fasta'):
                if 'count=' in seq.description:
                    count = int(seq.description.split('count=')[-1].strip())
                else:
                    count = 1
                new_seqid =seq.id.replace("'",'').replace('"','')
                all_counts[new_seqid]=count
                print(">"+new_seqid, file=wfh)
                print(str(seq.seq),file=wfh)

        # Work on this file now as the query
        args.query_fn = tmp_query_fn
        # 1) Run diamond
        if str(args.far_dmnd_fn)!='None':
            runDiamond(args.query_fn, args.far_dmnd_fn, final_far_diamond_fn, evalue=1e-4)
            fardb_recs = { rec.query:rec for rec in readDiamondTsv(final_far_diamond_fn) }
        else:
            fardb_recs = {}

        runDiamond(args.query_fn,  dmnd_ref_prot_fn, final_diamond_fn, evalue=MAX_DIAMOND_EVALUE)
        #shutil.copy(diamond_fn, final_diamond_fn)
        # 2) Filter results of DIAMOND, output summary message about rejections
        diamond_recs = { rec.query:rec for rec in readDiamondTsv(final_diamond_fn) }
        filt_seqids = set( rec.query for rec in diamond_recs.values() )
        # 2a) Cut short if nothing left
        if len(diamond_recs)==0:
            print('WARNING: no matches to reference database, check input', file=sys.stderr)
        else:
            all_seqids = {}

            # 3) Run alignment, with cut insertions, remember to quash standard output, translate
            with open(filt_query_fn, 'w') as wfh:
                for seq in SeqIO.parse(args.query_fn, 'fasta'):
                    all_seqids[seq.id]=len(seq)
                    #print(seq.id, seq.description)
                    if seq.id in filt_seqids:
                        print(">" + seq.description, file=wfh)
                        #if recs[seq.id].qframe<0:
                        #    print(str(seq.reverse_complement().seq),file=wfh)
                        #else:
                        print(str(seq.seq), file=wfh)
            runGenewiseDb(filt_query_fn, args.ref_hmm_fn, final_wise_fn)
            # CURRENT: Sum of exons is taken into account (hence readGenewiseStructsForSingle) ALTERNATIVE_TO_DO:Only first predicted gene/fragment is taken into account
            genewise = { rec.seqid:rec for rec in readWiseGeneStructsForSingle(final_wise_fn) }

            # For the ones that were identified by DIAMOND but seemingly have no match, rerun one by one against raw proteins (F. equiseti bug)
            counter = 0           
            for seq in SeqIO.parse(filt_query_fn, 'fasta'):
                if seq.id not in genewise:
                    counter+=1
                    addit_wise_fn = final_wise_fn + '.%00d' % counter
                    runGenewiseDbForSeqsUsingProts( { seq.id: str(seq.seq) }, args.ref_prot_raw_fn, addit_wise_fn, diamondr=diamond_recs )
                    gw = { rec.seqid:rec for rec in readWiseGeneStructsForSingle(addit_wise_fn) }
                    if seq.id in gw:
                        genewise[seq.id] = gw[seq.id]
                            
            wise2filt = {rec.seqid:rec.seqid.split(WISE2_MULTIGENE_SEP,1)[0] for rec in genewise.values()}
            filt2wise = defaultdict(list)
            for k in wise2filt:
                filt2wise[wise2filt[k]].append(k)
            
            with open(final_query_translated_fn, 'w') as wfh:
                for seq in SeqIO.parse(filt_query_fn, 'fasta'):
                    for seq_id in list(filt2wise[seq.id]):
                        if seq_id!=seq.id:
                            print(seq_id, seq.id, "PROBLEM", file=sys.stderr)
                            all_seqids[seq_id]=len(seq)
                            filt_seqids.add(seq_id)
                        try:
                            translated = translateExons(seq, genewise[seq_id].exons)
                            print(">%s" % seq_id, file=wfh)                            
                            print(translated,file=wfh)
                        except:
                            print(seq_id, genewise[seq_id], file=sys.stderr)
                            #filt_seqids.remove(seq_id)
                            pass
            # Rerun if any of the filtered seqids are missing from results
            
            runMafft(final_query_translated_fn, args.ref_prot_ali_fn, final_query_prot_ali_fn, seqids=filt_seqids)
            #print("DONE MAFFT")
            # Translate
            #translateMsaFn(final_query_ali_fn, final_query_prot_ali_fn)
            # TODO: 3a) Output error message if alignment file empty        
            # 4) Run EPA-ng
            runEpaNg(final_query_prot_ali_fn, args.ref_prot_ali_fn, args.ref_prot_tree_fn, args.out_dirn)
            #print("DONE EPA")
            #shutil_copy('epa_result.*')
            # 5) Run gappa examine graft to convert jplace into the final tree
            runGappaGraft(final_epa_jplace_fn, args.out_dirn)

            result_tree = PhyloTree(final_epa_newick_fn)
            ref_tree = PhyloTree(args.ref_prot_tree_fn)
            # This is not kept
            # # Run an additional check of final_epa_newick_fn to adjust branch lengths
            # result_tree = PhyloTree(final_epa_newick_fn)
            # ref_tree = PhyloTree(args.ref_prot_tree_fn)
            # to_change = set()
            # for node in ref_tree.traverse():
            #     if node.dist==0. and not node is ref_tree:
            #         #try:
            #         #    counterpart = result_tree&result_tree.get_common_ancestor([result_tree&leaf.name for leaf in node])
            #         #except:
            #         #    print([ leaf.name for leaf in node ])
            #         #    raise
            #         to_change.add(frozenset(leaf.name for leaf in node))
            #         #counterpart.dist=0.
            #         #print(counterpart.dist)
            #         #print(counterpart.up.dist)
            #         #print("SET for", [leaf.name for leaf in node])
            # for node in result_tree.traverse(strategy='postorder'):
            #     leaves = frozenset(leaf.name for leaf in node if leaf.name not in filt_seqids)
            #     if leaves in to_change:
            #         node.dist=0.
            #         to_change.remove(leaves)
            #     if len(to_change)==0:
            #         break
            with open(final_epa_newick_fn, 'w') as wfh:
                print(result_tree.write(format=0), file=wfh)

            # Initialise labelling of both trees
            ino2node = {}
            for ino, node in enumerate(ref_tree.traverse(strategy='postorder'), start=0):
                node.labels = set()
                ino2node[ino]=node
            for node in result_tree.traverse(strategy='postorder'):
                node.labels = set()
            labels=set(clade.clade_name for clade in clade_data.values())
            for clade in clade_data.values():
                ancestor_ref = ref_tree.get_common_ancestor(clade.tip1, clade.tip2)
                #print(ancestor_ref)
                for node in ancestor_ref.traverse():
                    node.labels.add(clade.clade_name)
                ancestor_res = result_tree.get_common_ancestor(clade.tip1, clade.tip2)
                for node in ancestor_res.traverse():
                    node.labels.add(clade.clade_name)
            # Find inclusive labels, so we can avoid putting parents if not needed
            label_parents = defaultdict(set) 
            for label1 in labels:
                for label2 in labels:
                    if label1!=label2:
                        guard=True
                        for node in ref_tree:
                            if label1 in node.labels and label2 not in node.labels:
                                guard=False
                                break
                        if guard:
                            label_parents[label1].add(label2)


            ## Output the final placements                        
            # 
            # - amplicon name (include ####2,3 etc if there are multiple genes annotated)
            # - DIAMOND prefilter results: seqid, evalue, percent identity
            # - classified labels, including OUTLIER and OTHER special labels 
            # - comma separated weights for classified labels - weight is 1.0 by default for outlier and/or NO_SIMILARITY
            # - DIAMOND fardb (default Sprot) results

            with open(summary_fn, 'w') as out_fh:
                fnames = ['amplicon_id', 'count', 'length', 'labels', 'weights', 'diamond_hit_id', 'diamond_evalue', 'diamond_idperc', 'fardb_hit_id', 'fardb_evalue', 'fardb_percid', 'fardb_description' ]
                jplace_r = { v.seqid:v for v in parseJplace(final_epa_jplace_fn) }
                print(*fnames, sep='\t', file=out_fh)
                #print(filt_seqids)
                for seqid in sorted(all_seqids):
                    orig_seq_id = seqid.split(WISE2_MULTIGENE_SEP,1)[0]
                    if orig_seq_id in filt_seqids:
                        r = defaultdict(float)
                        try:
                            prec = jplace_r[seqid]
                        except:
                            # OUTLIER2 - where there is no predicted coding sequence to place
                            prec = None

                        if prec is not None:
                            total = 0.
                            for placement in prec.placements:
                                node = ino2node[placement.edge_num]
                                weight = placement.like_weight_ratio
                                for label in node.labels:
                                    r[label]+=weight
                                if node.labels:
                                    total+=weight
                            # If there is not enough for 90% confidence prediction add OTHER classifier to signify unlabelled
                            if total<0.9:
                                r['OTHER']=1.0-total
                            # Eliminate any classifier, which is totally subsumed by another (i.e. child clade has same weight sum as parent)
                            rlabels = set(r)
                            for rlabel in rlabels:
                                for parent in label_parents[rlabel]:
                                    if r[parent]<=r[rlabel]:
                                        del r[parent]

                            #print(seqid, r)
                            #print(ino, node.labels)
                            rls = sorted(r.items(), key=lambda x: x[1], reverse=True)
                            #print(rls)
                            labels = ','.join([x[0] for x in rls])
                            weights = ','.join( ["%.1f" % (x[1]*100) for x in rls])

                            drec = diamond_recs[seqid]
                            frec = fardb_recs.get(seqid,None) or StructObject(hit=None, evalue=None, percid=None, stitle=None)
                            fields = [ seqid, all_counts[seqid], all_seqids[seqid], labels, weights, drec.hit, drec.evalue, "%.2f" % drec.percid, frec.hit, frec.evalue, "%.2f" % frec.percid if frec.percid else None, frec.stitle ]
                        else:
                            try:
                                drec = diamond_recs[seqid]
                                fields = [seqid, all_counts[seqid], all_seqids[seqid], 'OUTLIER_NOCDS', '*',  drec.hit, drec.evalue, "%.2f" % drec.percid, frec.hit, frec.evalue, "%.2f" % frec.percid if frec.percid else None, frec.stitle ]
                            except:
                                fields = [seqid, all_counts[seqid], all_seqids[seqid], 'OUTLIER_NOCDS', '*',  None, None, None, frec.hit, frec.evalue, "%.2f" % frec.percid if frec.percid else None, frec.stitle ]
                    else:
                        frec = fardb_recs.get(seqid,None) or StructObject(hit=None, evalue=None, percid=None, stitle=None)
                        fields = [seqid, all_counts[seqid], all_seqids[seqid], 'OUTLIER', '*', None, None, None, frec.hit, frec.evalue, "%.2f" % frec.percid if frec.percid else None, frec.stitle ]
                    print(*fields, sep='\t', file=out_fh)
            ## Generate output tree including annotations (use Phylotree)
            subprocess.call("python /scripts/ete3_visualizeAllTrees.py {out_dirn} {annots_fn} {clades_fn}".format(out_dirn=args.out_dirn, annots_fn=args.annots_fn, clades_fn=args.clades_fn), shell=True)
            #sys.exit(0)
            if 'CUSERID' in os.environ and 'CGROUPID' in os.environ:
                os.chown(args.out_dirn, int(os.environ['CUSERID']), int(os.environ['CGROUPID']))
            

            
