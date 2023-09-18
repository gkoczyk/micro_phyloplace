import ete3
import ete3.treeview
from querier3.utils import *

out_dir = sys.argv[1]
GTREE_FN = os.path.join(out_dir, 'epa_result.newick')
ANNOT_FN = sys.argv[2]
OUTPUT_ANNOTS_FN = os.path.join(out_dir, 'summary.tsv')
CLADES_FN = sys.argv[3]
if CLADES_FN=='None':
    CLADES_FN=None
WIDTH=283


gtree = ete3.PhyloTree(GTREE_FN)
gtree.ladderize()
# Read in annotations
annots = defaultdict(list)
for v in parseTsvList(ANNOT_FN):
    annots[v.master_seqid].append(v)

#annots = { v.seqid:v for v in parseTsvList(ANNOT_FN) }
output_annots = {v.amplicon_id: v for v in parseTsvList(OUTPUT_ANNOTS_FN) }
if CLADES_FN:
    clades = {v.clade_name:v for v in parseTsvList(CLADES_FN) }
else:
    clades = {}
ts = ete3.treeview.TreeStyle()
ts.scale = 2000 
ts.show_branch_support=False
ts.show_leaf_name=False
ts.min_leaf_separation=130
ts.draw_guiding_lines=True

MAJOR_SEP=8
MINOR_SEP=4
def nodeLayoutFunc(*args, **kwargs):
    default_position = 'branch-right'
    node = args[0]
    s = node.img_style = ete3.NodeStyle(hz_line_width=20, vt_line_width=20, hz_line_color='Black', vt_line_color='Black', size=19, shape='square', fgcolor='Black', bgcolor='white')
    #if gtree.get_common_ancestor('colletotrichum_graminicola__gene___GLRG_11770', 'ilyonectria_sp.__estExt_Genewise1.C_7_t50272') in node.get_ancestors():
    #        s['hz_line_color']=s['vt_line_color'] = s['fgcolor'] = 'SaddleBrown'
    if node.is_leaf():        
        if node.name in annots or node.name not in output_annots:
            annotvs = annots.get(node.name, [ StructObject(label='unannotated', taxon='-', klass='-', order='-', family='-') ] )
            if 11:   
                    colno=1
                    position = default_position
                    try:
                        main, domno = node.name.rsplit('__', 1)
                        domno=int(domno)
                    except:
                        main = node.name
                        domno = None
                    for i, annotv in enumerate(annotvs, start=1):
                        label = annotv.label
                        taxon = annotv.taxon
                        nameFace = ete3.TextFace(" %s" % label, fsize=35, fgcolor='black')
                        nameFace.margin_top=MAJOR_SEP
                        nameFace.margin_bottom=MINOR_SEP
                        ete3.faces.add_face_to_node(nameFace, node, colno, position=position)
                        # Taxonomic info
                        txs = [ str(v) for v in [annotv.klass, annotv.order, annotv.family, taxon] if v is not None ] 
                        taxstr = '   ' + '; '.join( str(x) for x in txs)
                        taxFace = ete3.TextFace(taxstr, fsize=30, fgcolor='black', fstyle='italic')
                        taxFace.margin_top=MINOR_SEP
                        taxFace.margin_bottom=MAJOR_SEP
                        ete3.faces.add_face_to_node(taxFace, node, colno, position=position)
                        if i==1:
                            nameFace.margin_top+=5
                        elif i>=len(annotvs):
                            taxFace.margin_bottom+=5
        else: #node.name in amplicon_seqids:
            nameFace = ete3.TextFace(" " + node.name, fsize=35, fgcolor='Red')
            annotv = output_annots[node.name]
            nameFace.margin_top=MAJOR_SEP
            nameFace.margin_bottom=MINOR_SEP
            ete3.faces.add_face_to_node(nameFace, node, 1, position=default_position)#), position="aligned")
            #desctxt = "count=%s; length=%s; placement:weight=%s:%s; sprot_hit:perc_ident=%s:%s;" % (annotv.count, annotv.length, annotv.labels, annotv.weights, annotv.fardb_hit_id, 
            #                                                                                          annotv.fardb_percid)
            if annotv is not None:
                desctxt = ("%s sequences; "% annotv.count if annotv.count>1 else '')
                desctxt += "" + "/".join( '%s(%s)' % (x,y) for x,y in zip(annotv.labels.split(','), str(annotv.weights).split(',')))
                descFace = ete3.TextFace( " " + desctxt, fsize=30, fgcolor='Red', fstyle='italic')
                descFace.margin_top=MINOR_SEP
                descFace.margin_bottom=MAJOR_SEP
            ete3.faces.add_face_to_node(descFace, node, 1, position=default_position)
            #ete3.faces.add_face_to_node(nameFace, node, 1, position=default_position)
            #s['bgcolor']='Pink'
            #node.dist=0.5
            #nameFace.background.color = 'Pink'
            #nameFace.margin_left=30
            #nameFace.margin_right=2
            #nameFace.margin_top=5
            #nameFace.margin_bottom=5
    if not any(leaf.name in annots for leaf in node): #all(leaf.name in amplicon_seqids for leaf in node):
        s['vt_line_type']=0
        s['hz_line_type']=0
        s['hz_line_color']='Red'
        s['vt_line_color']='Red'
        s['fgcolor']='Red'
    else:
        if not node.is_leaf():
            if node.support>0 and node.support<=1:
                node.support*=100
            if node.support>=100:
                supstr = '\n\n*'
            else:
                supstr = '\n\n%d' % node.support
            supportFace = ete3.TextFace(supstr, fsize=32, fgcolor='black')
            ete3.faces.add_face_to_node(supportFace, node, 0, position="float")
            
ts.layout_fn = layout_fn = nodeLayoutFunc

#print(annots['sp|Q03131|ERYA1_SACER__2'])

for clade_name in sorted(clades):
    clade = clades[clade_name]
    tip1 = clade.tip1
    tip2 = clade.tip2
    clade_tree = gtree.get_common_ancestor(tip1, tip2) #'colletotrichum_graminicola__gene___GLRG_11770', 'ilyonectria_sp.__estExt_Genewise1.C_7_t50272')
    nst = ete3.NodeStyle()
    nst['bgcolor']='White'
    clade_tree.set_style(nst)
    x = clade_tree.render( os.path.join(out_dir, clade_name+".svg"), w=WIDTH, units="mm", tree_style=ts)

nst = ete3.NodeStyle()
nst['bgcolor']='White'
gtree.set_style(nst)
x = gtree.render( os.path.join(out_dir, 'ALL_PLACEMENTS'+".svg"), w=WIDTH, units="mm", tree_style=ts)

