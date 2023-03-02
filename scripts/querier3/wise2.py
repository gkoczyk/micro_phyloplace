# Number of functions to help with Wise2 package output

from querier3.utils import *
# Separator for multigenes, beyond first returned result
WISE2_MULTIGENE_SEP='####'

# Reading GeneWise struct for single gene situation (i.e. only one gene expected per query)
@handlify(mode='r')
def readWiseGeneStructsForSingle(fh):
    WAITING=0
    ACCEPTING=1
    state = WAITING
    seen = set()
    for line in fh:
        if line.startswith('>Results'):
            state=ACCEPTING
            seqid = line.split()[4]
            if seqid in seen:
                STATE=WAITING
                continue
            seen.add(seqid)
            rec = StructObject(seqid=seqid, exons=[])
        if state==ACCEPTING:
            if line.startswith('//'):
                if rec:
                    if rec.exons and rec.exons[0].frame<0:
                        rec.exons.sort(key=lambda v:-v.xstart)
                    yield rec
                state=WAITING
            else:
                line = line.strip()
                if line.startswith('Exon'):
                    fields = line.split()
                    xstart=int(fields[1])
                    xend=int(fields[2])
                    frame=int(fields[4])+1
                    if xstart>xend:
                        xstart,xend=xend,xstart
                        frame*=-1
                    rec.exons.append(StructObject(xstart=xstart, xend=xend, frame=frame))

                    # Reading GeneWise struct for single gene situation (i.e. only one gene expected)
@handlify(mode='r')
def readWiseGeneStructsForMulti(fh):
    WAITING=0
    ACCEPTING=1
    state = WAITING
    for line in fh:
        if line.startswith('>Results'):
            state=ACCEPTING
            seqid = line.split()[4]
            rec = StructObject(seqid=seqid, exons=[])
        if state==ACCEPTING:
            if line.startswith('//'):
                if rec:
                    if rec.exons and rec.exons[0].frame<0:
                        rec.exons.sort(key=lambda v:-v.xstart)
                    yield rec
                state=WAITING
            if line.startswith('Gene'):
                fields = line.strip().split()
                if len(fields)==2:
                    gno = int(fields[1])
                    if gno>1:
                        if rec:
                            if rec.exons and rec.exons[0].frame<0:
                                rec.exons.sort(key=lambda v:-v.xstart)
                            yield rec
                        rec = StructObject(seqid=seqid+WISE2_MULTIGENE_SEP+"%d"%gno, exons=[])
            else:
                line = line.strip()
                if line.startswith('Exon'):
                    fields = line.split()
                    xstart=int(fields[1])
                    xend=int(fields[2])
                    frame=int(fields[4])+1
                    if xstart>xend:
                        xstart,xend=xend,xstart
                        frame*=-1
                    rec.exons.append(StructObject(xstart=xstart, xend=xend, frame=frame))
