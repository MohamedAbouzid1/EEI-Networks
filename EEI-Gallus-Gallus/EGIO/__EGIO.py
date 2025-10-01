import argparse
import string
import ctypes
import os
import copy
import math
import pandas as pd
import sys
from multiprocessing import Pool
##======================================================##
##                      basic function                  ##
##======================================================##
align= ctypes.cdll.LoadLibrary(os.getcwd() + '/___pairwisealign.so')
align.pairwisealign.argtypes = [ctypes.c_char_p]
align.pairwisealign.restype = ctypes.c_char_p

def pairwisealign(seqref,seqtst,method):
    ## pairalignment: all is global alignment, but with different algorithm
    ## 1. local: Smith-Waterman
    ## 2. global: Needleman-Wunsch
    ## 3. glocal: based on Smith-Waterman, and fill the the whole sequence with gap
    ## 3. gap optimise was added to correct gap location caused by penalty score matrix
    if not (method in ["local","global","glocal"]):
        raise SyntaxError("pairwisealign: no available method called " + method)
    if len(seqref) == 0 or len(seqtst) == 0:
        raise SyntaxError("pairwisealign: illegal sequence (no sequence at all)")
    
    if method == "global":
        gapopen = -10
    elif method == "glocal":
        gapopen = -10
    elif method == "local":
        gapopen = -10
    ##=====================================##
    ##====== local paiewise alignment =====##
    ##=====================================##
    if isinstance(seqref,list):
        seqref = "".join(seqref)
    if isinstance(seqtst,list):
        seqtst = "".join(seqtst)
    ##------------------------------------
    ##------------------------------------
    hsseqlo = ctypes.c_char_p(bytes(seqref, 'utf-8'))
    ttseqlo = ctypes.c_char_p(bytes(seqtst, 'utf-8'))
    sumalignlotmp = align.pairwisealign(hsseqlo,ttseqlo,method,gapopen)
    sumalignlo = sumalignlotmp.decode().split("|")
    lohs = sumalignlo[0]
    lott = sumalignlo[1]
    ##=====================================##
    ##========= correct alignment =========##
    ##=====================================##
    listsum = {"seq1":lohs,"seq2":lott}
    return listsum

def calIdentity(seqref, seqtst, method):
    ## calculate sequence identity
    if not (method in ["local","global"]):
        raise SyntaxError("calIdentity: no available method called " + method)
    if len(seqref) == 0 or len(seqtst) == 0:
        raise SyntaxError("calIdentity: illegal sequence (no sequence at all)")    
    
    lenseq1 = len(seqref)
    lenseq2 = len(seqtst)
    seqalign = pairwisealign(seqref, seqtst, method)
    seq1 = seqalign["seq1"]
    seq2 = seqalign["seq2"]
    ##============== count exact match ================
    matchcount = 0
    for i in range(0,len(seq1)):
        if seq1[i] == seq2[i]:
            matchcount = matchcount + 1
    ##=================== cauculate ===================
    identmp = round(matchcount/len(seq1),2) if len(seq1) > 0 else 0
    iden = [identmp,round(len(seq1.replace("-",""))/lenseq1,2),round(len(seq2.replace("-",""))/lenseq2,2)]
    
    return iden

def scorealignment(seqa,seqb):
    ## this is a score function, in which continuous gap
    ## will be less punished as a single gap, and the decrease 
    ## is caculated by a decrease factor: alpha = 0.8
    ## match: 10
    ## mismatch: -10
    ## gap: -8
    ## continuous gap[i]: gap[i-1]*0.8
    alpha = 0.8

    if isinstance(seqa,str):
        seqa = list(seqa)
    if isinstance(seqb,str):
        seqb = list(seqb)
    score = 0
    gapscore = -8
    if seqa[0] == seqb[0]:
        score = 10
    else:
        score = -10
    for i in range(1,len(seqa)):
        if seqa[i] == seqb[i]:
            score = score + 10
        else:
            if seqa[i] != "-" and seqb[i] != "-":
                score = score - 10    
            elif seqa[i] == "-" and seqb[i] != "-":
                if seqa[i-1] == "-":
                    gapscore = gapscore*alpha
                else:
                    gapscore = -8
                score = score + gapscore
            elif seqa[i] != "-" and seqb[i] == "-":
                if seqb[i-1] == "-":
                    gapscore = gapscore*alpha
                else:
                    gapscore = -8
                score = score + gapscore
    return round(score,3)

def listwhich(lsdata,deal,value):
    if deal == "==":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] == value]
    if deal == "!=":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] != value]
    if deal == ">":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] > value]
    if deal == ">=":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] >= value]
    if deal == "<":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] < value]
    if deal == "<=":
        judge = [x for x in range(0,len(lsdata)) if lsdata[x] <= value]         
    return(judge)

def intersect(ls1, ls2):
    store = [x for x in ls1 if x in ls2]
    return store

def unique(listdat):
    store = []
    if len(listdat) > 0:
        for i in range(0,len(listdat)):
            if not(listdat[i] in store):
                store.append(listdat[i])
    store.sort(reverse = False)
    return store

def union(ls1, ls2):
    if not(isinstance(ls1,list)) and not(isinstance(ls2,list)):
        unionstore = unqiue([ls1, ls2])
    elif isinstance(ls1,list) and not(isinstance(ls2,list)):
        unionstore = unique([ls2]+ls1)
    elif isinstance(ls2,list) and not(isinstance(ls1,list)):
        unionstore = unique([ls1]+ls2)
    elif isinstance(ls2,list) and isinstance(ls1,list):
        unionstore = unique(ls2+ls1)
    return unionstore

def typeinto(listd,colname,astype=int):
    coli = listd["coln"].index(colname)
    for i in range(0,len(listd["dat"])):
        listd["dat"][i][coli] = astype(listd["dat"][i][coli])
    return listd

def selectcol(listd,target):
    seli = listd["coln"].index(target)
    if not(target in listd["coln"]):
        raise SyntaxError("selectcol: wrong colname is given, please check")
    store =[listd["dat"][x][seli] for x in range(0,len(listd["dat"]))]
    return store

def listasint(listd,deal="min"):
    ## only be used in homo-colinearity test
    ## transform list in list to int
    if not(isinstance(listd,list)):
        store = [listd]
    else:
        store = []
        for i in range(0,len(listd)):
            if isinstance(listd[i],list):
                if deal == "min":
                    store.append(min(listd[i]))
                if deal == "max":
                    store.append(max(listd[i]))
                if deal == "all":
                    for j in range(0,len(listd[i])):
                        store.append(listd[i][j])
            else:
                store.append(listd[i])
    return store

def selectele(lsdata,targeti):
    store = [lsdata[x] for x in targeti]
    return store

def selectrow(listd,targeti):
    store = [listd["dat"][x] for x in targeti]
    return {"dat":store,"coln":listd["coln"]}

def selectlsls(listd,rowi,coln="all"):
    if coln == "all":
        coln = listd["coln"]
    coli = [listd["coln"].index(x) for x in coln]
    store = [[listd["dat"][y][x] for x in coli] for y in rowi]
    listsum = {"dat":store,"coln":coln}
    return listsum

def listcal(lsdata,calcu,value):
    if calcu == "+":
        listtmp = [x+value for x in lsdata]
    if calcu == "-":
        listtmp = [x-value for x in lsdata]
    if calcu == "*":
        listtmp = [x*value for x in lsdata]
    if calcu == "/":
        if value == 0:
            raise SyntaxError("listcal: 0 can not be divided")
        listtmp = [x/value for x in lsdata]
    if calcu == "%":
        if value == 0:
            raise SyntaxError("listcal: 0 can not be divided")
        listtmp = [x%value for x in lsdata]
    if calcu == "^":
        listtmp = [x^value for x in lsdata]
    return listtmp

def listlistcal(ls1,calcu,ls2):
    if len(ls1) != len(ls2):
        raise SyntaxError("listlistcal: the two lists do not have the same length")
    
    if calcu == "+":
        listtmp = [ls1[x]+ls2[x] for x in range(0,len(ls1))]
    if calcu == "-":
        listtmp = [ls1[x]-ls2[x] for x in range(0,len(ls1))]
    if calcu == "*":
        listtmp = [ls1[x]*ls2[x] for x in range(0,len(ls1))]
    if calcu == "/":
        if 0 in ls2:
            raise SyntaxError("listlistcal: 0 can not be divided")
        listtmp = [ls1[x]/ls2[x] for x in range(0,len(ls1))]
    if calcu == "%":
        if 0 in ls2:
            raise SyntaxError("listlistcal: 0 can not be divided")
        listtmp = [ls1[x]%ls2[x] for x in range(0,len(ls1))]
    if calcu == "^":
        listtmp = [ls1[x]^ls2[x] for x in range(0,len(ls1))]
    return listtmp

def listdrop(lsdata,dropi):
    store = [lsdata[x] for x in range(0,len(lsdata)) if not(x in dropi)]
    return store

def listddrop(listd,dropi):
    store = [listd['dat'][x] for x in range(0,len(listd['dat'])) if not(x in dropi)]
    return {"dat":store,"coln":listd["coln"]}

def listdinsert(listd,insi,element):
    if len(listd["dat"]) == 0:
        listd["dat"].insert(insi,element)
    else:
        if len(element) != len(listd["dat"][0]):
            raise SyntaxError("listdinsert: insert element do not map to the data dimension (column)")
        listd["dat"].insert(insi,element)
    return {"dat":listd["dat"],"coln":listd["coln"]}

def transpose(lsls):
    new_list = [[lsls[j][i] for j in range(0,len(lsls))] for i in range(0,len(lsls[0]))]
    return new_list

def listsort(listd,by,increase):
    if len(by) > 2:
        raise SyntaxError("listsort: the max element to sort is 2")
    listdat = listd["dat"]
    if len(by) == 1:
        by1 = listd["coln"].index(by[0])
        def takeone(elem):
            return elem[by1]
        if increase:
            listdat.sort(key =takeone, reverse = False)
        else:
            listdat.sort(key =takeone, reverse = True)
    
    else:
        by1 = listd["coln"].index(by[0])
        by2 = listd["coln"].index(by[1])
        def taketwo(elem):
            return elem[by1], elem[by2]
        if increase:
            listdat.sort(key =taketwo, reverse = False)
        else:
            listdat.sort(key =taketwo, reverse = True)    
    listsum = {"dat":listdat,"coln":listd["coln"]}
    return listsum

def setdiff(ls1,ls2):
    difflist=[ls1[x] for x in range(0,len(ls1)) if not(ls1[x] in ls2)]
    return difflist

def DFtoDict(dataframe, key_base, header=True):
    colname = list(dataframe.columns)

    if header:
        if not(key_base in colname):
            SyntaxError("DFtoDict: no column name called: " + key_base)

        keyi = colname.index(key_base)
    else:
        keyi = key_base
    
    print("transform data frame into dictionary")
    store = []
    storetmp = []
    keytmp = []
    for i in range(0, len(dataframe)):
        if not(dataframe.iloc[i].tolist()[keyi] in keytmp):
            if len(storetmp) > 0:
                store.append(storetmp)
            storetmp = []
            keytmp.append(dataframe.iloc[i].tolist()[keyi])
            storetmp.append(dataframe.iloc[i].tolist())
        else:
            storetmp.append(dataframe.iloc[i].tolist())
        
        if i % 50000 == 0 and i > 0:
            print(str(i) + " records are analyzed")
    print(str(i) + " records are analyzed")
    print("done!!")

    if len(storetmp) > 0:
        store.append(storetmp)        
    
    cdsDt = {}
    for k, v in zip(keytmp,store):
        cdsDt[k] = v

    return cdsDt

##======================================================##
##                  translate from cDNA                 ##
##======================================================##
def translate(cdna):
    nt_to_aa = {'AAA':'K','AAT':'N','AAG':'K','AAC':'N','ATA':'I','ATT':'I','ATG':'M','ATC':'I','AGA':'R','AGT':'S','AGG':'R','AGC':'S','ACA':'T','ACT':'T','ACG':'T','ACC':'T','TAA':'','TAT':'Y','TAG':'','TAC':'Y','TTA':'L','TTT':'F','TTG':'L','TTC':'F','TGA':'','TGT':'C','TGG':'W','TGC':'C','TCA':'S','TCT':'S','TCG':'S','TCC':'S','GAA':'E','GAT':'D','GAG':'E','GAC':'D','GTA':'V','GTT':'V','GTG':'V','GTC':'V','GGA':'G','GGT':'G','GGG':'G','GGC':'G','GCA':'A','GCT':'A','GCG':'A','GCC':'A','CAA':'Q','CAT':'H','CAG':'Q','CAC':'H','CTA':'L','CTT':'L','CTG':'L','CTC':'L','CGA':'R','CGT':'R','CGG':'R','CGC':'R','CCA':'P','CCT':'P','CCG':'P','CCC':'P'}
    cdnalen = len(cdna)
    proseq = ""
    for i in range(0,cdnalen // 3):
        try:
            aminoa = nt_to_aa[cdna[i*3:(i*3+3)]]
        except:
            aminoa = "X"
        proseq = proseq + aminoa
    return proseq

def extract_CDSMatrix_from_Exon(orfpos,exoncom,strtmp):
    ## extract CDS component from exon component
    orfst, orfen = orfpos.split("|")
    orfst, orfen = int(orfst), int(orfen)
    if orfst < 0:
        orfst = 1

    cds_pos = exoncom.split("|")
    cds_pos_int = [0 for i in range(0,len(cds_pos))]
    cdslen = [0 for i in range(0,len(cds_pos))]
    cdssum = [0 for i in range(0,len(cds_pos))]
    for i in range(0,len(cds_pos)):
        cds_pos_int[i] = [int(cds_pos[i].split(":")[0]), int(cds_pos[i].split(":")[1])]
        cdslen[i] = cds_pos_int[i][1] - cds_pos_int[i][0] +1
        cdssum[i] = sum(cdslen)
    ##==============
    stidx = min(listwhich(cdssum,">=",orfst))
    deltst = cdssum[stidx] - orfst
    ##--------------
    enidx = min(listwhich(cdssum,">=",orfen))
    delten = cdssum[enidx] - orfen
    ##--------------
    if strtmp == "+":
        cds_pos_int[stidx][0] = cds_pos_int[stidx][1] - deltst
        cds_pos_int[enidx][1] = cds_pos_int[enidx][1] - delten
    if strtmp == "-":
        cds_pos_int[stidx][1] = cds_pos_int[stidx][0] + deltst
        cds_pos_int[enidx][0] = cds_pos_int[enidx][0] + delten
    
    CDS = cds_pos_int[stidx:(enidx+1)]
    return CDS

def tranform_matrix_str(cdscom):
    ## transform matrix to string
    cdsstrtmp = ["." for _ in range(0,len(cdscom))]
    for i in range(0,len(cdscom)):
        cdsstrtmp[i] = str(cdscom[i][0])+":"+str(cdscom[i][1])
    return "|".join(cdsstrtmp)

def extract_pro_seq(cdscom,cdshscds,orfpos):
    ## translate CDS matrix into AA sequence
    orfst, orfen = orfpos.split("|")
    orfst, orfen = int(orfst), int(orfen)
    seq = ""
    if orfst < 0:
        seq = seq + "".join(["N" for _ in range(0,abs(orfst))])
    ##----------------------------    
    for i in range(0,len(cdscom)):
        #idx = intersect(cdshscds[(cdshscds.Start==cdscom[i][0])].index.tolist(),cdshscds[(cdshscds.End==cdscom[i][1])].index.tolist())
        idx = intersect(listwhich(selectcol(cdshscds,"start"),"==",cdscom[i][0]),listwhich(selectcol(cdshscds,"end"),"==",cdscom[i][1]))
        #seq = seq + cdshscds.loc[idx[0]]["Seq"]
        seq = seq + selectcol(cdshscds,"seq")[idx[0]]
    ##----------------------------
    aaseq = translate(seq)
    return aaseq

def translatecdna(genhscdsdat,cdshscds,strandhsstr):
    ## translate exon com matrix to cds com matrix
    exonnew = []
    for i in range(0,len(genhscdsdat["dat"])):
        exoncom = selectcol(genhscdsdat,"Exoncom")[i]
        orfpos = selectcol(genhscdsdat,"Orf")[i]
        cdscom = extract_CDSMatrix_from_Exon(orfpos,exoncom,strandhsstr)
        proseq = extract_pro_seq(cdscom,cdshscds,orfpos)
        exonnew.append([selectcol(genhscdsdat,"EnsemblT")[i],orfpos,tranform_matrix_str(cdscom),proseq])
    return {"dat":exonnew,"coln":["EnsemblT","Orf","CDS","AAseq"]}

##======================================================##
##                   group CDS-exons ***                ##
##======================================================##
def mergeOverlap(datad,grpstr):
    ## organize exon group based on exon region
    if grpstr == 1:
        datad = listsort(datad,by=["start","end"],increase = True)
        sti = listwhich(datad["coln"],"==","start")[0]
        eni = listwhich(datad["coln"],"==","end")[0]
        sqi = listwhich(datad["coln"],"==","seq")[0]
        stff = datad["dat"][0][sti]
        enff = datad["dat"][0][eni]
        seqff = datad["dat"][0][sqi]
        for i in range(1,len(datad["dat"])):
            sttmp = datad["dat"][i][sti]
            entmp = datad["dat"][i][eni]
            seqtmp = datad["dat"][i][sqi]
            if enff >= sttmp and enff < entmp:
                seq_ahead_tmp = list(seqff)[0:(sttmp - stff)]
                seqff = "".join(seq_ahead_tmp) + seqtmp
                stff = min(stff,sttmp)
                enff = max(enff, entmp)
    if grpstr == -1:
        datad = listsort(datad,by=["end","start"],increase = False)
        ## strand = -1, the large position is small
        sti = listwhich(datad["coln"],"==","start")[0]
        eni = listwhich(datad["coln"],"==","end")[0]
        sqi = listwhich(datad["coln"],"==","seq")[0]
        stff = datad["dat"][0][sti]
        enff = datad["dat"][0][eni]
        seqff = datad["dat"][0][sqi]
        for i in range(1,len(datad["dat"])):
            sttmp = datad["dat"][i][sti]
            entmp = datad["dat"][i][eni]
            seqtmp = datad["dat"][i][sqi]
            if stff > sttmp and stff <= entmp:
                seq_ahead_tmp = list(seqff)[0:(enff-entmp)]
                seqff = "".join(seq_ahead_tmp) + seqtmp
                stff = min(stff,sttmp)
                enff = max(enff, entmp)

    return seqff

def regRuler(regidxd,genstr):
    ## estimate exon ruler
    regidxd = typeinto(regidxd,"start",int)
    regidxd = typeinto(regidxd,"end",int)
    regname = regidxd["coln"]
    if genstr == 1:
        regcomd = listsort(regidxd,by=["start","end"],increase = True)
    if genstr == -1:
        regcomd = listsort(regidxd,by=["end","start"],increase = False)
    ## get group of each region
    group = []
    gtmp = 1
    sti = listwhich(regname,"==","start")[0]
    eni = listwhich(regname,"==","end")[0]
    stref = regcomd["dat"][0][sti]
    enref = regcomd["dat"][0][eni]
    for i in range(0,len(regcomd["dat"])):
        sttst = regcomd["dat"][i][sti]
        entst = regcomd["dat"][i][eni]
        if entst < stref or sttst > enref:
            gtmp = gtmp + 1
            group.append(gtmp)
            stref = sttst
            enref = entst
        else:
            group.append(gtmp)
            stref = min(stref, sttst)
            enref = max(enref,entst)
    ## summarize each group
    gid = unique(group)
    gsum = []
    starttmp = selectcol(regcomd,"start")
    endtmp = selectcol(regcomd,"end")
    #regseq = selectcol(regcomd,"seq")
    for i in range(0,len(gid)):
        rgri = listwhich(group,"==",gid[i])
        gst = min(selectele(starttmp,rgri))
        gen = max(selectele(endtmp,rgri))
        if len(rgri) == 1:
            gsum.append([gid[i],gst,gen,"nonunited"])
        else:
            gsum.append([gid[i],gst,gen,"united"])
    gsumname = ["group","start","end","etype"]
    gsumd = {"dat":gsum,"coln":gsumname}
    ## calculate the aligned position
    glentmp = listlistcal(selectcol(gsumd,"end"),"-",selectcol(gsumd,"start"))
    glen = listcal(glentmp,"+",1)
    ginst = [[] for _ in range(0,len(glen))]
    ginst[0] = 1
    ginen = [[] for _ in range(0,len(glen))]
    ginen[0]= glen[0]
    for i in range(1,len(glen)):
        ginst[i] = ginen[i-1] + 1
        ginen[i] = ginen[i-1] + glen[i]
    posinfo=[]
    for i in range(0,len(glen)):
        posinfo.append([glen[i],ginst[i],ginen[i]])
    posinfod = {"dat":posinfo,"coln":["len","instart","inend"]}    
    ## get represent region id and sequence
    gsti = listwhich(gsumname,"==","start")[0]
    geni = listwhich(gsumname,"==","end")[0]
    regi = [[] for _ in list(range(0,len(gsum)))]
    seqtmp = [[] for _ in list(range(0,len(gsum)))]
    regseqi = listwhich(regname,"==","seq")[0]
    #-----------------------
    for i in range(0,len(gsum)):
        j1 = listwhich(starttmp, ">=", gsum[i][gsti])
        j2 = listwhich(endtmp, "<=", gsum[i][geni])
        j = intersect(j1,j2)
        regi[i] = str(gsum[i][gsti])+"-"+str(gsum[i][geni])
        if len(j) == 1:
            seqtmp[i] = regcomd["dat"][j[0]][regseqi]
        else:
            # look for new start, end of the first region
            grpseqtmp = []
            for ji in range(0,len(j)):
                grpseqtmp.append(regcomd["dat"][j[ji]])
            grpseqtmpd = {"dat":grpseqtmp, "coln":["id","start","end","seq"]}
            seqmerge = mergeOverlap(datad = grpseqtmpd, grpstr = genstr)
            seqtmp[i] = seqmerge
    #-----------------------
    addinfo = []
    for i in range(0,len(regi)):
        addinfo.append([regi[i],seqtmp[i]])
    addinfod = {"dat":addinfo,"coln":["id","seq"]}
    ## get final region ruler
    gsum_t = transpose(gsumd["dat"])
    posinfo_t = transpose(posinfod["dat"])
    addinfo_t = transpose(addinfod["dat"])
    #-----------------------
    ruler_t = gsum_t + posinfo_t + addinfo_t
    ruler = transpose(ruler_t)
    rulername = gsumd["coln"] + posinfod["coln"] + addinfod["coln"]
    rulerd = {"dat":ruler,"coln":rulername}
    return rulerd

def findgrouppos(ex_len, seqlist, gapingroup = False):
    ## find group pos based on group length and aligned sequence
    if not(gapingroup): ## gap should not be included into neither of the groups
        ecrange = []
        count = 0
        posidx = -1
        pos_st = 0
        ex_idx = 0
        startgrp = True
        for i in range(0,len(seqlist)):
            posidx = posidx + 1
            if seqlist[i] != "-":
                if startgrp:
                    pos_st = posidx
                    startgrp = False
                count = count + 1
                if count == ex_len[ex_idx]:
                    pos_en = posidx
                    ecrange.append([pos_st,pos_en])
                    count = 0
                    startgrp = True
                    ex_idx = ex_idx + 1
        if len(ecrange) != len(ex_len) or ecrange[len(ex_len)-1][1]>=len(seqlist):
            raise SyntaxError("findgrouppos: nucleotide loss during group range screening (no gap)")
    else: ## inter-group gap should be included into the second group, which is required in [optimiseLocal] function
        ecrange = []
        count = 0
        posidx = -1
        pos_st = 0
        ex_idx = 0
        for i in range(0,len(seqlist)):
            posidx = posidx + 1
            if seqlist[i] != "-":
                count = count + 1
                if count == ex_len[ex_idx]:
                    pos_en = posidx
                    ecrange.append([pos_st,pos_en])
                    count = 0
                    pos_st = pos_en + 1
                    ex_idx = ex_idx + 1
                    #break
        if ecrange[len(ex_len)-1][1] > (len(seqlist)-1):
            raise SyntaxError("findgrouppos: nucleotide loss during group range screening (with gap)")
        else: ## confirm bottom gap is included into the group
            ecrange[len(ex_len)-1][1] = len(seqlist)-1
    ##=====================================##
    ##======= get ruler-exon region =======##
    ##=====================================##

    return ecrange

def matchhomoexon(homopair,regr,regt,identhres):
    seqi = listwhich(regr["coln"],"==","seq")[0]
    ## match homo exons
    if len(homopair) == 0:
        newmatch = []
    if len(homopair) == 1:
        newmatch = [[[homopair[0][0]],[homopair[0][1]]]]
    else:
        ## species1 as anchor
        newmatchhsa = []
        sttmp = [homopair[0][0]]
        entmp = [homopair[0][1]]
        for i in range(1,len(homopair)):
            if homopair[i][0] != homopair[i-1][0]:
                newmatchhsa.append([sttmp,entmp])
                sttmp = [homopair[i][0]]
                entmp = [homopair[i][1]]
            else:
                sttmp = union(sttmp,homopair[i][0])
                entmp = union(entmp,homopair[i][1])
        newmatchhsa.append([sttmp,entmp])
        ###===  add a filter to deal with multiple duplicated exons  ===##
        if len(newmatchhsa) == 1:
            newmatchtmp = newmatchhsa
        else:
            newmatchtmp = []
            maxttgrp = 0
            for i in range(0,len(newmatchhsa)):
                if len(newmatchhsa[i][1]) == 1:
                    newmatchtmp.append(newmatchhsa[i])
                    maxttgrp = newmatchhsa[i][1][0]
                else:
                    idxtmp = listwhich(newmatchhsa[i][1],">",maxttgrp)
                    grptttmp = selectele(newmatchhsa[i][1],idxtmp)
                    if len(grptttmp) == 1:
                        newmatchtmp.append([newmatchhsa[i][0],grptttmp])
                        maxttgrp = grptttmp[0]
                    if len(grptttmp) > 1:
                        grpdiftt = [1]
                        for j in range(1,len(grptttmp)):
                            grpdiftt.append(grptttmp[j]-grptttmp[j-1])
                        choosettidx = listwhich(grpdiftt,">",1)
                        if len(choosettidx) > 0:
                            stopttidx = min(choosettidx)
                        else:
                            stopttidx = len(grpdiftt)
                        newmatchtmp.append([newmatchhsa[i][0],grptttmp[0:stopttidx]])
                        maxttgrp = min(grptttmp[0:stopttidx])
        ## species2 as anchor
        if len(newmatchtmp) == 1:
            newmatch1 = newmatchtmp
        else:
            newmatchtmpd = {"dat":newmatchtmp,"coln":["hsag","ptrg"]}
            newmatch1 = []
            sttmp = newmatchtmp[0][0]
            entmp = newmatchtmp[0][1]
            record = True
            for i in range(1,len(newmatchtmp)):
                if newmatchtmp[i][1] != newmatchtmp[i-1][1]:
                    if record:
                        newmatch1.append([sttmp,entmp])
                        sttmp = newmatchtmp[i][0]
                        entmp = newmatchtmp[i][1]
                    if min(newmatchtmp[i][1]) > max(listasint(selectcol(newmatchtmpd,"ptrg")[0:i],"max")):
                        record = True
                        sttmp = newmatchtmp[i][0]
                        entmp = newmatchtmp[i][1]
                    else:
                        record = False
                else:
                    sttmp = union(sttmp,newmatchtmp[i][0])
                    entmp = union(entmp,newmatchtmp[i][1])
            if record:               
                newmatch1.append([sttmp,entmp])
        ###  add a filter to deal with non-continuous multiple match    ##
        newmatch = []
        for i in range(0,len(newmatch1)):
            if len(newmatch1[i][0]) == 1 and len(newmatch1[i][1]) == 1:
                newmatch.append(newmatch1[i])
            elif len(newmatch1[i][0]) == 1 and len(newmatch1[i][1]) > 1:
                grpdiftt = [1]
                for j in range(1,len(newmatch1[i][1])):
                    grpdiftt.append(newmatch1[i][1][j]-newmatch1[i][1][j-1])
                choosettidx = listwhich(grpdiftt,">",1)
                if len(choosettidx) > 0:
                    stopttidx = min(choosettidx)
                else:
                    stopttidx = len(grpdiftt)                
                newmatch.append([newmatch1[i][0],newmatch1[i][1][0:stopttidx]])
            elif len(newmatch1[i][0]) > 1 and len(newmatch1[i][1]) == 1:
                grpdifhs = [1]
                for m in range(1,len(newmatch1[i][0])):
                    grpdifhs.append(newmatch1[i][0][m]-newmatch1[i][0][m-1])
                choosehsidx = listwhich(grpdifhs,">",1)
                if len(choosehsidx) > 0:
                    stophsidx = min(choosehsidx)
                else:
                    stophsidx = len(grpdifhs)
                newmatch.append([newmatch1[i][0][0:stophsidx],newmatch1[i][1]])
            else:
                grpdiftt = [1]
                for j in range(1,len(newmatch1[i][1])):
                    grpdiftt.append(newmatch1[i][1][j]-newmatch1[i][1][j-1])
                choosettidx = listwhich(grpdiftt,">",1)
                if len(choosettidx) > 0:
                    stopttidx = min(choosettidx)
                else:
                    stopttidx = len(grpdiftt)
                ##--------------------------------##
                grpdifhs = [1]
                for m in range(1,len(newmatch1[i][0])):
                    grpdifhs.append(newmatch1[i][0][m]-newmatch1[i][0][m-1])
                choosehsidx = listwhich(grpdifhs,">",1)
                if len(choosehsidx) > 0:
                    stophsidx = min(choosehsidx)
                else:
                    stophsidx = len(grpdifhs)
                ##--------------------------------##
                newmatch.append([newmatch1[i][0][0:stophsidx],newmatch1[i][1][0:stopttidx]])
    ##------ optimise reverse transcription mediated exon fusion ------------##
    optimisereverse = []
    for i in range(0,len(newmatch)):
        if (len(newmatch[i][0]) == 1 and len(newmatch[i][1]) > 1) or (len(newmatch[i][0]) > 1 and len(newmatch[i][1]) == 1):
            optimisereverse.append([newmatch[i][0],newmatch[i][1]])
    ##----------- final optimise ------------## colinearity test
    if len(newmatch) > 0:
        for i in range(0,len(newmatch)):
            if len(newmatch[i][0]) > 1:
                new0 = [newmatch[i][0][0]]
                for j in range(1,len(newmatch[i][0])):
                    if newmatch[i][0][j] - newmatch[i][0][j-1] == 1:
                        new0.append(newmatch[i][0][j])
                    else:
                        newmatch[i][0][j]=newmatch[i][0][j-1]
                newmatch[i][0] = new0
            if len(newmatch[i][1]) > 1:
                new1 = [newmatch[i][1][0]]
                for j in range(1,len(newmatch[i][1])):
                    if newmatch[i][1][j] - newmatch[i][1][j-1] == 1:
                        new1.append(newmatch[i][1][j])
                    else:
                        newmatch[i][1][j]=newmatch[i][1][j-1]
                newmatch[i][1] = new1
    newmatchtmp = copy.deepcopy(newmatch)
    ##---------------------------------##
    dropi = []
    for i in range(0,len(newmatch)):
        if newmatch[i] in optimisereverse:
            continue
        if newmatch[i][0] != [0] and newmatch[i][1] != [0]:
            hsseqtmp = ""
            for j in range(0,len(newmatch[i][0])):
                hsseqtmp = hsseqtmp + regr["dat"][listwhich(selectcol(regr,"group"),"==",newmatch[i][0][j])[0]][seqi]
        ##-------------------------------------##
            ttseqtmp = ""
            for j in range(0,len(newmatch[i][1])):
                ttseqtmp = ttseqtmp + regt["dat"][listwhich(selectcol(regt,"group"),"==",newmatch[i][1][j])[0]][seqi]
        if calIdentity(hsseqtmp, ttseqtmp, method="local")[0] < identhres:
            dropi.append(i)
    if len(dropi) > 0:
        newmatch = listdrop(newmatch,dropi)

    return newmatch

def alignmentexon(regr, regt, identhres=0.8, coverthres=0.8, microexon=10, score=2, penalty=-2, gap=-1):
    ##-------------------------------------##
    seqi = listwhich(regr["coln"],"==","seq")[0]
    leni = listwhich(regr["coln"],"==","len")[0]
    ##==========================================================##
    ##  pairwise alignment among exons, this step is necessary  ##
    ##  to esitmate 1-N or N-1 homologous exons                 ##
    ##==========================================================##
    seq2 = selectcol(regr,"group")
    seq2.sort(reverse = True)
    seq1 = selectcol(regt,"group")
    seq1.sort(reverse = True)
    ##-------------------------------------##
    ##       build identity matrix         ##
    ##-------------------------------------##
    exnim = [[0.0]*(len(seq2)+1) for _ in range(0,len(seq1)+1)]
    for i in range(0,len(regt["dat"])):
        ri = len(regt["dat"]) - i -1
        hsgrpseq = regt["dat"][ri][seqi]
        for j in range(0,len(regr["dat"])):
            rj = len(regr["dat"]) - j -1
            ttgrpseq = regr["dat"][rj][seqi]
            identmp = calIdentity(hsgrpseq,ttgrpseq,method="local")
            if identmp[1] >= coverthres or identmp[2] >= coverthres:
                exnim[i+1][j+1] = identmp[0]
    ##------------------------------------------##
    ##     build score and arrow matrix         ##
    ## 1. match(identity>=0.8,coverage>=0.8): 6 ##
    ## 2. mismatch: -6                          ##
    ## 3. exon jump, exon split or fusion: -3   ##
    ##------------------------------------------##
    scorem = [[0]*(len(seq2)+1) for _ in range(0,len(seq1)+1)]
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if exnim[i][j] < identhres:
                scorei = penalty  ## mismatch penalty
            else:
                scorei = score   ## match score
            smat = [0,0,0]
            smat[0] = scorem[i-1][j] + scorei + gap ## position penalty
            smat[1] = scorem[i][j-1] + scorei + gap ## position penalty
            smat[2] = scorem[i-1][j-1] + scorei
            scorem[i][j] = max(smat)
    ##==========================================================##
    ##  test whether global or local alignment                  ##
    ##==========================================================##
    if exnim[len(exnim)-1][len(exnim[0])-1] >= identhres:
        ##-------------------------------------##
        ## global alignment  
        h=len(seq1)
        v=len(seq2)
    else:
        ## local alignment
        max_score = scorem[0][0]
        for i in range(0,len(seq1)+1):
            for j in range(1,len(seq2)+1):
                if max_score <= scorem[i][j]:
                    max_score = copy.deepcopy(scorem[i][j])
                    h = i
                    v = j
    ##==========================================================##
    ##  dynamic programming to find best exon matches           ##
    ##==========================================================##
    alnseq1 = []
    alnseq2 = []
    ##-------------------------------------##
    seq1_num = h-1
    seq2_num = v-1
    if exnim[h][v] >= identhres:
        alnseq1.append(seq1[seq1_num])
        alnseq2.append(seq2[seq2_num])

    while True:
        if v == 0 or h == 0:
            break

        score0 = scorem[h-1][v]
        score1 = scorem[h][v-1]
        score2 = scorem[h-1][v-1]

        maxscore = max([score0,score1,score2])

        if score2 == maxscore:
            seq1_num -= 1
            seq2_num -= 1
            v -= 1
            h -= 1
            if exnim[h][v] >= identhres:
                alnseq1.append(seq1[seq1_num])
                alnseq2.append(seq2[seq2_num])

        elif score0 == maxscore:
            seq1_num -= 1
            h -= 1
            if exnim[h][v] >= identhres:
                alnseq1.append(seq1[seq1_num])
                alnseq2.append(seq2[seq2_num])

        elif score1 == maxscore:
            seq2_num -= 1
            v -= 1
            if exnim[h][v] >= identhres:
                alnseq1.append(seq1[seq1_num])
                alnseq2.append(seq2[seq2_num])       

    ##==========================================================##
    ##               drop microexon to avoid mismatch           ##
    ##==========================================================##
    ##==========================================================##
    ##  find microexon                                          ##
    ##==========================================================##
    grpi = listwhich(regr["coln"],"==","group")[0]
    leni = listwhich(regr["coln"],"==","len")[0]
    hsreglen = []
    hsreggrp = []
    for i in range(0,len(regr["dat"])):
        hsreggrp.append(regr["dat"][i][grpi])
        hsreglen.append(regr["dat"][i][leni])
    dropgrphs = selectele(hsreggrp,listwhich(hsreglen,"<=",microexon))
    ##------------------------------------------##
    ttreglen = []
    ttreggrp = []
    for i in range(0,len(regt["dat"])):
        ttreggrp.append(regt["dat"][i][grpi])
        ttreglen.append(regt["dat"][i][leni])
    dropgrptt = selectele(ttreggrp,listwhich(ttreglen,"<=",microexon))
    ##==========================================================##
    ##  construct possible orthologous exon relationship        ##
    ##==========================================================##
    homoexonpair = []
    for i in range(0,len(alnseq1)):
        if alnseq2[i] in dropgrphs or alnseq1[i] in dropgrptt:
            continue
        else:
            homoexonpair.append([alnseq2[i],alnseq1[i]])
    ##==========================================================##
    ## optimise multiple to one or one to multiple relationships #
    ##==========================================================##
    # species1 as anchor
    newmatchhsa = []
    if len(homoexonpair) == 1:
        newmatchhsa = [[[homoexonpair[0][0]],[homoexonpair[0][1]]]]
    elif len(homoexonpair) > 1:
        sttmp = [homoexonpair[0][0]]
        entmp = [homoexonpair[0][1]]
        for i in range(1,len(homoexonpair)):
            if homoexonpair[i][0] != homoexonpair[i-1][0]:
                newmatchhsa.append([sttmp,entmp])
                sttmp = [homoexonpair[i][0]]
                entmp = [homoexonpair[i][1]]
            else:
                sttmp = union(sttmp,homoexonpair[i][0])
                entmp = union(entmp,homoexonpair[i][1])
        newmatchhsa.append([sttmp,entmp])
    #species2 as anchor
    possiblehomo = []
    if len(newmatchhsa) == 1:
        possiblehomo = copy.deepcopy(newmatchhsa)
    elif len(newmatchhsa) > 1:  
        sttmp = newmatchhsa[0][0]
        entmp = newmatchhsa[0][1]
        for i in range(1,len(newmatchhsa)):
            if newmatchhsa[i][1] != newmatchhsa[i-1][1]:
                possiblehomo.append([sttmp,entmp])
                sttmp = newmatchhsa[i][0]
                entmp = newmatchhsa[i][1]
            else:
                sttmp = union(sttmp,newmatchhsa[i][0])
                entmp = union(entmp,newmatchhsa[i][1])
        possiblehomo.append([sttmp,entmp])
    ##==========================================================##
    ##  got final orthologous exons                             ##
    ##==========================================================##
    ffhomopair = []
    for i in range(0,len(possiblehomo)):
        if len(possiblehomo[i][0]) == 1 and len(possiblehomo[i][1]) == 1:
            ffhomopair.append(possiblehomo[i])
        else:
            hsgrp = possiblehomo[i][0]
            ttgrp = possiblehomo[i][1]
            ##---------------------------##
            hsseqtmp = ""
            for hsg in hsgrp:
                hsgi = listwhich(selectcol(regr,"group"),"==",hsg)[0]
                hsseqtmp = hsseqtmp + selectcol(regr,"seq")[hsgi]
            ttseqtmp = ""
            for ttg in ttgrp:
                ttgi = listwhich(selectcol(regt,"group"),"==",ttg)[0]
                ttseqtmp = ttseqtmp + selectcol(regt,"seq")[ttgi]
            ##------------------
            wholealign = pairwisealign(hsseqtmp,ttseqtmp,method="global")
            scorewtest = scorealignment(wholealign["seq1"],wholealign["seq2"])                    
            ##---------------------------##
            scorelmax = []
            grptmp = []
            for hsg in hsgrp:
                hsgi = listwhich(selectcol(regr,"group"),"==",hsg)[0]
                hsseqtmp = selectcol(regr,"seq")[hsgi]
                for ttg in ttgrp:
                    ttgi = listwhich(selectcol(regt,"group"),"==",ttg)[0]
                    ttseqtmp = selectcol(regt,"seq")[ttgi]
                    grptmp.append([[hsg],[ttg]])
                    localalign = pairwisealign(hsseqtmp,ttseqtmp,method="global")
                    scoreltest = scorealignment(localalign["seq1"],localalign["seq2"])                                     
                    scorelmax.append(scoreltest)
                
            selecttari = listwhich(scorelmax,">=",scorewtest) ## judge whether the exon map is caused by repeat exon
            if len(selecttari) > 0:
                ffhomopair = ffhomopair + selectele(grptmp,[listwhich(scorelmax,"==",max(scorelmax))[0]])
            else:
                ffhomopair.append(possiblehomo[i])
    homoexongrp = {"dat":copy.deepcopy(ffhomopair),"coln":["hsag","ptrg"]}
    ##==========================================================##
    ##  retrive lost exons in multiple-one or one-multiple      ##
    ##  relationships because of microexon or lower identity    ##
    ##==========================================================##
    if len(homoexongrp["dat"]) > 0:
        for i in range(0,len(homoexongrp["dat"])):
            if len(homoexongrp["dat"][i][0]) > 1:
                homoexongrp["dat"][i][0] = list(range(min(homoexongrp["dat"][i][0]),max(homoexongrp["dat"][i][0])+1))
            if len(homoexongrp["dat"][i][1]) > 1:
                homoexongrp["dat"][i][1] = list(range(min(homoexongrp["dat"][i][1]),max(homoexongrp["dat"][i][1])+1))
    ##==========================================================##
    ##  return final results                                    ##
    ##==========================================================##   
    return homoexongrp

##======================================================##
##     orthology and colinearity test of exon           ##
##======================================================##
def sortgroup(grouppair, anchor):
    ## resort group pairs, this is important to export ordered homo-exon pairs
    ## re-organize group         
    newpaird = {"dat":grouppair,"coln":["hsa","ptr"]}
    if len(anchor) == 0:
        newpair = grouppair
    else:
        newpair = []
        for ii in range(0,len(anchor)):
            ## test whether list-list (1:N, N-1, N-N) exist
            if ii == 0:
                if isinstance(anchor[ii][0],list):
                    hsamax = min(anchor[ii][0])
                    ptrmax = min(anchor[ii][1])
                    hsamin = min(anchor[ii][0])
                    ptrmin = min(anchor[ii][1])
                else:
                    hsamax = anchor[ii][0]
                    ptrmax = anchor[ii][1]
                    hsamin = anchor[ii][0]
                    ptrmin = anchor[ii][1]                        
            else:
                if isinstance(anchor[ii][0],list):
                    hsamax = min(anchor[ii][0])
                    ptrmax = min(anchor[ii][1])
                else:
                    hsamax = anchor[ii][0]
                    ptrmax = anchor[ii][1]
                ##----------------------------------------
                if isinstance(anchor[ii-1][0],list):
                    hsamin = max(anchor[ii-1][0])
                    ptrmin = max(anchor[ii-1][1])
                else:
                    hsamin = anchor[ii-1][0]
                    ptrmin = anchor[ii-1][1]                
            ##----------------------------------------
            if ii == 0:
                hsii = union(listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"<",hsamax),listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"==",0))
                ttii = union(listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"<",ptrmax),listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"==",0))
                selecti = intersect(hsii, ttii)
                sttmp = selectrow(newpaird,selecti)
                newpair = newpair + listsort(sttmp,by=["hsa","ptr"],increase=True)["dat"]
                newpair = newpair + [anchor[ii]]
            else:
                hsii1 = intersect(listwhich(listasint(selectcol(newpaird,"hsa"),deal="min"),">",hsamin),listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"<",hsamax))
                ttii1 = intersect(listwhich(listasint(selectcol(newpaird,"ptr"),deal="min"),">",ptrmin),listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"<",ptrmax))
                hsii = union(hsii1,listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"==",0))
                ttii = union(ttii1,listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"==",0))
                selecti = intersect(hsii, ttii)
                sttmp = selectrow(newpaird,selecti)
                newpair = newpair + listsort(sttmp,by=["hsa","ptr"],increase=True)["dat"]
                newpair = newpair + [anchor[ii]]
            ##----------------------------------------
            if ii == len(anchor) - 1:
                hsii = union(listwhich(listasint(selectcol(newpaird,"hsa"),deal="min"),">",hsamax),listwhich(listasint(selectcol(newpaird,"hsa"),deal="max"),"==",0))
                ttii = union(listwhich(listasint(selectcol(newpaird,"ptr"),deal="min"),">",ptrmax),listwhich(listasint(selectcol(newpaird,"ptr"),deal="max"),"==",0))
                selecti = intersect(hsii, ttii)
                sttmp = selectrow(newpaird,selecti)
                newpair = newpair + listsort(sttmp,by=["hsa","ptr"],increase=True)["dat"]
    ##----------------------------------------                     
    return newpair

def findnonhomoregion(homogrp):
    ## find local groups that are required tooptimise
    if len(homogrp["dat"]) > 1:
        total = []
        localhs = []
        localpt = []
        if homogrp["dat"][0][0] == 0:
            localpt.append(homogrp["dat"][0][1])
        if homogrp["dat"][0][1] == 0:
            localhs.append(homogrp["dat"][0][0])
        for i in range(1,len(homogrp["dat"])):
            if homogrp["dat"][i][0] == 0:
                if homogrp["dat"][i-1][0] == 0 or homogrp["dat"][i-1][1] == 0:
                    localpt.append(homogrp["dat"][i][1])

            elif homogrp["dat"][i][1] == 0:
                if homogrp["dat"][i-1][1] == 0 or homogrp["dat"][i-1][0] == 0:
                    localhs.append(homogrp["dat"][i][0])
            
            elif homogrp["dat"][i][0] != 0 and homogrp["dat"][i][1] != 0:
                if len(localhs) > 0 and len(localpt) > 0:
                    total.append([localhs, localpt])
                localhs = []
                localpt = []

        if len(localhs) > 0 and len(localpt) > 0:
            total.append([localhs, localpt])               

    else:
        total = []

    return unique(total)

def findlocalmatch(regr,regt,cdshs,cdstt,identhres):
    hsgrp = selectcol(regr,"group")
    ttgrp = selectcol(regt,"group")
    ##------------------------------------------------------##
    idenml = []  ##
    #etyi = listwhich(regr["coln"],"==","etype")[0]
    epsti = listwhich(regr["coln"],"==","start")[0]
    epeni = listwhich(regr["coln"],"==","end")[0]
    cdsseqi = listwhich(cdshs["coln"],"==","seq")[0]
    for hsg in hsgrp:
        hsgi = listwhich(selectcol(regr,"group"),"==",hsg)[0]
        for ttg in ttgrp:
            ttgi = listwhich(selectcol(regt,"group"),"==",ttg)[0]
            hseidx = intersect(listwhich(selectcol(cdshs,"start"),">=",regr["dat"][hsgi][epsti]),listwhich(selectcol(cdshs,"end"),"<=",regr["dat"][hsgi][epeni]))
            tteidx = intersect(listwhich(selectcol(cdstt,"start"),">=",regt["dat"][ttgi][epsti]),listwhich(selectcol(cdstt,"end"),"<=",regt["dat"][ttgi][epeni]))
            idenme = []
            idenmmax = []            
            for jjj in hseidx:
                hsseqtmp = cdshs["dat"][jjj][cdsseqi]
                for iii in tteidx:
                    ttseqtmp = cdstt["dat"][iii][cdsseqi]
                    seqidenltmp = calIdentity(hsseqtmp,ttseqtmp,method="local")
                    idenme.append(seqidenltmp)
                    idenmmax.append(round(seqidenltmp[0]*seqidenltmp[1]*seqidenltmp[2],2))
            
            selecttmpi = listwhich(idenmmax,"==",max(idenmmax))[0]
            idenml.append([hsg,ttg] + idenme[selecttmpi])
    ##------------------------------------------------------##
    groupviol = []
    for i in range(0,len(idenml)):
        if idenml[i][2] >= identhres and (idenml[i][3] >= identhres or idenml[i][4] >= identhres):
            groupviol.append([idenml[i][0],idenml[i][1]])
    ##==============================================##
    ## filter
    if len(groupviol) > 1:
        dropi = []
        groupviold = {"dat":groupviol,"coln":["hsag","ptrg"]}
        for i in range(1,len(groupviol)):
            if not(groupviol[i][1] >= max(listasint(selectcol(groupviold,"ptrg")[0:i],"max"))):
                dropi = union(dropi,i)
        if len(dropi) > 0:
            groupviol = listdrop(groupviol,dropi)
    groupviol = unique(groupviol)
    ##==============================================##
    ## find idenm
    idenmld = {"dat":idenml,"coln":["hsag","ptrg","iden","hsacov","ptrcov"]}
    idenmlnew = []
    for i in range(0,len(groupviol)):
        selecti = intersect(listwhich(selectcol(idenmld,"hsag"),"==",groupviol[i][0]),listwhich(selectcol(idenmld,"ptrg"),"==",groupviol[i][1]))[0]
        idenmlnew.append(idenml[selecti])
    idenml = {"dat":copy.deepcopy(idenmlnew),"coln":["hsag","ptrg","iden","hsacov","ptrcov"]}
    ##======= detect possible 1-1 groups, which is the reference list if 1-N, N-1 even N-N failed
    if len(groupviol) > 0:
        seqi = listwhich(regr["coln"],"==","seq")[0]
        groupvio = matchhomoexon(groupviol,regr,regt,identhres)
        unigroupviol = []
        if len(groupvio) > 0:
            for i in range(0,len(groupvio)):
                if not(isinstance(groupvio[i][0],list)):
                    unigroupviol.append([[groupviol[i][0]],[groupviol[i][1]]])
                else:
                    idenmtestmax = []
                    testgrp = []
                    for hsgtmp in groupvio[i][0]:
                        for ttgtmp in groupvio[i][1]:
                            testgrp.append([hsgtmp,ttgtmp])
                            selectideni = intersect(listwhich(selectcol(idenml,"hsag"),"==",hsgtmp),listwhich(selectcol(idenml,"ptrg"),"==",ttgtmp))[0]
                            idenmtestmax.append(round(idenml["dat"][selectideni][2]*idenml["dat"][selectideni][3]*idenml["dat"][selectideni][4],2))    
                    targrp = testgrp[listwhich(idenmtestmax,"==",max(idenmtestmax))[0]]
                    unigroupviol.append([[targrp[0]],[targrp[1]]])
        unigroupviold = {"dat":unigroupviol,"coln":["hsag","ptrg"]}
        ##------------- reorganize ------------
        eletest = []
        for i in range(0,len(groupvio)):
            eletest.append(groupvio[i] in unigroupviol)
        if len(listwhich(eletest,"==",True)) == len(groupvio):
            groupvio = copy.deepcopy(unigroupviol)
        else:
            existst = []
            existen = []
            for j in range(0,len(groupvio)):
                existst = existst + groupvio[j][0]
                existen = existen + groupvio[j][1]
            selecti = []
            for j in range(0,len(unigroupviol)):
                if not(unigroupviol[j][0][0] in existst) and not(unigroupviol[j][1][0] in existen):
                    selecti.append(j)
            if len(selecti) > 0:
                groupvio = groupvio + selectele(unigroupviol,selecti)
                groupvio.sort()
        ##------------- retest ------------
        newhomopair = []
        for i in range(0,len(groupvio)):
            if len(groupvio[i][0]) == 1 and len(groupvio[i][1]) == 1:
                newhomopair.append(groupvio[i])
            else:
                ##===== calculate group to group identity
                identmpl = []
                for j in range(0,len(groupvio[i][0])):
                    for k in range(0,len(groupvio[i][1])):
                        identmp = calIdentity(regr["dat"][listwhich(selectcol(regr,"group"),"==",groupvio[i][0][j])[0]][seqi],regt["dat"][listwhich(selectcol(regt,"group"),"==",groupvio[i][1][k])[0]][seqi],method="local")
                        identmpl.append(identmp[0]*identmp[1]*identmp[2])
                ##----- calculate merged group identity
                hsseq = ""
                for j in range(0,len(groupvio[i][0])):
                    hsseq = hsseq + regr["dat"][listwhich(selectcol(regr,"group"),"==",groupvio[i][0][j])[0]][seqi]
                ttseq = ""
                for j in range(0,len(groupvio[i][1])):
                    ttseq = ttseq + regt["dat"][listwhich(selectcol(regt,"group"),"==",groupvio[i][1][j])[0]][seqi]
                identmp = calIdentity(hsseq,ttseq,method="local")
                ##------ judge
                dropidx = listwhich(identmpl,">",identmp[0]*identmp[1]*identmp[2]) ## judge whether the exon map is caused by repeat exon
                if len(dropidx) == 0:
                    newhomopair.append(groupvio[i])
                else: ## only keep the first one
                    ## find which unigroupviol contains target group, choose the group pair with maximum identity
                    selectunii = []
                    for hsgtmp in groupvio[i][0]:
                        for ttgtmp in groupvio[i][1]:
                            selectideni = intersect(listwhich(selectcol(unigroupviold,"hsag"),"==",[hsgtmp]),listwhich(selectcol(unigroupviold,"ptrg"),"==",[ttgtmp]))
                            selectunii = union(selectunii,selectideni)
                    if len(selectunii) > 0:
                        newhomopair.append(unigroupviold["dat"][selectunii[0]])
    else:
        newhomopair = []
    return newhomopair

def colinearity_test(homopair):
    if len(homopair) == 0:
        passhomo2 =  {"dat":[],"coln":["hsag","ptrg"]}
    else:
        homogdict = {"dat":unique(homopair),"coln":["hsag","ptrg"]}
        homogdict = listsort(homogdict,by=["hsag","ptrg"],increase=True)
        ##--------------- the start group: homopair[0] --------------##
        minttg = homogdict["dat"][0][1]
        minhsg = homogdict["dat"][0][0]
        passhomo = {"dat":[[minhsg,minttg]],"coln":["hsag","ptrg"]}
        ttgtmp = homogdict["dat"][0][1]
        ##-- find the hsa group that need to be optimised ---
        unighs = unique(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"hsag"),">=",minhsg)))
        for i in range(1,len(unighs)):
            ptrgtmp = selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i]))
            #minttg = max(selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i])))
            maxttgminus = max(selectele(selectcol(passhomo,"ptrg"),listwhich(selectcol(passhomo,"hsag"),"<",unighs[i])))
            ptrtari = listwhich(ptrgtmp,">",maxttgminus)
            if len(ptrtari) > 0:
                ttgtmptmp = ptrgtmp[min(ptrtari)]
                if len(ttgtmptmp) != 1 or len(ttgtmp) != 1:
                    ttgtmp = setdiff(ttgtmptmp,ttgtmp)
                else:
                    ttgtmp = copy.deepcopy(ttgtmptmp)
                
                if len(ttgtmp) > 0:
                    passhomo["dat"].append([unighs[i],ttgtmp])
        ##--------------- the start group: [minhsg,minttg] --------------##
        minttg = min(selectcol(homogdict,"ptrg"))
        minhsg = min(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"ptrg"),"==",minttg)))
        passhomo1 = {"dat":[[minhsg,minttg]],"coln":["hsag","ptrg"]}
        ttgtmp = copy.deepcopy(minttg)
        ##-- find the hsa group that need to be optimised ---
        unighs = unique(selectele(selectcol(homogdict,"hsag"),listwhich(selectcol(homogdict,"hsag"),">=",minhsg)))
        for i in range(1,len(unighs)):
            ptrgtmp = selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i]))
            #minttg = max(selectele(selectcol(homogdict,"ptrg"),listwhich(selectcol(homogdict,"hsag"),"==",unighs[i])))
            maxttgminus = max(selectele(selectcol(passhomo1,"ptrg"),listwhich(selectcol(passhomo1,"hsag"),"<",unighs[i])))
            ptrtari = listwhich(ptrgtmp,">",maxttgminus)
            if len(ptrtari) > 0:
                ttgtmptmp = ptrgtmp[min(ptrtari)]
                if len(ttgtmptmp) != 1 or len(ttgtmp) != 1:
                    ttgtmp = setdiff(ttgtmptmp,ttgtmp)
                else:
                    ttgtmp = copy.deepcopy(ttgtmptmp)
                
                if len(ttgtmp) > 0:
                    passhomo1["dat"].append([unighs[i],ttgtmp])
        ## choose the best matches
        if len(passhomo1["dat"]) > len(passhomo["dat"]):
            passhomo2 = passhomo1
        else:
            passhomo2 = passhomo
      
    return passhomo2

def mergehomogroup(homog1,homog2):
    if len(homog1["dat"]) == 0 and len(homog2["dat"]) == 0:
        homogrptmp = {"dat":[],"coln":["hsag","ptrg"]}
    elif len(homog1["dat"]) == 0:
        homogrptmp = copy.deepcopy(homog2)
    elif len(homog2["dat"]) == 0:
        homogrptmp = copy.deepcopy(homog1)
    else:
        selecti = [] ## for homog2
        dropi = []   ## for homog1
        for i in range(0,len(homog2["dat"])):
            rig = homog2["dat"][i][0]
            leg = homog2["dat"][i][1]
            idxtmp = []
            storenon = []
            storeexi = []
            for j in range(0,len(rig)):
                for m in range(0,len(leg)):
                    idxtmptmp = intersect(listwhich(selectcol(homog1,"hsag"),"==",[rig[j]]),listwhich(selectcol(homog1,"ptrg"),"==",[leg[m]]))
                    if len(idxtmptmp) > 0 :
                        idxtmp.append(idxtmptmp[0])
                    ## test if possible to replace
                    if not(([rig[j]] in selectcol(homog1,"hsag")) or ([leg[m]] in selectcol(homog1,"ptrg")) ):
                        storenon.append(True)
                    else:
                        storenon.append(False)
                    ## test if possible to add non-record
                    if not(([rig[j]] in selectcol(homog1,"hsag")) and ([leg[m]] in selectcol(homog1,"ptrg")) ):
                        storeexi.append(True)
                    else:
                        storeexi.append(False)
                    
            if len(idxtmp) > 0 and len(listwhich(storeexi,"==",True)) == len(storeexi)-1:
                dropi.append(idxtmp[0])
                selecti.append(i)
            elif len(idxtmp) == 0 and len(listwhich(storenon,"==",True)) == len(storenon):
                selecti.append(i)
        homog = []
        homog1["dat"] = listdrop(homog1["dat"],dropi)
        homog2["dat"] = selectele(homog2["dat"],selecti)
        homog = homog + homog1["dat"] + homog2["dat"]
        homogrptmp = colinearity_test(homog)

    return homogrptmp

def homo_colinearity_test(species1, species2,regrhs,regrtt,cdshs,cdstt,homocdsdat,identhres=0.8, coverthres=0.8, minexon=10, mapscore=2, misscore=-2, gapscore=-1):
    ##=====================================##
    ##============= find index ============##
    ##=====================================##
    hsai = listwhich(homocdsdat["coln"],"==",species1)[0]
    ptri = listwhich(homocdsdat["coln"],"==",species2)[0]
    idsti = listwhich(cdshs["coln"],"==","start")[0]
    ideni = listwhich(cdshs["coln"],"==","end")[0]
    grpi = listwhich(regrhs["coln"],"==","group")[0]
    ##==================================================================##
    ##===== construct orthologous exon matrix based on exon BLASTN =====##
    ##==================================================================##
    if len(homocdsdat["dat"]) > 0:
        homogroup = []
        for hci in range(0,len(homocdsdat["dat"])):
            hsae = homocdsdat["dat"][hci][hsai]
            ptre = homocdsdat["dat"][hci][ptri]
            ##----------------------------------------
            hsast = cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hsae)[0]][idsti]
            hsaen = cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hsae)[0]][ideni]
            idxhs = intersect(listwhich(selectcol(regrhs,"start"),"<=",hsast),listwhich(selectcol(regrhs,"end"),">=",hsaen))[0]
            hsgrp = regrhs["dat"][idxhs][grpi]
            ##----------------------------------------
            ptrst = cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",ptre)[0]][idsti]
            ptren = cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",ptre)[0]][ideni]
            idxtt = intersect(listwhich(selectcol(regrtt,"start"),"<=",ptrst),listwhich(selectcol(regrtt,"end"),">=",ptren))[0]
            ttgrp = regrtt["dat"][idxtt][grpi]
            ##----------------------------------------
            homogroup.append([[hsgrp],[ttgrp]])
        homogrp1 = {"dat":unique(homogroup),"coln":["hsag","ptrg"]}
        homogrp1 = listsort(homogrp1,by=["hsag","ptrg"],increase=True)
        ##== colinearity_test of conserved exons that are called by BLASTN =##
        homogrp1 = colinearity_test(homogrp1["dat"])
        ##=====================================================================##
    else:
        homogrp1 = {"dat":[],"coln":["hsag","ptrg"]}
    ##==================================================================##
    ## construct orthologous exon matrix based on exon based alignment  ##
    ##==================================================================##
    homogrp2 = alignmentexon(regr=regrhs, regt=regrtt, microexon=minexon, identhres=identhres, coverthres=coverthres, score=mapscore, penalty=misscore, gap=gapscore)
    homogrp2tmp = copy.deepcopy(homogrp2)
    homogrp2tmptmp = colinearity_test(homogrp2["dat"])
    while homogrp2tmptmp != homogrp2tmp:
        homogrp2tmp = copy.deepcopy(homogrp2tmptmp)
        homogrp2tmptmp = colinearity_test(homogrp2tmptmp["dat"])
    homogrp2 = copy.deepcopy(homogrp2tmptmp)
    ##==================================================================##
    ## correct the dynamic programming method using BLASTN based results##
    ##==================================================================##
    homogrp = mergehomogroup(homogrp1,homogrp2)    
    ##----------- organize orthologous exon group -------------##
    for i in range(0,len(homogrp["dat"])):
        if len(homogrp["dat"][i][0]) == 1 and len(homogrp["dat"][i][1]) == 1:
            homogrp["dat"][i][0] = homogrp["dat"][i][0][0]
            homogrp["dat"][i][1] = homogrp["dat"][i][1][0]
    anchor = copy.deepcopy(homogrp["dat"])
    ##==========================================================##
    ##============ add non record exonsinto homo list ==========##
    ##==========================================================##
    hsgroup = selectcol(regrhs,"group")
    ttgroup = selectcol(regrtt,"group")
    for hsg in hsgroup:
        if not(hsg in listasint(selectcol(homogrp,"hsag"),deal="all")):
            if len(homogrp["dat"])>0 and hsg < max(listasint(selectcol(homogrp,"hsag"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"hsag"),deal="min"),">",hsg)),[hsg,0])
            else:
                homogrp["dat"].append([hsg,0])
    for ttg in ttgroup:
        if not(ttg in listasint(selectcol(homogrp,"ptrg"),deal="all")):
            if len(homogrp["dat"])>0 and ttg < max(listasint(selectcol(homogrp,"ptrg"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"ptrg"),deal="min"),">",ttg)),[0,ttg])
            else:
                homogrp["dat"].append([0,ttg])
    ##==========================================================##
    ## sort group in standard order so that findnonhomoregion   ##
    ## can estimate correct position range                      ##
    ##==========================================================##
    homogrp["dat"] = sortgroup(grouppair=homogrp["dat"], anchor=anchor)
    ##=============================================================##
    ## find non-orthologous exon that are needed violent alignment ##
    ##=============================================================##
    total = findnonhomoregion(homogrp)
    ##-------------------------------------##
    localmatch = []
    if len(total) > 0:
        for ii in range(0,len(total)):
            hsgi = []
            for j in total[ii][0]:
                hsgi.append(listwhich(selectcol(regrhs,"group"),"==",j)[0])
            ptgi = []
            for j in total[ii][1]:
                ptgi.append(listwhich(selectcol(regrtt,"group"),"==",j)[0])
            violentalign = findlocalmatch(regr=selectrow(regrhs,hsgi), regt=selectrow(regrtt,ptgi),cdshs=cdshs,cdstt=cdstt,identhres=identhres)
            localmatch = localmatch + violentalign
    ##---------- drop and loop to insert ---------## 
    if len(localmatch) > 0:
        for i in range(0,len(localmatch)):
            if len(localmatch[i][0]) == 1 and len(localmatch[i][1]) == 1:
                localmatch[i][0] = localmatch[i][0][0]
                localmatch[i][1] = localmatch[i][1][0]
        ## drop all the non-orthologous exon that are needed violent alignment
        dropi = []
        for i in range(0,len(total)):
            if isinstance(total[i][0],list):
                for j in total[i][0]:
                    dropi.append(listwhich(selectcol(homogrp,"hsag"),"==",j)[0])
                for j in total[i][1]:
                    dropi.append(listwhich(selectcol(homogrp,"ptrg"),"==",j)[0])
            else:      
                dropi.append(listwhich(selectcol(homogrp,"hsag"),"==",total[i][0])[0])
                dropi.append(listwhich(selectcol(homogrp,"ptrg"),"==",total[i][1])[0])
        if len(dropi) > 0:
            homogrp["dat"] = listdrop(homogrp["dat"],dropi)
        ### homogrp might be empty after dropping, for example retro-gene
        if len(homogrp["dat"]) == 0:
            homogrp["dat"] = copy.deepcopy(localmatch)
            ### insert localmatch to orthologous group
        else:
            for i in range(0,len(localmatch)):
                if isinstance(localmatch[i][0],list):
                    if listasint(localmatch[i][0],"max")[0] > max(listasint(selectcol(homogrp,"hsag"),"max")):
                        homogrp["dat"].append(localmatch[i])
                    else:
                        inserti = min(listwhich(listasint(selectcol(homogrp,"hsag"),"min"),">",listasint(localmatch[i][0],"max")[0]))
                        homogrp["dat"].insert(inserti,localmatch[i])
                else:
                    if localmatch[i][0] > max(listasint(selectcol(homogrp,"hsag"),"max")):
                        homogrp["dat"].append(localmatch[i])
                    else:
                        inserti = min(listwhich(listasint(selectcol(homogrp,"hsag"),"min"),">",localmatch[i][0]))
                        homogrp["dat"].insert(inserti,localmatch[i])
        ##----------- organize orthologous exon group -------------##
        for i in range(0,len(homogrp["dat"])):
            if isinstance(homogrp["dat"][i][0],list) and len(homogrp["dat"][i][0]) == 1 and len(homogrp["dat"][i][1]) == 1:
                homogrp["dat"][i][0] = homogrp["dat"][i][0][0]
                homogrp["dat"][i][1] = homogrp["dat"][i][1][0]
    ##=====================================================
    hsgroup = selectcol(regrhs,"group")
    ttgroup = selectcol(regrtt,"group")
    for hsg in hsgroup:
        if not(hsg in listasint(selectcol(homogrp,"hsag"),deal="all")):
            if hsg < max(listasint(selectcol(homogrp,"hsag"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"hsag"),deal="min"),">",hsg)),[hsg,0])
            else:
                homogrp["dat"].append([hsg,0])
    for ttg in ttgroup:
        if not(ttg in listasint(selectcol(homogrp,"ptrg"),deal="all")):
            if ttg < max(listasint(selectcol(homogrp,"ptrg"),deal="all")):
                homogrp["dat"].insert(min(listwhich(listasint(selectcol(homogrp,"ptrg"),deal="min"),">",ttg)),[0,ttg])
            else:
                homogrp["dat"].append([0,ttg])
    ##=====================================================##
    ##====== test if lost group during optimise       =====##
    ##=====================================================##
    refhsa = selectcol(regrhs,"group")
    refptr = selectcol(regrtt,"group")
    tsthsa = []
    tstptr = []
    for i in range(0,len(homogrp["dat"])):
        if isinstance(homogrp["dat"][i][0],list):
            tsthsa = tsthsa + homogrp["dat"][i][0]
            tstptr = tstptr + homogrp["dat"][i][1]
        else:
            if homogrp["dat"][i][0] != 0:
                tsthsa.append(homogrp["dat"][i][0])
            if homogrp["dat"][i][1] != 0:
                tstptr.append(homogrp["dat"][i][1])  
    ##----------------------------------------
    diffhsa = setdiff(refhsa,tsthsa)
    diffptr = setdiff(refptr,tstptr)
    if len(diffhsa) != 0 or len(diffptr) != 0:
        raise SyntaxError("homo_colinearity_test: group loss after homo-corlinearity optimise")
    ##----------------------------------------
    if len(homogrp["dat"])>1:
        for i in range(1,len(homogrp["dat"])):
            if homogrp["dat"][i][0] != 0:
                test1 = min(listasint(homogrp["dat"][i][0])) > max(listasint(selectcol(homogrp,"hsag")[0:i],"max"))
            else:
                test1 = True
            if homogrp["dat"][i][1] != 0:
                test2 = min(listasint(homogrp["dat"][i][1])) > max(listasint(selectcol(homogrp,"ptrg")[0:i],"max"))
            else:
                test2 = True           
            if not(test1) or not(test2):
                raise SyntaxError("homo_colinearity_test: corlinearity test fail"+str()+":"+str())

    return homogrp

def orgexongroup(genhs,regrhs,strandhs,chromhs,gentt,regrtt,strandtt,chromtt,refhomo):
    ## organize ruler group information
    store = []
    newgroup = 0
    sti = listwhich(regrhs["coln"],"==","start")[0]
    eni = listwhich(regrhs["coln"],"==","end")[0]
    seqi = listwhich(regrhs["coln"],"==","seq")[0]
    for i in range(0,len(refhomo["dat"])):
        newgroup = newgroup + 1
        if not(isinstance(refhomo["dat"][i][0],list)):
            hsg = refhomo["dat"][i][0]
            ttg = refhomo["dat"][i][1]
            if hsg == 0 and ttg != 0:
                hspos = "None"
                ttstart = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][sti]
                ttend = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][eni]
                ttpos = "chr"+str(chromtt)+":"+str(ttstart)+":"+str(ttend)+":"+str(strandtt)
                iden = [0]
                grptype = "0-1"
            elif hsg != 0 and ttg == 0:
                hsstart = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][sti]
                hsend = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][eni]
                hspos = "chr"+str(chromhs)+":"+str(hsstart)+":"+str(hsend)+":"+str(strandhs)
                ttpos = "None"
                iden = [0]
                grptype = "1-0" 
            else:
                hsstart = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][sti]
                hsend = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][eni]
                hspos = "chr"+str(chromhs)+":"+str(hsstart)+":"+str(hsend)+":"+str(strandhs)
                ##----------------------------------------
                ttstart = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][sti]
                ttend = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][eni]
                ttpos = "chr"+str(chromtt)+":"+str(ttstart)+":"+str(ttend)+":"+str(strandtt)
                ##----------------------------------------
                hsseqtmp = regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsg)[0]][seqi]
                ttseqtmp = regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttg)[0]][seqi]
                ##----------------------------------------
                iden = calIdentity(hsseqtmp, ttseqtmp, method="global")
                grptype = "1-1"
        else:
            hsgtmp = []
            hssttmp = []
            hsentmp = []
            hsseqtmp = ""
            for hsii in refhomo["dat"][i][0]:
                hsseqtmp = hsseqtmp + regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsii)[0]][seqi]
                hssttmp.append(regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsii)[0]][sti])
                hsentmp.append(regrhs["dat"][listwhich(selectcol(regrhs,"group"),"==",hsii)[0]][eni])
                hsgtmp.append(hsii)
            ##----------------------------------------
            ttgtmp = []
            ttsttmp = []
            ttentmp = []
            ttseqtmp = ""
            for ttii in refhomo["dat"][i][1]:
                ttseqtmp = ttseqtmp + regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttii)[0]][seqi]
                ttsttmp.append(regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttii)[0]][sti])
                ttentmp.append(regrtt["dat"][listwhich(selectcol(regrtt,"group"),"==",ttii)[0]][eni])                
                ttgtmp.append(ttii)
            ##----------------------------------------
            hspos = "chr"+str(chromhs)+":"+str(min(hssttmp))+":"+str(max(hsentmp))+":"+str(strandhs)
            ttpos = "chr"+str(chromtt)+":"+str(min(ttsttmp))+":"+str(max(ttentmp))+":"+str(strandtt)
            #### it is important to use global alignment         
            iden = calIdentity(hsseqtmp, ttseqtmp, method="global")
            if len(hsgtmp) == 1 and len(ttgtmp) > 1:
                grptype = "1-N"
            elif len(hsgtmp) > 1 and len(ttgtmp) == 1:
                grptype = "N-1"
            else:
                grptype = "N-N"
        
        store.append([newgroup,genhs,hspos,gentt,ttpos,iden[0],grptype])

    return {"dat":store, "coln":["group","hsagen","hsapos","ptrgen","ptrpos","identity","grouptype"]}

##======================================================##
##          detection of orthologous isoforms           ##
##======================================================##
def regComSum(genestd,regidxd,dattype="cds/utr"):
    genename = genestd["coln"]
    genest = genestd["dat"]
    regidxname = regidxd["coln"]
    regidx = regidxd["dat"]

    regnum = len(regidx)
    samnum = len(genest)
    coli = listwhich(genename,"==","ensemblt")[0]
    isonam = ["id","start"]
    for n in list(range(0,len(genest))):
        isonam.append(genest[n][coli])
    
    commatrix = []
    starti = listwhich(regidxname,"==","start")[0]
    endi = listwhich(regidxname,"==","end")[0]
    idi = listwhich(regidxname,"==","id")[0]
    typei = listwhich(genename,"==",dattype)[0]
    for i in list(range(0,regnum)):
        regposabb = ":".join([str(regidx[i][starti]),str(regidx[i][endi])]) # compare the string rather pos num
        regstatus = []
        regstatus.append(regidx[i][idi])                  # the first column stores cds region id
        regstatus.append(int(regidx[i][starti]))  
        for j in list(range(0,samnum)):
            regpostmp = genest[j][typei].split("|")
            regstatmp = 0
            locnum = len(regpostmp)
            for m in list(range(0,locnum)):
                if regpostmp[m] == regposabb:
                    regstatmp = 1
                    break
            regstatus.append(regstatmp)
        commatrix.append(regstatus)
    ## drop void cds region
    commatd = {"dat":commatrix,"coln":isonam}
    countcom = [0 for _ in range(0,len(commatrix))]
    for i in range(0,len(commatrix)):
        countcom[i] = sum(commatrix[i][1:])
    seli = listwhich(countcom,"!=",0)
    regsum = selectlsls(commatd,seli,"all")
    return regsum

def consertiveIso(exoncon,exonttcon,cdshs,cdstt,homoexondat,identhres):
    ##========================================================##
    ##-------------- estimate exon group list  ---------------##
    ##========================================================##
    sti = listwhich(cdshs["coln"],"==","start")[0]
    eni = listwhich(cdshs["coln"],"==","end")[0]
    grpi = listwhich(homoexondat["coln"],"==","group")[0]
    ##---------- hsa -------------
    hsaexongroup = []
    for hse in exoncon:
        st = int(cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hse)[0]][sti])
        en = int(cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",hse)[0]][eni])
        grptmpi = intersect(listwhich(selectcol(homoexondat,"hsast"),"<=",st),listwhich(selectcol(homoexondat,"hsaen"),">=",en))[0]
        hsaexongroup.append(homoexondat["dat"][grptmpi][grpi])
    hsaegrp = unique(hsaexongroup)
    ##---------- hsa -------------
    ptrexongroup = []
    for pte in exonttcon:
        st = int(cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",pte)[0]][sti])
        en = int(cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",pte)[0]][eni])
        grptmpi = intersect(listwhich(selectcol(homoexondat,"ptrst"),"<=",st),listwhich(selectcol(homoexondat,"ptren"),">=",en))[0]
        ptrexongroup.append(homoexondat["dat"][grptmpi][grpi])
    ptregrp = unique(ptrexongroup)
    ##--------- compare ----------
    if hsaegrp != ptregrp:
        passtest = 0
    else:
        passtest = 1
    ##========================================================##
    ##----------- need further sequence test ??? -------------##
    ##========================================================##
    infostr = ""
    if passtest == 1:
        seqi = listwhich(cdshs["coln"],"==","seq")[0]
        idenlist = []
        infotmp = []
        for i in range(0,len(hsaegrp)):
            hsae = selectele(exoncon,listwhich(hsaexongroup,"==",hsaegrp[i]))
            ptre = selectele(exonttcon,listwhich(ptrexongroup,"==",hsaegrp[i]))

            hsaseq = ""
            ptrseq = ""
            for he in hsae:
                hsaseq = hsaseq + cdshs["dat"][listwhich(selectcol(cdshs,"id"),"==",he)[0]][seqi]
            for pe in ptre:
                ptrseq = ptrseq + cdstt["dat"][listwhich(selectcol(cdstt,"id"),"==",pe)[0]][seqi]

            identmp = calIdentity(hsaseq, ptrseq, method="global")[0]
            idenlist.append(identmp)
            infotmp.append(hsae[0]+"-"+ptre[0]+"|"+str(identmp))    
        
        if len(listwhich(idenlist,"<",identhres)) > 0:
            passtest = 0

        infostr = ";".join(infotmp)
    
    suminfo = {"passtest":passtest, "info":infostr}
    
    return suminfo

def orthoiso(species1, species2, genhs, geninfohs, strandhs, cdshs, gentt, geninfott, strandtt, cdstt, conexondat, identhres):
    ##--------------- reorganize conexondat  ----------------##
    homoexon = []
    grpi = listwhich(conexondat["coln"],"==","Group")[0]
    hsapi = listwhich(conexondat["coln"],"==",species1+"Pos")[0]
    ptrpi = listwhich(conexondat["coln"],"==",species2+"Pos")[0]
    ideni = listwhich(conexondat["coln"],"==","Iden")[0]
    for i in range(0,len(conexondat["dat"])):
        ##---------- hsa -------------
        if conexondat["dat"][i][hsapi] == "None":
            hsast = 0
            hsaen = 0
        else:
            hsagrouppos = conexondat["dat"][i][hsapi].split(":")
            hsast = int(hsagrouppos[1])
            hsaen = int(hsagrouppos[2])
        ##---------- ptr -------------
        if conexondat["dat"][i][ptrpi] == "None":
            ptrst = 0
            ptren = 0
        else:
            ptrgrouppos = conexondat["dat"][i][ptrpi].split(":")
            ptrst = int(ptrgrouppos[1])
            ptren = int(ptrgrouppos[2])
        ##---------- organize ------------
        homoexon.append([conexondat["dat"][i][grpi],hsast,hsaen,ptrst,ptren,conexondat["dat"][i][ideni]])
        
    homoexondat = {"dat":homoexon,"coln":["group","hsast","hsaen","ptrst","ptren","iden"]}
    ##-------------     estimate exon group     --------------##
    sti = listwhich(cdshs["coln"],"==","start")[0]
    eni = listwhich(cdshs["coln"],"==","end")[0]   
    ##---------- hsa -------------
    for i in range(0,len(cdshs["dat"])):
        grptmpi = intersect(listwhich(selectcol(homoexondat,"hsast"),"<=",int(cdshs["dat"][i][sti])),listwhich(selectcol(homoexondat,"hsaen"),">=",int(cdshs["dat"][i][eni])))[0]
        cdshs["dat"][i].insert(0,homoexondat["dat"][grptmpi][grpi])
    cdshs["coln"].insert(0,"group")
    ##---------- ptr -------------
    for i in range(0,len(cdstt["dat"])):
        grptmpi = intersect(listwhich(selectcol(homoexondat,"ptrst"),"<=",int(cdstt["dat"][i][sti])),listwhich(selectcol(homoexondat,"ptren"),">=",int(cdstt["dat"][i][eni])))[0]
        cdstt["dat"][i].insert(0,homoexondat["dat"][grptmpi][grpi])
    cdstt["coln"].insert(0,"group")

    ##============ estimate cds ruler and gap =================
    cdsmtt = regComSum(genestd=geninfott,regidxd=cdstt,dattype="cds")
    if strandtt == 1:
        isott = listsort(cdsmtt,by=["start"],increase = True)
    if strandtt == -1:
        isott = listsort(cdsmtt,by=["start"],increase = False)
    isott = selectlsls(isott,list(range(0,len(isott["dat"]))),selectele(isott["coln"],[0]+list(range(2,len(isott["coln"])))))
    ##--------------------------------
    cdsmhs = regComSum(genestd=geninfohs,regidxd=cdshs,dattype="cds")
    if strandhs == 1:
        isohs = listsort(cdsmhs,by=["start"],increase = True)
    if strandhs == -1:
        isohs = listsort(cdsmhs,by=["start"],increase = False)
    isohs = selectlsls(isohs,list(range(0,len(isohs["dat"]))),selectele(isohs["coln"],[0]+list(range(2,len(isohs["coln"])))))
    ##------------ estimate homologous isoform  --------------##
    isoformhs = selectele(isohs["coln"],list(range(1,len(isohs["coln"]))))
    conseriso = []
    conhsgen = []
    conttgen = []
    isohsa = []
    exoniden = []
    for j in range(1,len(isohs["coln"])):
        exoncon = selectele(selectcol(isohs,"id"),listwhich(selectcol(isohs,isohs["coln"][j]),"==",1))
        tariso = "."
        for e in range(1,len(isott["coln"])):
            exonttcon = selectele(selectcol(isott,"id"),listwhich(selectcol(isott,isott["coln"][e]),"==",1))
            ## test consertive
            conser = consertiveIso(exoncon,exonttcon,cdshs,cdstt,homoexondat,identhres=identhres)
            if conser["passtest"] == 1:
                ##================================================##
                ### test ORF shift using R language!!!!!!!!!!!!!! ##
                ##================================================##
                tariso = isott["coln"][e]
                isohsa.append(isoformhs[j-1])
                conseriso.append(isott["coln"][e])
                conhsgen.append(genhs)
                conttgen.append(gentt)
                exoniden.append(conser["info"])

    summ = {"dat":transpose([conhsgen,conttgen,isohsa,conseriso,exoniden]),"coln":[species1,species2,species1+"iso",species2+"iso","exonmatch"]}
    return summ

##======================================================##
##                   analysis process                   ##
##======================================================##
def EGIO(data_dict):
    genhs = data_dict["gene1"]

    print("analysis of gene " + genhs + " ")

    chromhs = data_dict["cds1"]["dat"][0][2]
    strandhsstr = data_dict["cds1"]["dat"][0][5]
    if strandhsstr == "-":
        strandhs = -1
    else:
        strandhs = 1
    
    genhscdstmp1 = selectlsls(data_dict["iso1"],list(range(0,len(data_dict["iso1"]["dat"]))),coln=["EnsemblT","Orf","Exoncom"])
    genhscdsdat = listddrop(genhscdstmp1,listwhich(selectcol(genhscdstmp1,"Orf"),"==",".|."))
    
    cdshs = selectlsls(data_dict["cds1"],list(range(0,len(data_dict["cds1"]["dat"]))),coln=["ID","Start","End","Seq"])
    cdshs["coln"] = ["id","start","end","seq"]

    geninfohs = translatecdna(genhscdsdat,cdshs,strandhsstr)
    geninfohs["coln"] = ["ensemblt","orf","cds","aaseq"]
    #####==============================================#####
    chromtt = data_dict["cds2"]["dat"][0][2]
    gentt = data_dict["gene2"]
    strandttstr = data_dict["cds2"]["dat"][0][5]
    if strandttstr == "-":
        strandtt = -1
    else:
        strandtt = 1

    genttcdstmp = selectlsls(data_dict["iso2"],list(range(0,len(data_dict["iso2"]["dat"]))),coln=["EnsemblT","Orf","Exoncom"])
    genttcdsdat = listddrop(genttcdstmp,listwhich(selectcol(genttcdstmp,"Orf"),"==",".|."))
    
    cdstt = selectlsls(data_dict["cds2"],list(range(0,len(data_dict["cds2"]["dat"]))),coln=["ID","Start","End","Seq"])
    cdstt["coln"] = ["id","start","end","seq"]

    geninfott = translatecdna(genttcdsdat,cdstt,strandttstr)
    geninfott["coln"] = ["ensemblt","orf","cds","aaseq"]
    #####=================== blastexon  ===============#####
    if len(data_dict["blastn"]["dat"]) > 0:
        homocdsdat = selectlsls(data_dict["blastn"],listwhich(selectcol(data_dict["blastn"],"EnsemblG2"),"==",gentt),coln=[data_dict["species1"],data_dict["species2"]])
    else:
        homocdsdat = {"dat":[], "coln":[data_dict["species1"],data_dict["species2"]]}
    ##------------------ estimate cds ruler ------------------
    cdsrulerhs = regRuler(regidxd=cdshs,genstr=strandhs)
    cdsrulertt = regRuler(regidxd=cdstt,genstr=strandtt)
    ##------------ estimate CDS homologous group  ------------
    refhomo = homo_colinearity_test(data_dict["species1"], data_dict["species2"], cdsrulerhs,cdsrulertt,cdshs,cdstt,homocdsdat, identhres=data_dict["identity"], coverthres=data_dict["coverage"], minexon=data_dict["mexon"], mapscore=data_dict["match"], misscore=data_dict["mismatch"], gapscore=data_dict["gap"])
    ##=======================================================================##
    ##                         website back end scripts                      ##
    ##=======================================================================##
    sumexon = orgexongroup(genhs,cdsrulerhs,strandhs,chromhs,gentt,cdsrulertt,strandtt,chromtt,refhomo)
    ##=======================================================================##
    orthoexon = pd.DataFrame(sumexon["dat"])
    orthoexon.columns = ["Group",data_dict["species1"]+"EnsemblG",data_dict["species1"]+"Pos",data_dict["species2"]+"EnsemblG",data_dict["species2"]+"Pos","Iden","Type"]
    #####==============================================#####
    #####============= read blastexon data ============#####
    #####==============================================#####
    homoexon = []
    if len(orthoexon) > 0:
        for hei in range(0,len(orthoexon)):
            homoexon.append([orthoexon.iloc[hei]["Group"],orthoexon.iloc[hei][data_dict["species1"]+"Pos"],orthoexon.iloc[hei][data_dict["species2"]+"Pos"],orthoexon.iloc[hei]["Iden"],orthoexon.iloc[hei]["Type"]])
    conexondat = {"dat":homoexon, "coln":["Group",data_dict["species1"]+"Pos",data_dict["species2"]+"Pos","Iden","Type"]}
    #####==============================================#####
    #####============= detect homo isoofrm ============#####
    #####==============================================#####
    sumiso = orthoiso(data_dict["species1"], data_dict["species2"],genhs, geninfohs, strandhs, cdshs, gentt, geninfott, strandtt, cdstt, conexondat,identhres=data_dict["identity"])    
    #####==============================================#####
    return {"oe":sumexon, "ot":sumiso}

##****************************************************************************************************************************##
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'detection of orthologous exons and isoforms using EGIO')
    parser.add_argument('--orthog',
        nargs = '?',
        help = 'path of file contains orthologous gene pairs')
    parser.add_argument('--blastn',
        nargs = '?',
        help = 'path of file contains homologous exon pairs')
    parser.add_argument('--species1',
        nargs = '?',
        help = 'query species')
    parser.add_argument('--species2',
        nargs = '?',
        help = 'subject species')
    parser.add_argument('--isocom1',
        nargs = '?',
        default = None,
        help = 'path of file contains isoform information of species1')
    parser.add_argument('--cdscom1',
        nargs = '?',
        default = None,
        help = 'path of file contains cds information of species1')
    parser.add_argument('--isocom2',
        nargs = '?',
        default = None,
        help = 'path of file contains isoform information of species2')
    parser.add_argument('--cdscom2',
        nargs='?',
        default = None,
        help = 'path of file contains cds information of species2')
    parser.add_argument('--identhres',
        nargs = '?',
        default = 0.8,
        type = float,
        help = 'identity threshold of pairwise alignment')
    parser.add_argument('--coverthres',
        nargs = '?',
        default = 0.8,
        type = float,
        help = 'recoprocal coverage threshold of pairwise alignment')
    parser.add_argument('--match',
        nargs = '?',
        default = 2,
        type = int,
        help = 'match score of dynamic programming')
    parser.add_argument('--mismatch',
        nargs = '?',
        default = -2,
        type = int,
        help = 'mismatch score of dynamic programming')
    parser.add_argument('--gap',
        nargs = '?',
        default = -1,
        type = int,
        help = 'gapscore score of dynamic programming')
    parser.add_argument('--minexon',
        nargs = '?',
        default = 0,
        type = int,
        help = 'exons less than minexon size to be optimised after dynamic programming')

    parser.add_argument('--pnum',
        nargs = '?',
        default = 1,
        type = int,
        help = 'cpu core numbers to run EGIO')
    
    args = parser.parse_args()
    #================================
    orthogene = pd.read_table(str(args.orthog),header=0,sep='\t')
    blastexon = pd.read_table(str(args.blastn),header=0,sep='\t')

    isocom1 = pd.read_table(str(args.isocom1),header=0,sep='\t')
    isocom_proc1 = isocom1[isocom1['Orf'] != ".|."]
    cdscom1 = pd.read_table(str(args.cdscom1),header=0,sep='\t')

    isocom2 = pd.read_table(str(args.isocom2),header=0,sep='\t')
    isocom_proc2 = isocom2[isocom2['Orf'] != ".|."]
    cdscom2 = pd.read_table(str(args.cdscom2),header=0,sep='\t')
    #---------------------------------
    blastexon1 = pd.merge(blastexon,cdscom1[['EnsemblG',"ID"]],left_on=args.species1,right_on="ID")
    
    cdscomtmp2 = cdscom2[['EnsemblG',"ID"]]
    cdscomtmp2.columns = ['EnsemblG2',"ID2"]
    blastexon2 = pd.merge(blastexon1,cdscomtmp2,left_on=args.species2,right_on="ID2")
    
    blaste = blastexon2[[args.species1,args.species2,"EnsemblG","EnsemblG2"]]
    blaste = blaste.sort_values("EnsemblG")
    #---------------------------------
    print("organize blastn data")
    bedict = DFtoDict(blaste, key_base = "EnsemblG", header=True)
    print("organize isoform data of " + args.species1)
    isodict1 = DFtoDict(isocom_proc1, key_base = "EnsemblG", header=True)
    print("organize cds data of " + args.species1)
    cdsdict1 = DFtoDict(cdscom1, key_base = "EnsemblG", header=True)
    print("organize isoform data of " + args.species2)
    isodict2 = DFtoDict(isocom_proc2, key_base = "EnsemblG", header=True)
    print("organize cds data of " + args.species2)
    cdsdict2 = DFtoDict(cdscom2, key_base = "EnsemblG", header=True)
    ##------------------------------
    p = Pool(args.pnum)
    res_l = []
    for gi in range(0,len(orthogene)):
        gen1 = list(orthogene[args.species1])[gi]
        gen2 = list(orthogene[args.species2])[gi]

        if not(gen1 in isodict1):
            print("Error in analyais of " + gen1 + ": "+ gen1 + " is not recorded in the transcriptome or it is no longer a protein coding gene")
            continue
        
        if not(gen2 in isodict2):
            print("Error in analyais of " + gen1 + ": "+ gen2 + " is not recorded in the transcriptome or it is no longer a protein coding gene")
            continue

        if not(gen1 in bedict):
            bedata = {"dat":[],"coln":list(blaste.columns)}
        else:
            bedata = {"dat":bedict[gen1],"coln":list(blaste.columns)}
        
        data_in = { "blastn":bedata,
                    "iso1":{"dat":isodict1[gen1],"coln":list(isocom_proc1.columns)},
                    "cds1":{"dat":cdsdict1[gen1],"coln":list(cdscom1.columns)},
                    "iso2":{"dat":isodict2[gen2],"coln":list(isocom_proc2.columns)},
                    "cds2":{"dat":cdsdict2[gen2],"coln":list(cdscom2.columns)},
                    "match":args.match,
                    "mismatch":args.mismatch,
                    "gap":args.gap,
                    "identity":args.identhres,
                    "coverage":args.coverthres,
                    "mexon":args.minexon,
                    "species1":args.species1,
                    "species2":args.species2,
                    "gene1":gen1,
                    "gene2":gen2
                    }
        
        res = p.apply_async(EGIO,args=(data_in,))
        res_l.append(res)

    p.close()
    p.join()

    exoniden = [['Group',str(args.species1)+'EnsemblG',str(args.species1)+'Pos',str(args.species2)+'EnsemblG',str(args.species2)+'Pos','Iden','Type']]
    isoiden = [[str(args.species1),str(args.species2),'iso'+str(args.species1),'iso'+str(args.species2),'exoniden']]
    for res in res_l:
        datatmp = res.get()
        exoniden = exoniden + datatmp["oe"]["dat"]
        isoiden = isoiden + datatmp["ot"]["dat"]

    exonidendf = pd.DataFrame(exoniden)
    isoidendf = pd.DataFrame(isoiden)

    exonidendf.to_csv(os.getcwd()+'/ExonGroup_testpro_'+str(args.species1)+'_'+str(args.species2)+'.txt', sep='\t', header=False, index=False)
    isoidendf.to_csv(os.getcwd()+'/OrthoIso_testpro_'+str(args.species1)+'_'+str(args.species2)+'.txt', sep='\t', header=False, index=False)
  
