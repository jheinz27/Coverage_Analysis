import pybedtools
import argparse
import pysam

parser = argparse.ArgumentParser(description='Check coverages in given region')
parser.add_argument('-r', metavar='<input.bam>', required=True, help='reference alignment file')
parser.add_argument('-t', metavar='<input.bam>', required=True, help='target alignment file' )
parser.add_argument('-reg', metavar='<roi.txt>', required=False, help='Copy number variation region')
parser.add_argument('-o', metavar='<out.txt>', required=True, help='outfile name')
parser.add_argument('-rreads',metavar= '<input.bam>', required= False, default = 'na', help = 'reference reads')
parser.add_argument('--index', action= "store_true", help= 'will index bam files if index files have not been created')


args = parser.parse_args()
refFile = args.r
targetFile = args.t
regions = args.reg
outname = args.o
rreads = args.rreads
index = args.index


def updateCov(dictToUse, chrm, pos,cover):
    for k in dictToUse.keys():
        keysInfo = k.split(';')

        if keysInfo[0] == chrm:

            if int(keysInfo[1]) <= pos < int(keysInfo[2]):
                v = dictToUse[k]
                temp = int(v) + int(cover)
                dictToUse[k] = temp


def getCovOutput(tool, dicttoUpdate):
    totBases = 0
    totCov = 0
    sample = []
    i = 1
    xcovs = 0
    xBases = 0
    dict_chrmCovs = {}
    for t in tool.fn:

        l = t.strip()
        s = l.split('\t')

        loc = int(s[1])
        totBases += 1
        cov = int(s[2])

        totCov += cov
        curChr = s[0]
        updateCov(dicttoUpdate, curChr, loc, cov)

        if curChr in dict_chrmCovs.keys():
            vals = dict_chrmCovs[curChr]
            bases = vals[0] + 1
            coverages = vals[1] + cov
            dict_chrmCovs[curChr] = [bases,coverages]
        else:
            dict_chrmCovs[curChr] = [1, cov]
    #    if curChr == 'chrX':
     #       xcovs += cov
      #      xBases += 1

        i += 1

    #aveCov = totCov / totBases
    #aveXCov = xcovs / xBases
    for key in dict_chrmCovs.keys():
        v = dict_chrmCovs[key]
        dict_chrmCovs[key] = v[1] / v[0]


    return dict_chrmCovs

#mapping Qualities
def getAveMQ(inputSamFile) :
    counts = 0
    totmq = 0
    for read in inputSamFile.fetch():
        counts += 1
        mq = read.mapping_quality
        totmq += int(mq)

    aveMQ = totmq / counts
    return aveMQ


def getRegionAveMQ(chrm, start, end, file):
    regionCounts = 0
    regionSum = 0
    for read in file.fetch(chrm, start, end):
        regionCounts += 1
        regionSum += int(read.mapping_quality)

    regionAve = regionSum / regionCounts
    return regionAve


def getStrings(info, cov_dict, dict, MQ, samfile):

    inf = info.split(';')
    chrm = inf[0]
    chrm_cov = cov_dict[chrm]
    start = int(inf[1])
    end = int(inf[2])
    cov = dict[info]
    regCov = cov / (int(end) - int(start))
    norm = regCov / chrm_cov
    #if isRef:
    roiMQ = getRegionAveMQ(inf[0],start,end, samfile)
    #else:
     #   roiMQ = getRegionAveMQ(inf[0], start, end, tarSamfile)
    normMQ = roiMQ / MQ
    loc = inf[0] + ':'+ str(inf[1]) + '-' + str(inf[2])
    outString =  loc + '\t' + str(regCov) + '\t' + str(round(norm,4)) + '\t' + str(round(normMQ,4))
    return outString

#function calls

if index:
    pysam.index(refFile)
    pysam.index(targetFile)

a = pybedtools.BedTool(refFile)
ref_covs = a.genome_coverage(d=True, stream= True, ibam= refFile)

b = pybedtools.BedTool(targetFile)
tar_covs = a.genome_coverage(d=True, stream= True, ibam= targetFile)

dict_roi_ref = {}
dict_roi_tar = {}
dict_tar2ref = {}

with open(regions,'r') as g:
    for line in g:
        if not line.startswith('#'):
            l = line.strip()
            s = l.split('\t')
            chr = s[0]
            refkey = chr + ';'+ s[1] + ';' + s[2]
            targetKey = s[3] + ';'+ s[4] + ';' + s[5]
            dict_roi_ref[refkey] = 0
            dict_roi_tar[targetKey] = 0
            dict_tar2ref[targetKey] = refkey

dict_roi_ref2ref = dict_roi_ref.copy()

refSamfile = pysam.AlignmentFile(refFile, "rb")
tarSamfile =  pysam.AlignmentFile(targetFile, "rb")

refAveCov = getCovOutput(ref_covs, dict_roi_ref)
print('Average Coverge Reference: ')
print(dict_roi_ref)

ref2refAveCov ={}

rr = False
if rreads != 'na':
    if index:
        pysam.index(rreads)
    ref2refSamFile = pysam.AlignmentFile(rreads, "rb")
    c = pybedtools.BedTool(rreads)
    ref2ref_covs = c.genome_coverage(d=True, stream=True, ibam= rreads)
    ref2refAveCov = getCovOutput(ref2ref_covs, dict_roi_ref2ref)
    print('ref2ref covs')
    print(dict_roi_ref2ref)
    ref2refAveMQ = getAveMQ(ref2refSamFile)
    rr = True


refAveMQ = getAveMQ(refSamfile)
print('Average Whole Genome Mapping Quality Reference: ' + str(refAveMQ))

tarAveCov= getCovOutput(tar_covs,dict_roi_tar)
print('Average Coverge ref')
print(refAveCov)


tarAveMQ = getAveMQ(tarSamfile)
print('Average Whole Genome Mapping Quality Reference: ' + str(tarAveMQ))

with open(outname, 'w') as out:
    print('dict_infos')
    for k in dict_tar2ref.keys():
        ref_key = dict_tar2ref[k]

        refInfo = getStrings(ref_key, refAveCov, dict_roi_ref, refAveMQ, refSamfile)
        tarInfo = getStrings(k, tarAveCov, dict_roi_tar, tarAveMQ, tarSamfile)
        if rr:
            ref2refInfo = getStrings(ref_key, ref2refAveCov, dict_roi_ref2ref, ref2refAveMQ, ref2refSamFile)
            out.write(ref2refInfo + '\t' + refInfo + '\t' + tarInfo + '\n')
        else:
            out.write(refInfo + '\t' + tarInfo + '\n')