import argparse, os, sys, logging
import operator
import random
import shutil
import pysam
import urllib
import requests
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser(description='')

parser.add_argument('data', help='The data.compareNovel file output by compareNovel.py')
parser.add_argument('-g', '--gtex_data', default='/projects/dmacmillanprj2/polya/gtex/gtex_tissue_sample_ids', help='Gtex data file for gtex samples')
#parser.add_argument('',help='')
#parser.add_argument("-l", "--log", dest="logLevel", default='WARNING', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level. Default = "WARNING"')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))

args = parser.parse_args()

# Globally used variables
gtex_tissues = [
'Bladder',
'Brain',
'Breast',
'Colon',
'Esophagus',
'Kidney',
'Liver',
'Lung',
'Ovary',
'Pancreas',
'Prostate',
'Skin',
'Stomach',
'Thyroid',
'Uterus',
'Whole'
]

ccle_tissues = [
'autonomic_ganglia',
'biliary_tract',
'bone',
'breast',
'central_nervous_system',
'endometrium',
'haematopoietic_and_lymphoid_tissue',
'kidney',
'large_intestine',
'liver',
'lung',
'oesophagus',
'ovary',
'pancreas',
'pleura',
'prostate',
'salivary_gland',
'skin',
'small_intestine',
'soft_tissue',
'stomach',
'thyroid',
'upper_aerodigestive_tract',
'urinary_tract'
]

gtex_to_ccle_tissues = {
'Bladder': None,
'Brain': 'central_nervous_system',
'Breast': 'breast',
'Colon': 'large_intestine',
'Esophagus': 'oesophagus',
'Kidney': 'kidney',
'Liver': 'liver',
'Lung': 'lung',
'Ovary': 'ovary',
'Pancreas': 'pancreas',
'Prostate': 'prostate',
'Skin': 'skin',
'Stomach': 'stomach',
'Thyroid': 'thyroid',
'Uterus': None,
'Whole': None
}

ccle_to_gtex_tissues = {v:k for k,v in gtex_to_ccle_tissues.items()}

operators = {
             'lt': operator.lt,
             'gt': operator.gt,
             'ge': operator.ge,
             'le': operator.le,
             'eq': operator.eq,
             'ne': operator.ne,
             }
webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
ucsc_url = 'http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hgt.customText=http://bcgsc.ca/downloads/dmacmillan/'

class Gtex:

    def __init__(self, tissue, gender, _id, gtex, path):
        self.tissue = tissue
        self.gender = gender
        self._id = _id
        self.gtex = gtex
        self.path = path

class Data:

    def __init__(self, chrom, gene, tissue, dist, pas, _id, cline, cs, score, med_left, med_right, med_diff, utr_attr):
        self.chrom = chrom
        self.gene = gene
        self.tissue = tissue
        self.dist = dist
        self.pas = pas
        self._id = _id
        self.cline = cline
        self.cs = cs
        self.score = score
        self.med_left = med_left
        self.med_right = med_right
        self.med_diff = med_diff
        self.utr_attr = utr_attr

output_columns = ['CHROM', 'GENE', 'TISSUE', 'DIST_FROM_ANNOT',
                  'PAS', 'ID', 'CELL_LINE', 'CLEAVAGE_SITE_CENTROID',
                  'SCORE', 'MEDIAN_LEFT', 'MEDIAN_RIGHT',
                  'MED_DIFF', 'CLOSEST_UTR3_ATTR']

if not os.path.isdir(args.outdir):
    try:
        os.makedirs(args.outdir)
    except OSError:
        pass

# Logging
#logging.basicConfig(filename=os.path.join(args.outdir, 'logfile'), level=getattr(logging, args.logLevel))

def readData(data_path, header=True):
    data = []
    with open(data_path, 'r') as f:
        if header:
            f.readline()
        for line in f:
            line = Data(*line.strip().split('\t'))
            data.append(line)
    return data

# kwargs = attribute: [value, operator]
def filterData(data, operators, **kwargs):
    result = []
    for dat in data:
        skip = False
        for filt in kwargs:
            if ( operators[kwargs[filt][0]]( float(getattr(dat, filt)), kwargs[filt][1] ) ):
                skip = True
                break
        if skip:
            continue
        result.append(dat)
    return result

def groupData(data):
    result = {}
    for dat in data:
        if dat.cs not in result:
            result[dat.cs] = []
        result[dat.cs].append(dat)
    return result

def removeMultipleSameSamples(gdata):
    result = {}
    for cs in gdata:
        num_unique_samples = len(set([x._id for x in gdata[cs]]))
        num_sites = len(gdata[cs])
        if (num_sites > 1) and (num_unique_samples == 1):
            continue
        result[cs] = gdata[cs]
    return result

def generateBrowserTrack(chrom, cs, window = 2000):
    start = cs - window
    end = cs + window
    browser = 'browser position {}:{}-{}\n'.format(chrom, start, end)
    browser += 'browser hide all\n'
    browser += 'browser pack knownGene\n'
    browser += 'browser pack refGene\n'
    browser += 'browser pack acembly\n'
    browser += 'browser pack ensGene\n'
    return browser

def generateKleatTrack(ccle_id, strand):
    webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
    if strand == '+':
        path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added/{}.plus.bg'.format(ccle_id)
    else:
        path = '/projects/dmacmillanprj2/polya/ccle/filteredKleats/kleats_plus_added/{}.minus.bg'.format(ccle_id)
    #if not os.path.isfile(os.path.join(webdir, os.path.basename(path))):    
        #sprint('\rCopying KLEAT track for {} ...'.format(ccle_id))
        #shutil.copy(path, webdir)
    f = open(path, 'r')
    data = ('').join(f.readlines()[1:])
    f.close()
    return 'track type=bedGraph name={}_novel_site description="Cleavage Sites for {} transcripts in {}" color="0,0,255" visibility=full\n{}\n'.format(ccle_id, strand, ccle_id, data)

def generateCcleTracks(*ccles):
    result = ''
    webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
    for ccle in ccles:
        colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
        ccle_bw_path = '/projects/btl/polya/ccle_bigwigs/{}.bw'.format(ccle)
        if not os.path.isfile('{}{}.bw'.format(webdir, ccle)):    
            #sprint('\rCopying bigwig track for {} ...'.format(ccle))
            try:
                shutil.copy(ccle_bw_path, webdir)
            except IOError:
                continue
        big_data_url = 'http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(ccle)
        track_line = 'track type=bigWig name="{}" description="genomecov bigwig track for {}" color="{}" visibility=full bigDataUrl={}\n'.format(ccle, ccle, colour, big_data_url)
        result += track_line
    return result

def generateGtexTracks(*gtexs):
    result = ''
    webdir = '/gsc/www/bcgsc.ca/downloads/dmacmillan/'
    for gtex in gtexs:
        colour = '{},{},{}'.format(random.randint(0,255),random.randint(0,255),random.randint(0,255))
        gtex_bw_path = '/projects/btl/polya/gtex_bigwigs/{}.bw'.format(gtex)
        if not os.path.isfile('{}{}.bw'.format(webdir, gtex)):    
            #sprint('\rCopying bigwig track for {} ...'.format(gtex))
            try:
                shutil.copy(gtex_bw_path, webdir)
            except IOError:
                continue
        big_data_url = 'http://bcgsc.ca/downloads/dmacmillan/{}.bw'.format(gtex)
        track_line = 'track type=bigWig name="{}" description="genomecov bigwig track for {}" color="{}" visibility=full bigDataUrl={}\n'.format(gtex, gtex, colour, big_data_url)
        result += track_line
    return result

def createUcscUrl(name, to_write, webdir, ucsc_url):
    outpath = os.path.join(webdir, name)
    with open(outpath, 'w') as f:
        f.write(to_write)
    return ucsc_url + name + '&hgt.psOutput=on'

def readGtex(gtex):
    result = []
    with open(gtex, 'r') as f:
        for line in f:
            line = Gtex(*line.strip().split('\t'))
            result.append(line)
    return result

def generateNovelTrack(single_data):
    result = 'track name="Novel Site" description="Novel site {}" color="255,50,50" visibility=full\n'.format(single_data.cs)
    result += '{}\t{}\t{}\t{}\n'.format(single_data.chrom, int(single_data.cs)-1, single_data.cs, single_data.score)
    return result

def createImage(single_data, ccle_id_out, gtex_id, webdir, ucsc_url, outdir):
    browser = generateBrowserTrack(single_data.chrom, int(single_data.cs))
    ccle_tracks = generateCcleTracks(single_data._id, ccle_id_out)
    gtex_tracks = generateGtexTracks(gtex_id)
    novel_track = generateNovelTrack(single_data)
    to_write = browser + ccle_tracks + novel_track + gtex_tracks
    name ='{}_{}'.format(single_data._id, single_data.cs)
    url = createUcscUrl(name, to_write, webdir, ucsc_url)
    page = requests.get(url)
    soup = BeautifulSoup(page.content, 'html.parser')
    pdf = soup.find_all('a')[-7]['href']
    pdf = 'http://genome.ucsc.edu/cgi-bin/' + pdf
    pdf_out = os.path.join(outdir, name + '.pdf')
    urllib.urlretrieve(pdf, pdf_out)
    return pdf_out

def getCcleNotMatching(single_data, agdata):
    for cs in set(agdata).difference(single_data.cs):
        for dat in agdata[cs]:
            if (dat.tissue != single_data.tissue):
                return dat._id

# kwargs = attribute: [operator, value]
def getGtex(gtex, operators, **kwargs):
    for g in gtex:
        for k in kwargs:
#            print getattr(g,k)
#            print kwargs[k][0]
#            print kwargs[k][1]
            if k == 'tissue':
                try:
                    tissue = gtex_to_ccle_tissues[g.tissue]
                except KeyError:
                    print 'no gtex tissue: {}'.format(tissue)
                    continue
                if (operators[kwargs[k][0]](tissue, kwargs[k][1])):
                    return g._id
            if (operators[kwargs[k][0]](getattr(g, k), kwargs[k][1])):
                return g._id

all_gtex = readGtex(args.gtex_data)
data = readData(args.data)
data = data[1:]
#print len(data)
#data = filterData(data, operators, med_diff=['lt', 80])
#print len(data)
gdata = groupData(data)
agdata = removeMultipleSameSamples(gdata)

for cs in agdata:
    num_samples = len(agdata[cs])
    num_tissues = len(set([x.tissue for x in agdata[cs]]))
    if (num_samples > 1) and (num_tissues == 1):
        print '*'*20
        print cs
        single_data = agdata[cs][0]
        gtex_id = getGtex(all_gtex, operators, tissue=['eq', single_data.tissue])
        if not gtex_id:
            print 'No gtex could be found'
            continue
        ccle_out = getCcleNotMatching(single_data, agdata)
        print createImage(single_data, ccle_out, gtex_id, webdir, ucsc_url, args.outdir)
