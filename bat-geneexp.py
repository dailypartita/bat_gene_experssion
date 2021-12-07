#!/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/software/.conda/anaconda-3-2021./bin/python
#%% initial import
import os
import sys
import argparse
import pandas as pd

os.system('conda init bash')
os.system('source /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/software/.conda/anaconda-3-2021./etc/profile.d/conda.sh')
os.system('conda activate /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/software/.conda/anaconda-3-2021./envs/rna-seq')

def file_finder(_root_dir, _pattern, _result_list):
    # search <file name>._pattern in _root_dir
    _items = os.listdir(_root_dir)
    for _item in _items:
        _path = os.path.join(_root_dir, _item)
        if os.path.isdir(_path):
            file_finder(_path, _pattern, _result_list)
        elif _path.split('/')[-1].split('.')[-1] == _pattern:
            _result_list.append(_path)
    return _result_list

parser = argparse.ArgumentParser()
parser.add_argument('-mode', help='<build> build hisat2 index to species in -db ||<main> filter & mapping & htseq stat || <summary> make summary for htseq stat result by samples')
parser.add_argument('-db', help='path of database root', default='/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/database')
parser.add_argument('-sp', help='which species as referance sequence', default='Rousettus_aegyptiacus')
parser.add_argument('-gtf', help='path of gtf (for which species)', default='/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/database/Rousettus_aegyptiacus/GCF_014176215.1_mRouAeg1.p_genomic.gtf')
parser.add_argument('-sample_name', help='uniq sample name for fq1&2 input')
parser.add_argument('-r1', help='raw reads1')
parser.add_argument('-r2', help='raw reads2')
parser.add_argument('-htc_result_path', help='path of all .htc file root dir', default='/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/output')
parser.add_argument('-tmp', help='temp file stroge path', default='/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/temp')
parser.add_argument('-o', help='path to storge merged FPKM <merged.htc>', default='/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/output')
parser.add_argument('-p', help='processor to use', default=8)
args = parser.parse_args()

samtools = '/ldfssz1/ST_INFECTION/P20Z10200N0206_pathogendb/yangkaixin1/software/miniconda/bin/samtools'

#%% make index for each species
def mk_index():
    if not args.db:
        print('parameter <db> should be defined for step-0: build index')
        exp_handel()
    if args.sp:
        print('make index for %s'%args.sp)
        _sp = args.db + '/' + args.sp
        for j in os.listdir(_sp):
            j = _sp + '/' + j
            if j[-4:] == '.fna':
                reffa = j
            if j[-4:] == '.gtf':
                _gtf = j
        print('python /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/database/shell/makeindex.py -reffa %s -gtf %s -o %s'%(reffa, _gtf, _sp))
        os.system('python /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/database/shell/makeindex.py -reffa %s -gtf %s -o %s'%(reffa, _gtf, _sp))
    else:
        for _sp in os.listdir(args.db):
            _sp = os.path.join(args.db, _sp)
            if os.path.isdir(_sp) and _sp.split('/')[-1] != 'shell':
                for j in os.listdir(_sp):
                    j = os.path.join(_sp, j)
                    if j[-4:] == '.fna':
                        reffa = j
                    if j[-4:] == '.gtf':
                        _gtf = j
                print('python /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/database/shell/makeindex.py -reffa %s -gtf %s -o %s'%(reffa, _gtf, _sp))
                os.system('python /jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/database/shell/makeindex.py -reffa %s -gtf %s -o %s'%(reffa, _gtf, _sp))
    return 

#%% filter and mapping reads(PE) to ref genome
def fliter():
    print('fastp filter...')
    filter_path = args.tmp + '/1_filter/' + args.sample_name
    try:
        os.makedirs(filter_path)
    except:
        pass
    print('fastp -w %s -i %s -I %s --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --detect_adapter_for_pe -q 20 -u 20 --length_required 50 -n 2 -y -c -p --disable_trim_poly_g -R %s  -o %s/r1_filted.fq -O %s/r2_filted.fq -h %s/fastp_report.html -j %s/fastp_report.json'%(args.p, args.r1, args.r2, args.sample_name, filter_path, filter_path, filter_path, filter_path ))
    os.system('fastp -w %s -i %s -I %s --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --detect_adapter_for_pe -q 20 -u 20 --length_required 50 -n 2 -y -c -p --disable_trim_poly_g -R %s  -o %s/r1_filted.fq -O %s/r2_filted.fq -h %s/fastp_report.html -j %s/fastp_report.json'%(args.p, args.r1, args.r2, args.sample_name, filter_path, filter_path, filter_path, filter_path ))
    return filter_path + '/r1_filted.fq', filter_path + '/r2_filted.fq'

def mapping(_fr1, _fr2):
    print('hisat2 mapping...')
    bam_dir = args.tmp + '/2_alligment/' + args.sample_name
    try:
        os.makedirs(bam_dir)
    except:
        pass
    _index = args.db + '/' + args.sp + '/hisat2_index'
    sampath = bam_dir + '/' + args.sample_name
    print('alliging %s ...'%args.sample_name)
    print('hisat2 -p %s --dta -x %s -1 %s -2 %s -S %s.sam'%(args.p, _index, _fr1, _fr2, sampath))
    print('%s view -@ %s -S -b %s.sam > %s.bam'%(samtools, args.p, sampath, sampath))
    print('%s sort -@ %s %s.bam  %s_sorted'%(samtools, args.p, sampath, sampath))
    print('%s index -@ %s %s_sorted.bam'%(samtools, args.p, sampath))
    os.system('hisat2 -p %s --dta -x %s -1 %s -2 %s -S %s.sam'%(args.p, _index,  _fr1, _fr2, sampath))
    os.system('%s view -@ %s -S -b %s.sam > %s.bam'%(samtools, args.p, sampath, sampath))
    os.system('%s sort %s.bam -o %s_sorted.bam'%(samtools, sampath, sampath))
    os.system('%s index %s_sorted.bam'%(samtools, sampath))
    return sampath + '_sorted.bam'

#%% stat FPKM of mapping result
def htseq_count(_bam_path):
    _htc_path = args.o + '/' + args.sample_name
    print('htseq counting...')
    print('htseq-count -s no -r pos -f bam -n 20 %s %s > %s.htc'%(_bam_path, args.gtf, _htc_path))
    os.system('htseq-count -s no -r pos -f bam -n 20 %s %s > %s.htc'%(_bam_path, args.gtf, _htc_path))

#%% summary FPKM by samples
def sum_htc(_htc_root_path):
    _temp, _raw_data = [], []
    _htc_path_list = file_finder(_htc_root_path, 'htc', _temp)
    for htc in _htc_path_list:
        _df = pd.read_csv(htc, sep='\t', header=None, index_col = 0, names=['fpkm-%s'%htc.split('/')[-1].split('.')[0]])
        _raw_data.append(_df[:])
    _demo = _raw_data[0]
    for _data in _raw_data:
        _result = pd.merge(_demo, _data, how='inner', left_index=True, right_index=True)
    _result.to_csv('%s/merged_htc.csv'%args.o, index_label='gene')

#%%else chunks
def remove_temp():
    os.system('rm %s/1_filter/%s/*.fq -rf'%(args.tmp, args.sample_name))
    os.system('rm %s/2_alligment/%s/*.sam -rf'%(args.tmp, args.sample_name))

def exp_handel():
    os.system('python %s -h'%sys.argv[0])
    exit()

#%% process part
if args.mode == 'build':
    if not args.db:
        print('in build mode, -db mast be defined')
        exit()
    if args.sp:
        mk_index()
    else:
        mk_index()
    print('%s finish'%args.mode)

elif args.mode == 'main':
    if not args.tmp or not args.sample_name or not args.db or not args.gtf or not args.sp or not args.r1 or not args.r2:
        print('parameter <sample_name>, <r1>, <r2> should be defined.')
        exp_handel()
    #fr1, fr2 = fliter()
    bam_path = mapping(args.r1, args.r1)
    htseq_count(bam_path)
    remove_temp()
    print('%s finish'%args.mode)

elif args.mode == 'summary':
    if not args.htc_result_path:
        print('parameter <htc_result_path> should be defined.')
        exp_handel()
    sum_htc(args.o)
    print('%s finish'%args.mode)

else:
    exp_handel()
