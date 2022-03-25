python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py preproc -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py mapping -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg -jcs 2
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py collapse -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py refine -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg -refine
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py find_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py visual_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg -g Zm00001d050245
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py rank_iso -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py allele_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg
#python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py go -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg -tg test_data/targetGene1.txt,test_data/targetGene2.txt,test_data/targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3 -filterby pvalue -cutoff 0.1
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py report -basic -asp -geneStruc -asas -html -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo_hybrid.cfg
