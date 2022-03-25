python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py preproc -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py mapping -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -c -jcs 2
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py collapse -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py refine -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -refine
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py find_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py visual_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -g Zm00001d050245
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py rank_iso -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py palen_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py go -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py report -basic -asp -geneStruc -palen -html -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_ont.cfg
