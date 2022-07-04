run_dir="/data/CrazyHsu_data/Projects/isoseq/iFLAS_20220704"
python ${run_dir}/iflas.py preproc -cfg ${run_dir}/test_data/test_ont.cfg
python ${run_dir}/iflas.py mapping -cfg ${run_dir}/test_data/test_ont.cfg -c -jcs 5
python ${run_dir}/iflas.py collapse -cfg ${run_dir}/test_data/test_ont.cfg
python ${run_dir}/iflas.py refine -cfg ${run_dir}/test_data/test_ont.cfg -refine
python ${run_dir}/iflas.py find_as -cfg ${run_dir}/test_data/test_ont.cfg
python ${run_dir}/iflas.py visual_as -cfg ${run_dir}/test_data/test_ont.cfg -g Zm00001d050245
python ${run_dir}/iflas.py rank_iso -cfg ${run_dir}/test_data/test_ont.cfg
python ${run_dir}/iflas.py palen_as -cfg ${run_dir}/test_data/test_ont.cfg
python ${run_dir}/iflas.py go -cfg ${run_dir}/test_data/test_ont.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3
python ${run_dir}/iflas.py report -basic -asp -geneStruc -palen -html -cfg ${run_dir}/test_data/test_ont.cfg
