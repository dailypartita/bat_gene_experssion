# python /home/kaixinyang/projects/AFRICA-GAM/expression/shell/bat-geneexp.py \
# -mode main \
# -r1 210607_SEQ022_FP100002544BL_L01_SP2104130637/FP100002544BL_L01_1048-2048_1.fq.gz \
# -r2 210607_SEQ022_FP100002544BL_L01_SP2104130637/FP100002544BL_L01_1048-2048_2.fq.gz \
# -sample_name R2103020425_1_2021-06-09 #for Rousettus_aegyptiacus only, other sp need add <-gtf> and <-sp>
python /home/kaixinyang/projects/AFRICA-GAM/expression/shell/bat-geneexp.py \
-mode main \
-r1 /210607_SEQ022_FP100002544BL_L01_SP2104130622/FP100002544BL_L01_1033-2033_1.fq.gz \
-r2 /210607_SEQ022_FP100002544BL_L01_SP2104130622/FP100002544BL_L01_1033-2033_2.fq.gz \
-sample_name R2103020410_1_2021-06-09
# -tmp /home/kaixinyang/projects/AFRICA-GAM/expression/temp-test
