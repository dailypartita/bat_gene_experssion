### 蝙蝠基因表达矩阵计算

脚本根目录 `/jdfssz1/ST_HEALTH/P20Z10200N0206/kaixinyang/projects/AFRICA-GAM/expression/shell`

1. 建hisat2 INDEX

   * 参考 `1-work-build_all_index.sh`
2. 过滤reads和比对，计算每个样本的FPKM

   * 参考 `2-work-mapping_each_sample.sh`
3. 统计所有样本的FPKM

   * 参考 `3-work-make_summary.sh`
4. 计算每个样本的FPKM需要将 `2-work-mapping_each_sample.sh`按行拆开投递

   * 参考 `ignition.sh`

   > 以上说明只针对埃及果蝠(Rousettus_aegyptiacus)

5. 操作流程
   1. 首先用所有样本的sample_name，reads1路径，reads2路径（支持压缩的fq文件）生成 `2-work-mapping_each_sample.sh`
   2. 用 `ignition.sh`给定的参数（8核，2G）按行拆开投递`2-work-mapping_each_sample.sh`，该脚本自动清空所有temp中产生的中间大文件。
   3. 第二步跑完之后在测试节点上运行 `3-work-make_summary.sh`统计所有结果，生成结果在 `output/merged_htc.csv`目录下
   4. 如果有样本跑失败了直接重新投这个样本在 `2-work-mapping_each_sample.sh`中对应的那一行脚本即可，所有样本运行成功后执行第三步统计并输出结果
6. 权限
   * `temp` 和 `output` 777，其余可读可执行
7. 耗时
   * 每个样本耗时20min以内 
