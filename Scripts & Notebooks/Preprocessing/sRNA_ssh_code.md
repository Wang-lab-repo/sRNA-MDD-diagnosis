# small RNA pipeline in ssh
# 0. 传输检查
```bash
cd /mnt/d/depression/nanjingqinghsaonian/Patients
find -type f -name "*_1.fq.gz" > name.list
cat name.list | while read id
do  
  echo "<=== ${id} ===>"
  md5sum ${id}
done > md5.log
```
# 1. qc_trim
```bash
cd /public3/home/a2s000171/xrz/depression
mkdir -p sequence trim fastqc clean

base_dir=/public3/home/a2s000171/xrz/depression
raw_data=$base_dir/sequence
find $raw_data -type f -name "*_1.fq.gz" > $base_dir/name.list
find $raw_data -type f -name "*_1.fq.gz"  -exec basename {} \; > $base_dir/short_name.list

cat name.list | while read id
do  
  echo "<=== ${id} ===>"
  md5sum ${id}
done > md5.log

# 多线程跑trim
# cat name.list | parallel -k -j 8 --no-run-if-empty --linebuffer "
#   fastqc -t 8 -o ./fastqc {}
#   trim_galore -q 25 --phred33 --fastqc --length 18 --stringency 3 -o ./trim {}
# "
# cat name.list | while read id
# do
#   file_name=$(basename "${id}" | cut -d '_' -f 1-2)
#   echo "${file_name%%.*}"
#   sed "s|{}|${id}|g" qc_trim.sh > "./trim/${file_name%%.*}_qc_trim.sh"
# done




cd /public3/home/a2s000171/xrz/depression
################################ 
vim qc_trim.sh
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

source /public3/soft/module/module.sh
module load gcc/8.3.0 mpi/intel/17.0.7-thc
export PATH=/public3/home/a2s000171/xrz/biosoft/parallel-20231122/install/bin:$PATH
export PATH=/public3/home/a2s000171/hc/package/anaconda3/bin:$PATH
export PATH=/public3/home/a2s000171/hc/package/TrimGalore-master:$PATH

base_dir=/public3/home/a2s000171/xrz/depression

cat name.list | parallel --no-run-if-empty --linebuffer -k -j 61 "
  file_name=\$(basename {})
  fastqc -o \"$base_dir/fastqc\" {}
  trim_galore -q 25 --phred33 --fastqc --length 18 --stringency 3 -o \"$base_dir/trim\" {} > \"$base_dir/trim/\${file_name%%.*}_trim.log\" 2>&1
"
################################ 

sbatch qc_trim.sh -p amd_256 -N 1 -n 64 -J srna_trim -o ./trim/trim.log  2>&1
```


# 2. SPORTS


```bash
# trim中fq.gz转fq 
cd /public3/home/a2s000171/xrz/depression/trim
for i in *_trimmed.fq.gz
do
  gzip -dc $i > ../clean/${i%%.*}.fq
done


# 写成脚本
cd /public3/home/a2s000171/xrz/depression
##############################
vim sports2.sh
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

source /public1/soft/modules/module.sh
module load R/3.6.1
module load anaconda/3-Python3.8.8-zsq
source activate py35

export PATH=/public3/home/a2s000171/hc/package/sports1.1-master/source:$PATH

base_dir=/public3/home/a2s000171/xrz/depression
clean_data=$base_dir/clean
genome_dir=/public3/home/a2s000171/hc/rna_seq/reference/Homo_sapiens
output_dir=$base_dir/output2

find $clean_data -name '*_trimmed.fq' > $output_dir/sports_dataset.txt

sports.pl -i $output_dir/sports_dataset.txt -p 62 \
-g $genome_dir/genome/hg38/genome -m $genome_dir/miRBase/miRBase -r $genome_dir/rRNAdb/human_rRNA \
-t $genome_dir/GtRNAdb/hg19/hg19-tRNAs -w $genome_dir/piRBase/piR_human \
-e $genome_dir/Ensembl/release-89/Homo_sapiens.GRCh38.ncrna \
-o $output_dir

####################################
sbatch -p amd_256 -N 1 -n 64 -J srna sports2.sh -o ./output2/sports.log 2>&1
# 第一次结果储存在output文件夹中，因为没有调用py35重新跑一遍，结果储存在output2文件夹中。  
```


# 3. 单独miRNA

```bash
cd /public3/home/a2s000171/xrz/depression
mkdir -p mirna
base_dir=/public3/home/a2s000171/xrz/depression
clean_data=$base_dir/clean
mirna_data=$base_dir/mirna_data
find $clean_data -name '*_trimmed.fq' > mirna_name.list
find $clean_data -name '*_trimmed.fq' -exec basename {} \;| sed 's/_trimmed\.fq$//'  > ./mirna/short_mirna_name.list


# 写成脚本
cd /public3/home/a2s000171/xrz/depression
##############################
vim mirna.sh
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

export PATH=/public3/home/a2s000171/hc/package/anaconda3/bin:$PATH

base_dir=/public3/home/a2s000171/xrz/depression
genome_dir=/public3/home/a2s000171/hc/rna_seq/reference/Homo_sapiens
mirdeep_dir=/public3/home/a2s000171/hc/package/anaconda3/pkgs/mirdeep2-2.0.1.3-hdfd78af_1/bin
clean_data=$base_dir/clean
mirna_data=$base_dir/mirna

cat short_mirna_name.list | parallel --no-run-if-empty --linebuffer -k -j 61 "
    perl $mirdeep_dir/mapper.pl $clean_data/{}_trimmed.fq -h -e -j -k AGATCGGAAGAGC \
    -l 18 -m -p $genome_dir/genome/hg38/genome -s $mirna_data/{}_mapped.fa -t $mirna_data/{}_vs_genome.arf -v
    perl $mirdeep_dir/quantifier.pl -p $genome_dir/miRBase/hairpin.fa -m $genome_dir/miRBase/mature.fa \
    -r $mirna_data/{}_mapped.fa -t hsa -y {}_mirna"
# -h 如果不是fasta，用该参数处理成fasta
# -e fastq: 表示输入文件是fastq
# -j 移除ATCGUNatcgun以外的字符
# -k 表示去除接头序列,Illumina universal adapter
# -l 18 剔除长度在18 bp以下的序列
# -m 合并相同的reads
# -p bowite索引;-s 处理后的read;-t 处理后比对文件

####################################
sbatch -p amd_256 -N 1 -n 64 mirna.sh -o mirna.log 2>&1
# 第一次结果储存在output文件夹中，因为没有调用py35重新跑一遍，结果储存在output2文件夹中。 
```

 
# 4. 单独piRNA
## 本地太慢了
1. 下载1.0.0版本的bowtie，与超算保持一致  
* conda
```bash
# 配置conda
cd /mnt/d/biosoft
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/bin/activate

conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
# conda config --show channels
conda create -n py310 python=3.10.13
# 3.10.13

https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip/download
```
* 直接下载
```bash
cd /mnt/d/biosoft
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip
unzip bowtie-1.0.0-linux-x86_64.zip
cd bowtie-1.0.0
chmod 777 bowtie
```

2. piRNA代码

```bash
# 安装软件
cd /mnt/d/biosoft
wget https://github.com/deweylab/RSEM/archive/v1.2.28.tar.gz && tar -zxf v1.2.28.tar.gz
cd RSEM-1.2.28
chmod +x rsem-calculate-expression

cd /mnt/d/depression
mkdir -p ./align ./pi_exp
find $clean_data -name '*_trimmed.fq' -exec basename {} \; | awk -F '_' '{print $1}' > ./piRNA_dataset.txt

vim pirna.sh
################################
#!/bin/bash

base_dir=/mnt/d/depression
clean_data=$base_dir/clean
align_out=$base_dir/align
pi_exp=$base_dir/pi_exp
ref_dir=/mnt/d/sRNA/reference/Homo_sapiens/piRBase/piR_human

cat $base_dir/piRNA_dataset.txt | parallel -k -j 6 --no-run-if-empty --linebuffer "
  bowtie "${ref_dir}" -v 1 -m 3 -q -S "${clean_data}/{}_1_trimmed.fq" \
  | samtools view -b -@ 4 > "${align_out}/{}.bam" 
  rsem-calculate-expression --no-bam-output --alignments \
  -p 4 --bam "${align_out}/{}.bam" "${ref_dir}" "${pi_exp}/{}.R1rsem" 
"

bash pirna.sh > ./pi_exp/log.txt 2>&1
```

## 超算代码

```bash
cd /public3/home/a2s000171/xrz/depression
mkdir -p ./pirna/align ./pirna/pi_exp

base_dir=/public3/home/a2s000171/xrz/depression
clean_data=$base_dir/clean
pirna_alin=$base_dir/pirna/align
find $clean_data -name '*_trimmed.fq' -exec basename {} \; | awk -F '_' '{print $1}' | sort > $pirna_alin/piRNA_dataset.txt


# 写成脚本
cd $base_dir
##############################
vim pirna.sh
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

export PATH=/public3/home/a2s000171/hc/package/anaconda3/bin:$PATH
source activate ts

base_dir=/public3/home/a2s000171/xrz/depression
clean_data=$base_dir/clean
pirna_alin=$base_dir/pirna/align
pi_exp=$base_dir/pirna/pi_exp
genome_dir=/public3/home/a2s000171/hc/rna_seq/reference/Homo_sapiens
ref_dir=$genome_dir/piRBase/piR_human

cat $pirna_alin/piRNA_dataset.txt | parallel --no-run-if-empty --linebuffer -k -j 16 "
  bowtie "${ref_dir}" -v 1 -m 3 -q -S "${clean_data}/{}_1_trimmed.fq" \
  | samtools view -b -@ 4 > "${pirna_alin}/{}.bam" 
  rsem-calculate-expression --no-bam-output --alignments \
  -p 4 --bam "${pirna_alin}/{}.bam" "${ref_dir}" "${pi_exp}/{}.R1rsem" 
  rm "${clean_data}/{}_1_trimmed.fq"
"
####################################
sbatch -p amd_256 -N 1 -n 64 pirna.sh -o ./pirna/pirna.log 2>&1
```



