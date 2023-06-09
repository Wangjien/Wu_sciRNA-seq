# 整体分析

##  1 建立分析环境

### 推荐使用conda安装软件，没有安装conda可以在网络中寻找资源学习，并且配置好需要的下载channel。
* conda channel (~/.condarc文件)
```bash
channels:
  - bioconda
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch-lts: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
```

### 分析需要的软件

* python 3.0
* R 4.0
* samtools
* STAR
* trim_galore

### 分析需要的数据

* 参考基因组(*.fa)和gene信息(*.gtf)
* 原始测序数据（双端测序，*.fastq.gz）

### 安装需要的软件

创建分析需要的虚拟环境，并且在虚拟环境中安装需要软件

* python3 samtools Seurat(R包)

```bash
conda create -n sciRNA python=3.0
conda activate sciRNA
conda install -y STAR samtools cutadapter r-seurat
```

* TrimGalore [TrimGalore]: https://github.com/FelixKrueger/TrimGalore/releases

```bash
mkdir software
cd ./software
wget -O TrimGalore.tar.gz https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz 
tar -zxvf TrimGalore.tar.gz && rm -rf TrimGalore.tar.gz 
# 将软件添加到PATH中
 echo export PATH=\"`pwd`:\$PATH\" >> ~/.bashrc && source ~/.bashrc
```
* python 需要的模块
```bash 
python3 -m pip install numpy
python3 -m pip install pandas
python3 -m pip install pyfaxis
python3 -m pip install matplotlib
python3 -m pip install gzip
python3 -m pip install pysam
```

### 下载需要的参考文件
#### GeneCode

* Human（https://www.gencodegenes.org/human/）

[![p9BCgnH.md.png](https://s1.ax1x.com/2023/05/09/p9BCgnH.md.png)](https://imgse.com/i/p9BCgnH)

[![p9BChNt.md.png](https://s1.ax1x.com/2023/05/09/p9BChNt.md.png)](https://imgse.com/i/p9BChNt)

```bash
mkdir -p ~/Reference/{GeneCode,UCSC}/{Human,Mouse}/{Anno,Fasta}
## 下载GeneCode中的人类注释信息
wget -P ~/Reference/GeneCode/Human/Anno/ \
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
## 下载GeneCode中的人类基因组信息
wget -P ~/Reference/GeneCode/Human/Fasta/ \
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz
```

* Mouse（[GENCODE - Mouse Release M32 (gencodegenes.org)](https://www.gencodegenes.org/mouse/)）

[![p9BPuvD.md.png](https://s1.ax1x.com/2023/05/09/p9BPuvD.md.png)](https://imgse.com/i/p9BPuvD)

```bash
## 下载GeneCode中的小鼠注释信息
wget -P ~/Reference/GeneCode/Mouse/Anno/ \
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz
## 下载GeneCode中的小鼠基因组信息
wget -P ~/Reference/GeneCode/Mouse/Fasta/ \
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.genome.fa.gz
```
#### UCSC
[![p9DeoSs.md.png](https://s1.ax1x.com/2023/05/10/p9DeoSs.md.png)](https://imgse.com/i/p9DeoSs)

### 生成STAR比对需要的Index文件

```bash
#!!!! 注意：--sjdboverhang 的值是根据测序Reads的长度-1得到的
## 人类数据
mkdir -p ~/Reference/Star_index/{STAR_Human,STAR_Mouse}
STAR --runMode genomeGenerate \
     --runThreadN 10 \
     --genomeDir ~/Reference/Star_index/STAR_Human/ \
     --genomeFastaFiles  \
     --sjdbGTFfile   \
     --sjdbOverhang 149 # 双端150bp
     
## 小鼠数据
STAR --runMode genomeGenerate \
     --runThreadN 2 \
     --genomeDir ~/Reference/Star_index/STAR_Mouse/ \
     --genomeFastaFiles ~/Reference/GeneCode/Mouse/Fasta/GRCm39.genome.fa \
     --sjdbGTFfile ~/Reference/GeneCode/Mouse/Anno/gencode.vM32.annotation.gtf \
     --sjdbOverhang 149 # 双端150bp
```

[![p9BiyeH.md.png](https://s1.ax1x.com/2023/05/09/p9BiyeH.md.png)](https://imgse.com/i/p9BiyeH)

