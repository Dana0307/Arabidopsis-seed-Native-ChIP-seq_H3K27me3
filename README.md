# Step 1 install the following software
 -[fastp](https://github.com/OpenGene/fastp) v0.20.0 # 去接头以及数据过滤

 -[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.3.4.3  # 基因组比对

 -[samtools](https://www.htslib.org) v1.21  # 比对文件操作

 -[sambamba](https://lomereiter.github.io/sambamba/) v1.0.1  # 去重合排序 samtools也可以做 但是这个比较快

 -[deepTools](https://deeptools.readthedocs.io/en/latest/) v3.5.6 # 数据可视化，区间计数等

 -[MACS3](https://macs3-project.github.io/MACS/) #用来call peak

 -[Bedtools](https://plink.readthedocs.io/en/latest/bedtools_int/) #peak 注释


# Step 2 prepare the Arabidopsis genome file
  义访已经创建好，放在服务器共用地址 (/data/genome_data/TAIR10_Tan/bowtie2_index/TAIR10)，可以直接用。如果是共用没有的可以自己建，详见`bowtie2`教程 (如何创建的步骤待更新)

# Step 3 Run the following pipeline 

## 3.1 fastp 
   实验室Tn5建库所用接头序列都是一样，所以可以统一如下设置。当然fastp对PE测序可以自动检测接头也可以不设置
   file 指代文件名
   -w 是指定线程数 最大可指定16个，一般我用16个
   输出文件一般看html里面的质控

  ```
  fastp -w 16 \
--adapter_sequence=CTGTCTCTTATACACATCT --adapter_sequence_r2=CTGTCTCTTATACACATCT \
-i ./file_R1.fastq.gz \
-I ./file_R2.fastq.gz \
-o ./file_R1.fastp.fastq.gz \
-O ./file_R2.fastp.fastq.gz \
-j ./file.json \
-h ./file.html
```

## 3.2 Bowtie2
-p 是指定线程数, 一般我用64个
-x 指定索引路径. 义访以建立好。放在服务器共用地址了，可以直接用。
```
 bowtie2 -p 64 -x /data/genome_data/TAIR10_Tan/bowtie2_index/TAIR10 \
-S ./file.sam \
-1 ./file_R1.fastp.fastq.gz \
-2 ./file_R2.fastp.fastq.gz
```

## 3.3 samtools
   将生成的sam格式文件转化为bam文件，节省空间。同时过滤比对质量大于20
   记得加-h，不然会丢失头文件信息
```
view -@ 64 -q 20 -bh ./file.sam -o ./file.bam
```

## 3.4 sambamba 
   义访推荐这个软件，他觉得比samtools排序和去重都做的比较快
```
sambamba sort -m 30G -t 64 \
--tmpdir /data/tan/tan_seed/tmp/ \
./file.bam \
-o ./file_sorted.bam

sambamba markdup -r -t 64 \
--tmpdir /data/tan/tan_seed/tmp/ \
./file_sorted.bam \
./file_dedup.bam

samtools index ./file_dedup.bam
```

## 3.5 bamCoverage `deeptools`
   用deeptools工具包里面的bamCoverage子命令可以将bam文件转化为bw文件在IGV里面可视化
   里面有好几种不同的归一化方式常用的是RPKM和RPGC，他们都根据总测序量做了归一化，但是RPGC多了一步归一化到基因组大小
   RPKM比较方便
   bs 是指定bin的大小
   ChIP-seq一般指定--extendReads，但是RNA-seq不能指定这个参数，会导致失去剪切信息
   所以需要指定基因组大小，拟南芥基因组大小是119481543(官方推荐的)
```
bamCoverage -p 64 -bs 10 --extendReads --normalizeUsing RPGC \
--effectiveGenomeSize 119481543 \
-b ./file_dedup.bam \
-o ./file_dedup.bw
```
## 3.6 MACS3
  call peak 现在用的最多的还是MACS2或者MACS3，义访用的是MACS3
  H3K27me3一般用broad模式，需要指定--broad-cutoff 0.05
  --max-gap 50 是一个经验值，需要根据不同样本K27的特征去设定
```
macs3 callpeak \
--broad --max-gap 50 \
--broad-cutoff 0.05 \
--nolambda \
-f BAMPE \
-g 119481543 \
-q 0.05 \
-n file \
-c ./file_input.bam \
-t ./file_dedup.bam \
--outdir ./callpeak/
```

## 3.7 peak 筛选 （以H3K27me3 为例）
   peak筛选，对于H3K27me3我一般设置三个条件即：peak的宽度 > 200, foldchange > 3, q-value < 0.01，来筛选更可靠的peak
   其实这也是一个经验值可以根据自己的样本调整
   一般去掉叶绿体和线粒体上的信息
```
cat ./callpeak/file_peaks.broadPeak |awk '{if((($3-$2) > 200) && ($7 > 3) && ($9 > 2)) print $0}' OFS="\t" \
|grep -v ChrC |grep -v ChrM > ./callpeak/file_peaks_fli.bed
```

## 3.8 peak 注释 （以H3K27me3 为例）
 peak注释，可以用`Bedtools`注释，根据overlap信息简单粗暴非常适合K27这种定位在gene body区域的信号
 可以自己设置overlap的条件，可以调整的地方也很多
 这里是我种子里面设置的条件：peak和gene overlap的长度 > 100bp
 或者peak和gene overlap的长度站peak或者基因的50%

```
bedtools intersect -wo -a ./callpeak/file_peaks_fli.bed -b /data/genome_data/TAIR10_Tan/TAIR10.bed \
|awk '{if($NF>=100) print $10}' > ./PCG/file_peak_100bp.PCG

bedtools intersect -wo -f 0.5 -F 0.5 -e -a ./peakfile/file_peaks_fli.bed -b /data/genome_data/TAIR10_Tan/TAIR10.bed \
|cut -f 10 > ./PCG/file_peak_50%.PCG

sort ./PCG/file_peak_50%.PCG ./PCG/file_peak_50%.PCG |uniq > ./PCG/file_peak.PCG

grep -f ./PCG/file_peak.PCG /data/genome_data/TAIR10_Tan/TAIR10.bed \
> ./peakfile/file_peak_PCG.bed

rm ./PCG/file_peak_100bp.PCG ./PCG/file_peak_50%.PCG
```

## 3.9 ComputeMatrix plot heatmap
 注释到基因上之后就可以根据注释到的基因去画一个heatmap判断修饰的pattern，用到的也是deeptools工具包里面的子命令
 K27在基因上的分布一般用scale-regions模式
```
computeMatrix scale-regions -b 1000 -a 1000 -p 32 -m 2000 \
-R ./peakfile/file_peak_PCG.bed \
-S ./file_dedup.bam \
--skipZeros --missingDataAsZero \
-o ./file.gz

plotHeatmap \
-m ./file.gz \
-out ./file.pdf \
--averageType "mean" \
--startLabel "TSS" \
--endLabel "TES" \
--regionsLabel file_peak_PCG \
--yAxisLabel "RPGC" \
--colorMap Reds \
--heatmapWidth 6 \
--perGroup \
--heatmapHeight 15 \
--dpi 300
```

## multiBwigSummary plot Correlation
  数据质控，计算样本相关性也可以用deeptools工具包
  用multiBwigSummary计算基因组上每个bin上每个样本的值，然后用plotCorrelation计算样本间pearson相关系数
  一般跑完这个结合igv上track信息和heatmap的pattern，就大概知道数据质量咋样了
  其他具体分析可能需要学习R或者Python
```
multiBwigSummary bins -p 16 \
-b ./*.bw \
-out All_files_bin.npz --outRawCounts All_files_bin.tab

plotCorrelation \
--corData All_files_bin.npz \
--corMethod pearson \
--whatToPlot heatmap \
--plotFile Corrolation_of_heatmap.pdf
```

