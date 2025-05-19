#!/bin/bash

#The default working directory should be "Glial_Chimera_scRNA_2020" folder if not otherwise mentioned

module load samtools/1.9
module load deeptools/3.5.1
module load bedtools/2.30.0
module load ucsc/b1
module load kentutils/302.1.0


#folder with gene and enhancer annotation files
inFolder="data_for_import/CUT&Tag"

#bigwid folder
bigwigFolder="output/CUT&Tag/bigwig"

#output folder with plotting files
outFolder="output/CUT&Tag/genomicCoverage"

suffix="merged_average.bw"
bedList=("down" "up" "control")

###
# generate BED files of control/up/down-regulated genes and their putative enhancers, for plotting
###

for bed in "${bedList[@]}"
do
	geneList=$outFolder/"GPC.vs.PSC.GeneList.filter_PCT0.3_expr0.3_log2FC0.25_${bed}.txt"
	
	#generate bed file for the genes
	awk -F'\t' 'FNR==NR {lines[$1]; next} $4 in lines {print $0}' $geneList $inFolder/Homo_sapiens.GRCh38.106.transcript_proteinCoding_geneName.bed | sort -u > $outFolder/"${bed}_gene_fc0.25.bed"
	
	#extract enhancers for the genes, print GHid
	awk -F'\t' 'FNR==NR {lines[$1]; next} $2 in lines {print $0}' $geneList $inFolder/GeneHancer_AnnotSV_gene_association_scores_v5.18_elite.txt | awk '{print $1}' | sort -u > $outFolder/"${bed}_gene_fc0.25_enhancer.txt"
	#get coordinate BED file for the GHid
	awk -F'\t' 'FNR==NR {lines[$1]; next} $4 in lines {print $0}' $outFolder/"${bed}_gene_fc0.25_enhancer.txt" $inFolder/GeneHancer_AnnotSV_elements_v5.18_elite_enhancer.txt > $outFolder/"${bed}_gene_fc0.25_enhancer.bed"
	rm -f $outFolder/"${bed}_gene_fc0.25_enhancer.txt"
done


###
# generate bigwig files from bam files
# To do this step, please download the bam files from GEO database to "data_for_import/CUT&Tag/bam/" folder,
# please also make sure the file names match the pattern "WA09_CTd0_K27me3_1.bam" 
###

#loop through the bam files
for ifile in data_for_import/CUT&Tag/bam/*.bam
do
	sampleName=$(echo "$ifile" | awk -F'/' '{print $NF}' | sed -r 's/.bam//g')
	echo "Working on sample $sampleName, time point $ifolder"
		
	outputBigwig=$bigwigFolder/$sampleName
	
	#make bigwid files	
	samtools index $ifile
	bamCoverage -p 6 --binSize 10 --normalizeUsing RPKM --extendReads --bam $ifile -o "${outputBigwig}.bw"
	
done	


###
# merge bigwig files from replicates by average scores
#m(work inside output/CUT&Tag/bigwig folder)
###

#change working directory into the bigwig folder
cd $bigwigFolder

celltype="both"
timeList=("CTd0" "CTd120" "CTd180")
histoneList=("K27ac" "K27me3" "K4me3")

#loop through $timeList
for timepoint in "${timeList[@]}"
do	
	for histone in "${histoneList[@]}"
	do
		echo "Working on: $timepoint, $histone"
		
		#get the number of replicates
		nrep=$(ls -l *_${timepoint}_${histone}_?.bw | wc -l)
		
		if (($nrep>1))
		then
			echo "number of replicate: $nrep, merging multiple bigwig files"
			
			#merge bigwig (Merge together multiple bigWigs into a single output bedGraph, sum)
			bigWigMerge *_${timepoint}_${histone}_?.bw temp.bdg
			
			##bedClip (Remove lines from bed file that refer to off-chromosome places)
			bedClip temp.bdg $genomeSize temp.clip.bdg		
			##sort bed by chromosome and then by start position
			bedtools sort -i temp.clip.bdg > temp.sort.bdg			
			##merge bed (combines overlapping (by at least 1 bp) and/or bookended intervals into a single, “flattened” or “merged” interval)
			bedtools merge -i temp.sort.bdg -d -1 -c 4 -o mean > temp.merge.bdg
			##get mean value
			cat temp.merge.bdg | awk -v var="$nrep" '{OFS="\t"; print $1, $2, $3, $4/var}' > temp.norm.bdg
		
			#bedGraphToBigWig
			bedGraphToBigWig temp.norm.bdg $genomeSize ${celltype}_${timepoint}_${histone}_merged_average.bw

			rm temp*
		else
			echo "number of replicate: $nrep, no need to merge"
			cp *_${timepoint}_${histone}_?.bw ${celltype}_${timepoint}_${histone}_merged_average.bw
		fi
	done
done

#change working directory back out 
cd ../../../



###
# plot genomic coverage of the chromatin data for up/down/control genes with deeptools
# (work inside output/CUT&Tag/genomicCoverage folder)
###

#change working directory into the genomicCoverage folder for plotting
cd $outFolder

#plot H3K4me3 signal around TSS 

histone="K4me3"

for bed in "${bedList[@]}"
do
	geneBED="${bed}_gene_fc0.25.bed"
	computeMatrix reference-point --referencePoint TSS -p 12 --skipZeros --missingDataAsZero -b 2000 -a 2000 -R $geneBED -S ../bigwig/"${celltype}_CTd180_${histone}_${suffix}" ../bigwig/"${celltype}_CTd120_${histone}_${suffix}" ../bigwig/"${celltype}_CTd0_${histone}_${suffix}" --samplesLabel "${celltype}_CTd180_${histone}" "${celltype}_CTd120_${histone}" "${celltype}_CTd0_${histone}" -o "${celltype}_0_180_${histone}_matrix_TSS_${bed}.mat.gz" 
	#plot
	plotHeatmap --sortUsing sum -m $outFolder/"${celltype}_0_180_${histone}_matrix_TSS_${bed}.mat.gz" -out "${celltype}_0_180_${histone}_TSS_${bed}.pdf" --colorMap Reds --plotTitle ${histone}_${bed} --yMin 0 --yMax 600
	plotHeatmap --sortUsing sum -m $outFolder/"${celltype}_0_180_${histone}_matrix_TSS_${bed}.mat.gz" -out "${celltype}_0_180_${histone}_TSS_${bed}_overlap.pdf" --perGroup --colorMap Reds --plotTitle ${histone}_${bed} --yMin 0 --yMax 600
done


#plot H3K27me3 signal around TSS

histone="K27me3"

for bed in "${bedList[@]}"
do
	geneBED="${bed}_gene_fc0.25.bed"
	computeMatrix reference-point --referencePoint TSS -p 12 --skipZeros --missingDataAsZero -b 2000 -a 2000 -R $geneBED -S ../bigwig/"${celltype}_CTd180_${histone}_${suffix}" ../bigwig/"${celltype}_CTd120_${histone}_${suffix}" ../bigwig/"${celltype}_CTd0_${histone}_${suffix}" --samplesLabel "${celltype}_CTd180_${histone}" "${celltype}_CTd120_${histone}" "${celltype}_CTd0_${histone}" -o "${celltype}_0_180_${histone}_matrix_TSS_${bed}.mat.gz" 
	#plot
	plotHeatmap --sortUsing sum -m "${celltype}_0_180_${histone}_matrix_TSS_${bed}.mat.gz" -out "${celltype}_0_180_${histone}_TSS_${bed}.pdf" --colorMap Blues --plotTitle ${histone}_${bed} --yMin 0 --yMax 45
	plotHeatmap --sortUsing sum -m "${celltype}_0_180_${histone}_matrix_TSS_${bed}.mat.gz" -out "${celltype}_0_180_${histone}_TSS_${bed}_overlap.pdf" --perGroup --colorMap Blues --plotTitle ${histone}_${bed} --yMin 0 --yMax 45
done


##plot H3K27ac signal around enhancer regions

histone="K27ac"

for bed in "${bedList[@]}"
do 
	geneBED="${bed}_gene_fc0.25_enhancer.bed"
	computeMatrix scale-regions -p 12 --skipZeros --missingDataAsZero -b 1000 -a 1000 --regionBodyLength 2000 -R $geneBED -S ../bigwig/"${celltype}_CTd180_${histone}_${suffix}" ../bigwig/"${celltype}_CTd120_${histone}_${suffix}" ../bigwig/"${celltype}_CTd0_${histone}_${suffix}" --samplesLabel "${celltype}_CTd180_${histone}" "${celltype}_CTd120_${histone}" "${celltype}_CTd0_${histone}" -o "${celltype}_0_180_${histone}_matrix_ScaleGene_${bed}.mat.gz" 
	#plot
	plotHeatmap --sortUsing sum -m "${celltype}_0_180_${histone}_matrix_ScaleGene_${bed}.mat.gz" -out "${celltype}_0_180_${histone}_ScaleGene_${bed}.pdf" --colorMap Oranges --plotTitle ${histone}_${bed} --regionsLabel "Enhancer" --startLabel Start --endLabel End --xAxisLabel "distance (bp)" --yMin 0 --yMax 30
	plotHeatmap --sortUsing sum -m "${celltype}_0_180_${histone}_matrix_ScaleGene_${bed}.mat.gz" -out "${celltype}_0_180_${histone}_ScaleGene_${bed}_overlap.pdf" --perGroup --colorMap Oranges --plotTitle ${histone}_${bed} --regionsLabel "Enhancer" --startLabel Start --endLabel End --xAxisLabel "distance (bp)" --yMin 0 --yMax 30
done

#change working directory back out 
cd ../../../

