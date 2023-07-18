#!/bin/bash

### Hardcoded paths to tools (picard, gatk) and files, adjust for personal use
### Hardcoded Spark options for resource allocation at remove duplicate step
### Hardcoded call to "bam_list.txt" to reference bams

#For MJL only, but ensure your java env is correct
#MJL tested: Java 1.7, picard 3.0, gatk 4.4
#sudo archlinux-java set java-17-openjdk

mkdir final_bams
mkdir dup_metrics

while read x; do
	echo $x

	#Convert sam output to bam
	java -jar /usr/share/picard-tools/picard.jar SortSam \
		-I ${x}_align.bam -O ${x}_sorted.bam -SO coordinate

	#Add read group header to bam file
	java -jar /usr/share/picard-tools/picard.jar AddOrReplaceReadGroups \
		-RGLB Lane1 -RGPL Illumina -RGPU TTAGGC -RGSM $x -I ${x}_sorted.bam -O ${x}_header.bam
	rm ${x}_sorted.bam #int file

	# Mark and Remove duplicates (Spark configurations can vary depending on your computational resources)
	# Also produces metrics file with info about coverage, etc.
	java -Xmx16793998336 -jar /usr/share/java/gatk/GenomeAnalysisTK.jar MarkDuplicatesSpark \
		-I ${x}_header.bam -O ${x}_rmdup.bam -M ${x}_dup_metrics.txt \
		## Set here for a personal computer with >4 CPU cores and >16 GB RAM (4 processes @ 4GB / process)
		--conf 'spark.executor.instances=4' --conf 'spark.executor.cores=4' --conf 'spark.executor.memory=4G'

	#Index new bam file
	java -jar /usr/share/picard-tools/picard.jar BuildBamIndex \
		INPUT=${x}_rmdup.bam

	mv ${x}_dupmet.txt dup_metrics/
	mv ${x}_rmdup.* final_bams/
done < bam_list.csv
