#!/bin/bash

#$ -q rcc-30d

# S.c. Gene Conversion
date
cd /N/dc2/scratch/megbehri/SAM_MURI/140808/plate_1
AR=(118)

for i in "${AR[@]}"
do
	cd Sample_${i}
	gunzip *
	cat *R1_001.fastq > Sample${i}_R1.fastq
	cat *R2_001.fastq > Sample${i}_R2.fastq
        rm -rf Pindel*
	
	# clean and trim
	echo "#!/bin/bash" > SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "#$ -l vmem=50gb walltime=48:00:00" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "cd /N/dc2/scratch/megbehri/SAM_MURI/140808/plate_1/Sample_${i}" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "time cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -o Sample${i}_R1_rmadapter.fastq -p Sample${i}_R2_rmadapter.fastq Sample${i}_R1.fastq Sample${i}_R2.fastq" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "time cutadapt -q 15,10  -o Sample${i}_R1_filtered.fastq -p Sample${i}_R2_filtered.fastq Sample${i}_R1_rmadapter.fastq Sample${i}_R2_rmadapter.fastq" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "time cutadapt -u 15 -o Sample${i}_R1_trimmed.fastq -p Sample${i}_R2_trimmed.fastq Sample${i}_R1_filtered.fastq Sample${i}_R2_filtered.fastq" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "qsub -l walltime=20:00:00,vmem=64gb,nodes=1:ppn=4 SamMuri${i}_Assemble.sh" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "mkdir fastqc" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "fastqc -o fastqc/ Sample${i}_R1_trimmed.fastq" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "fastqc -o fastqc/ Sample${i}_R2_trimmed.fastq" >> SamMuri${i}_qc_clean.sh
	echo "" >> SamMuri${i}_qc_clean.sh
	echo "exit" >> SamMuri${i}_qc_clean.sh

#### Assembly and Pindel ####

	echo "#!/bin/bash" > SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "cd /N/dc2/scratch/megbehri/SAM_MURI/140808/plate_1/Sample_${i}/" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "echo \"BWA\" >&2" >> SamMuri${i}_Assemble.sh
	echo "time bwa mem /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna Sample${i}_R1_trimmed.fastq Sample${i}_R2_trimmed.fastq > Sample${i}.sam" >> SamMuri${i}_Assemble.sh
   	echo "" >> SamMuri${i}_Assemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb SamMuri${i}_compress_fastq.sh" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "echo \"SAMTOOLS\" >&2" >> SamMuri${i}_Assemble.sh
	echo "samtools view -bS -T /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna Sample${i}.sam > Sample${i}.bam" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "samtools sort Sample${i}.bam Sample${i}.sorted" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "samtools index Sample${i}.sorted.bam" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "samtools sort -n -m 4000000000 Sample${i}.sorted.bam Sample${i}.sorted.pindel" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb SamMuri${i}_RetroSeq.sh" >> SamMuri${i}_Assemble.sh
	echo "samtools index Sample${i}.sorted.pindel.bam" >> SamMuri${i}_Assemble.sh
	###Below Added to test if GATK helps clean up reads for Pindel	
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/AddOrReplaceReadGroups.jar I=Sample${i}.sorted.bam O=Sample${i}.sorted.fixed.bam SORT_ORDER=coordinate RGID=GeneConv RGLB=bar RGPL=illumina RGSM=Sample${i} RGPU=6 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> SamMuri${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=Sample${i}.sorted.fixed.bam O=Sample${i}.sorted.fixed.marked.bam M=Sample${i}.metrics CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT" >> SamMuri${i}_Assemble.sh
	echo "GATK" >&2
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fa -I Sample${i}.sorted.fixed.marked.bam -T RealignerTargetCreator -o Sample${i}.intervals" >> SamMuri${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fa -I Sample${i}.sorted.fixed.marked.bam -T IndelRealigner -targetIntervals Sample${i}.intervals --filter_bases_not_stored -o Sample${i}.sorted.fixed.marked.realigned.bam" >> SamMuri${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fa -I Sample${i}.sorted.fixed.marked.realigned.bam -glm BOTH -rf BadCigar -o Sample${i}.sorted.fixed.marked.realigned.vcf" >> SamMuri${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -T BaseRecalibrator -I Sample${i}.sorted.fixed.marked.realigned.bam -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fa -rf BadCigar --filter_bases_not_stored -knownSites Sample${i}.sorted.fixed.marked.realigned.vcf -o Sample${i}.recal_data.grp" >> SamMuri${i}_Assemble.sh
	echo "java -Xmx2g -classpath "/N/soft/rhel6/picard/picard-tools-1.107/" -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar  -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fa -I Sample${i}.sorted.fixed.marked.realigned.bam -T PrintReads -rf BadCigar  -o Sample${i}.mapped.bam -BQSR Sample${i}.recal_data.grp" >> SamMuri${i}_Assemble.sh
	#java -Xmx2g -classpath /N/soft/rhel6/picard/picard-tools-1.107/ -jar /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar -nt 4 -T UnifiedGenotyper -R ../../Ecoli_RefGenome/Ecoli_K12_MG1655.fa -I Sample${i}.sorted.fixed.marked.realigned.bam -rf BadCigar -o Sample${i}.sorted.fixed.marked.realigned.vcf" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "echo \"PINDEL\" >&2" >> SamMuri${i}_Assemble.sh
	#Below changed to use output from GATK
	echo "Sample${i}.sorted.fixed.marked.realigned.bam	300	Sample${i}" > Sample${i}_Pindel.config.txt
	echo "" >> SamMuri${i}_Assemble.sh
	echo "mkdir Pindel4" >> SamMuri${i}_Assemble.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel -f /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -i Sample${i}_Pindel.config.txt -l -k -M 4 -o Sample${i}.pindel4.results" >> SamMuri${i}_Assemble.sh
	echo " mv *.pindel4.* Pindel4/" >> SamMuri${i}_Assemble.sh
	echo "" >> SamMuri${i}_Assemble.sh
	echo "qsub -l walltime=2:00:00,vmem=20gb SamMuri${i}_convertPindel.sh" >> SamMuri${i}_Assemble.sh
	echo "exit" >> SamMuri${i}_Assemble.sh
	
 	# compress fastq again
	echo "#!/bin/bash" > SamMuri${i}_compress_fastq.sh
	echo "" >> SamMuri${i}_compress_fastq.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> SamMuri${i}_compress_fastq.sh
	echo "" >> SamMuri${i}_compress_fastq.sh
	echo "cd /N/dc2/scratch/megbehri/SAM_MURI/140808/plate_1/Sample_${i}/" >> SamMuri${i}_compress_fastq.sh
	echo "bzip2 Sample${i}_R*.fastq" >> SamMuri${i}_compress_fastq.sh
	echo "rm Sample${i}_R*_filtered.fastq" >> SamMuri${i}_compress_fastq.sh
	echo "rm Sample${i}_R*_rmadapter.fastq" >> SamMuri${i}_compress_fastq.sh
	echo "bzip2 Sample${i}_R*_trimmed.fastq" >> SamMuri${i}_compress_fastq.sh
	echo "" >> SamMuri${i}_compress_fastq.sh
	echo "exit" >> SamMuri${i}_compress_fastq.sh
	
	##### call novel insertion sequences######
	echo "#!/bin/bash" > SamMuri${i}_RetroSeq.sh
	echo "" >> SamMuri${i}_RetroSeq.sh
	echo "#$ -l vmem=50gb walltime=2:00:00 " >> SamMuri${i}_RetroSeq.sh
	echo "" >> SamMuri${i}_RetroSeq.sh
	echo "cd /N/dc2/scratch/megbehri/SAM_MURI/140808/plate_1/Sample_${i}/" >> SamMuri${i}_RetroSeq.sh
	echo "mkdir InsSeq" >> SamMuri${i}_RetroSeq.sh
	echo "cd InsSeq" >> SamMuri${i}_RetroSeq.sh
	echo "perl /N/dc2/scratch/megbehri/SAM_MURI/Tools/RetroSeq/bin/retroseq.pl -discover -eref /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/InsSeq/Ecoli_IS_RefFile.txt -bam ../Sample${i}.sorted.bam -output Sample${i}.IS.Reads -align" >> SamMuri${i}_RetroSeq.sh
	echo "perl /N/dc2/scratch/megbehri/SAM_MURI/Tools/RetroSeq/bin/retroseq.pl -call -bam ../Sample${i}.sorted.bam -ref /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -output Sample${i}.IS -input Sample${i}.IS.Reads -hets" >> SamMuri${i}_RetroSeq.sh
	echo "" >> SamMuri${i}_RetroSeq.sh
	echo "exit" >> SamMuri${i}_RetroSeq.sh
	
	##### convert pindel files to VCF #####
	echo "#!/bin/bash" > SamMuri${i}_convertPindel.sh
	echo "" >> SamMuri${i}_convertPindel.sh
	echo "#$ -l vmem=50gb walltime=48:00:00 " >> SamMuri${i}_convertPindel.sh
	echo "" >> SamMuri${i}_convertPindel.sh
	echo "cd /N/dc2/scratch/megbehri/SAM_MURI/140808/plate_1/Sample_${i}/Pindel4" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_INV -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.INV.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_BP -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.BP.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_CloseEndMapped -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.CEM.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_LI -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.LI.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_RP -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.RP.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_SI -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.SI.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_TD -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.TD.vcf" >> SamMuri${i}_convertPindel.sh
	echo "/N/dc2/scratch/megbehri/SAM_MURI/Tools/pindel/pindel2vcf -p Sample${i}.pindel4.results_D -r /N/dc2/scratch/megbehri/SAM_MURI/140808/Ecoli_RefGenome/Ecoli_K12_MG1655.fna -R Ecoli_MG1655 -d 06182015 -v Sample${i}.pindel4.D.vcf" >> SamMuri${i}_convertPindel.sh

	echo "" >> SamMuri${i}_convertPindel.sh
	echo "exit" >> SamMuri${i}_convertPindel.sh
	
	chmod u+x SamMuri${i}_qc_clean.sh
	chmod u+x SamMuri${i}_Assemble.sh
	chmod u+x SamMuri${i}_compress_fastq.sh
	chmod u+x SamMuri${i}_convertPindel.sh
	chmod u+x SamMuri${i}_RetroSeq.sh
	
	qsub -l walltime=2:00:00,vmem=20gb SamMuri${i}_qc_clean.sh
	cd ..
	
	done

