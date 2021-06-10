#!/usr/bin/env bash

#BSUB -P insertionSite         # project code
#BSUB -J combine      # job name
#BSUB -W 10:00                # wall-clock time (hrs:mins)
#BSUB -n 1      # number of cpu
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -e errors.%J     # error file name in which %J is replaced by the job ID
#BSUB -o output.%J     # output file name in which %J is replaced by the job ID


parameterFile=$1

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


#fastq format
#  *_R1_*fastq.gz
#  *_R2_*fastq.gz

module load python/3.7.0
module load fastqc/0.11.5
module load bwa/0.7.17
module load samtools/1.10
module load bedtools/2.29.2
module load R/3.5.1
#module load cutadapt/1.2.1
#cutadapt runningn with 3.7.0


targetPredictionsR=$SCRIPT_DIR/target_gene_prediction.R
myjoinBin=$SCRIPT_DIR/myjoin


temporary="temporary"
mkdir -p $temporary
cutPrimerDir="$temporary/cutPrimer"

bwaDir="$temporary/bwa"
fastqcDir="$temporary/fastqc"
fastqcBeforeDir=$fastqcDir"/beforeCutAdapt"
newFastqDir="$temporary/newFastq"
bam2BedDir="$temporary/bam2Bed"
bed2peakDir="$temporary/bed2peak"
finalDir="finalResults"



mkdir -p $cutPrimerDir
mkdir -p $fastqcBeforeDir
mkdir -p $newFastqDir
mkdir -p $bwaDir
mkdir -p $bam2BedDir
mkdir -p $bed2peakDir
mkdir -p $finalDir
################################################
################################################



declare -A PARAMS

while IFS=$'\t' read -r -a varVal
do
 var=${varVal[0]}
 val=${varVal[1]}
 PARAMS[$var]=$val
done < $parameterFile



inputData=${PARAMS[inputFASTQ]}
 if [ -v $inputData ];then
     echo "inputData directory not defined"
     exit -1
 fi


bsubQueue=${PARAMS[bsubQueue]}
 if [ -v $bsubQueue ];then
     echo "bsubQueue not defined in parameter file"
     exit -1
 fi




scratch=${PARAMS[scratch]}
 if [ -v $scratch ];then
     echo "scratch not defined in parameter file"
     exit -1
 fi


bwaIndex=${PARAMS[bwaIndex]}
 if [ -v $bwaIndex ];then
     echo "bwaIndex not defined in parameter file"
     exit -1
 fi

 overlapDistances=${PARAMS[overlapDistances]}
 if [ -v $overlapDistances ];then
     echo "overlapDistances not defined in parameter file"
     exit -1
 fi


 inputDataDir=${PARAMS[inputDataDir]}
 if [ -v $inputDataDir ];then
     echo "inputDataDir not defined in parameter file"
     exit -1
 fi



 knownSites=${PARAMS[knownSites]}
 if [ -v $knownSites ];then
     echo "knownSites not defined in parameter file"
     exit -1
 fi


logDir="$scratch/logs"
mkdir -p $logDir

echo "Logs in $logDir"
#################################################
 ##############################################


runCutadapt(){

    r1Base=$1
    r2Base=$2

    echo $r1Base $r2Base

    bsub -K -q $bsubQueue  -J "cutP R1"   -M 8000 -R rusage[mem=8000] -oo $logDir/cutprimer.$r1Base.stdout.txt -eo $logDir/cutprimer.$r1Base.stderr.txt \
	"zcat $inputData/$r1Base.fastq.gz | cutadapt -e 0.1 -g ATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTC --info-file $cutPrimerDir/$r1Base.R1 - > $cutPrimerDir/$r1Base.cut.fastq" &

	bsub -K -q $bsubQueue  -J "cutP R2"    -M 8000 -R rusage[mem=8000] -oo $logDir/cutprimer.$r1Base.stdout.txt -eo $logDir/cutprimer.$r1Base.stderr.txt \
	"zcat $inputData/$r2Base.fastq.gz | cutadapt -e 0.1 -g GACTGCGTATCAGT --info-file $cutPrimerDir/$r1Base.R2 - > $cutPrimerDir/$r2Base.cut.fastq" &
    wait

}
runFastqc(){


    bsub -q $bsubQueue  -J "fqc R1"   -M 8000 -R rusage[mem=8000] -oo $logDir/fastqc.$r1Base.stdout.txt -eo $logDir/fastqc.$r1Base.stderr.txt \
    "fastqc -o $fastqcBeforeDir  $inputData/$1"

    bsub -q $bsubQueue  -J "fqc R1"   -M 8000 -R rusage[mem=8000] -oo $logDir/fastqc.$r2Base.stdout.txt -eo $logDir/fastqc.$r2Base.stderr.txt \
    "fastqc -o $fastqcBeforeDir  $inputData/$2"
}

 newFastqJoin1 (){
    info1=$1
    cutPrimerDir=$2
    tmpDir=$3
    myjoinBin=$4
    paste $tmpDir/temp.fastq $cutPrimerDir/$info1 | cut -f 1,2,3,4,5,6,7,8,11,12 \
    | awk -F "\t"  '{if($10-$9>30){print}}' | cut -f 1-4 | sed 's/^/@/'| sed 's/\t/\n/g'

}

 newFastqJoin2 (){
    info1=$1
    cutPrimerDir=$2
    tmpDir=$3
    myjoinBin=$4
      paste $tmpDir/temp.fastq $cutPrimerDir/$info1 | cut -f 1,2,3,4,5,6,7,8,11,12 \
      | awk -F "\t"  '{if($10-$9>30){print}}' \
      | cut -f 5-8 | sed 's/^/@/'| sed 's/\t/\n/g'
}

newFastqPaste(){
    set +o posix  #to get to redirect into paste working
    fq_R1=$1
    fq_R2=$2
    cutPrimerDir=$3
    tmpDir=$4
    paste <(cat $cutPrimerDir/$fq_R1 | paste - - - -| sed 's/^.//') <(cat $cutPrimerDir/$fq_R2 | paste - - - -|sed 's/^.//') > $tmpDir/temp.fastq

}


export -f newFastqJoin1
export -f newFastqPaste
export -f newFastqJoin2

runNewFastq(){
    fq_R1=$1.cut.fastq
    fq_R2=$2.cut.fastq
    info1=$1.R1
    info2=$2.R2
    tmpDir=$3

    newFastq1Out=$1.new.fastq
    newFastq2Out=$2.new.fastq
    bsub -K -q $bsubQueue  -J "nfpaste"   -M 8000 -R rusage[mem=8000] -oo $logDir/newfastqPaste.$r1Base.stdout.txt -eo $logDir/newfastqPaste.$r1Base.stderr.txt \
    newFastqPaste $fq_R1  $fq_R2 $cutPrimerDir  $tmpDir
    bsub -K -q $bsubQueue  -J "join 1  "   -M 8000  -R "rusage[mem=8000]" -oo $logDir/newFastqJoin1.$1.stdout.txt -eo $logDir/newFastqJoin1.$1.stderr.txt \
    "newFastqJoin1 $info1 $cutPrimerDir  $tmpDir $myjoinBin> $newFastqDir/$newFastq1Out" &
     bsub -K -q $bsubQueue  -J "join 2 "   -M 8000 -R "rusage[mem=8000]" -oo $logDir/newfastqJoin2.$1.stdout.txt -eo $logDir/newFastqJoin2.$1.stderr.txt \
     "newFastqJoin2 $info1 $cutPrimerDir  $tmpDir $myjoinBin> $newFastqDir/$newFastq2Out" &
    wait
}


runBWA(){
    fq_R1=$1.new.fastq
    fq_R2=$2.new.fastq
    bwa_out=$1.bam
    bwaIndex=$3
    newFastqDir=$4
    bwaDir=$5

    echo "BWA $bwaDir"

    echo " bwa mem -t 4  $bwaIndex $newFastqDir/$fq_R1 $newFastqDir/$fq_R2 | samtools view -hbS - > $bwaDir/$bwa_out"
    bwa mem -t 4  $bwaIndex $newFastqDir/$fq_R1 $newFastqDir/$fq_R2 | samtools view -hbS - > $bwaDir/$bwa_out

    ##this was added by me later on
   samtools sort $bwaDir/$bwa_out -o $bwaDir/$bwa_out.sorted
    samtools index $bwaDir/$bwa_out.sorted $bwaDir/$bwa_out.sorted.bai
   samtools flagstat $bwaDir/$bwa_out.sorted > $bwaDir/$bwa_out.stat
}

export -f runBWA



runBWABsub(){

    bsub -K -q $bsubQueue  -J "bwa"   -M 8000 -R rusage[mem=8000] -oo $logDir/bwa.$1.stdout.txt -eo $logDir/bwa.$1.stderr.txt \
    "runBWA $1 $2 $bwaIndex  $newFastqDir  $bwaDir"
}



runBamToBed(){
    set +o posix  #to get to redirect into paste working
    baseFile=$1
    bwaDir=$2
    bam2BedDir=$3
    bwa_out=$baseFile.bam
    scratch=$4

    samtools view -H $bwaDir/$bwa_out > $bam2BedDir/$baseFile.sam
    samtools view  -F 4  $bwaDir/$bwa_out | perl -ne \
   '$line=$_; @rec=split("\t",$line); $rec[5]=~ s/[SH].*//; if( !($rec[5] =~ /^[0-9]+$/ && $rec[5] > 6)){print $line;}' >> $bam2BedDir/$baseFile.sam
	samtools sort -T  $scratch/$baseFile.sort.tmp -n <(samtools view -hbS $bam2BedDir/$baseFile.sam ) -o $bam2BedDir/$baseFile.sorted.bam

    echo "basefile $baseFile"
    echo "bam2beddir $bam2BedDir"
	bedtools bamtobed -i $bam2BedDir/$baseFile.sorted.bam -bedpe -mate1 > $bam2BedDir/$baseFile.temp
    echo "bedtools bamtobed -i $bam2BedDir/$baseFile.sorted.bam -bedpe -mate1 > $bam2BedDir/$baseFile.temp"
	cat $bam2BedDir/$baseFile.temp | perl -e \
	'while(<STDIN>){$line=$_; @rec=split("\t", $line); if($rec[0] eq $rec[3]){ $start = ($rec[1], $rec[4])[$rec[1] > $rec[4]]; $end = ($rec[2], $rec[5])[$rec[2] < $rec[5]]; print "$rec[0]\t$start\t$end\t$rec[6]\t$rec[7]\t$rec[8]\n";}}' | awk -F "\t" '{if(($3-$2)<1000){print}}'\
	 > $bam2BedDir/$baseFile.bed

	cut -f 1,2,3,6 $bam2BedDir/$baseFile.bed | sort -u | awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$4}' > $bam2BedDir/$baseFile.rmdup.bed
	#rm $bam2BedDir/$baseFile.temp $bam2BedDir/$baseFile.sam

}
export -f runBamToBed
 runBamToBedBsub(){
     bsub -K -q $bsubQueue  -J "bamtobed"   -M 4000 -R rusage[mem=4000] -oo $logDir/bamtobed.$1.stdout.txt -eo $logDir/bamtobed.$1.stderr.txt \
     "runBamToBed $1 $bwaDir $bam2BedDir $scratch"
}



runDistance(){
set +o posix
  fileBase=$1
    bam2BedDir=$2
    bed2peakDir=$3
    finalDir=$4
    scratch=$5
    distance=$6
    knownSites=$7
    targetPredictionsR=$8
    module unload betools
    module load bedtools/2.25.0
    outDir1=$bed2peakDir/peaks_$distance
    outDir2=$finalDir/peaks_$distance

    mkdir -p $outDir1
    mkdir -p $outDir2
    echo "fileBase: $fileBase"
    echo "bam2BedDir $bam2BedDir"

    #bam2BedDir="bam2Bed"
    #fileBase=2137707_MZ-Sp5_S2_L001_R1_001


	bedtools --version

    awk -F "\t" '{if($6=="+"){print $1"\t"$2"\t"$6} else {print $1"\t"$3"\t"$6}}' $bam2BedDir/$fileBase.rmdup.bed \
    | grep -v vector | sort | uniq -c | awk '{print $2"\t"$3"\t"$3+1"\t.\t"$1"\t"$4}' | sort -k1,1 -k2,2n | \
    bedtools merge -c 5 -s -o sum -d $distance -i - | awk '{print $1"\t"$2"\t"$3"\t.\t"$5"\t"$4}' | grep -v vector | \
    bedtools intersect -b <(grep -v vector $bam2BedDir/$fileBase.bed| awk '{if($6=="+"){print $1"\t"$2"\t"$2+1"\t.\t.\t"$6} else {print $1"\t"$3"\t"$3+1"\t.\t.\t"$6}}') -c -s -a - | \
    awk '{print $1"\t"$2"\t"$3"\t""'$fileBase'""\t"$5"\t"$6"\t"$7}' > $outDir1/$fileBase.peak.merge.txt


    echo "bedtools window -w $distance -a $outDir1/$fileBase.peak.merge.txt -b $knownSites -sm -c > $outDir1/$fileBase.2.txt"
	bedtools window -w $distance -a $outDir1/$fileBase.peak.merge.txt -b $knownSites -sm -c > $outDir1/$fileBase.2.txt

    echo "TP: $targetPredictionsRDir"


R --slave <<EOF

source("$targetPredictionsR")
inputBed <- read.table("$outDir1/$fileBase.2.txt", sep="\t", check.names =FALSE )
colnames(inputBed) <- c("seqnames", "start", "end", "name", "nUniqueReads", "strand", "nReads", "nOverlapWithSpike")
inputBed[["name"]] <- gsub("-", ".", paste("X", inputBed[["name"]], 1:nrow(inputBed), sep="_"))
d <- target_gene_prediction(inputBed)
inputBed=inputBed[order(inputBed$name),]


rownames(d) <- d[["peak_id"]]
inputBed[["gene"]] <- d[inputBed[["name"]], "nearest_gene_symbol"]
inputBed[["tss_distances"]] <- d[inputBed[["name"]], "tss_distances"]
inputBed[["gene_region"]] <- d[inputBed[["name"]], "gene_region"]
inputBed[["name"]] <- sub("^X_", "", inputBed[["name"]])
inputBed <- inputBed[order(inputBed[["nUniqueReads"]], decreasing=TRUE),]
write.table(inputBed, file="$outDir1/$fileBase.peak.merge.xls", sep="\t", quote=FALSE, row.names=FALSE)
EOF




R --slave <<EOF
options(java.parameters = "-Xmx8000m")
library("xlsx")
options(stringsAsFactors = FALSE)
inputBed <- read.table("$outDir1/$fileBase.peak.merge.xls", sep="\t", header=TRUE, check.names =FALSE )
inputBed <- subset(inputBed, nUniqueReads > 1 | nReads > 5 )
inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)

write.xlsx(inputBed, file="$outDir1/$fileBase.peak.merge.xlsx")
#inputBed <- subset(inputBed, nOverlapWithSpike == 0)
inputBed[["percent"]] <- inputBed[,"nUniqueReads"] / sum(inputBed[,"nUniqueReads"])
inputBed[["percent"]] <- round(100*inputBed[["percent"]], 2)

inputBed[["percentAllReads"]] <- inputBed[,"nReads"] / sum(inputBed[,"nReads"])
inputBed[["percentAllReads"]] <- round(100*inputBed[["percentAllReads"]], 2)
write.xlsx(inputBed, file="$outDir2/$fileBase.peak.merge.xls")
EOF
#write.xlsx(head(inputBed, n=20), file="$outDir2/$fileBase.top20.peak.merge.xlsx")

}
export -f runDistance
runBed2Peak(){
    fileBase=$1
    bam2BedDir=$2
    bed2peakDir=$3
    finalDir=$4
    scratch=$5
    distances=$6
    knownSites=$7
    targetPredictionsR=$8

    IFS=','
    read -ra DISTANCES <<<"$distances"

    for i in "${DISTANCES[@]}"; #accessing each element of array
    do

        runDistance $fileBase $bam2BedDir $bed2peakDir $finalDir $scratch $i $knownSites $targetPredictionsR




    done





}
export -f runBed2Peak



runBed2PeakBsub(){
 bsub -K -q $bsubQueue  -J "bed2peak"   -M 4000 -R rusage[mem=4000] -oo $logDir/bed2peak.$1.stdout.txt -eo $logDir/bed2peak.$1.stderr.txt \
     "runBed2Peak $1 $bam2BedDir $bed2peakDir $finalDir $scratch $overlapDistances $knownSites $targetPredictionsR"

}


runbed2bwBsub(){
     bsub -K -q $bsubQueue  -J "bed2bw"   -M 4000 -R rusage[mem=4000] -oo $logDir/bed2peak.$1.stdout.txt -eo $logDir/bed2peak.$1.stderr.txt \
     "runBed2Peak $1 $bam2BedDir $bed2peakDir $finalDir $scratch $overlapDistances $knownSites $targetPredictionsR"
}


 runSample (){
    r1Base=$1
    r2Base=`echo $1| sed -e 's/_R1_/_R2_/g'`
    tmpDir=$scratch/$r1Base
    mkdir -p $tmpDir
  runFastqc   $r1Base $r2Base &
  runCutadapt $r1Base $r2Base
   runNewFastq $r1Base $r2Base $tmpDir
   runBWABsub $r1Base $r2Base
   runBamToBedBsub $r1Base
   runBed2PeakBsub $r1Base


}


for f in `ls -1 $inputDataDir | grep -P "fastq.gz$" | grep '_R1_' `
do
    sampleBase=`echo $f | sed -e 's/.fastq.gz//g'`
    sleep 5
     runSample $sampleBase &
done



wait

echo "FINISHED"
