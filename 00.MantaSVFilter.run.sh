for i in ../BAM_directory/*.sorted.bam;do
        name=${i##*/}
        name=${name%%.*}
	#echo $name
	python 00.filter.py $name
done
