FLAG=$1
INPUT=$2
OUTPUT=$3

if [[ $FLAG == "test" ]]; then
	/home/ubuntu/tools/nextflow \
		main.nf \
		-profile docker,test \
		-resume \
		-entry OhioTestPrep \
		-work-dir ~/output/tbtest/work \
		--outdir ~/output/tbtest
fi

if [[ $FLAG == "run" ]]; then
	/home/ubuntu/tools/nextflow \
		main.nf \
		-profile docker \
		-resume \
		-entry OhioTBAnalyzer \
		-work-dir $OUTPUT/work \
		--outdir $OUTPUT \
		--input $INPUT
fi
