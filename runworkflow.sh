if [[ $1 == "prep" ]]; then
	/home/ubuntu/tools/nextflow \
		main.nf \
		-profile docker,test \
		-entry OhioTestPrep
fi

if [[ $1 == "run" ]]; then
	/home/ubuntu/tools/nextflow \
		main.nf \
		-profile docker,test \
		-entry OhioTBAnalyzer \
		-work-dir ~/output/tbtest
fi
