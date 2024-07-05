if [[ $1 == "prep" ]]; then
	/home/ubuntu/tools/nextflow \
		main.nf \
		-profile docker,test \
		-resume \
		-entry OhioTestPrep \
		-work-dir ~/output/tbtest/work \
		--outdir ~/output/tbtest
fi

if [[ $1 == "run" ]]; then
	/home/ubuntu/tools/nextflow \
		main.nf \
		-profile docker,test \
		-resume \
		-entry OhioTBAnalyzer \
		-work-dir ~/output/tbtest/work \
		--outdir ~/output/tbtest
fi
