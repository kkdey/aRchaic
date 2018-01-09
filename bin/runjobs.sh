# run the pipeline !!!
snakemake --snakefile Snakefile  \
	  --cluster-config cluster.json \
	  -c "sbatch -t {cluster.time} --mem {cluster.mem} -o {cluster.out} -e {cluster.err} --partition=sandyb --account=pi-mstephens" \
	  -j 30 \
	  -k \
	  -p \
	    $*
