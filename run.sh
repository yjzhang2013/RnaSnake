snakemake --cluster-config cluster.json --cluster "qsub -V -d ./ -S /bin/bash -q batch" --jobs 200
multiqc -o report Preprocess Mapping Quantification
