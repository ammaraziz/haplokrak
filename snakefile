import snakemake
from pathlib import Path
import pandas
from Bio import SeqIO

input_directory = "data/"
output = "analysis/"

SAMPLE_NAME, SAMPLE_NUMBER, PAIR = glob_wildcards(input_directory + "/{sample_name}_{sample_number}_{pair}_001.fastq.gz")
SAMPLES = [i + "_" + x for i, x in zip(SAMPLE_NAME, SAMPLE_NUMBER)]
print(SAMPLE_NAME)

rule all:
    input:
    	# kraken
        expand(output + "kraken/{sample}_kraken.status", sample = SAMPLES),
        expand(output + "kraken/{sample}_extract.status", sample = SAMPLES),
        # haploflow
        expand(output + "haploflow/{sample}.haploflow.status", sample = SAMPLES),
        expand(output + "aggregated/{sample}.fasta", sample = SAMPLES),
        # abricate
        expand(output + "abricate/{sample}.abricate.results", sample = SAMPLES),
        expand(output + "filtered/{sample}.filtered.fasta", sample = SAMPLES),
        # bowtie2
        expand(output + "bowtie/{sample}.bam", sample = SAMPLES)

rule kraken:
    input:
        r1 = expand(input_directory + "{{sample}}_{pair}_001.fastq.gz", pair = ["R1"]),
        r2 = expand(input_directory + "{{sample}}_{pair}_001.fastq.gz", pair = ["R2"])
    output:
        status = output + "kraken/{sample}_kraken.status",
        kraken_report = output + "kraken/{sample}.report",
        kraken_file = output + "kraken/{sample}.kraken"
    params:
        krakdb = "db/minikraken2/",
        outdir = "analysis/{sample}/",
        classified_reads = output + "kraken/{sample}_#_classified.fastq"
    threads: 20
    shell:"""
	touch {output.status}
    kraken2 \
    --db {params.krakdb} \
    --threads {threads} \
    --report {output.kraken_report} \
    --classified-out {params.classified_reads} \
    --output {output.kraken_file} \
    --use-names \
    --paired \
    --gzip-compressed \
    {input.r1} \
    {input.r2}
    """

rule extract:
    input:
        kraken_output = output + "kraken/{sample}.kraken",
        kraken_report = output + "kraken/{sample}.report"   
    output:
        entero_r1 = output + "kraken/{sample}_entero_1.fastq",
        entero_r2 = output + "kraken/{sample}_entero_2.fastq",
        status = output + "kraken/{sample}_extract.status"
    params:
        classified_r1 = output + "kraken/{sample}__1_classified.fastq",
        classified_r2 = output + "kraken/{sample}__2_classified.fastq",
        taxid = 464095 # picornavirus taxon id
    shell:"""
    touch {output.status}
    python extract_kraken_reads.py \
    -1 {params.classified_r1} \
    -2 {params.classified_r2} \
    -k {input.kraken_output} \
    -r {input.kraken_report} \
    -t {params.taxid} \
    -o {output.entero_r1} \
    -o2 {output.entero_r2} \
    --include-children
    """

checkpoint haploflow:
	threads: 10
	input:
             r1 = rules.extract.output.entero_r1,
             r2 = rules.extract.output.entero_r2
	output:
		status = output + "haploflow/{sample}.haploflow.status"
	params:
		output_dir = output + "haploflow/{sample}"
	shell:"""
		touch {output.status}
		haploflow --read-file \
		{input.r1} \
		{input.r2} \
		--out {params.output_dir}
		"""

def get_haplo_data(wildcards):
  checkpoint_output = checkpoints.haploflow.get(**wildcards).output[0]
  
  # split paths
  path = Path(checkpoint_output)
  name = path.name.split(".")[0]
  folder = path.parent

  # return 
  return expand(folder / "contigs.fa",
              sample=glob_wildcards(os.path.join(folder, "contigs.fa")))

rule aggregate_haploflow:
    input:
        get_haplo_data
    output:
        output + "aggregated/{sample}.fasta"
    shell:"""
        cp {input} {output}
        """

rule abricate:
	input:
		rules.aggregate_haploflow.output
	output:
		output + "abricate/{sample}.results"
	params:
		minid = 10
	threads: 5
	shell:"""
	abricate {input} \
	--db entero_wgfull \
	--minid {params.minid} \
	--threads {threads} > {output} 
	"""

# rule needs reworking
# must combine the same taxid together
rule filter_fasta:
	input:
		abricate = rules.abricate.output,
		assembly = rules.aggregate_haploflow.output
	output:
		output + "filtered/{sample}.filtered.fasta"
	run:
		df = pandas.read_csv(input[0], sep="\t")
		wanted = {df['SEQUENCE'][0] : df['RESISTANCE'][0]}
		records = (r for r in SeqIO.parse(input[1], "fasta") if r.id in wanted.keys())	

		with open(output[0],'a') as filtered:
			for record in records:
				record.id = f"{record.id}|{wanted[record.id]}"		
				record.description = ""
				SeqIO.write(record, filtered, "fasta")

rule build_index:
	input:
		rules.filter_fasta.output
	output:
		output + "filtered/{sample}.status.index"
	params:
		basename = output + "filtered/{sample}.filtered.fasta.index"
	shell:"""
	bowtie2-build -q {input} {params.basename}
	touch {output}
	"""

rule align_reads:
	input:
		assembly = rules.filter_fasta.output,
		r1 = rules.haploflow.input.r1,
		r2 = rules.haploflow.input.r2
	output:
		bam = output + "aligned/{sample}.bam",
		stats = output + "aligned/{sample}.statistics"
	threads: 20
	shell:"""
	bowtie2 -q \
	--no-unal \
	-p {threads} \
	-x {input.assembly}.index \
	-1 {input.r1} -2 {input.r2} \
	2> {output.stats} \
	| samtools view -bS - \
	| samtools sort -@ {threads} - 1> {output.bam} 2> /dev/null
	samtools index {output.bam}
	"""
