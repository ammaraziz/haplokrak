import snakemake
from pathlib import Path
import pandas
from Bio import SeqIO

IFQ = "../data/"
output = "analysis/"

SAMPLE_NAME, SAMPLE_NUMBER, PAIR = glob_wildcards(IFQ + "/{sample_name}_{sample_number}_{pair}_001.fastq.gz")
SAMPLES = [i + "_" + x for i, x in zip(SAMPLE_NAME, SAMPLE_NUMBER)]
print(SAMPLES)

rule all:
    input:
        expand(output + "{sample}/{sample}_kraken.status", sample = SAMPLES),
        expand(output + "{sample}/{sample}_extract.status", sample = SAMPLES),
        
        expand(output + "{sample}/{sample}/{sample}.haploflow.status", sample = SAMPLES),
        expand(output + "aggregated/{sample}/{sample}.haploflow.fasta", sample = SAMPLES),
        expand(output + "abricate/{sample}.abricate.results", sample = SAMPLES),
        expand(output + "filtered/{sample}.filtered.fasta", sample = SAMPLES),
        expand(output + "aligned/{sample}.bam", sample = SAMPLES)

rule kraken:
    threads: 2
    input:
        r1 = expand(IFQ + "{{sample}}_{pair}_001.fastq.gz", pair = ["R1"]),
        r2 = expand(IFQ + "{{sample}}_{pair}_001.fastq.gz", pair = ["R2"])
    output:
        status = output + "{sample}/{sample}_kraken.status",
        report = output + "{sample}/{sample}.report",
        kraken_file = output + "{sample}/{sample}.kraken"
    params:
        krakdb = "db/minikraken2/",
        outdir = "analysis/{sample}/",
        classified_reads = output + "{sample}/{sample}_#_classified.fastq",
        unclassified_reads = output + "{sample}/{sample}_#_unclassified.fastq"
    threads: 20
    shell:"""
	touch {output.status}
        kraken2 \
        --db {params.krakdb} \
        --threads {threads} \
        --report {output.report} \
        --unclassified-out {params.unclassified_reads} \
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
        classified_r1 = output + "{sample}/{sample}__1_classified.fastq",
        classified_r2 = output + "{sample}/{sample}__2_classified.fastq",
        kraken_output = output + "{sample}/{sample}.kraken",
        kraken_report = output + "{sample}/{sample}.report"   
    output:
        entero_r1 = output + "{sample}/{sample}_entero_1.fastq",
        entero_r2 = output + "{sample}/{sample}_entero_2.fastq",
        status = output + "{sample}/{sample}_extract.status"
    params:
        taxid = 464095
    shell:"""
    touch {output.status}
    python extract_kraken_reads.py \
    -1 {input.classified_r1} \
    -2 {input.classified_r2} \
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
             r1 = rules.extract.output.entero_1,
             r2 = rules.extract.output.entero_2
	output:
		status = output + "{sample}/{sample}.status"
	params:
		output_dir = output + "{sample}"
	shell:"""
		haploflow --read-file {input.r1} {input.r2} --out {params.output_dir}
		touch {output.status}
		"""

def aggregate_input(wildcards):
  checkpoint_output = checkpoints.haploflow.get(**wildcards).output[0]
  
  # split paths
  path = Path(checkpoint_output)
  name = path.name.split(".")[0]
  folder = path.parent

  # return 
  return expand(folder / "contigs.fa",
              sample=glob_wildcards(os.path.join(folder, "contigs.fa")))

rule aggregate:
    input:
        aggregate_input
    output:
        output + "aggregated/{sample}.fasta"
    shell:"""
        cp {input} {output}
        """

rule abricate:
	input:
		rules.aggregate.output
	output:
		output + "abricate/{sample}.results"
	threads: 10
	shell:"""
	abricate {input} --db entero_wgfull --minid 70 --threads {threads} > {output} 
	"""

rule filter_fasta:
	input:
		abricate = rules.abricate.output,
		assembly = rules.aggregate.output
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
