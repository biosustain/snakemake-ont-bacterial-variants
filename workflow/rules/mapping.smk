# --------------------------------------------------------------------------- #
# Read filtering                                                              #
# --------------------------------------------------------------------------- #
if config["filtering"] == "filtlong":
    rule filtlong:
        """Filtering of nanopore reads with filtlong"""
        message:
            "--- Read filtering with filtlong"
        input:
            get_fastq,
        output:
            os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
        log:
            os.path.join(outdir, "filtered_reads/logs/{sample}.log"),
        params:
            default=" ".join(config["filtlong"]["params"]),
        conda:
            "../envs/filtlong.yml"
        shell:
            "filtlong "
            "{params.default} "
            "{input} "
            "2> {log} | "
            "gzip > {output}"

elif config["filtering"] == "chopper":
    rule chopper:
        """Filtering of nanopore reads with chopper"""
        message:
            "--- Read filtering with chopper"
        input:
            reads=get_fastq,
        output:
            os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
        log:
            os.path.join(outdir, "filtered_reads/logs/{sample}.log"),
        params:
            default=lambda wildcards: " ".join(config["chopper"]["params"])
        threads:
            config["chopper"]["threads"]
        conda:
            "../envs/chopper.yml"
        shell:
            "chopper "
            "{params.default} "
            "--threads {threads} "
            "-i {input.reads} "
            "2> {log} | "
            "gzip > {output}"


# --------------------------------------------------------------------------- #
# Copy and index refernce files                                               #
# --------------------------------------------------------------------------- #
rule copy_index:
    """Copy and index refernce files"""
    message:
        "--- Copy and index refernce files"
    input:
        reference=lambda wildcards: GENOMES[wildcards.genome]
    output:
        reference=os.path.join(outdir, "inputs", "{genome}.fa"),
        index=os.path.join(outdir, "inputs", "{genome}.fa.fai"),
    log:
        stderr=os.path.join(outdir, "inputs", "logs", "{genome}.copy.stderr"),
        stdout=os.path.join(outdir, "inputs", "logs", "{genome}.copy.stdout"),
    conda:
        "../envs/minimap2.yml"
    shell:
        "cp {input.reference} {output.reference} && "
        "samtools faidx "
        "{output.reference} "
        "-o {output.index} "
        "1> {log.stdout} "
        "2> {log.stderr}"


# --------------------------------------------------------------------------- #
# Mapping against reference genome                                            #
# --------------------------------------------------------------------------- #
if config["aligner"] == "ngmlr":
    rule mapping:
        """Mapping of ONT data against reference genome with ngmlr"""
        message:
            "--- Mapping against reference genome with ngmlr"
        input:
            reads=os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
            reference=get_reference,
        output:
            sam=temp(os.path.join(outdir, "mapping/{sample}.sam")),
            bam=os.path.join(outdir, "mapping/{sample}.bam"),
            bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
            stats=temp(os.path.join(outdir, "mapping/{sample}.stats")),
        log:
            ngmlr=os.path.join(outdir, "mapping/logs/{sample}.ngmlr.log"),
            samtools=os.path.join(outdir, "mapping/logs/{sample}.samtools.log"),
        params:
            tmpdir=os.path.join(outdir, "mapping/{sample}"),
            other=" ".join(config["nglmr"]["params"])
        threads: config["ngmlr"]["threads"]
        conda:
            "../envs/ngmlr.yml"
        shell:
            "ngmlr "
            "-t {threads} "
            "-x ont "
            "{params.other} "
            "-r {input.reference} "
            "-q {input.reads} "
            "-o {output.sam} "
            "2> {log.ngmlr} && "
            "samtools sort "
            "{output.sam} "
            "-@ {threads} "
            "-T {params.tmpdir} "
            "-O bam "
            "> {output.bam} "
            "2> {log.samtools} && "
            "samtools index {output.bam} 2>> {log.samtools} && "
            "samtools stats {output.bam} > {output.stats} 2>> {log.samtools}"

elif config["aligner"] == "minimap2":
    rule mapping:
        """Mapping of ONT data against reference genome with minimap2"""
        message:
            "--- Mapping against reference genome with minimap2"
        input:
            reads=os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
            reference=get_reference,
        output:
            sam=temp(os.path.join(outdir, "mapping/{sample}.sam")),
            bam=os.path.join(outdir, "mapping/{sample}.bam"),
            bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
            stats=temp(os.path.join(outdir, "mapping/{sample}.stats")),
        log:
            minimap2=os.path.join(outdir, "mapping/logs/{sample}.minimap2.log"),
            samtools=os.path.join(outdir, "mapping/logs/{sample}.samtools.log"),
        params:
            tmpdir=os.path.join(outdir, "mapping/{sample}"),
            other=" ".join(config["minimap2"]["params"])
        threads: config["minimap2"]["threads"]
        conda:
            "../envs/minimap2.yml"
        shell:
            "minimap2 "
            "-a "
            "-t {threads} "
            "-x map-ont "
            "{params.other} "
            "-o {output.sam} "
            "{input.reference} "
            "{input.reads} "
            "2> {log.minimap2} && "
            "samtools sort "
            "{output.sam} "
            "-@ {threads} "
            "-T {params.tmpdir} "
            "-O bam "
            "> {output.bam} "
            "2> {log.samtools} && "
            "samtools index {output.bam} 2>> {log.samtools} && "
            "samtools stats {output.bam} > {output.stats} 2>> {log.samtools}"


# --------------------------------------------------------------------------- #
# Genome coverage                                                             #
# --------------------------------------------------------------------------- #
rule genomecoverage:
    """Generation of nucleotide-precise genome coverage depths"""
    message:
        "--- Genome coverage"
    input:
        os.path.join(outdir, "mapping/{sample}.bam"),
    output:
        temp(os.path.join(outdir, "mapping/{sample}.coverage")),
    log:
        os.path.join(outdir, "mapping/logs/{sample}.genomecoverage.log"),
    conda:
        # "../envs/ngmlr.yml"
        "../envs/minimap2.yml"
    shell:
        "bedtools genomecov "
        "-d "
        "-ibam {input} "
        "> {output} "
        "2> {log}"


# --------------------------------------------------------------------------- #
# Alignment ends                                                              #
# --------------------------------------------------------------------------- #
rule alignmentends:
    """Generation of 5' and 3' alignment end distributions by strand"""
    message:
        "--- Alignment end distribution"
    input:
        os.path.join(outdir, "mapping/{sample}.bam"),
    output:
        five_plus=temp(os.path.join(outdir, "mapping/{sample}.5ends.plus.counts")),
        five_minus=temp(os.path.join(outdir, "mapping/{sample}.5ends.minus.counts")),
        three_plus=temp(os.path.join(outdir, "mapping/{sample}.3ends.plus.counts")),
        three_minus=temp(os.path.join(outdir, "mapping/{sample}.3ends.minus.counts")),
    log:
        os.path.join(outdir, "mapping/logs/{sample}.alignmentends.log"),
    conda:
        # "../envs/ngmlr.yml"
        "../envs/minimap2.yml"
    shell:
        "bedtools genomecov -d -ibam {input} "
        "-5 -strand + > {output.five_plus} 2> {log} && "
        "bedtools genomecov -d -ibam {input} "
        "-5 -strand - > {output.five_minus} 2>> {log} && "
        "bedtools genomecov -d -ibam {input} "
        "-3 -strand + > {output.three_plus} 2>> {log} && "
        "bedtools genomecov -d -ibam {input} "
        "-3 -strand - > {output.three_minus} 2>> {log}"
