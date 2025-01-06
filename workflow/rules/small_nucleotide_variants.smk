# --------------------------------------------------------------------------- #
# Identification of small nucleotide variants using Medaka                    #
# --------------------------------------------------------------------------- #
rule snv_medaka:
    """Identification of small nucleotide variants using Medaka"""
    message:
        "--- SNV calling using Medaka"
    input:
        reads=os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
        reference=get_reference,
    output:
        medakadir=directory(os.path.join(outdir, "SNV/medaka/{sample}")),
        vcf=os.path.join(outdir, "SNV/medaka/{sample}.vcf"),
    log:
        stdout=os.path.join(outdir, "SNV/medaka/logs/{sample}.stdout"),
        stderr=os.path.join(outdir, "SNV/medaka/logs/{sample}.stderr"),
    params:
        model=get_medaka_model,
    threads: config["medaka"]["threads"]
    conda:
        "../envs/medaka.yml"
    shell:
        "medaka_variant "
        "-m {params.model} "
        "-t {threads} "
        "-o {output.medakadir} "
        "-i {input.reads} "
        "-r {input.reference} "
        "> {log.stdout} "
        "2> {log.stderr} && "
        "cp {output.medakadir}/medaka.annotated.vcf {output.vcf}"


# --------------------------------------------------------------------------- #
# Download models for Clair3                                                  #
# --------------------------------------------------------------------------- #
rule models_clair3:
    """Downloading basecalling models for Clair3"""
    message:
        "--- Downloading models for Clair3"
    output:
        directory(os.path.join(outdir, "SNV/clair3/model")),
    log:
        stdout=os.path.join(outdir, "SNV/clair3/logs/model.stdout"),
        stderr=os.path.join(outdir, "SNV/clair3/logs/model.stderr"),
    params:
        download_model_for_clair3,
    conda:
        "../envs/basic.yml"
    shell:
        "wget "
        "--directory-prefix {output} "
        "{params[0][0]} "
        "> {log.stdout} "
        "2> {log.stderr} && "
        "tar -xf {output}/{params[0][1]}.tar.gz -C {output} -v >> {log.stdout} && "
        "rm -rf {output}/{params[0][1]}.tar.gz && "
        "mv {output}/{params[0][1]}/* {output} && "
        "rm -rf {output}/{params[0][1]}"


# --------------------------------------------------------------------------- #
# Identification of small nucleotide variants using Clair3                    #
# --------------------------------------------------------------------------- #
rule snv_clair3:
    """Identification of small nucleotide variants using Clair3"""
    message:
        "--- SNV calling using Clair3"
    input:
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
        bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
        reference=get_reference,
        model=rules.models_clair3.output,
    output:
        clair3dir=directory(os.path.join(outdir, "SNV/clair3/{sample}")),
        vcf=os.path.join(outdir, "SNV/clair3/{sample}.vcf"),
    log:
        stdout=os.path.join(outdir, "SNV/clair3/logs/{sample}.stdout"),
        stderr=os.path.join(outdir, "SNV/clair3/logs/{sample}.stderr"),
    params:
        " ".join(config["clair3"]["params"]),
    threads: config["clair3"]["threads"]
    conda:
        "../envs/clair3.yml"
    shell:
        "run_clair3.sh "
        "--threads {threads} "
        "{params} "
        "--bam_fn {input.bam} "
        "--ref_fn {input.reference} "
        "--model_path {input.model} "
        "--output {output.clair3dir} "
        "> {log.stdout} "
        "2> {log.stderr} && "
        "cp {output.clair3dir}/merge_output.vcf.gz {output.vcf}.gz && "
        "gunzip {output.vcf}.gz"


# --------------------------------------------------------------------------- #
# Identification of small nucleotide variants using NanoCaller                #
# --------------------------------------------------------------------------- #
rule snv_NanoCaller:
    """Identification of small nucleotide variants using NanoCaller"""
    message:
        "--- SNV calling using NanoCaller"
    input:
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
        bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
        reference=get_reference,
        # model=rules.models_clair3.output,
    output:
        NanoCallerdir=directory(os.path.join(outdir, "SNV/NanoCaller/{sample}")),
        vcf=os.path.join(outdir, "SNV/NanoCaller/{sample}.vcf"),
    log:
        stdout=os.path.join(outdir, "SNV/NanoCaller/logs/{sample}.stdout"),
        stderr=os.path.join(outdir, "SNV/NanoCaller/logs/{sample}.stderr"),
    params:
        " ".join(config["NanoCaller"]["params"]),
    threads: config["NanoCaller"]["threads"]
    conda:
        "../envs/NanoCaller.yml"
    shell:
        "NanoCaller "
        "--cpu {threads} "
        "{params} "
        "--bam {input.bam} "
        "--ref {input.reference} "
        "--output {output.NanoCallerdir} "
        "--prefix {wildcards.sample} "
        "--sample {wildcards.sample} "
        "> {log.stdout} "
        "2> {log.stderr} && "
        "cp {output.NanoCallerdir}/{wildcards.sample}.vcf.gz {output.vcf}.gz && "
        "gunzip {output.vcf}.gz"


# --------------------------------------------------------------------------- #
# Identification of small nucleotide variants using DeepVariant (Docker)      #
# --------------------------------------------------------------------------- #
rule snv_DeepVariant:
    """Identification of small nucleotide variants using DeepVariant (Docker)"""
    message:
        "--- SNV calling using DeepVariant (Docker)"
    input:
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
        bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
        reference=get_reference,
        reference_fai=get_reference_fai,
    output:
        vcf=os.path.join(outdir, "SNV/DeepVariant/{sample}.vcf"),
        gvcf=os.path.join(outdir, "SNV/DeepVariant/{sample}.g.vcf"),
    log:
        stdout=os.path.join(outdir, "SNV/DeepVariant/logs/{sample}.stdout"),
        stderr=os.path.join(outdir, "SNV/DeepVariant/logs/{sample}.stderr"),
    params:
        contigs=get_contig_ids_from_reference,
        config=" ".join(config["DeepVariant"]["params"]),
    container:
        "docker://google/deepvariant:1.6.1"
    threads: config["DeepVariant"]["threads"]
    shell:
        "/opt/deepvariant/bin/run_deepvariant "
        "--ref={input.reference} "
        "--output_vcf={output.vcf} "
        "--output_gvcf={output.gvcf} "
        "--num_shards {threads} "
        "--haploid_contigs={params.contigs} "
        "{params.config} "
        "--reads={input.bam} "
        "> {log.stdout} "
        "2> {log.stderr}"


# --------------------------------------------------------------------------- #
# Identification of variants using Breseq                                     #
# --------------------------------------------------------------------------- #
# At least for now, it doesn't look like using previous alignments is an option:
# https://github.com/barricklab/breseq/issues/378
rule breseq:
    """Identification of variants using Breseq"""
    message:
        "--- Identification of variants using Breseq"
    input:
        reads=os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
        reference=get_reference,
        gff=get_annotation,
    output:
        breseqdir=directory(os.path.join(outdir, "SNV/breseq/{sample}")),
        annotated=os.path.join(outdir, "SNV/breseq/{sample}/data/annotated.gd"),
        out_vcf=os.path.join(outdir, "SNV/breseq/{sample}/output/output.vcf"),
        final_vcf=os.path.join(outdir, "SNV/breseq/{sample}.vcf"),
    log:
        stdout=os.path.join(outdir, "SNV/breseq/logs/{sample}.stdout"),
        stderr=os.path.join(outdir, "SNV/breseq/logs/{sample}.stderr"),
    params:
        " ".join(config["breseq"]["params"]),
    threads: config["breseq"]["threads"]
    conda:
        "../envs/breseq.yml"
    shell:
        "breseq "
        "-r {input.reference} "
        "-r {input.gff} "
        "-j {threads} "
        "--nanopore "
        "{params} "
        "-o {output.breseqdir} "
        "{input.reads} "
        "> {log.stdout} "
        "2> {log.stderr} && "

        "cp {output.out_vcf} {output.final_vcf}"
        # "cp {output.annotated} {output.gd_file} && "
