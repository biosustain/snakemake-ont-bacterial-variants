# --------------------------------------------------------------------------- #
# Convert VCF to GD and annotate                                              #
# --------------------------------------------------------------------------- #
rule gdtools_annotate:
    """Convert vcf to gd format and annotate"""
    message:
        "--- Convert vcf to gd format and annotate"
    input:
        vcf_file=os.path.join(outdir, "variant_reports/{sample}/{sample}.{tool}.vcf"),
        reference=get_reference,
        gff=get_annotation,
    output:
        gd_file=temp(os.path.join(outdir, "variant_reports/{sample}/{sample}.{tool}.gd")),
        annotated=os.path.join(outdir, "variant_reports/{sample}/{sample}.{tool}.annotated.gd"),
    log:
        stdout=os.path.join(outdir, "variant_reports/logs/{sample}.{tool}.annotate.stdout"),
        stderr=os.path.join(outdir, "variant_reports/logs/{sample}.{tool}.annotate.stderr"),
    conda:
        "../envs/breseq.yml"
    shell:
        "gdtools VCF2GD "
        "{input.vcf_file} "
        "-o {output.gd_file} "
        "1> {log.stdout} "
        "2> {log.stderr} && "
        "gdtools ANNOTATE "
        "-r {input.gff} "
        "-r {input.reference} "
        "-f GD "
        "-o {output.annotated} "
        "{output.gd_file} "
        "1> {log.stdout} "
        "2> {log.stderr}"


# --------------------------------------------------------------------------- #
# Compare with gdtools                                                       #
# --------------------------------------------------------------------------- #
rule gdtools_compare:
    """Compare variants with gdtools"""
    message:
        "--- Compare variants with gdtools"
    input:
        # gd_files=expand(os.path.join(outdir, "variant_reports/{sample}/{sample}.{tool}.gd"), sample=SAMPLES, tool=["breseq", "clair3", "medaka", "NanoCaller", "cutesv", "sniffles2"]),
        # gd_files=expand(os.path.join(outdir, "variant_reports/{sample}/{sample}.{tool}.gd"), sample=SAMPLES, tool=["breseq", "clair3", "NanoCaller", "DeepVariant",  "cutesv", "sniffles2"]),
        # gd_files=expand(os.path.join(outdir, "variant_reports/{sample}/{sample}.{tool}.gd"), sample=SAMPLES, tool=["breseq", "clair3", "DeepVariant",  "cutesv", "sniffles2"]),
        gd_files=get_group_gd_files,
        reference=get_group_reference,
        gff=get_group_annotation,
    output:
        html=os.path.join(outdir, "variant_reports/group_{group}/comparison.html"),
        tsv=os.path.join(outdir, "variant_reports/group_{group}/comparison.tsv"),
    log:
        stdout=os.path.join(outdir, "variant_reports/logs/group_{group}/compare.stdout"),
        stderr=os.path.join(outdir, "variant_reports/logs/group_{group}/compare.stderr"),
    conda:
        "../envs/breseq.yml"
    shell:
        "gdtools COMPARE "
        "-r {input.gff} "
        "-r {input.reference} "
        "-f HTML "
        "-o {output.html} "
        "{input.gd_files} "
        "1> {log.stdout} "
        "2> {log.stderr} && "
        "gdtools COMPARE "
        "-r {input.gff} "
        "-r {input.reference} "
        "-f table "
        "-o {output.csv} "
        "{input.gd_files} "
        "1> {log.stdout} "
        "2> {log.stderr}"