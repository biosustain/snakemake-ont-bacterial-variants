import os
import pandas as pd


# --------------------------------------------------------------------------- #
# General helper functions                                                    #
# --------------------------------------------------------------------------- #
def read_sample_sheet():
    samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
    # Create temporary columns to fill NaN only for grouping purposes
    samples["_annotation"] = samples["annotation"].fillna("none")
    samples["_masked_regions"] = samples["masked_regions"].fillna("none")
    # assign groups for annotated comparison for each unique set of reference,
    # annotation and masked_regions combinations
    samples["group"] = (
        samples.groupby(["reference", "_annotation", "_masked_regions"]).ngroup() + 1
    )
    samples["group"] = samples["group"].astype(str)
    # Drop the temporary columns to restore original NaN values
    samples = samples.drop(columns=["_annotation", "_masked_regions"])
    return samples


def list_reference_genomes():
    paths = SAMPLEINFO["reference"].unique()
    genomes = {}
    for path in paths:
        ident = os.path.basename(path)
        ident = os.path.splitext(ident)[0]
        if not ident in genomes.keys():
            genomes[ident] = path
    return genomes


def generate_results(wildcards):
    return [
        os.path.join(outdir, f"variant_reports/{sample}/{sample}_overview.html")
        for sample in SAMPLES
    ] + [
        os.path.join(outdir, f"variant_reports/group_{group}/comparison.html")
        for group in SAMPLEINFO.group.unique()
    ]


# --------------------------------------------------------------------------- #
# Helper functions using wildcards                                            #
# --------------------------------------------------------------------------- #
def get_fastq(wildcards):
    return SAMPLEINFO.loc[wildcards.sample, "fastq"]


def get_original_reference(wildcards):
    return SAMPLEINFO.loc[wildcards.sample, "reference"]


def get_reference(wildcards):
    return os.path.join(
        outdir,
        "inputs",
        os.path.splitext(
            os.path.basename(SAMPLEINFO.loc[wildcards.sample, "reference"])
        )[0]
        + ".fa",
    )


def get_reference_fai(wildcards):
    return get_reference(wildcards) + ".fai"


def get_contig_ids_from_reference(wildcards):
    """Specific input for deepvariant"""
    contig_ids = []
    with open(get_reference(wildcards), "r") as f:
        for line in f:
            if line.startswith(">"):
                contig_id = line[1:].split()[
                    0
                ]  # Extract ID after '>' and split at whitespace
                contig_ids.append(contig_id)
    return '"' + ",".join(contig_ids) + '"'


def get_annotation(wildcards):
    if not pd.isna(SAMPLEINFO.loc[wildcards.sample, "annotation"]):
        if os.path.exists(SAMPLEINFO.loc[wildcards.sample, "annotation"]):
            return SAMPLEINFO.loc[wildcards.sample, "annotation"]
    else:
        return ""


def download_model_for_clair3(wildcards):
    path2model = {
        "r1041_e82_400bps_sup_v4.2.0": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v420.tar.gz",
        "r1041_e82_400bps_sup_v4.3.0": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v430.tar.gz",
        "r1041_e82_400bps_sup_v5.0.0": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v500.tar.gz",
    }
    model = path2model[config["basecalling_model"]]
    model_name = os.path.basename(model)
    model_name = model_name.replace(".tar.gz", "")
    return (model, model_name)


def get_medaka_model(wildcards):
    modeldetails = config["basecalling_model"].split("_")
    modeldetails.insert(-1, "variant")
    return "_".join(modeldetails)


def get_vcf(wildcards):
    assigntype = {
        "medaka": "SNV",
        "clair3": "SNV",
        "NanoCaller": "SNV",
        "cutesv": "SV",
        "sniffles2": "SV",
    }
    return os.path.join(outdir, assigntype[wildcards.tool], "{tool}/{sample}.vcf")


def get_bed_file_for_filtering(wildcards):
    return str(SAMPLEINFO.loc[wildcards.sample, "masked_regions"])


def get_shared_variants(wildcards):
    if config["remove_common_variants"] and len(SAMPLES) > 1:
        return [
            str(os.path.join(outdir, "SNV/clair3/common_variants.vcf")),
            # str(os.path.join(outdir, "SNV/medaka/common_variants.vcf")),
            # str(os.path.join(outdir, "SNV/NanoCaller/common_variants.vcf")),
            str(os.path.join(outdir, "SNV/DeepVariant/common_variants.vcf")),
            str(os.path.join(outdir, "SV/cutesv/common_variants.vcf")),
            str(os.path.join(outdir, "SV/sniffles2/common_variants.vcf")),
            str(os.path.join(outdir, "SNV/breseq/common_variants.vcf")),
        ]
    else:
        return ()


def get_frequency_filtering_parameters(wildcards):
    if wildcards.tool == "clair3":
        return f"-i 'FORMAT/AF >= {config[" frequency_threshold "][wildcards.tool]}'"
    elif wildcards.tool == "DeepVariant":
        return f"-i 'FORMAT/VAF >= {config[" frequency_threshold "][wildcards.tool]}'"
    elif wildcards.tool == "sniffles2":
        return f"-i 'INFO/AF >= {config[" frequency_threshold "][wildcards.tool]}'"
    elif wildcards.tool == "cutesv":
        return f"-i 'INFO/AF >= {config[" frequency_threshold "][wildcards.tool]}'"
    elif wildcards.tool == "breseq":
        return f"-i 'INFO/AF >= {config[" frequency_threshold "][wildcards.tool]}'"


def get_group_reference(wildcards):
    reference = SAMPLEINFO[SAMPLEINFO["group"] == wildcards.group].reference.unique()
    if len(reference) == 1:
        return reference[0]
    else:
        print(reference)
        raise ValueError(f"Multiple or no references found for group {wildcards.group}")


def get_group_annotation(wildcards):
    gff = SAMPLEINFO[SAMPLEINFO["group"] == wildcards.group].annotation.unique()
    if len(gff) == 1:
        return gff[0]
    else:
        raise ValueError(
            f"Multiple or no annotations found for group {wildcards.group}"
        )


def get_group_gd_files(wildcards):
    group_samples = SAMPLEINFO[SAMPLEINFO["group"] == wildcards.group][
        "sample"
    ].tolist()
    gd_files = [
        os.path.join(outdir, f"variant_reports/{sample}/{sample}.{tool}.gd")
        for sample in group_samples
        for tool in ["breseq", "clair3", "DeepVariant", "cutesv", "sniffles2"]
    ]
    return gd_files
