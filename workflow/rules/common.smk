__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022, Martin Rippin"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas
import yaml

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from snakemake.utils import min_version
from snakemake.utils import validate

min_version("6.10")


### Set and validate config file
configfile: "config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


### Load and validate resources file
config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file
samples = pandas.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


### Read and validate units file
units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    sample="|".join(get_samples(samples)),
    unit="N|T|R",
    read="fastq[1|2]",


### Functions ###
def get_in_fastq(units, wildcards):
    return expand(
        "prealignment/merged/{{sample}}_{{type}}_{read}.fastq.gz",
        read=["fastq1", "fastq2"],
    )


def get_in_fq(wildcards):
    input_list = []
    for unit in get_units(units, wildcards, wildcards.type):
        prefix = "prealignment/merged/{}_{}".format(unit.sample, unit.type)
        input_unit = "{}_fastq1.fastq.gz {}_fastq2.fastq.gz {}".format(
            prefix,
            prefix,
            "'@RG\\tID:{}\\tSM:{}\\tPL:{}\\tPU:{}\\tLB:{}'".format(
                "{}_{}.{}.{}".format(unit.sample, unit.type, unit.lane, unit.barcode),
                "{}_{}".format(unit.sample, unit.type),
                unit.platform,
                "{}.{}.{}".format(unit.flowcell, unit.lane, unit.barcode),
                "{}_{}".format(unit.sample, unit.type),
            ),
        )
        input_list.append(input_unit)
    return " --in-fq ".join(input_list)


def compile_output_list(wildcards: snakemake.io.Wildcards):
    output_list = ["qc/multiqc/MultiQC.html"]
    output_list.append(
        ["alignment/bwa_mem/%s_%s.bam" % (sample, type)
        for sample in get_samples(samples)
        for type in get_unit_types(units, sample)]
        )
    return output_list


#    output_list.append(["cnv_sv/cnvkit_vcf/%s_T.vcf" % (sample) for sample in get_samples(samples)])
#    output_list.append(["cnv_sv/cnvkit_diagram/%s_T.png" % (sample) for sample in get_samples(samples)])
#    output_list.append(["cnv_sv/manta/%s.ssa.%s.vcf" % (sample, diagnosis) for sample in get_samples(samples) for diagnosis in ["aml", "all"]])
#    output_list.append(
#        ["parabricks/mutectcaller/%s.vep.%s.vcf" % (sample, diagnosis)
#            for sample in get_samples(samples)
#            for diagnosis in ["aml", "all"]])
