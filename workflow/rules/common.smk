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

min_version("7.8.0")


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
units = pandas.read_table(config["units"], dtype=str).set_index(["sample", "type", "flowcell", "barcode", "lane"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

#Msamples = pandas.read_table("M.txt", dtype=str).set_index("samples", drop=False)
#Fsamples = pandas.read_table("F.txt", dtype=str).set_index("samples", drop=False)


### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    type="N|T|R",
    read="fastq[1|2]",
    sex="M|F",


### Functions

def compile_output_list(wildcards):
    output_list = ["qc/multiqc/multiqc_DNA.html"]
    output_list.append([
        "vcf_final/%s.vcf.gz" % sample
        for sample in get_samples(samples)
    ])
    output_list.append("conifer/calls/calls_F_svd6.txt")
    output_list.append("conifer/calls/calls_M_svd6.txt")
    output_list.append([
        "cnv_sv/exomedepth_F/%s_N.SV.filter.txt" % sample
        for sample in get_samples(samples)
    ])
    output_list.append([
        "cnv_sv/exomedepth_M/%s_N.SV.filter.txt" % sample
        for sample in get_samples(samples)
    ])
    output_list.append([
        "prealignment/merged/{}_{}_{}.fastq.gz".format(sample, t, read)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for read in ["fastq1", "fastq2"]
    ])
    output_list.append([
        "compression/spring/%s_%s_%s_%s_%s.spring" % (sample, flowcell, lane, barcode, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for flowcell in set([u.flowcell for u in units.loc[(sample,t,)].dropna().itertuples()])
        for barcode in set([u.barcode for u in units.loc[(sample,t,)].dropna().itertuples()])
        for lane in set([u.lane for u in units.loc[(sample,t,)].dropna().itertuples()])
    ])
    output_list.append([
            "compression/crumble/%s_N.crumble.cram" % sample
            for sample in get_samples(samples)
    ])
    return output_list
