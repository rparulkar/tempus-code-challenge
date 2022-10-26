import argparse
import asyncio
import json
import math
import requests
import copy

from timeit import default_timer
from concurrent.futures import ThreadPoolExecutor

START_TIME = default_timer()

out_header = ["Chromosome", "Position", "Reference", "Alternate", "Total Depth",
              "Alternate Depth", "Alternate AF", "MAF", "Gene Name", "ENSEMBL Gene ID",
              "ENSEMBL Transcript ID", "Canonical", "Most Severe Consequence", "Variant Type"]


# Class to represent variant records
class VarRecord:
    def __init__(self, vcf_row_dict, sample="sample"):
        self.chrom = vcf_row_dict["CHROM"]
        self.pos = vcf_row_dict["POS"]
        self.ref = vcf_row_dict["REF"]
        self.alt = vcf_row_dict["ALT"]
        self.filter = vcf_row_dict["FILTER"]
        self.format = vcf_row_dict["FORMAT"]
        self.sample = dict(zip(self.format.split(":"), vcf_row_dict[sample].split(":")))
        self.hgvs = self.construct_hgvs()
        self.annotation = None

    # Gets basic coverage information from variant record
    def read_coverage(self):
        depth = float(self.sample["NR"])
        alt_depth = float(self.sample["NV"])
        return str(depth), str(alt_depth), str(alt_depth / depth)

    # Creates an object for each allele, for multi allelic sites
    def parse_multi_allelic_sites(self):
        alleles = [self.ref] + self.alt.split(',')
        multi_allele_dict = dict(zip(alleles, list(range(len(alleles)))))
        nr_dict = dict(zip(self.sample["GT"].split("/"), self.sample["NR"].split(",")))
        nv_dict = dict(zip(self.sample["GT"].split("/"), self.sample["NV"].split(",")))

        new_vars = []
        for allele in self.alt.split(","):
            num = multi_allele_dict[allele]
            nr, nv = nr_dict[str(num)], nv_dict[str(num)]
            allele_var = copy.deepcopy(self)
            allele_var.sample["NR"] = nr
            allele_var.sample["NV"] = nv
            allele_var.alt = allele
            allele_var.hgvs = allele_var.construct_hgvs()
            new_vars.append(allele_var)
        return new_vars

    # Constructs hgvs string
    # I think this is the most basic construction, not totally confident it's considered
    # "proper" hgvs notation.
    def construct_hgvs(self):
        if len(self.ref) == 1:  # SNV
            hgvs = "{}:g.{}{}>{}".format(self.chrom, self.pos, self.ref, self.alt)
        else:  # All other variant types
            end = str(int(self.pos) + len(self.ref) - 1)
            hgvs = "{}:g.{}_{}{}>{}".format(self.chrom, self.pos, end, self.ref, self.alt)
        return hgvs

    # Adds annotation as attribute after checking if annotation and record position matches
    def add_annotation(self, annotation):
        if int(self.pos) != annotation["start"]:
            record = "{}-{}".format(self.chrom, self.pos)
            ann_record = "{}-{}".format(annotation["seq_region_name"], annotation["start"])
            print("Annotation position does not match variant record position, {} vs {}".format(record, ann_record))
        self.annotation = annotation

    # Parses annotation associated with Var Record
    def parse_annotation_consequence(self, most_severe_consequence):
        # Sets the consequence dict key based on some assumptions
        if most_severe_consequence == "regulatory_region_variant":
            consequence_key = "regulatory_feature_consequences"
        elif most_severe_consequence == "TF_binding_site_variant":
            consequence_key = "motif_feature_consequences"
        else:
            consequence_key = "transcript_consequences" if "transcript_consequences" in self.annotation.keys() else "intergenic_consequences"

        # Currently this is looking for consequences that match the most severe consequence
        # as well as the canonical transcripts. Only picks canonical transcripts if applicable.
        # For intergenic variants, this should always be a 1:1 mapping.
        severe_consequences = [consequence for consequence in self.annotation[consequence_key] if
                               most_severe_consequence in consequence["consequence_terms"]]
        canonical_consequences = [consequence for consequence in severe_consequences if
                                  consequence.get("canonical", 0) == 1]
        consequences = canonical_consequences if len(canonical_consequences) > 0 else severe_consequences
        return [(consequence.get("gene_symbol", "None"), consequence.get("gene_id", "None"),
                 consequence.get("transcript_id", "None"), consequence.get("canonical", 0)) for consequence in
                consequences]

    # Parses minor allele frequency, if available
    def parse_minor_allele_frequency(self, key="colocated_variants"):
        if key in self.annotation:
            for x in self.annotation["colocated_variants"]:
                if x.get("minor_allele_freq", None) and self.alt == x["minor_allele"]:
                    return str(x["minor_allele_freq"])
        return "None"

    # Write output file with appropriate columns
    def get_records_with_annotations(self):
        depth, alt_depth, vaf = self.read_coverage()
        maf = self.parse_minor_allele_frequency()
        most_severe_consequence = self.annotation["most_severe_consequence"]
        consequences = self.parse_annotation_consequence(most_severe_consequence)

        records = []
        for consequence in consequences:
            gene, gene_id, transcript_id, canonical = consequence
            records.append([self.chrom, self.pos, self.ref, self.alt, self.filter, depth, alt_depth, vaf, maf,
                            gene, gene_id, str(transcript_id), str(canonical), most_severe_consequence,
                            self.annotation["variant_class"]])

        return records


# Groups and sorts annotation, to give us a 1:1 mapping with variant records
def sort_annotations(annotations):
    """
    :param annotations: List of annotations, represented as JSON
    :return: Sorted, grouped by chromosome dictionary of annotations
    """
    annotation_dict = {}

    # First assign per chromosome
    for annotation in annotations:
        chrom = annotation["seq_region_name"]
        if chrom not in annotation_dict:
            annotation_dict[chrom] = []
        annotation_dict[chrom].append(annotation)

    # Sort by start position
    for chrom in annotation_dict:
        annotation_dict[chrom] = sorted(annotation_dict[chrom], key=lambda x: int(x["start"]))

    return annotation_dict


# Function to query VEP API
def vep_api(server, session, hgvs_list, task_num, ext="/vep/human/hgvs/?canonical=1&variant_class=1"):
    """
    :param server: VEP server to use, depends on genome assembly
    :param session: Requests session
    :param hgvs_list: Subset of hgvs annotations
    :param task_num: Task number, can be used to identify the subset of hgvs notations
    :param ext: URL endpoint, in this case the hgvs endpoint with some optional parameters included
    :return: Return annotation as JSON
    """
    url = server + ext
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    with session.post(url, headers=headers, data=json.dumps({"hgvs_notations": hgvs_list})) as response:
        data = response.json()
        if response.status_code != 200:  # If not OK status code, print failure
            print("FAILURE::{0}".format(task_num))

        # Tracks time of task
        elapsed_time = default_timer() - START_TIME
        completed_at = "{:5.2f}s".format(elapsed_time)
        print("{0:<30} {1:>20}".format(task_num, completed_at))

        return data


# Function to generate tasks, based on VEP API POST parameters
async def get_annotations(server, max_post_size, hgvs_list):
    """
    :param server: VEP server to use, depends on genome assembly
    :param max_post_size: Maximum number of hgvs to include in a single POST call
    :param hgvs_list: List of hgvs annotations
    :return: List of annotations from VEP API
    """
    # Number of posts needed to get annotations for each variant record
    num_posts = math.floor(len(hgvs_list) / max_post_size)
    num_posts = num_posts + 1 if (len(hgvs_list) % max_post_size) != 0 else num_posts

    print("{0:<30} {1:>20}".format("No", "Completed at"))

    results = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        with requests.Session() as session:
            loop = asyncio.get_event_loop()

            # Generates tasks per hgvs subset, based on the number of POST calls needed
            tasks = []
            for num in range(num_posts):
                start = num * max_post_size
                end = (num + 1) * max_post_size
                hgvs_subset = hgvs_list[start:end]
                task = loop.run_in_executor(executor, vep_api, *(server, session, hgvs_subset, num))
                tasks.append(task)

            # Gather tasks and add response to results
            for response in await asyncio.gather(*tasks):
                results.extend(response)

    return results


# Main entrypoint
def main(vcf_fh, output_file, server, max_post_size, vep_annotations=False):
    """
    :param vcf_fh: Open file handle for input VCF data
    :param output_file: File path for output TSV
    :param server: VEP server to use, depends on genome assembly
    :param max_post_size: Maximum number of loci to include in a single POST call
    :param vep_annotations: Optional parameter to include a pre-parsed annotations JSON file
    :return: TSV output file as well as a annotations JSON file, if annotations not provided
    """
    header = None
    var_records = []
    for line in vcf_fh:
        if line.startswith("##"):  # Skip metadata fields
            continue
        if line.startswith("#"):  # Keep header
            header = line.lstrip("#").rstrip("\n").split("\t")
            continue
        var_row = dict(zip(header, line.strip("\n").split("\t")))
        if len(var_row["ALT"].split(",")) > 1:  # Tease apart multi allelic sites
            var = VarRecord(var_row)
            var_records += var.parse_multi_allelic_sites()
        else:
            var_records.append(VarRecord(var_row))

    hgvs_list = [var.hgvs for var in var_records]  # Get hgvs annotations per variant record

    # Sets up asynchronous calling
    loop = asyncio.get_event_loop()
    future = asyncio.ensure_future(get_annotations(server, max_post_size, hgvs_list))
    annotations = loop.run_until_complete(future)

    # Annotations should map exactly 1:1 with variant records
    # Just need to group by chromosome, sort and collate to a single list
    annotation_dict = sort_annotations(annotations)
    sorted_annotations = []
    for key in annotation_dict.keys():
        sorted_annotations.extend(annotation_dict[key])

    # Map var record with annotation and write to output
    out_fh = open(output_file, 'w')
    out_fh.write("\t".join(out_header))
    out_fh.write("\n")
    count = 0
    for var, annotation in zip(var_records, sorted_annotations):
        var.add_annotation(annotation)
        records = var.get_records_with_annotations()
        for record in records:
            out_fh.write("\t".join(record))
            out_fh.write("\n")

        if len(records) > 0:  # If these records are non-zero, annotations were added successfully
            count += 1
    out_fh.close()

    # Print to confirm that each variant record is accounted for
    # Expectation is these 2 numbers should be the same, equal to the number of variant rows in original vcf
    print(count)
    print(len(var_records))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate variants with VEP')
    parser.add_argument('input_vcf', type=open,
                        help='Input VCF File')
    parser.add_argument('output_file',
                        help='Output TSV File')
    parser.add_argument('--hg19', action='store_true', help='If genome is hg19 and not hg38')
    args = parser.parse_args()

    vep_server = "http://grch37.rest.ensembl.org" if args.hg19 else "https://rest.ensembl.org"  # Sets VEP server
    vep_post_size = 300 if args.hg19 else 200

    main(args.input_vcf, args. output_file, vep_server, vep_post_size)
