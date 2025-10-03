import csv
import argparse
import os

def parse_gtf_for_exon_mappings(gtf_file):
    """Parse GTF file to create exon ID to coordinate mapping"""
    id_to_coord = {}
    coord_to_id = {}

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if fields[2] != 'exon':
                continue

            chrom = fields[0]
            start = fields[3]
            end = fields[4]
            strand = "-1" if fields[6] == "-" else "1"

            attr = fields[8]
            exon_id = None
            for attr_pair in attr.split(';'):
                if 'exon_id' in attr_pair:
                    exon_id = attr_pair.split('"')[1]
                    break

            if exon_id:
                coord = f"chr{chrom}:{start}:{end}:{strand}"
                id_to_coord[exon_id] = coord
                coord_to_id[coord] = exon_id

    return {'id_to_coord': id_to_coord, 'coord_to_id': coord_to_id}


def save_combined_mapping_to_tsv(mapping_dict, output_file):
    """Save combined exon_id, coord, and reverse mapping to a TSV file"""
    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["exon_id", "coord", "reverse_lookup"])

        for exon_id, coord in mapping_dict['id_to_coord'].items():
            reverse_lookup = mapping_dict['coord_to_id'].get(coord, "N/A")
            writer.writerow([exon_id, coord, reverse_lookup])


def main():
    parser = argparse.ArgumentParser(description="Parse a GTF file to create exon ID to coordinate mappings and save to TSV.")
    parser.add_argument("--gtf", required=True, help="Path to the input GTF file for the species.")
    parser.add_argument(
        "--output",
        required=False,
        help="Path to the output TSV file. If not provided, will create alongside the GTF with suffix .exon_coord_mapping.tsv",
    )

    args = parser.parse_args()

    gtf_file_path = args.gtf
    if args.output:
        output_file_path = args.output
    else:
        base_name = os.path.splitext(os.path.basename(gtf_file_path))[0]
        default_name = f"{base_name}.exon_coord_mapping.tsv"
        output_file_path = os.path.join(os.path.dirname(gtf_file_path), default_name)

    mapping = parse_gtf_for_exon_mappings(gtf_file_path)
    save_combined_mapping_to_tsv(mapping, output_file_path)
    print(f"Mapping saved to {output_file_path}")


if __name__ == "__main__":
    main()
