import argparse

def parse_arguments():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Predict novel human EEIs based on mouse EEIs and orthology data')
    parser.add_argument('--egio', required=True, help='EGIO output file with orthologous exon mappings')
    parser.add_argument('--eei', required=True, help='Mouse EEI network file')
    parser.add_argument('--human_exon_map', required=True, help='Human exon ID to coordinate mapping file')
    parser.add_argument('--mouse_exon_map', required=True, help='Mouse exon ID to coordinate mapping file')
    parser.add_argument('--human_eei', required=False, help='Optional existing human EEI file to filter out known interactions')
    parser.add_argument('--identity_threshold', type=float, default=0.8, 
                        help='Minimum sequence identity to consider orthology (default: 0.8)')
    parser.add_argument('--output', required=True, help='Output file for predicted EEIs')
    return parser.parse_args()
