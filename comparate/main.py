
import os
import argparse
import pandas as pd

from comparate.core.HandleFile import HandleFile
from comparate.core.Comparative import Comparate


def obtain_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--route_csv_truth', required=True, type=str, help='Path to CSV file.')
    parser.add_argument('--route_csv_other', required=True, type=str, help='Path to CSV file.')
    parser.add_argument('--out', type='str', required=True, help='Output path.')
    parser.add_argument('--route_te', type=str, help='Path to GFF file.')
    parser.add_argument('--threshold_complete_match', type=float, default=.5, help='threshold_complete_match')
    parser.add_argument('--threshold_overlap_match', type=float, default=.2, help='threshold_overlap_match')
    return parser.parse_args()

def execute_main_program():
    
    args = obtain_arguments()

    route_out = args.out+'/' if not args.out.endswith('/') else args.out

    if not os.path.exists(route_out):
        os.mkdir(route_out)

    ROUTE_SUMMARY_TRUTH = route_out+'summary_truth.csv'
    ROUTE_SUMMARY_OTHER = route_out+'summary_other.csv'
    ROUTE_COMPARATE = route_out+'comparate.csv'
        
    instance_HandleFile_truth = HandleFile(args.route_csv_truth)
    instance_HandleFile_other = HandleFile(args.route_csv_other)
    instance_comparate = Comparate(args.route_csv_truth, args.route_csv_other, threshold_complete_match=args.threshold_complete_match, threshold_overlap_match=args.threshold_overlap_match)

    summary_truth = []
    summary_truth.append(instance_HandleFile_truth.summary_file())
    summary_truth.extend(instance_HandleFile_truth.summary_chr())

    summary_other = []
    summary_other.append(instance_HandleFile_other.summary_file())
    summary_other.extend(instance_HandleFile_other.summary_chr())

    pd.DataFrame(summary_truth).to_csv(ROUTE_SUMMARY_TRUTH)
    pd.DataFrame(summary_other).to_csv(ROUTE_SUMMARY_OTHER)

    comparate_files = []
    complete, overlap = instance_comparate.comparate_files()
    comparate_files.append(complete)
    comparate_files.append(overlap)

    comparate_files.extend(instance_comparate.comparate_chrs())

    pd.DataFrame(comparate_files).to_csv(ROUTE_COMPARATE)

