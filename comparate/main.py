
import os
import argparse
import pandas as pd

from comparate.core.HandleFile import HandleFile
from comparate.core.Comparative import Comparate
from comparate.core.Comparative_TE import Comparative_TE
from comparate.core.Comparative_Pseudogene import Comparative_Pseudogene


def obtain_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--route_csv_truth', required=True, type=str, help='Path to CSV file.')
    parser.add_argument('--out', type=str, required=True, help='Output path.')
    parser.add_argument('--route_csv_other', type=str, help='Path to CSV file.')
    parser.add_argument('--route_TE_truth', type=str, help='Path to TE (truth CSV) in GFF file.')
    parser.add_argument('--route_TE_other', type=str, help='Path to TE (other CSV) in GFF file.')
    parser.add_argument('--compare', action='store_true', help='compare route_csv_truth and route_csv_other.')
    parser.add_argument('--route_pseudogene_truth', type=str, help='Path to Pseudogene tr file.')
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
    ROUTE_GRAPH_TE_IR = route_out+'te_ir_truth.png'
    ROUTE_GRAPH_TE_GEN = route_out+'te_gen_truth.png'
    ROUTE_GRAPH_TE_IR_OTHER = route_out+'te_ir_other.png'
    ROUTE_GRAPH_TE_GEN_OTHER = route_out+'te_gen_other.png'
    ROUTE_GRAPH_PSEUDO_IR = route_out+'pseudogene_ir.png'
    ROUTE_GRAPH_PSEUDO_GENES = route_out+'pseudogene_gen.png'
        
    instance_HandleFile_truth = HandleFile(args.route_csv_truth)
    summary_truth = []
    summary_truth.append(instance_HandleFile_truth.summary_file())
    summary_truth.extend(instance_HandleFile_truth.summary_chr())
    pd.DataFrame(summary_truth).to_csv(ROUTE_SUMMARY_TRUTH)

    if args.route_pseudogene_truth:
        instance_comparative_Pseudogene = Comparative_Pseudogene(args.route_csv_truth, args.route_pseudogene_truth)
        result_gen_pseudo, result_ir_pseudo = instance_comparative_Pseudogene.comparate_pseudogenes(threshold=.5)
        instance_comparative_Pseudogene.generate_graph(result_gen_pseudo, "Conteo de los Pseudogenes en genes", f"Pseudogenes en genes ({instance_comparative_Pseudogene.name_truth})", ROUTE_GRAPH_PSEUDO_IR)
        instance_comparative_Pseudogene.generate_graph(result_ir_pseudo, "Conteo de los Pseudogenes en regiones intergénicas", f"Pseudogenes en regiones intergénicas ({instance_comparative_Pseudogene.name_truth})", ROUTE_GRAPH_PSEUDO_GENES)

        
    if args.route_TE_truth:
        instance_comparate_TE = Comparative_TE(args.route_csv_truth, args.route_TE_truth)
        result_ir_te_unicos, result_gen_te_unicos = instance_comparate_TE.comparate_te()
        instance_comparate_TE.generate_graph(result_ir_te_unicos, "Conteo de los TE en regiones intergénicas clasificados como gen", f"Elementos transponibles en regiones intergénicas ({instance_comparate_TE.name_truth})", ROUTE_GRAPH_TE_IR)
        instance_comparate_TE.generate_graph(result_gen_te_unicos, "Conteo de los TE en genes clasificados como gen", f"Elementos transponibles en genes ({instance_comparate_TE.name_truth})", ROUTE_GRAPH_TE_GEN)

    if args.route_csv_other:
        instance_HandleFile_other = HandleFile(args.route_csv_other)
        summary_other = []
        summary_other.append(instance_HandleFile_other.summary_file())
        summary_other.extend(instance_HandleFile_other.summary_chr())
        pd.DataFrame(summary_other).to_csv(ROUTE_SUMMARY_OTHER)

    if args.route_csv_other and args.compare:
        instance_comparate = Comparate(args.route_csv_truth, args.route_csv_other, threshold_complete_match=args.threshold_complete_match, threshold_overlap_match=args.threshold_overlap_match)
        comparate_files = []
        complete, overlap = instance_comparate.comparate_files()
        comparate_files.append(complete)
        comparate_files.append(overlap)

        comparate_files.extend(instance_comparate.comparate_chrs())

        pd.DataFrame(comparate_files).to_csv(ROUTE_COMPARATE)

    if args.route_csv_other and args.route_TE_other:
        instance_comparate_TE = Comparative_TE(args.route_csv_other, args.route_TE_other)
        result_ir_te_unicos, result_gen_te_unicos = instance_comparate_TE.comparate_te()
        instance_comparate_TE.generate_graph(result_ir_te_unicos, "Conteo de los TE en regiones intergénicas clasificados como gen", f"Elementos transponibles en regiones intergénicas ({instance_comparate_TE.name_truth})", ROUTE_GRAPH_TE_IR_OTHER)
        instance_comparate_TE.generate_graph(result_gen_te_unicos, "Conteo de los TE en genes clasificados como gen", f"Elementos transponibles en genes ({instance_comparate_TE.name_truth})", ROUTE_GRAPH_TE_GEN_OTHER)
