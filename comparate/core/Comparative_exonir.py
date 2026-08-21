
import pandas as pd
import numpy as np
from collections import defaultdict

from .Criterion import Criterion


class Comparate:

    def __init__(self, route_truth: str, route_other: str, threshold_complete_match, threshold_overlap_match):
        self.name_truth = route_truth.split('/')[-1]
        self.route_truth = route_truth
        self.name_other = route_other.split('/')[-1]
        self.route_other = route_other

        self.truth_file_read = self.obtain_gff(route_truth)
        self.other_file_read = self.obtain_gff(route_other)

        self.instance_criterion = Criterion(threshold_complete_match, threshold_overlap_match)


    def obtain_gff(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''
        - Read the GFF3 located in 'route' with the encoding specified in 'encoding'. 
          Additionally, add the old_idx column, indicating the location of each sample 
          in the original dataframe (for future modifications).
        '''
        data = pd.read_csv(route, index_col=0, encoding= encoding)
        data['old_idx'] = data.index
        data['strand'] = data['strand'].replace('.', '+')
        return data
    

    def arranged_dataFrame_function(self, file: pd.DataFrame):
        arranged = defaultdict(lambda: defaultdict(list))
        entire_dataFrame = file[file['type'].isin(['intergenic_region', 'gene'])]
        for record in entire_dataFrame.to_dict(orient='records'):
            if pd.isna(record['prob_gene']):
                continue
            arranged[record['chr']][record['strand']].append(record)
        return arranged

    def capture_complete_and_overlap(self, arranged_dataFrame, other_file_dataFrame):
        list_complete_match = []
        list_overlap_match = []
        entireOther_dataFrame = other_file_dataFrame[other_file_dataFrame['type'].isin(['intergenic_region', 'gene'])]
        for record in entireOther_dataFrame.to_dict(orient='records'):
            if pd.isna(record['prob_gene']):
                continue
            list_indicated = arranged_dataFrame[record['chr']][record['strand']]
            list_complete = []
            list_overlap = []
            for perharps_overlap in list_indicated:
                if record['start'] < perharps_overlap['start'] and record['end'] > perharps_overlap['end']:
                    list_complete.append(perharps_overlap)
                elif record['end'] > perharps_overlap['start'] and perharps_overlap['end'] > record['start']:
                    list_overlap.append(perharps_overlap)
            list_complete_match.append((record, list_complete))
            list_overlap_match.append((record, list_overlap))

        return list_complete_match, list_overlap_match
    
    def comparate_complete_match(self, list_complete_match):
        result_more_than_one = []
        result_different_gene = []

        for record, complete_match in list_complete_match:
            if len(complete_match) == 1 and complete_match[0]['type'] != record['type']:
                if record['type'] == 'intergenic_region' and complete_match[0]['type'] == 'gene':
                    result_different_gene.append(self.instance_criterion.single_complete_match_gene(record))
            elif len(complete_match) > 1:
                result_multiple_match = self.instance_criterion.multiple_complete_match(record, complete_match)
                if result_multiple_match != -1:
                    result_more_than_one.append(result_multiple_match)

        return result_different_gene, result_more_than_one

    
    def comparate_overlap_match(self, list_overlap_match):
        result_more_than_one = []
        result_different_gene = []

        for record, complete_match in list_overlap_match:
            if len(complete_match) == 1 and complete_match[0]['type'] != record['type']:
                if record['type'] == 'intergenic_region' and complete_match[0]['type'] == 'gene':
                    result_different_gene.append(self.instance_criterion.single_overlap_match_gene(record))
            elif len(complete_match) > 1:
                result_multiple_match = self.instance_criterion.multiple_overlap_match(record, complete_match)
                if result_multiple_match != -1:
                    result_more_than_one.append(result_multiple_match)
        
        return result_different_gene, result_more_than_one

    def comparate_files(self):
        new_record_complete = {'Mode': 'complete_match', 'chr': 'all',  'gene_truth': 0, 'gen_within_ir': 0}
        new_record_overlap = {'Mode': 'overlap_match', 'chr': 'all', 'gene_truth': 0, 'gen_within_ir': 0}

        arranged_dataFrame = self.arranged_dataFrame_function(self.truth_file_read)
        list_complete_match, list_overlap_match = self.capture_complete_and_overlap(arranged_dataFrame, self.other_file_read)
        result_complete_different_gene, result_complete_more_than_one = self.comparate_complete_match(list_complete_match)
        result_overlap_different_gene, result_overlap_more_than_one = self.comparate_overlap_match(list_overlap_match)
        
        if result_complete_different_gene:
            new_record_complete['gene_truth'] = sum(result_complete_different_gene) / len(result_complete_different_gene)
        if result_complete_more_than_one:
            new_record_complete['gen_within_ir'] = sum(result_complete_more_than_one) / len(result_complete_more_than_one)

        if result_overlap_different_gene:
            new_record_overlap['gene_truth'] = sum(result_overlap_different_gene) / len(result_overlap_different_gene)
        if result_overlap_more_than_one:
            new_record_overlap['gen_within_ir'] = sum(result_overlap_more_than_one) / len(result_overlap_more_than_one)

        return new_record_complete, new_record_overlap
    
    def comparate_chrs(self):
        chrs_in_file = np.unique(self.other_file_read['chr'])
        results_chrs = []
        arranged_dataFrame = self.arranged_dataFrame_function(self.truth_file_read)
        otherFile_dataFrame = self.other_file_read
        for key_chr in chrs_in_file:
            new_record_complete = {'Mode': 'complete_match', 'chr': key_chr,  'gene_truth': 0, 'gen_within_ir': 0}
            new_record_overlap = {'Mode': 'overlap_match', 'chr': key_chr, 'gene_truth': 0, 'gen_within_ir': 0}

            dataFrame_chr_otherFile = otherFile_dataFrame[otherFile_dataFrame['chr'] == key_chr]

            list_complete_match, list_overlap_match = self.capture_complete_and_overlap(arranged_dataFrame, dataFrame_chr_otherFile)
            result_complete_different_gene, result_complete_more_than_one = self.comparate_complete_match(list_complete_match)
            result_overlap_different_gene, result_overlap_more_than_one = self.comparate_overlap_match(list_overlap_match)

            if result_complete_different_gene:
                new_record_complete['gene_truth'] = sum(result_complete_different_gene) / len(result_complete_different_gene)
            if result_complete_more_than_one:
                new_record_complete['gen_within_ir'] = sum(result_complete_more_than_one) / len(result_complete_more_than_one)

            if result_overlap_different_gene:
                new_record_overlap['gene_truth'] = sum(result_overlap_different_gene) / len(result_overlap_different_gene)
            if result_overlap_more_than_one:
                new_record_overlap['gen_within_ir'] = sum(result_overlap_more_than_one) / len(result_overlap_more_than_one)

            results_chrs.append(new_record_complete)
            results_chrs.append(new_record_overlap)

        return results_chrs