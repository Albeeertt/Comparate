import numpy as np


class Criterion:

    def __init__(self, threshold_complete_match, threshold_overlap_match):
        self.threshold_complete_match = threshold_complete_match
        self.threshold_overlap_match = threshold_overlap_match

    def multiple_complete_match(self, record, list_match):
        types_match = []
        for match in list_match:
            types_match.append(match['type'])
        unique_values = np.unique(types_match)
        if len(unique_values) == 1 and unique_values[0] == 'intergenic_region' and unique_values[0] != record['type']:
            if record['prob_intergenic_region'] >= self.threshold_complete_match:
                return True, 'intergenic_region'
            return False, 'intergenic_region'
        elif record['type'] == 'intergenic_region' and (len(unique_values) > 1 or (len(unique_values) == 1 and unique_values[0] == 'gene' and unique_values[0] != record['type'])):
            if record['prob_gene'] >= self.threshold_complete_match:
                return True, 'gene'
            return False, 'gene'
        else:
            return -1, -1
        
    def multiple_overlap_match(self, record, list_match):
        types_match = []
        for match in list_match:
            types_match.append(match['type'])
        unique_values = np.unique(types_match)
        if len(unique_values) == 1 and unique_values[0] == 'intergenic_region':
            return -1
        elif len(unique_values) > 1 or (len(unique_values) == 1 and unique_values[0] == 'gene'):
            if record['prob_gene'] >= self.threshold_overlap_match:
                return True
            return False
        
    def single_complete_match_gene(self, record):
        if record['prob_gene'] >= self.threshold_complete_match:
            return True
        return False
    
    def single_overlap_match_gene(self, record):
        if record['prob_gene'] >= self.threshold_overlap_match:
            return True
        return False 

    def single_complete_match_ir(self, record):
        if record['prob_intergenic_region'] >= self.threshold_complete_match:
            return True
        return False

    def single_overlap_match_ir(self, record):
        if record['prob_intergenic_region'] >= self.threshold_overlap_match:
            return True
        return False