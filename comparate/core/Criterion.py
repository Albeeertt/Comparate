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
        if len(unique_values) == 1 and unique_values[0] == 'intergenic_region':
            return -1
        elif len(unique_values) > 1 or (len(unique_values) == 1 and unique_values[0] == 'gene'):
            if record['prob_gene'] >= self.threshold_complete_match:
                return True
            return False
        
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