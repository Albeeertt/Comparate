import pandas as pd
import numpy as np

class HandleFile:

    def __init__(self, route_file: str):
        self.name = route_file.split('/')[-1]
        self.route_file = route_file
        self.file_read = self.obtain_gff(route_file)
        self.HARD_LIMITS_PROB = (.5, .6)
        self.SOFT_LIMITS_PROB = (.6, .7)

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
    
    def summary_file(self):
        new_record = {'Annotation': self.name, 'chr': 'all', 'total_gene': 0, 'total_ir': 0, 'genes_contaminated': 0, 'genes_contaminated_permissive': 0, 'correct_gene': 0, 'ir_contaminated': 0, 'ir_contaminated_permissive': 0, 'correct_ir': 0} 
        for record in self.file_read.to_dict(orient='records'):
            if pd.isna(record['prob_gene']):
                continue
            if record['type'] == 'gene':
                new_record['total_gene'] += 1
                if self.HARD_LIMITS_PROB[0] <= record['prob_gene'] <= self.HARD_LIMITS_PROB[1]:
                    new_record['genes_contaminated'] += 1
                if self.SOFT_LIMITS_PROB[0] < record['prob_gene'] <= self.SOFT_LIMITS_PROB[1]:
                    new_record['genes_contaminated_permissive'] += 1
                if record['Result'] == 'gene':
                    new_record['correct_gene'] += 1
            elif record['type'] == 'intergenic_region':
                new_record['total_ir'] += 1
                if self.HARD_LIMITS_PROB[0] <= record['prob_intergenic_region'] <= self.HARD_LIMITS_PROB[1]:
                    new_record['ir_contaminated'] += 1
                if self.SOFT_LIMITS_PROB[0] < record['prob_intergenic_region'] <= self.SOFT_LIMITS_PROB[1]:
                    new_record['ir_contaminated_permissive'] += 1
                if record['Result'] == 'intergenic_region':
                    new_record['correct_ir'] += 1
            
        new_record['Recall-genes'] = new_record['correct_gene'] / new_record['total_gene']
        new_record['Recall-ir'] = new_record['correct_ir'] / new_record['total_ir']

        new_record['accuracy'] = (new_record['correct_gene']+new_record['correct_ir'])/(new_record['total_gene']+new_record['total_ir'])
        return new_record
    
    def summary_chr(self):
        chrs_in_file = np.unique(self.file_read['chr'])
        entire_dataFrame = self.file_read
        list_result_chrs = []
        for key_chr in chrs_in_file:
            new_record = {'Annotation': self.name, 'chr': key_chr, 'total_gene': 0, 'total_ir': 0, 'genes_contaminated': 0, 'genes_contaminated_permissive': 0, 'correct_gene': 0, 'ir_contaminated': 0, 'ir_contaminated_permissive': 0, 'correct_ir': 0} 
            only__chr_dataFrame = entire_dataFrame[entire_dataFrame['chr'] == key_chr]
            for record in only__chr_dataFrame.to_dict(orient='records'):
                if pd.isna(record['prob_gene']):
                    continue
                if record['type'] == 'gene':
                    new_record['total_gene'] += 1
                    if self.HARD_LIMITS_PROB[0] <= record['prob_gene'] <= self.HARD_LIMITS_PROB[1]:
                        new_record['genes_contaminated'] += 1
                    if self.SOFT_LIMITS_PROB[0] < record['prob_gene'] <= self.SOFT_LIMITS_PROB[1]:
                        new_record['genes_contaminated_permissive'] += 1
                    if record['Result'] == 'gene':
                        new_record['correct_gene'] += 1
                elif record['type'] == 'intergenic_region':
                    new_record['total_ir'] += 1
                    if self.HARD_LIMITS_PROB[0] <= record['prob_intergenic_region'] <= self.HARD_LIMITS_PROB[1]:
                        new_record['ir_contaminated'] += 1
                    if self.SOFT_LIMITS_PROB[0] < record['prob_intergenic_region'] <= self.SOFT_LIMITS_PROB[1]:
                        new_record['ir_contaminated_permissive'] += 1
                    if record['Result'] == 'intergenic_region':
                        new_record['correct_ir'] += 1
                    
            new_record['Recall-genes'] = new_record['correct_gene'] / new_record['total_gene']
            new_record['Recall-ir'] = new_record['correct_ir'] / new_record['total_ir']
            new_record['accuracy'] = (new_record['correct_gene']+new_record['correct_ir'])/(new_record['total_gene']+new_record['total_ir'])
            list_result_chrs.append(new_record)
        return list_result_chrs