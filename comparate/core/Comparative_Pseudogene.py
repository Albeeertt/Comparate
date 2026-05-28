
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from typing import List, Dict, Tuple

class Comparative_Pseudogene:

    def __init__(self, route_truth: str, route_pseudogene: str):
        self.name_truth = route_truth.split('/')[-1]
        self.route_truth = route_truth
        self.route_pseudogene = route_pseudogene

        self.truth_file_read = self.obtain_csv(route_truth)
        self.pseudogene_file_read = self.obtain_gff(route_pseudogene)


    def obtain_gff(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:

        data = pd.read_csv(route, comment='#', sep='\t', header=None, encoding=encoding)
        data.columns = ['chr',	'db',	'type',	'start',	'end',	'score',	'strand',	'phase',	'attributes']
        data['old_idx'] = data.index
        data['strand'] = data['strand'].replace('.', '+')
        return data


    def obtain_csv(self, route: str, encoding: str = 'utf-8') -> pd.DataFrame:
        '''
        - Read the GFF3 located in 'route' with the encoding specified in 'encoding'. 
          Additionally, add the old_idx column, indicating the location of each sample 
          in the original dataframe (for future modifications).
        '''
        data = pd.read_csv(route, index_col=0, encoding= encoding)
        data['old_idx'] = data.index
        data['strand'] = data['strand'].replace('.', '+')
        return data
    

    def obtain_structure_pseudogenes(self) -> Dict:
        structure_pseudogenes = defaultdict(lambda: defaultdict(list))
        for record in self.pseudogene_file_read.to_dict(orient='records'):
            if record['type'] == 'pseudogene':
                structure_pseudogenes[record['chr']][record['strand']].append(record)
        return structure_pseudogenes
    
    def comparate_pseudogenes(self, threshold: float = .5) -> Tuple[Dict, Dict]:
        result_ir_pseudo = {
            f'prob_mayor_{threshold}': 0,
            f'prob_menor_{threshold}': 0
        }
        result_gen_pseudo = {
            f'prob_mayor_{threshold}': 0,
            f'prob_menor_{threshold}': 0
        }
        structure_pseudogenes = self.obtain_structure_pseudogenes()
        for record in self.truth_file_read.to_dict(orient='records'):
            if pd.isna(record['prob_gene']):
                continue
        
            list_record_chr: List[Dict] = structure_pseudogenes[record['chr']][record['strand']]
            posee_pseudo: bool = False
            for record_pseudo in list_record_chr:
                if record['start'] <= record_pseudo['start'] and record['end'] >= record_pseudo['end']:
                    posee_pseudo = True
                    break
            
            if posee_pseudo and record['type'] == 'gene':
                if record['prob_gene'] >= threshold:
                    result_gen_pseudo[f'prob_mayor_{threshold}'] += 1
                else:
                    result_gen_pseudo[f'prob_menor_{threshold}'] += 1
            elif posee_pseudo and record['type'] == 'intergenic_region':
                if record['prob_gene'] >= threshold:
                    result_ir_pseudo[f'prob_mayor_{threshold}'] += 1 
                else:
                    result_ir_pseudo[f'prob_menor_{threshold}'] += 1
        return result_gen_pseudo, result_ir_pseudo
    

    def generate_graph(self, results_to_graph: Dict, y_label: str, title: str, out: str):
        labels = ['pseudogenes']
        mayor = [results_to_graph['prob_mayor_0.5']]
        menor = [results_to_graph['prob_menor_0.5']]

        x = np.arange(len(labels))
        width = 0.4
        
        plt.figure(figsize=(12,6))
        plt.bar(x - width/2, mayor, width, label='prob > 0.5')
        plt.bar(x + width/2, menor, width, label='prob < 0.5')
        plt.xticks(x, labels, rotation=45, ha='right')
        plt.ylabel(y_label)
        plt.title(title)
        plt.legend()

        plt.tight_layout()

        plt.savefig(out, dpi=300, bbox_inches='tight')