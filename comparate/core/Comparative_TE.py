
import pandas as pd
import numpy as np
from collections import defaultdict

import matplotlib.pyplot as plt

from typing import Dict, List

class Comparative_TE: 

    def __init__(self,  route_truth: str, te_route_truth: str):
        self.name_truth = route_truth.split('/')[-1]
        self.route_truth = route_truth
        self.route_te = te_route_truth

        self.truth_file_read = self.obtain_csv(route_truth)
        self.te_file_read = self.obtain_gff(te_route_truth)

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
    
    def obtain_structure_TE(self, value_kimura: float = 1.) -> Dict:
        structure_te = defaultdict(lambda: defaultdict(list))
        for record in self.te_file_read.to_dict(orient='records'):
            dict_attributes = { part.split('=')[0]:part.split('=')[1] for part in record['attributes'].split(';')}
            if dict_attributes.get('KIMURA80', -1) != -1:
                record['KIMURA80'] = float(dict_attributes['KIMURA80'])
                if record['KIMURA80'] <= value_kimura:
                    structure_te[record['chr']][record['strand']].append(record)
            else:
                structure_te[record['chr']][record['strand']].append(record)
        return structure_te
    
    def comparate_te(self, threshold: float = .5):
        result_ir_te_unicos = defaultdict(lambda: {
            f'prob_mayor_{threshold}': 0,
            f'prob_menor_{threshold}': 0
        })
       
        result_gen_te_unicos = defaultdict(lambda: {
            f'prob_mayor_{threshold}': 0,
            f'prob_menor_{threshold}': 0
        }) 
        structure_te = self.obtain_structure_TE()
        for record in self.truth_file_read.to_dict(orient='records'):
            if pd.isna(record['prob_gene']):
                continue

            tipo_te_lista = []

            list_records_chr: List[Dict] = structure_te[record['chr']][record['strand']]
            posee_te: bool = False
            for record_te in list_records_chr:
                # if record['end'] > record_te['start'] and record_te['end'] > record['start']:
                if record['start'] <= record_te['start'] and record['end'] >= record_te['end']:
                    posee_te = True
                    tipo_te_lista.append(record_te['type'])

            tipo_te_lista = np.unique(tipo_te_lista)
                    
            if posee_te and record['type'] == 'gene':
                if record['prob_gene'] >= threshold:
                    for tipoo in tipo_te_lista:
                        result_gen_te_unicos[tipoo][f'prob_mayor_{threshold}'] += 1
                else:
                    for tipoo in tipo_te_lista:
                        result_gen_te_unicos[tipoo][f'prob_menor_{threshold}'] += 1
            elif posee_te and record['type'] == 'intergenic_region':
                if record['prob_gene'] >= threshold:
                    for tipoo in tipo_te_lista:
                        result_ir_te_unicos[tipoo][f'prob_mayor_{threshold}'] += 1            
                else:
                    for tipoo in tipo_te_lista:
                        result_ir_te_unicos[tipoo][f'prob_menor_{threshold}'] += 1      

        return result_ir_te_unicos, result_gen_te_unicos

    def generate_graph(self, results_to_graph: Dict, y_label: str, title: str, out: str):
        labels = list(results_to_graph.keys())
        mayor = [results_to_graph[k]['prob_mayor_0.5'] for k in labels]
        menor = [results_to_graph[k]['prob_menor_0.5'] for k in labels]

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

