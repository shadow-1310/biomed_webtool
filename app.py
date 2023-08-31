import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# import streamlit as st
import os

def make_df(path):
    all_files = []
    ids = []
    chains = []
    units = []
    all_species = []
    sequences = []
    tot = []

    def process_lines(lines, file):
        if not lines:
            return

        lines = [l.strip('\n') for l in lines]
        meta_data = lines[0].split('|')
        sequence = lines[1]
        
        pdb_id = meta_data[0].strip('>')
        chain = meta_data[1]
        unit_details = meta_data[2]
        species = meta_data[3]

        ids.append(pdb_id)
        chains.append(chain)
        units.append(unit_details)
        all_species.append(species)
        sequences.append(''.join(sequence))
        tot.append(len(''.join(sequence)))

        print(len(sequences[-1]))
        all_files.append(file)

    def read_files(BASE):
        files = os.listdir(BASE)
        for file in files:
            file_path = os.path.join(BASE, file)
            print(file_path)
            with open(file_path, 'r') as f:
                file_name = os.path.basename(file_path)
                curr = []
                while True:
                    line = f.readline()
                    if not line:
                        process_lines(curr, file_name)
                        break
                    if line[0] == '>':
                        process_lines(curr, file_name)
                        curr = []

                    curr.append(line)
    read_files(path)

    dict = {
            'file_name' : all_files,
            'PDB_ID' : ids,
            'chain_number' : chains,
            'unit' : units,
            'species_name' : all_species,
            'AA_sequence' : sequences,
            'total_lenth' : tot,
            }
    df = pd.DataFrame(dict)

    return df

def calc_percent(string, aa, tot):
    """
    A function to calculate the percent of a given Amino Acid within that whole sequence
    It takes 3 parameters, 
    string: amino-acid sequence
    aa : amino acid whose percent is to be calculated
    tot : total length of the given amino acid sequence
    """
    count = {aa : 0}

    for a in string:
        if a == aa:
            count[aa] += 1

    percent = np.round((count[aa] / tot)*100, 2)
    return percent 

class aa_data:
    def __init__(self, path):
        self.df = make_df(path)

    def add_percent(self, amino_acids):
        for aa in amino_acids:
            self.df[aa] = self.df.apply(lambda x: calc_percent(x['AA_sequence'], aa, x['total_lenth']), axis=1)
    
    def make_exclusion(self, pdb_ids):
        self.df = self.df[~self.df['PDB_ID'].isin(pdb_ids)]


    def plot_box(self):
        self.df_plot = self.df.iloc[:, 7:]
        plt.style.use("ggplot")
        plt.figure(figsize=(15, 15))
        sns.boxplot(self.df_plot, orient="v")
        plt.xticks(rotation='vertical', size = 14)
        plt.yticks(size = 14)
        plt.savefig("boxplot.eps",format = 'eps', bbox_inches='tight', dpi=1200)
        plt.show()

if __name__ == '__main__':
    BASE = "data"
    a = aa_data(BASE)
    a.add_percent(['A', 'C', 'D', 'E', 'H','I', 'J', 'K','F'])
    a.make_exclusion(['7QRU_4', '7QRU_5', '7QRU_7'])
    print(a.df)
    print(a.df.columns)
    a.plot_box()


# aa_sequence = 

