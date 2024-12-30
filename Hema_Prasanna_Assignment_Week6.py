import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.colors as mcolors

plt.style.use('BME163.mplstyle') 
def parse_args():
    parser = argparse.ArgumentParser(description="Generate a heatmap from gene expression data.")
    parser.add_argument('-e', '--exp_file', required=True, help='Expression data file path')
    parser.add_argument('-p', '--phase_file', required=True, help='Phase data file path')
    parser.add_argument('-g', '--genes', help='Comma-separated list of genes for tick marks', default='')
    return parser.parse_args()

def normalize_data(values):
    values = np.array(values, dtype=float)
    return ((values - np.min(values)) / (np.max(values) - np.min(values))) * 100

def load_data(filepath, delimiter='\t', header=True):
    with open(filepath, 'r') as file:
        if header:
            headers = file.readline().strip().split(delimiter)
        data = [line.strip().split(delimiter) for line in file]
    return headers, data

def interpolate_data(data, factor=2):
    new_length = data.shape[0] * factor
    x_original = np.arange(data.shape[0])
    x_new = np.linspace(0, data.shape[0] - 1, new_length)
    data_interpolated = np.array([np.interp(x_new, x_original, data[:, i]) for i in range(data.shape[1])]).T
    return data_interpolated

def main():
    args = parse_args()
    exp_headers, exp_data = load_data(args.exp_file)
    phase_headers, phase_data = load_data(args.phase_file)

    gene_index = exp_headers.index('Ensembl_ID')
    FPKM_indices = [i for i, header in enumerate(exp_headers) if 'FPKM_CT' in header]
    phase_dict = {row[0].strip(): float(row[phase_headers.index('Peak_phase(CT)')]) for row in phase_data}

    exp_dict = {row[gene_index].strip(): np.array([float(row[i]) for i in FPKM_indices]) for row in exp_data if row[gene_index].strip() in phase_dict}
    sorted_genes = sorted(exp_dict.keys(), key=lambda x: phase_dict[x])
    normalized_matrix = np.array([normalize_data(exp_dict[gene]).astype(int) for gene in sorted_genes])

    normalized_matrix = interpolate_data(normalized_matrix, factor=5) 

    viridian = mcolors.LinearSegmentedColormap.from_list("viridian", [
        (253/255, 231/255, 37/255), (94/255, 201/255, 98/255),
        (33/255, 145/255, 140/255), (59/255, 82/255, 139/255),
        (68/255, 1/255, 84/255)]).reversed()

    fig = plt.figure(figsize=(5, 3))
    panel1 = plt.axes([0.7/5 , 0.3/3 , 0.75/5 , 2.5/3],frameon=True)
    panel2 = plt.axes([1.5/5, 1.45/3, 0.1/5, 0.2/3], frameon=True) 

    heatmap = panel1.imshow(normalized_matrix, aspect='auto', cmap=viridian, interpolation='nearest')
    tick_positions = np.arange(len(FPKM_indices))
    panel1.set_xticks(tick_positions)
    tick_labels = [f'{3 * i}' if i % 2 == 0 else '' for i in range(len(FPKM_indices))]
    panel1.set_xticklabels(tick_labels, fontsize=8)
    panel1.set_xlabel('CT',fontsize=8)

    specific_genes = ['Clock', 'Insig2', 'Nr1d1', 'Dbp', 'Per3', 'Per1', 'Per2', 'Cry1', 'Arntl']
    gene_positions = [5 ,441, 1890, 3225, 3696, 3990, 4569, 5098, 6300] 
    panel1.set_yticks(gene_positions)
    panel1.set_yticklabels(specific_genes, fontsize=8,)

    panel1.tick_params(axis='x', length=1.5, pad=2) 
    panel1.tick_params(axis='y', length=1.5, pad=2) 


    cbar = plt.colorbar(heatmap, cax=panel2, orientation='vertical')
    cbar.ax.set_yticklabels(['Min', 'Max'], fontsize=8)

    plt.savefig("/Users/hemap/Downloads/Hema_Prasanna_Assignment_Week6_corrected.png", dpi=600)
    plt.show()

if __name__ == "__main__":
    main()
