
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate sequence logos for splice sites.')
    parser.add_argument('-s', '--splice_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-A', '--A_image', required=True)
    parser.add_argument('-T', '--T_image', required=True)
    parser.add_argument('-G', '--G_image', required=True)
    parser.add_argument('-C', '--C_image', required=True)
    return parser.parse_args()

def read_sequences(file_path):
    sequences = {'5SS': [], '3SS': []}
    with open(file_path, 'r') as file:
        seq = ''
        category = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences[category].append(seq)
                category = '5SS' if "5'" in line else '3SS'
                seq = ''
            else:
                seq += line
        if seq:
            sequences[category].append(seq)
    return sequences

def calculate_frequencies(sequences):
    length = 20 
    nucleotides = 'ATGC'
    freqs = {'5SS': np.zeros((length, 4)), '3SS': np.zeros((length, 4))}
    idx = {n: i for i, n in enumerate(nucleotides)}
    
    for key, seqs in sequences.items():
        for seq in seqs:
            for i, nuc in enumerate(seq):
                freqs[key][i, idx[nuc]] += 1
        freqs[key] /= len(seqs)

    return freqs

def plot_sequence_logos(freqs, args):
    figureWidth, figureHeight = 5, 2
    plt.figure(figsize=(figureWidth, figureHeight))
    panelWidth, panelHeight = 1.5, 0.5
    panel1 = plt.axes([0.5 / figureWidth, 0.3, panelWidth / figureWidth, panelHeight / figureHeight], frameon=True)
    panel2 = plt.axes([2.2 / figureWidth, 0.3, panelWidth / figureWidth, panelHeight / figureHeight], frameon=True)
    
    panel1.tick_params(axis='y', labelsize=4)


    for panel, key in zip([panel1, panel2], ['5SS', '3SS']):
        for pos in range(20):
            base_freqs = [(freqs[key][pos, i], i) for i in range(4)]
            base_freqs.sort(reverse=True, key=lambda x: x[0])  
            total_info = -np.sum([bf[0] * np.log2(bf[0] + 1e-10) for bf in base_freqs])
            info_norm = 2 - total_info
            bottom = 0
            for freq, idx in reversed(base_freqs[:3]): 
                height = freq * info_norm
                img = mpimg.imread([args.A_image, args.T_image, args.G_image, args.C_image][idx])
                panel.imshow(img, extent=(pos-10, pos-9, bottom, bottom+height), aspect='auto')
                bottom += height
        
        panel.set_xlim(-10, 10)
        panel.set_ylim(0, 2)
        panel.set_xticks([-10, -5, 0, 5, 10])
        panel.set_xticklabels(['-10', '-5', '0', '5', '10'],fontsize=7.4)
        panel.set_xlabel('Distance to\nSplice Site',fontsize=7.4)
        panel.set_title(f"{key.replace('SS', '\'SS')}",fontsize=7.4)

       
        panel.axvline(0, color='black', linewidth=0.5)

   
    panel1.set_ylabel('Bits',fontsize=7.4)
    panel1.set_yticks([0, 1, 2],)
    
   
    panel2.set_yticklabels([])
    panel2.set_yticks([])

    plt.tight_layout()
    plt.savefig(args.output_file, dpi=300)

def main():
    args = parse_arguments()
    sequences = read_sequences(args.splice_file)
    frequencies = calculate_frequencies(sequences)
    plot_sequence_logos(frequencies, args)

if __name__ == '__main__':
    main()