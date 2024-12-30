

import argparse
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Plot cell positions and densities.')
    parser.add_argument('-p', '--position', required=True, help='Path to the position.tsv file')
    parser.add_argument('-c', '--celltype', required=True, help='Path to the celltype.tsv file')
    parser.add_argument('-o', '--output', required=True, help='Output PNG file path')
    args = parser.parse_args()

    position_data, celltype_data = read_data(args.position, args.celltype)
    data = merge_data(position_data, celltype_data)

    plot_data(data, args.output)

def read_data(position_file, celltype_file):
    position_data = {}
    with open(position_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            position_data[parts[0]] = (float(parts[1]), float(parts[2]))

    celltype_data = {}
    with open(celltype_file, 'r') as file:
        next(file) 
        for line in file:
            parts = line.strip().split()
            celltype_data[parts[2]] = parts[1]

    return position_data, celltype_data

def merge_data(position_data, celltype_data):
    data = {}
    for barcode, position in position_data.items():
        if barcode in celltype_data:
            cell_type = celltype_data[barcode]
            data[barcode] = {'position': position, 'celltype': cell_type}
    return data

def plot_cells(panel, data, colors):
    celltype_positions = {}
    for barcode, info in data.items():
        cell_type = info['celltype']
        positions = celltype_positions.get(cell_type, [])
        positions.append(info['position'])
        celltype_positions[cell_type] = positions
        

        color = colors.get(cell_type, 'gray') 
        panel.plot(info['position'][0], info['position'][1], 'o', 
                   markeredgecolor='black', markerfacecolor=color, markersize=4, mew=0.5)
    

    for cell_type, positions in celltype_positions.items():
        positions = np.array(positions)
        median_x = np.median(positions[:, 0])
        median_y = np.median(positions[:, 1])

    #outline effect for labels
        panel.text(median_x, median_y, cell_type, fontsize=9, ha='center', va='center',
                   color='white', weight='bold')

        panel.text(median_x + 0.05, median_y, cell_type, fontsize=8, ha='center', va='center',
                   color='black', weight='bold')

def calculate_density(data, minimum_distance=8/72, panel_width=1.5, panel_height=1.5):
    points = [info['position'] for info in data.values()]
    xrange = max(p[0] for p in points) - min(p[0] for p in points)
    yrange = max(p[1] for p in points) - min(p[1] for p in points)
    density = []
    for i, (x1, y1) in enumerate(points):
        count = sum(1 for j, (x2, y2) in enumerate(points) if i != j and np.sqrt(
            ((x1 - x2) / xrange * panel_width)**2 + ((y1 - y2) / yrange * panel_height)**2
        ) < minimum_distance)
        density.append(min(count, 100))
    return density

def plot_density(panel, data, density):
    viridis_cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(vmin=0, vmax=100)
    points = [info['position'] for info in data.values()]
    colors = viridis_cmap(norm(density))
    panel.scatter([p[0] for p in points], [p[1] for p in points], 
                  color=colors, edgecolor='none', s=4**2)

def plot_data(data, output_file):
    figure_width, figure_height = 5, 3
    plt.figure(figsize=(figure_width, figure_height))
    plt.style.use('BME163')

    panel_width, panel_height = 1.5, 1.5
    relative_panel_width = panel_width / figure_width
    relative_panel_height = panel_height / figure_height

    panel1 = plt.axes([0.5/figure_width,0.5/figure_height,relative_panel_width, relative_panel_height])
    panel2 = plt.axes([2.5/figure_width,0.5/figure_height, relative_panel_width, relative_panel_height])
    panel3 = plt.axes([4.0/figure_width,1.1/figure_height,0.1/figure_width,0.3/figure_height])
    
    panel1.set_xlabel('tSNE 2')
    panel1.set_ylabel('tSNE 1')
    panel2.set_xlabel('tSNE 2')
    panel2.set_ylabel('tSNE 1')

    colors = {'monocyte': 'red', 'neuron': 'blue', 'glia': 'green', 'tCell': 'purple', 'bCell': 'cyan'}

    plot_cells(panel1, data, colors)
    density = calculate_density(data)
    scatter = plot_density(panel2, data, density)

   
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(density), vmax=max(density)))
    sm.set_array([])  
    cbar = plt.colorbar(sm, cax=panel3)
    cbar.set_ticks([min(density), max(density)])
    cbar.set_ticklabels(['Min', 'Max'])

    panel2.text(0.05, 0.05, 'Density', transform=panel2.transAxes, fontsize=8, ha='left', va='bottom')
    plt.savefig(output_file, dpi=600)

if __name__ == '__main__':
    main()
