
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse

# Argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Plot gene locus data.")
    parser.add_argument('-p5', '--psl5', required=True, help='Path to BME163_Input_Data_5.psl file')
    parser.add_argument('-p6', '--psl6', required=True, help='Path to BME163_Input_Data_6.psl file')
    parser.add_argument('-g', '--gtf', required=True, help='Path to gencode.vM12.annotation.gtf file')
    parser.add_argument('-c', '--coordinates', required=True, help='Coordinates in the format chr:start-end')
    parser.add_argument('-o', '--output', required=True, help='Output file name')
    return parser.parse_args()

# Panel configuration
def panel_edit(panel):
    panel.tick_params(bottom=False, labelbottom=False,
                      left=False, labelleft=False,
                      right=False, labelright=False,
                      top=False, labeltop=False)
    panel.axes.get_yaxis().set_ticks([])
    panel.axes.get_xaxis().set_ticks([])
    panel.set_xlim(45232000, 45241000)

# GTF file parsing
def parse_gtf(gtf_file, chromosome, start, end):
    transcripts = []
    gtfdict = {}
    for line in open(gtf_file):
        if line.startswith("#"):
            continue
        split_list = line.strip().split('\t')
        if split_list[0] == chromosome and int(split_list[3]) >= start and int(split_list[4]) <= end:
            feature_type = split_list[2]
            if feature_type == "exon" or feature_type == "CDS":
                transcript = split_list[8].split('transcript_id "')[1].split('"')[0]
                if transcript not in gtfdict:
                    gtfdict[transcript] = []
                gtfdict[transcript].append([split_list[0], int(split_list[3]), int(split_list[4]), feature_type])
    
    for transcript, parts in gtfdict.items():
        starts = []
        ends = []
        block_starts = []
        block_widths = []
        types = []
        for part in parts:
            starts.append(part[1])
            ends.append(part[2])
            block_starts.append(part[1])
            block_widths.append(part[2] - part[1])
            types.append(part[3])
        transcripts.append([parts[0][0], min(starts), max(ends), block_starts, block_widths, False, types])
    return sorted(transcripts, key=lambda x: x[2])

# PSL file parsing
def parse_psl(psl_file, chromosome, start, end):
    reads = []
    with open(psl_file) as file:
        for line in file:
            if line.startswith("start"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 21:
                continue
            try:
                block_sizes = list(map(int, fields[18].rstrip(',').split(',')))
                block_starts = list(map(int, fields[20].rstrip(',').split(',')))
                read = [fields[13], int(fields[15]), int(fields[16]), block_starts, block_sizes, False]
                if fields[13] == chromosome and (start < int(fields[15]) < end or start < int(fields[16]) < end):
                    reads.append(read)
            except ValueError:
                continue
    return reads

# function to stack reads
def stack_reads(reads):
    stacked_reads = []
    for read in reads:
        placed = False
        for stack in stacked_reads:
            if read[1] > stack[-1][2]: 
                stack.append(read)
                placed = True
                break
        if not placed:
            stacked_reads.append([read])
    return stacked_reads

# flot transcripts
def plot_transcripts(panel, transcripts):
    y_pos = 0
    y_increment = 1  
    for transcript in transcripts:
       
        rectangle = mplpatches.Rectangle((transcript[1], y_pos + 0.23), transcript[2] - transcript[1], 0.05, facecolor='grey', edgecolor='black', linewidth=0.25)
        panel.add_patch(rectangle)
        
        for i, block_start in enumerate(transcript[3]):
            height = 0.5 if transcript[6][i] == "CDS" else 0.25
            rectangle = mplpatches.Rectangle((block_start, y_pos), transcript[4][i], height, facecolor='grey', edgecolor='black', linewidth=0.25)
            panel.add_patch(rectangle)
        y_pos += y_increment
    panel.set_ylim(0, y_pos + 1)

# plot reads for the middle panel (sorted by end)
def plot_reads_sorted_by_end(panel, reads, color):
    y_pos = 0
    y_increment = 1  
    stacked_reads = stack_reads(sorted(reads, key=lambda x: x[2]))
    for stack in stacked_reads:
        for read in stack:
           
            rectangle = mplpatches.Rectangle((read[1], y_pos + 0.18), read[2] - read[1], 0.05, facecolor=color, edgecolor=color, linewidth=0)
            panel.add_patch(rectangle)
            
            for block_start, block_width in zip(read[3], read[4]):
                rectangle = mplpatches.Rectangle((block_start, y_pos), block_width, 0.5, facecolor=color, edgecolor=color, linewidth=0)
                panel.add_patch(rectangle)
            y_pos += y_increment
    panel.set_ylim(0, y_pos + 1)

# plot reads for the bottom panel (sorted by start)
def plot_reads_sorted_by_start(panel, reads, color):
    y_pos = 0
    y_increment = 1  
    stacked_reads = stack_reads(sorted(reads, key=lambda x: x[1]))
    for stack in stacked_reads:
        for read in stack:
            
            rectangle = mplpatches.Rectangle((read[1], y_pos + 0.18), read[2] - read[1], 0.05, facecolor=color, edgecolor=color, linewidth=0)
            panel.add_patch(rectangle)
            
            for block_start, block_width in zip(read[3], read[4]):
                rectangle = mplpatches.Rectangle((block_start, y_pos), block_width, 0.5, facecolor=color, edgecolor=color, linewidth=0)
                panel.add_patch(rectangle)
            y_pos += y_increment
    panel.set_ylim(0, y_pos + 1)

# plot reads for the condensed panel (sorted by start)
def plot_reads_sorted_by_start_condensed(panel, reads, color):
    y_pos = 0
    y_increment = 0.2  
    stacked_reads = stack_reads(sorted(reads, key=lambda x: x[1]))
    for stack in stacked_reads:
        for read in stack:
            
            rectangle = mplpatches.Rectangle((read[1], y_pos + 0.18), read[2] - read[1], 0.05, facecolor=color, edgecolor=color, linewidth=0)
            panel.add_patch(rectangle)
            
            for block_start, block_width in zip(read[3], read[4]):
                rectangle = mplpatches.Rectangle((block_start, y_pos), block_width, 0.5, facecolor=color, edgecolor=color, linewidth=0)
                panel.add_patch(rectangle)
            y_pos += y_increment
    panel.set_ylim(0, y_pos * 1.10)  

# histogram for the bottom panel
def plot_histogram(panel, reads, start, end):
    coverage = np.zeros(end - start)
    for read in reads:
        for block_start, block_width in zip(read[3], read[4]):
            coverage[block_start - start:block_start - start + block_width] += 1
    panel.bar(range(start, end), coverage, width=1.0, color=(88/255, 85/255, 120/255), edgecolor='none')
    panel.set_ylim(max(coverage) + 1, 0)  
    panel.set_xlim(start, end)

   
    panel.set_yticks([])
    panel.set_xticks([])

# main function
def main():
    args = parse_args()
    chromosome, region = args.coordinates.split(':')
    start, end = map(int, region.split('-'))

    # parse files
    gtf_data = parse_gtf(args.gtf, chromosome, start, end)
    psl5_data = parse_psl(args.psl5, chromosome, start, end)
    psl6_data = parse_psl(args.psl6, chromosome, start, end)

    # filter reads
    filtered_psl5_data = [read for read in psl5_data if read[0] == chromosome and (start < read[1] < end or start < read[2] < end)]
    filtered_psl6_data = [read for read in psl6_data if read[0] == chromosome and (start < read[1] < end or start < read[2] < end)]

    # create figure
    plt.figure(figsize=(5, 6))

    # create panels
    panel0 = plt.axes([0.1/5, 3.9/6, 4/5, 1.5/6])
    panel1 = plt.axes([0.1/5, 2.2/6, 4/5, 1.5/6])
    panel2 = plt.axes([0.1/5, 0.5/6, 4/5, 1.5/6])
    panel3 = plt.axes([0.1/5, 0.1/6, 4/5, 0.4/6])

    # edit panels
    panel_edit(panel0)
    panel_edit(panel1)
    panel_edit(panel2)
    panel_edit(panel3)

    # plot data
    plot_transcripts(panel0, gtf_data)  # grey plot
    plot_reads_sorted_by_end(panel1, filtered_psl5_data, (230/255, 87/255, 43/255))  # orange plot
    plot_reads_sorted_by_start_condensed(panel2, filtered_psl6_data, (88/255, 85/255, 120/255))  # blue plot
    plot_histogram(panel3, filtered_psl6_data, start, end)  # coverage histogram

    
    plt.savefig(args.output, dpi=2400)
    plt.close()

if __name__ == "__main__":
    main()


