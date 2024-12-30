'''
import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib.patches as mplpatches

# Function to calculate log2(values + 1)
def log_transform(values):
    return np.log2(np.array(values) + 1)

# Argument parser for input and output files
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputFile', type=str, help='Input file path')
parser.add_argument('-o', '--outputFile', type=str, help='Output file path')
args = parser.parse_args()

# Set the style for the plot
plt.style.use('BME163')

# Read the data from the input file and transform
data = np.loadtxt(args.inputFile, usecols=(1, 2), skiprows=1)
x_values = log_transform(data[:, 0])
y_values = log_transform(data[:, 1])

# Define the colors
iBlue = (88/255, 85/255, 120/255)
Grey = 'grey'
iGreen = (120/255, 172/255, 145/255)

# Figure and panel dimensions
figureWidth = 3
figureHeight = 3
panelWidth = 1.5
panelHeight = 1.5
leftPanelWidth = 0.25
topPanelHeight = 0.25

# set up the figure and axes
figure = plt.figure(figsize=(figureWidth, figureHeight))

main_panel = figure.add_axes([leftPanelWidth/figureWidth, 0.1, panelWidth/figureWidth, panelHeight/figureHeight])
left_panel = figure.add_axes([0.05, 0.1, leftPanelWidth/figureWidth, panelHeight/figureHeight], sharey=main_panel)
top_panel = figure.add_axes([leftPanelWidth/figureWidth, 0.85, panelWidth/figureWidth, topPanelHeight/figureHeight], sharex=main_panel)

# plot the scatterplot
main_panel.scatter(x_values, y_values, s=10, color=iBlue, alpha=0.1, edgecolor='none')



# save the figure
plt.savefig(args.outputFile, dpi=600)
'''



import matplotlib.pyplot as plt
import numpy as np
import argparse
import matplotlib.patches as mplpatches

# Function to calculate log2(values + 1)
def log_transform(values):
    return np.log2(np.array(values) + 1)

# Argument parser for input and output files
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputFile', type=str, help='Input file path')
parser.add_argument('-o', '--outputFile', type=str, help='Output file path')
args = parser.parse_args()

# Set the style for the plot
plt.style.use('BME163')

# Read the data from the input file and transform
data = np.loadtxt(args.inputFile, usecols=(1, 2), skiprows=1)
x_values = log_transform(data[:, 0])
y_values = log_transform(data[:, 1])

# Define the colors
iBlue = (88/255, 85/255, 120/255)
Grey = 'grey'
iGreen = (120/255, 172/255, 145/255)

# Create the figure and axes
figure = plt.figure(figsize=(3, 3))
main_panel = figure.add_axes([0.2, 0.2, 1.5/3, 1.5/3])  # Main scatter plot

# Histograms (transformed using log)
x_hist, x_bins = np.histogram(x_values, bins=50)
y_hist, y_bins = np.histogram(y_values, bins=50)
x_hist = log_transform(x_hist)
y_hist = log_transform(y_hist)
x_hist_norm = x_hist / max(x_hist) * (1.5/3)
y_hist_norm = y_hist / max(y_hist) * (1.5/3)

# Histogram axes
left_panel = figure.add_axes([0.1, 0.2, 0.08, 1.5/3], sharey=main_panel)  # Left histogram
top_panel = figure.add_axes([0.2, 0.72, 1.5/3, 0.08], sharex=main_panel)  # Top histogram


# plot the scatterplot
main_panel.scatter(x_values, y_values, s=10, color=iBlue, alpha=0.1, edgecolor='none')

# Plot histograms using bar patches
for i in range(len(x_bins)-1):
    height=x_hist_norm[i] * 0.5
    top_panel.add_patch(mplpatches.Rectangle((x_bins[i], 0), x_bins[i+1] - x_bins[i], x_hist_norm[i], height,
                                facecolor=iGreen, edgecolor= 'black',linewidth=0.2))
   #top_panel.add_patch(plt.Rectangle((x_bins[i], 0), x_bins[i+1] - x_bins[i], x_hist[i], facecolor=Grey, edgecolor=Grey))
for i in range(len(y_bins)-1):
    left_panel.add_patch(mplpatches.Rectangle((0, y_bins[i]), y_hist_norm[i], y_bins[i+1] - y_bins[i], height,
                                facecolor=Grey, edgecolor='black',linewidth=0.2))

# Set limits for the panels
main_panel.set_xlim(0, max(x_values))
main_panel.set_ylim(0, max(y_values))
top_panel.set_xlim(main_panel.get_xlim())
left_panel.set_ylim(main_panel.get_ylim())

# save the figure
plt.savefig(args.outputFile, dpi=600)
