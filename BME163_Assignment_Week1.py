import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse

# setup the CL argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--outFile', type=str, required=True, help='Output file name')
args = parser.parse_args()

# set BME163 stylesheet
plt.style.use('BME163')

# create figure
figure = plt.figure(figsize=(5, 2)) 

# define the left panel for circles
left_panel = figure.add_axes([0.04, 0.1, 1.0/5.0, 1/2])  
left_panel.set_xlim(0, 15)
left_panel.set_ylim(0, 16)
left_panel.set_aspect('equal', adjustable='datalim')  
left_panel.axis('off')  

# define the colors for the circles 
circle_colors = [
    (225/255, 13/255, 50/255),  # RB1
    (242/255, 50/255, 54/255),  # RB2
    (239/255, 99/255, 59/255),  # RB3
    (244/255, 138/255, 30/255), # RB4
    (248/255, 177/255, 61/255), # RB5
    (143/255, 138/255, 86/255), # RB6
    (32/255, 100/255, 113/255), # RB7
    (42/255, 88/255, 132/255),  # RB8
    (56/255, 66/255, 156/255),  # RB9
    (84/255, 60/255, 135/255),  # RB10
    (110/255, 57/255, 115/255), # RB11
    (155/255, 42/255, 90/255),  # RB12
]

# plot the circles in the left panel
for i in range(len(circle_colors)):
    circle = mplpatches.Circle((2 + 1*i, 8), 1, facecolor='none', edgecolor=circle_colors[i], linewidth=1)
    left_panel.add_patch(circle)

# define the right panel for the heatmap
right_panel = figure.add_axes([0.35,0.1, 2/5, 1/2])  
right_panel.set_xlim(0, 2)
right_panel.set_ylim(0, 1)
right_panel.axis('off')  

# create horizontal gradients
num_rectangles = 300  
colors_viridis = plt.get_cmap('viridis')(np.linspace(0, 1, num_rectangles))
colors_plasma = plt.get_cmap('plasma')(np.linspace(0, 1, num_rectangles))

# draw gradient using multiple thin rectangles
rect_height = 0.5 
for i in range(num_rectangles):
    rect_x = i / num_rectangles * 2 
    color_v = colors_viridis[i]
    color_p = colors_plasma[i]
    viridis_rect = mplpatches.Rectangle((rect_x, 0), 2/num_rectangles, rect_height, color=color_v, edgecolor='black', linewidth=1.5)
    plasma_rect = mplpatches.Rectangle((rect_x, 0.5), 2/num_rectangles, rect_height, color=color_p, edgecolor='black', linewidth=1.5)
    right_panel.add_patch(viridis_rect)
    right_panel.add_patch(plasma_rect)

# outline panels and rectangles black
left_outline = mplpatches.Rectangle((0, 0), 1, 1, transform=left_panel.transAxes, fill=False, edgecolor='black', linewidth=1.5)
right_outline = mplpatches.Rectangle((0, 0), 1, 1, transform=right_panel.transAxes, fill=False, edgecolor='black', linewidth=1.5)
left_panel.add_patch(left_outline)
right_panel.add_patch(right_outline)


right_panel.axhline(y=0.49, color='black', linewidth=0.7)

#save figure
figure.savefig(args.outFile, dpi=600)
