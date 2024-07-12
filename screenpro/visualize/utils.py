import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd

# variables
almost_black = '#111111'
dark2 = ['#1b9e77',
         '#d95f02',
         '#7570b3',
         '#e7298a',
         '#66a61e',
         '#e6ab02',
         '#a6761d',
         '#666666']
blue_yellow = matplotlib.colors.LinearSegmentedColormap.from_list(
    'BuYl', [(0, '#ffff00'), (.49, '#000000'), (.51, '#000000'), (1, '#0000ff')])
blue_yellow.set_bad('#999999', 1)
yellow_blue = matplotlib.colors.LinearSegmentedColormap.from_list(
    'YlBu', [(0, '#0000ff'), (.49, '#000000'), (.51, '#000000'), (1, '#ffff00')])
yellow_blue.set_bad('#999999', 1)

# plt.rcParams['font.sans-serif'] = [
#     'Helvetica', 'Arial', 'Verdana', 'Bitstream Vera Sans'
# ]
# plt.rcParams['font.size'] = 8
# plt.rcParams['font.weight'] = 'regular'
# plt.rcParams['text.color'] = almost_black
#
# axisLineWidth = .5
# plt.rcParams['axes.linewidth'] = axisLineWidth
# plt.rcParams['lines.linewidth'] = 1.5
#
# plt.rcParams['axes.facecolor'] = 'white'
# plt.rcParams['axes.edgecolor'] = almost_black
# plt.rcParams['axes.labelcolor'] = almost_black
# # plt.rcParams['axes.color_cycle'] = dark2_all
#
# plt.rcParams['patch.edgecolor'] = 'none'
# plt.rcParams['patch.linewidth'] = .25
# # plt.rcParams['patch.facecolor'] = dark2_all[0]
#
# plt.rcParams['savefig.dpi'] = 1000
# plt.rcParams['savefig.format'] = 'svg'
#
# plt.rcParams['legend.frameon'] = False
# plt.rcParams['legend.handletextpad'] = .25
# plt.rcParams['legend.fontsize'] = 8
# plt.rcParams['legend.numpoints'] = 1
# plt.rcParams['legend.scatterpoints'] = 1
#
# plt.rcParams['ytick.direction'] = 'out'
# plt.rcParams['ytick.color'] = almost_black
# plt.rcParams['ytick.major.width'] = axisLineWidth
# plt.rcParams['xtick.direction'] = 'out'
# plt.rcParams['xtick.color'] = almost_black
# plt.rcParams['xtick.major.width'] = axisLineWidth


def cleanAxes(ax, top=False, right=False, bottom=True, left=True):
    ax.grid('off')
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)

    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')

    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()

