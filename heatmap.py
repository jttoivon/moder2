#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
from matplotlib import ticker
import sys
import getopt
import os

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def readarray(lines):
    result=[]
    for line in lines:
        line = line.rstrip('\n')
        tmp=line.split("\t")
        result.append(tmp)
    return np.array(result)

def float_to_string(f):
    if -0.0002 <= f <= 0.0:
        return "-"
    else:
        return '%.0f' % (f*1000)

def myfloat(s):
    "Converts string to float. Allows use of unicode minus sign (u'\\u2212')."
    r = "".join(["-" if x == '\u2212' else x for x in s])
    return float(r)

def make_heatmap(data, xlabels, ylabels, fmt, title="", outputfile="", fontsize=32.0, cell_labels=False):
#    plt.style.use('ggplot')
    
    linewidth=1.0
#    fontsize=32.0
    labelfontsize=fontsize*0.6
    tickfontsize=fontsize*0.6
#    rcParams['axes.titlepad'] = 20     # Padding between the title and the plot, requires recent version of matplotlib
    
    width=data.shape[1]
    height=data.shape[0]  # number of orientations
    fig = plt.figure()
    ax = plt.subplot(111)

    ##################
    #
    # Numbers in cells
    #
    ##################
    
    if cell_labels:
        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
        #        plt.text(x + 0.5, y + 0.5, '%.4f' % data[y, x],
                plt.text(x, y, float_to_string(data[y, x]),
                         horizontalalignment='center',
                         verticalalignment='center',
                         )

                    
    cmap = plt.get_cmap('YlOrRd')
    subcmap = truncate_colormap(cmap, 0.0, 0.8)
    subcmap.set_under(color='white', alpha=None)
    
#    rcParams['lines.solid_joinstyle'] = "round"
    plt.imshow(data, vmin=0.0, cmap=subcmap, interpolation='nearest', aspect='equal')
    divider = make_axes_locatable(ax)
    if title:
        plt.title(title, fontsize=fontsize)

    #########################
    #
    # Set x and y tick labels
    #
    #########################

    plt.yticks(fontsize=tickfontsize)
    ax.yaxis.set_ticks(np.arange(0,height,1))
    if isinstance(ylabels[0], int):
        ax.set_yticklabels(["%i" % i for i in ylabels])
    else:
        ax.set_yticklabels(["%s" % s for s in ylabels])
        
    plt.xticks(fontsize=tickfontsize)
    ax.xaxis.set_ticks(np.arange(0,width,1))
    if isinstance(xlabels[0], int):
        ax.set_xticklabels(["%i" % i for i in xlabels])
    else:
        ax.set_xticklabels(["%s" % s for s in xlabels])
        
    ###############################
    #
    # No major nor minor tick marks
    #
    ###############################

    plt.tick_params(axis='both', which='both', bottom=False, top=False, right=False, left=False)
    
    ###################
    #
    # Make grid
    #
    ###################

    # These lines will create grid on minor ticks, that is, make cell borders
    ax.set_xticks(np.arange(-0.5, width, 1), minor=True);
    ax.set_yticks(np.arange(-0.5, height, 1), minor=True);
    # Gridlines based on minor ticks
    ax.grid(which='minor', color='black', linestyle='-', linewidth=1, solid_joinstyle='round')

    ##################
    #
    # Make colorbar
    #
    ##################

    cax = divider.append_axes("right", size="5%", pad=0.05)
    try:
        cb=plt.colorbar(cax=cax)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        ##cb.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        cb.update_ticks()
        temp=cax.get_yticklabels()
        for i,t in enumerate(temp):
            #print temp[i].get_text()
#            temp[i].set_text("%.0f" % (float(temp[i].get_text())*1000))   # Multiply values in colorbar by 1000
            temp[i] = "%.0f" % (float(temp[i].get_text())*1000)   # Multiply values in colorbar by 1000
            #print temp[i].get_text()
        
        #cb.update_ticks()
        cax.set_yticklabels(temp, fontsize=tickfontsize)
    except ValueError:   # If labels contain unicode minus, then something went wrong and better not show colorbar
        cb.remove()
        #print "Unicode error!"
        pass
#    cax.yaxis.set_tick_params(labelright=False)   # No tick labels in colorbar

    if outputfile:
        plt.savefig(outputfile, format=fmt, bbox_inches="tight")
    else:
        plt.show()


        
def main():        
    endings = ["png", "pdf", "ps", "eps", "svg"]
    prog=os.path.basename(sys.argv[0])
    title=""



    usage="""Usage:
    %s [ options ] cobfile [ imagefile ]

    -t, --title arg\t\tUses the given title in the figure


    Visualizes a cob file (.cob). The extension of the imagefile
    (%s) chooses the format of the image.

    """ % (prog, ", ".join(endings))



    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'ht:', ["help", "title="])
    except getopt.GetoptError as e:
        print(e)
        sys.stderr.write(usage)
        sys.exit(1)
    optdict = dict(optlist)
    args = [sys.argv[0]]+ args

    for o, arg in optlist:
        if o in ["-h", "--help"]:
            print(usage)
            sys.exit(0)
        elif o in ["-t", "--title"]:
            title=arg


    if len(args) < 2:
        sys.stderr.write("Give at least one parameter\n")
        sys.stderr.write(usage)
        sys.exit(1)




    input=args[1]
    try:
        outputfile=args[2]
        parts = outputfile.split(".")
        ending=parts[-1]
        if len(parts) >= 2 and ending in endings:
            fmt=ending
        else:
            sys.stderr.write("The extension of the outputfile should be one of the following: %s\n" % (", ".join(endings)))
            sys.stderr.write("Exiting!\n")
            sys.exit(1)
    except IndexError:
        fmt=""
        outputfile=""

    try:
        with open(input, "r") as f:
            lines=f.readlines()
    except IOError:
        sys.stderr.write("Could not read file %s\n" % input)
        sys.exit(1)

    cob = readarray(lines)
    drange = list(map(int, cob[0,1:]))

    #print drange
    #print cob

    #data = np.random.rand(3, 9)
    data=cob[1:,1:].astype(float)
    #print data.shape
    if data.shape[0] in [1,2]:
        use_rna=True
    else:
        use_rna=False

    if use_rna:
        orients=["HT", "TH"]
    else:
        orients=["HT", "HH", "TT", "TH"]

    number_of_orientations=data.shape[0]
    orients = orients[0:number_of_orientations]

    vfunc = np.vectorize(lambda x: x if x > 0.0 else -0.0002)   # Modify zero values to -0.0002 in order for the colormap to work better.
    data=vfunc(data)
    make_heatmap(data, drange, orients, fmt, title=title, outputfile=outputfile, cell_labels=True)

if __name__ == '__main__':
    main()
