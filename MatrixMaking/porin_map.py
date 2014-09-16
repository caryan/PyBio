from Bio import SeqIO

import numpy as np

from bokeh.plotting import *
from bokeh.objects import HoverTool
from collections import OrderedDict

import glob

fileNames = glob.glob("more_sequences/op*")


for fileName in fileNames:

    recs = [r for r in SeqIO.parse(fileName, "fasta")]

    recs[0].description = "PAO1_" + recs[0].description

    goodRecs = []

    goodRecs.append(recs[0])

    strains = [r.description.split("_")[0] for r in recs]
    strainsDone = [strains[0]]

    for r,strain in zip(recs[1:], strains[1:]):
        if strain in strainsDone:
            print("Duplicate strain {}... dropping".format(strain))
            if "*" in str(r):
                print("Uhoh! It had a * too!")
            continue
        strainsDone.append(strain)
        goodRecs.append(r)

    recs = goodRecs
    strains = strainsDone

    numStrains = len(strains)
    numProtiens = len(recs[0].seq)

    pa = np.zeros((numStrains, numProtiens), dtype=np.uint8)

    refSeqList = list(str(recs[0].seq))

    for ct,r in enumerate(recs):
        for ct2, c in enumerate(str(r.seq)):
            if c == refSeqList[ct2]:
                pa[ct,ct2] = 0
            else:
                if c == "*":
                    pa[ct,ct2] = 2
                else:
                    pa[ct,ct2] = 1

    colours = ["#1D8281", "#FBD258", "#DB634F"]

    diffCounts = np.sum(pa>0, axis=0)

    source = ColumnDataSource(data = {"strain":np.tile(strains, (numProtiens, 1)).flatten(order="F"),
                                        "pos":np.tile(np.arange(numProtiens), (numStrains,1)).flatten(order="C"),
                                        "refProtien":np.tile(refSeqList, (pa.shape[0],1)).flatten(order="C"),
                                        "curProtien":[c for r in recs for c in str(r.seq)],
                                        "diffCounts":np.tile(diffCounts, (numStrains, 1)).flatten(order="C"),
                                        "colour":[colours[_] for _ in pa.flatten()]})


    reset_output()

    porin = fileName.split("/")[1].split("_")[0]
    output_file(porin + ".html", title=porin)

    figure()

    rect("pos", "strain", 0.95, 0.95, source=source,
         x_axis_location="above",
         y_range=strains,
         x_range=[-0.5, pa.shape[1]],
         color="colour", line_color=None,
         title=porin + " POA1 Differences",
         plot_width=1500, plot_height=6000,
         tools="resize,previewsave,hover")


    #   - remove the axis and grid lines
    #   - remove the major ticks
    #   - make the tick labels smaller
    #   - set the x-axis orientation to vertical, or angled
    grid().grid_line_color = None
    axis().axis_line_color = None
    axis().major_tick_line_color = None
    axis().major_label_text_font_size = "5pt"
    axis().major_label_standoff = 0
    xaxis().major_label_orientation = np.pi/3

    hover = curplot().select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("strain", "@strain"),
        ("site", "@pos"),
        ("protien", "@curProtien"),
        ("PAO1", "@refProtien"),
        ("tot. diff.", "@diffCounts"),
    ])


    save()
