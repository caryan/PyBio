from Bio import SeqIO

import time
import numpy as np

import sympy as sy

from bokeh.browserlib import view
from bokeh.document import Document
from bokeh.session import Session
from bokeh.objects import Plot, HoverTool, PanTool, WheelZoomTool, Range1d, FactorRange, ColumnDataSource
from bokeh.glyphs import Rect, Line
from bokeh.widgets import Select, HBox, VBox, Dialog, Slider, TextInput
from requests.exceptions import ConnectionError

from collections import OrderedDict

import glob

fileName = "more_sequences/oprD_fasta.fa"

document = Document()
session = Session()
session.use_doc('porin_map_server')
session.load_document(document)

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

porin = fileName.split("/")[1].split("_")[0]

plot = Plot(title=porin + " POA1 Differences", 
                x_range=Range1d(start=-0.5, end=pa.shape[1]), y_range=FactorRange(factors=strains),
                plot_width=1000, plot_height=3000)

plot.add_glyph(source, Rect(x="pos", y="strain", width=0.95, height=0.95, fill_color="colour", line_color=None))

plot.add_tools(PanTool(), WheelZoomTool(), HoverTool(tooltips = OrderedDict([
    ("strain", "@strain"),
    ("site", "@pos"),
    ("protien", "@curProtien"),
    ("PAO1", "@refProtien"),
    ("tot. diff.", "@diffCounts")
    ])
    ) )


document.add(plot)
session.store_document(document)

if __name__ == "__main__":
    link = session.object_link(document.context)
    print("Please visit %s to see the plots" % link)
    view (link)
    print("\npress ctrl-C to exit")
    try:
        while True:
            session.load_document(document)
            time.sleep(5)
    except KeyboardInterrupt:
        print()
    except ConnectionError:
        print("Connection to bokeh-server was terminated")
