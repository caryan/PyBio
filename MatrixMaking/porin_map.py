from Bio import SeqIO

import time
import numpy as np
import pandas as pd

from bokeh.browserlib import view
from bokeh.document import Document
from bokeh.session import Session
from bokeh.objects import Plot, HoverTool, PanTool, WheelZoomTool, ResetTool, \
                            Range1d, FactorRange, ColumnDataSource
from bokeh.glyphs import Rect, Line
from bokeh.widgets import Select, HBox, VBox, TableColumn, HandsonTable
from requests.exceptions import ConnectionError

from collections import OrderedDict

import glob

porins = [_.split("/")[1].split("_")[0] for _ in glob.glob("more_sequences/op*_fasta.fa") ]

MICcols = ["MIC meropenem (ug/ml)", "MIC doripenem (ug/ml)", "MIC ceftazidime (ug/ml)", "MIC piperacillin/tazobactam (ug/ml)",
            "MIC cefepime (ug/ml)", "MIC aztreonam (ug/ml)"]

sortBy = MICcols[0]
porin = "oprD"

document = Document()
session = Session(load_from_config=False)
session.use_doc('porin_map_server')
session.load_document(document)

source = ColumnDataSource(data=dict())
source2 = ColumnDataSource(data=dict())

#Load the meta information
df = pd.read_csv("mastersheet_porin_pics_sorting_data.txt", sep="\t")

allowedStrains = ["ARC" + str(_) for _ in df["Isolate"]]

strains = []
seqLength = 0

def update_data():

    global strains
    global seqLength

    fileName = "more_sequences/" + porin + "_fasta.fa"
    recs = [r for r in SeqIO.parse(fileName, "fasta")]

    recs[0].description = "PAO1_" + recs[0].description

    goodRecs = []

    goodRecs.append(recs[0])

    strains = [r.description.split("_")[0] for r in recs]
    strainsDone = [strains[0]]

    for r,strain in zip(recs[1:], strains[1:]):
        if strain in strainsDone:
            # print("Duplicate strain {}... dropping".format(strain))
            # if "*" in str(r):
            #     print("Uhoh! Dropped  It had a * too!")
            continue
        if strain not in allowedStrains:
            continue
        strainsDone.append(strain)
        goodRecs.append(r)

    recs = goodRecs
    strains = strainsDone


    #Sort by the asked for 
    dfSorted = df.sort(sortBy)

    sortedStrains = ["PAO1"] + ["ARC" + str(_) for _ in dfSorted["Isolate"]]
    strainIdx = [sortedStrains.index(s) for s in strains]

    sortedTup = sorted(zip(strainIdx, strains, recs), key=lambda t: t[0])

    _  = zip(*sorted(zip(strainIdx, strains, recs), key=lambda t: t[0]))

    strains = list(_[1])
    recs = list(_[2])

    numStrains = len(strains)
    numProtiens = len(recs[0].seq)

    pa = np.zeros((numStrains, numProtiens), dtype=np.uint8)

    refSeqList = list(str(recs[0].seq))
    seqLength = len(refSeqList)

    for ct,r in enumerate(recs):
        for ct2, c in enumerate(str(r.seq)):
            if c == refSeqList[ct2]:
                pa[ct,ct2] = 0
            else:
                if c == "-":
                    pa[ct,ct2] = 1
                elif c == "*":
                    pa[ct,ct2] = 2
                else:
                    pa[ct,ct2] = 3

    colours = ["#7fc97f", "#beaed4", "#fdc086", "#ffff99"]

    diffCounts = np.sum(pa>0, axis=0)

    source.data = {"strain":np.tile(strains, (numProtiens, 1)).flatten(order="F"),
                                        "pos":np.tile(np.arange(numProtiens), (numStrains,1)).flatten(order="C"),
                                        "refProtien":np.tile(refSeqList, (pa.shape[0],1)).flatten(order="C"),
                                        "curProtien":[c for r in recs for c in str(r.seq)],
                                        "diffCounts":np.tile(diffCounts, (numStrains, 1)).flatten(order="C"),
                                        "colour":[colours[_] for _ in pa.flatten()],
                                        "MIC":np.tile(dfSorted[sortBy].values, (numProtiens,1)).flatten(order="F")
                                        }

    #Load in the A,B,D data
    posBis = []
    strainBis = []
    coloursBis = []
    for ct,name in enumerate(["class A", "class B", "class D"]):
        posBis.extend([seqLength + 2.5 + 5*ct for _ in range(numStrains)])
        strainBis.extend(strains)
        coloursBis.extend(["#000000" if _ else "#ffffff" for _ in (dfSorted[name].values != '0')])

    source2.data = {
                    "strain": strainBis,
                    "pos": posBis,
                    "colour": coloursBis
                    }

    print("Loaded data for {} into source".format(porin))


update_data()


plot = Plot(title="oprD vs. PAO1 Differences", plot_width=1000, plot_height=3000,
            x_range=Range1d(start=-0.5, end=seqLength + 15),
            y_range=FactorRange(factors=strains))

plot.add_glyph(source, Rect(x="pos", y="strain", width=0.95, height=0.95, fill_color="colour", line_color=None))

plot.add_glyph(source2, Rect(x="pos", y="strain", width=5, height=0.95, fill_color="colour", line_color=None))

plot.add_tools(PanTool(), WheelZoomTool(), ResetTool(), HoverTool(tooltips = OrderedDict([
    ("strain", "@strain"),
    ("site", "@pos"),
    ("protien", "@curProtien"),
    ("PAO1", "@refProtien"),
    ("tot. diff.", "@diffCounts"),
    ("MIC", "@MIC")
    ])
    ) )

def update_plot():
    print("Updating plot ranges with new seqLength = {}".format(seqLength))
    plot.title = porin + " vs. PAO1 Differences"
    plot.x_range.end = seqLength + 15
    plot.y_range.factors = strains

def on_porin_change(obj, attr, old, new):
    "Load the new porin data"
    print("Updating to porin: " + new)
    global porin
    porin = new
    update_data()
    update_plot()
    session.store_document(document)

porinSelect = Select(title="Porin", value="oprD", options=porins)
porinSelect.on_change("value", on_porin_change)


def on_sort_change(obj, attr, old, new):
    global sortBy
    sortBy = new
    print("Now sorting by: {}".format(new))
    update_data()
    update_plot()
    session.store_document(document)

sortSelect = Select(title="Sort by", value=MICcols[0], options=MICcols)
sortSelect.on_change("value", on_sort_change)

controls = HBox(children=[porinSelect, sortSelect])

# columns = [TableColumn(field="Strain", header="Strain")]
# for m in MICcols:
#     columns.append(TableColumn(field=m, type="numeric", header=m))
# table = HandsonTable(source=source2, columns=columns)
# import ipdb; ipdb.set_trace()

layout = VBox(children=[controls, plot])
document.add(layout)

session.store_document(document)

if __name__ == "__main__":
    link = session.object_link(document.context)
    print("Please visit %s to see the plots" % link)
    view (link)
    print("\npress ctrl-C to exit")
    try:
        while True:
            session.load_document(document)
            time.sleep(2)
    except KeyboardInterrupt:
        print()
    except ConnectionError:
        print("Connection to bokeh-server was terminated")
