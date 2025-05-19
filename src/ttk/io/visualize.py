import collections
import os
import tempfile

import numpy as np
import py3Dmol
from Bio import AlignIO, SeqIO
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, Grid, Plot, Range1d
from bokeh.models.glyphs import Rect, Text
from bokeh.plotting import figure

from ttk.io import PDBFile, align

# import matplotlib.colors as mcolors
colors = [
    "azure",
    "coral",
    "deepskyblue",
    "yellow",
    "peachpuff",
    "violet",
    "lavenderblush",
    "bisque",
    "lawngreen",
    "pink",
    "goldenrod",
    "orangered",
    "olive",
    "grey",
    "lightsteelblue",
    "lightseagreen",
    "blue",
    "cyan",
    "lightslategrey",
    "indigo",
    "darkgray",
    "hotpink",
    "lightsalmon",
    "indianred",
    "mediumorchid",
    "cornsilk",
    "linen",
    "skyblue",
    "darkcyan",
    "darkgrey",
    "darkturquoise",
    "slategrey",
    "honeydew",
    "dodgerblue",
    "blanchedalmond",
    "chocolate",
    "lightcoral",
    "darkslategrey",
    "tomato",
    "forestgreen",
    "mediumseagreen",
    "gold",
    "burlywood",
    "lemonchiffon",
    "olivedrab",
    "darkseagreen",
    "lightslategray",
    "lightgoldenrodyellow",
    "sandybrown",
    "navajowhite",
    "mistyrose",
    "darkmagenta",
    "yellowgreen",
    "lightgray",
    "oldlace",
    "mediumturquoise",
    "papayawhip",
    "rebeccapurple",
    "dimgray",
    "beige",
    "steelblue",
    "lime",
    "darkorange",
    "teal",
    "darkslateblue",
    "sienna",
    "brown",
    "turquoise",
    "seashell",
    "darksalmon",
    "palegoldenrod",
    "chartreuse",
    "plum",
    "darkorchid",
    "wheat",
    "peru",
    "springgreen",
    "rosybrown",
    "black",
    "darkblue",
    "midnightblue",
    "magenta",
    "thistle",
    "saddlebrown",
    "mediumslateblue",
    "ivory",
    "mediumaquamarine",
    "lightgrey",
    "palevioletred",
    "khaki",
    "darkkhaki",
    "slategray",
    "gainsboro",
    "cornflowerblue",
    "lightblue",
    "floralwhite",
    "antiquewhite",
    "dimgrey",
    "fuchsia",
    "orange",
    "lavender",
    "snow",
    "deeppink",
    "aliceblue",
    "mintcream",
    "darkviolet",
    "palegreen",
    "lightpink",
    "mediumspringgreen",
    "limegreen",
    "paleturquoise",
    "lightyellow",
    "royalblue",
    "aqua",
    "mediumpurple",
    "darkred",
    "lightskyblue",
    "firebrick",
    "mediumblue",
    "gray",
    "orchid",
    "aquamarine",
    "white",
    "darkgoldenrod",
    "purple",
    "red",
    "darkslategray",
    "navy",
    "lightcyan",
    "mediumvioletred",
    "maroon",
    "tan",
    "salmon",
    "whitesmoke",
    "darkgreen",
    "cadetblue",
    "greenyellow",
    "silver",
    "powderblue",
    "seagreen",
    "ghostwhite",
    "green",
    "crimson",
    "moccasin",
    "darkolivegreen",
    "lightgreen",
    "slateblue",
    "blueviolet",
]

alphabet = list("ARNDCQEGHILKMFPSTWYV")
clr_dict = {}
for i, each in enumerate(alphabet):
    clr_dict[each] = colors[i]
clr_dict["-"] = "white"


def visualise_topology(topology, show_interface=False):
    if show_interface:
        interfaces = collections.defaultdict(set)
        num_chains = len(topology.chains)
        for i in range(num_chains):
            for j in range(i + 1, num_chains):
                c1, c2 = topology.chains[i], topology.chains[j]
                c1_res, c2_res = find_interface(c1, c2)
                c1_res_idx, c2_res_idx = [x.index for x in c1_res], [
                    x.index for x in c2_res
                ]
                interfaces[i] = interfaces[i].union(c1_res_idx)
                interfaces[j] = interfaces[j].union(c2_res_idx)

    pdb_content = PDBFile(use_res_index=True).to_content(topology)
    view = py3Dmol.view(width=500, height=500)
    view.addModelsAsFrames(pdb_content)

    view.setHoverable(
        {},
        True,
        """function(atom,viewer,event,container) {
                       if(!atom.label) {
                        atom.label = viewer.addLabel(atom.resi +" "+ atom.resn + ":" + atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                       }}""",
        """function(atom,viewer) {
                       if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                       }
                    }""",
    )

    for idx, chain in enumerate(topology.chains):
        view.setStyle({"chain": chain.id}, {"cartoon": {"color": colors[idx]}})

    if show_interface:
        for chain_idx, res_idx in interfaces.items():
            view.addStyle(
                {"chain": topology.chains[chain_idx].id, "resi": list(res_idx)},
                {"stick": {"colorscheme": "purpleCarbon"}},
            )

    view.setStyle({"hetflag": True}, {"stick": {}})
    #  view.addStyle({
    #  'chain': 'A',
    #  'resi': ["45", "20", "22", "156", "15"]
    #  }, {'stick': {
    #  'colorscheme': 'purpleCarbon'
    #  }})
    # view.setStyle({'model': -1, 'resi':["246"]}, {"cartoon": {'color': "yellow"}})
    # view.setStyle({'hetflag':False}, {"cartoon": {'color': 'spectrum'}})
    # view.rotate(180,'y',4000)
    view.zoomTo()
    view.show()


def visualise_pdb(pdb_fpath):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(open(pdb_fpath, "r").read(), "pdb")

    view.setHoverable(
        {},
        True,
        """function(atom,viewer,event,container) {
                       if(!atom.label) {
                        atom.label = viewer.addLabel(atom.resi +" "+ atom.resn + ":" + atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                       }}""",
        """function(atom,viewer) {
                       if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                       }
                    }""",
    )

    view.setStyle({"model": -1}, {"cartoon": {"color": "spectrum"}})
    view.setStyle({"hetflag": True}, {"stick": {}})

    view.zoomTo()
    view.show()


def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs = clr_dict
    colors = [clrs[i] for i in text]
    return colors


def align_n_view_seqs(seqs_dict, muscle_exe="muscle", fontsize="9pt", plot_width=800):
    seq_record_ls = align.build_seq_record_from_seqs(seqs_dict)
    aln = align.muscle_alignment(seq_record_ls, muscle_exe=muscle_exe)
    p = view_alignment(aln, fontsize=fontsize, plot_width=plot_width)
    return p, aln


def view_alignment(aln, fontsize="9pt", plot_width=800, starting_idx=0):
    """Bokeh sequence alignment view"""

    # make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)
    N = len(seqs[0])
    S = len(seqs)
    width = 0.4

    x = np.arange(1, N + 1)
    y = np.arange(0, S, 1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy + 0.5
    h = 1 / S
    # now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs) * 15 + 50
    x_range = Range1d(0, N + 1, bounds="auto")
    if N > 100:
        viewlen = 100
    else:
        viewlen = N
    # view_range is for the close up view
    view_range = (starting_idx, starting_idx + viewlen)
    tools = "xpan, xwheel_zoom, reset, save"

    # entire sequence view (no text, with zoom)
    p = figure(
        title=None,
        plot_width=plot_width,
        plot_height=200,
        x_range=x_range,
        y_range=(0, S),
        tools=tools,
        min_border=0,
        toolbar_location="below",
    )
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.6,
    )
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False

    # sequence text view with ability to scroll along x axis
    p1 = figure(
        title=None,
        plot_width=plot_width,
        plot_height=plot_height,
        x_range=view_range,
        y_range=ids,
        tools="xpan,reset",
        min_border=0,
        toolbar_location="below",
    )  # , lod_factor=1)
    glyph = Text(
        x="x",
        y="y",
        text="text",
        text_align="center",
        text_color="black",
        text_font="monospace",
        text_font_size=fontsize,
    )
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color=None,
        fill_alpha=0.4,
    )
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p], [p1]], toolbar_location="below")
    return p
