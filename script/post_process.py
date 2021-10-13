import tecplot
import os

tecplot.new_layout()
frame = tecplot.active_frame()
files = os.listdir()
files = [f for f in files if ".szplt" in f]
data = tecplot.data.load_tecplot_szl(files)
plot = frame.plot()
for fieldmap in plot.active_fieldmaps:
    for zone in fieldmap.zones:
        fm = plot.fieldmap(zone)
        if "qpoints" in zone.name:
            fm.scatter.size = 1.
            fm.scatter.symbol_type = tecplot.constant.SymbolType.Geometry
            fm.scatter.symbol().shape = tecplot.constant.GeomShape.Diamond
        else:
            fm.scatter.show = False
        if not ("edges" in zone.name):
            fm.mesh.show = False
        if not ("interior" in zone.name):
            fm.contour.show = False
tecplot.save_layout("foo.lay")
