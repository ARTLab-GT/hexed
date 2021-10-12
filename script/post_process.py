import tecplot
import os

tecplot.new_layout()
frame = tecplot.active_frame()
files = os.listdir()
files = [f for f in files if ".szplt" in f]
files = [f for f in files if not ("qpoint" in f)]
files = [f for f in files if not ("interior" in f)]
data = tecplot.data.load_tecplot_szl(files)
plot = frame.plot()
for fieldmap in plot.active_fieldmaps:
    for zone in fieldmap.zones:
        name = zone.name
    if not ("edges" in zone.name):
        plot.fieldmap(zone).mesh.show = False
#plot.fieldmap(0).mesh.show = False
tecplot.save_layout("foo.lay")
