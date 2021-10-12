import tecplot
import os

tecplot.new_layout()
frame = tecplot.active_frame()
files = os.listdir()
files = [f for f in files if "_interior.szplt" in f]
data = tecplot.data.load_tecplot_szl(files)
tecplot.save_layout("foo.lay")
