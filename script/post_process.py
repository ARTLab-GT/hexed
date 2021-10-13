import tecplot
import os

# basic plot setup
tecplot.new_layout()
frame = tecplot.active_frame()
files = os.listdir()
files = [f for f in files if ".szplt" in f]
data = tecplot.data.load_tecplot_szl(files)
variable_names = frame.dataset.variable_names
n_dim = len([var for var in variable_names if "pos" in var])
if n_dim == 3:
    plot_type = tecplot.constant.PlotType.Cartesian3D
elif n_dim == 2:
    plot_type = tecplot.constant.PlotType.Cartesian2D
else:
    plot_type = tecplot.constant.PlotType.XYLine
plot = frame.plot(plot_type)
plot.activate()

# configure Zone Style
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
            fm.edge.show = False

# compute new variables
heat_rat = 1.4
gas_const = 287.058
equations = ""
for i_dim in range(n_dim):
    equations += f"{{velocity{i_dim}}} = {{state{i_dim}}}/{{state{n_dim}}}\n"
equations += "{velocity_magnitude} = ("
for i_dim in range(n_dim):
    equations += f"{{velocity{i_dim}}}*{{velocity{i_dim}}} + "
equations = equations[:-3] + ")**0.5\n"
equations += f"""
{{density}} = {{state{n_dim}}}
{{pressure}} = ({{state{n_dim + 1}}} - 0.5*{{density}}*{{velocity_magnitude}}**2)*{heat_rat - 1}
"""[1:]
tecplot.data.operate.execute_equation(equations, ignore_divide_by_zero=True)

tecplot.save_layout("all.lay")
