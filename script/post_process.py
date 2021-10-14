import tecplot
import os
import sys

# basic plot setup
tecplot.new_layout()
frame = tecplot.active_frame()
if len(sys.argv) <= 1:
    work_dir = "."
else:
    work_dir = sys.argv[1]
files = os.listdir(work_dir)
files = [f"{work_dir}/{f}" for f in files if ".szplt" in f]
if not files:
    raise Exception("No `.szplt` files found in specified directory.")
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

# new variables
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
equations += f"""
{{temperature}} = {{pressure}}/({{density}}*{gas_const})
{{sound_speed}} = ({heat_rat}*{{pressure}}*{{density}})**0.5
"""[1:]
equations += "{mach} = {velocity_magnitude}/{sound_speed}\n"
tecplot.data.operate.execute_equation(equations, ignore_divide_by_zero=True)
if n_dim >= 2:
    plot.contour(0).variable = data.variable("mach")
    plot.vector.u_variable = data.variable("velocity0")
    plot.vector.v_variable = data.variable("velocity1")
if n_dim >= 3:
    plot.vector.w_variable = data.variable("velocity2")

# Zone Style
if n_dim > 1:
    for fieldmap in plot.active_fieldmaps:
        for zone in fieldmap.zones:
            if "qpoints" in zone.name:
                fieldmap.scatter.size = 1.
                fieldmap.scatter.symbol_type = tecplot.constant.SymbolType.Geometry
                fieldmap.scatter.symbol().shape = tecplot.constant.GeomShape.Diamond
            else:
                fieldmap.scatter.show = False
            if not ("edges" in zone.name):
                fieldmap.mesh.show = False
            if not ("interior" in zone.name):
                fieldmap.contour.show = False
                fieldmap.edge.show = False
else:
    plot.delete_linemaps()
    n_zones = data.num_zones
    for i_zone in range(n_zones):
        zone = data.zone(i_zone)
        linemap = plot.add_linemap(zone=zone, x=data.variable("position0"), y=data.variable("mach"))
        if "interior" in zone.name:
            linemap.line.show = True
        else:
            linemap.line.show = False
        if "qpoints" in zone.name:
            linemap.symbols.show = True
            linemap.symbols.symbol_type = tecplot.constant.SymbolType.Geometry
            linemap.symbols.symbol().shape = tecplot.constant.GeomShape.Diamond
            linemap.symbols.color = tecplot.constant.Color.Black
        else:
            linemap.symbols.show = False
    plot.view.fit()

tecplot.save_layout(f"{work_dir}/all.lay")
