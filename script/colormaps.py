## \file colormaps.py
# write matplotlib [perceptually uniform](https://matplotlib.org/stable/tutorials/colors/colormaps.html#sequential)
# colormaps in a format that can be imported by tecplot

## \cond

from matplotlib.pyplot import get_cmap
import os

names = ["magma", "inferno", "plasma", "viridis", "cividis", "twilight", "turbo"]
n_point = 50 # this seems to be roughly the maximum that tecplot can handle

dir_name = "colormaps"
if dir_name not in os.listdir():
    os.mkdir(dir_name)

for name in names:
    try:
        cmap = get_cmap(name)
        text = f"""
#!MC 1410
$!CreateColorMap
  Name = 'hexed_{name}'
  NumControlPoints = {n_point}
"""[1:]
        for i_point in range(n_point):
            frac = i_point/(n_point - 1)
            color = cmap(frac)
            rgb = "\n      ".join([f"{'RGB'[i]} = {int(color[i]*255)}" for i in range(3)])
            text += f"""
  ControlPoint {i_point + 1}
  {{
    ColorMapFraction = {frac}
    LeadRGB
    {{
      {rgb}
    }}
    TrailRGB
    {{
      {rgb}
    }}
  }}
"""[1:]
        with open(f"{dir_name}/{name}.map", "w") as output_file:
            output_file.write(text)
    except:
        pass

## \endcond
