# ROSS Report
[![Build Status](https://travis-ci.com/ross-rotordynamics/ross-report.svg?branch=master)](https://travis-ci.com/ross-rotordynamics/ross-report)
[![Build status](https://ci.appveyor.com/api/projects/status/5s4tdwsm7nk519go/branch/master?svg=true)](https://ci.appveyor.com/project/GabrielBachmannArimCarneirodeAlbuquerque/ross-report/branch/master)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

This package aims to build a report for [ROSS](https://github.com/ross-rotordynamics/ross), a rotordynamics python package.

## Quick start to graphics submodule
The ross-report's structure is very simple, in the highest level we have an object called Layout, which constructs the whole html page.  
Layout is composed of Pages which is composed of Content objects (like Text or PlotlyFigure), to arrange all the components on the html version you can simply put them in order.
For a static PDF version of the report you can use `CTRL + P` on a browser like [chrome](https://www.google.com/intl/pt-BR/chrome/).

In this first example we're going to see how to construct a simple `hello_world.html` with some text and a Plotly Figure.
```Python
from report.graphics import *
import plotly.graph_objects as go

simple_text = """Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.
"""
another_simple_text = simple_text[:-2][::-1].capitalize() + '.'   # actually the same

fig = go.FigureWidget(data=go.Bar(y=[2, 3, 1]))
fig.layout["title"]["text"] = "Boring bars"
fig.layout["height"] = 500
fig.layout["width"] = 750

page = Page(content=[
                     Text(simple_text),
                     PlotlyFigure(fig),
                     Text(another_simple_text),
                     ]
            )

layout = Layout(pages=page, main_title="Hello World")
html_str = layout.render_html_str()

with open("hello_world.html", "w") as f:
    f.write(html_str)

```
The `hello_world.html` should look like [this](https://rawcdn.githack.com/ross-rotordynamics/ross-report/eb0d73c4462cd584f0f2ec4cc40047a91e952918/hello_world.html).
## Complete example
In this example we first instantiate a ross.Report, run some analysis with the analysis subpackage and then export an html file containing the graphical report. 
```python
from report.analysis import *
from report.graphics import *
import ross as rs
import numpy as np

i_d = 0
o_d = 0.05
n = 6
L = [0.25 for _ in range(n)]

shaft_elem = [
    rs.ShaftElement(
        l,
        i_d,
        o_d,
        material=rs.steel,
        shear_effects=True,
        rotary_inertia=True,
        gyroscopic=True,
    )
    for l in L
]

disk0 = rs.DiskElement.from_geometry(
    n=2, material=rs.steel, width=0.07, i_d=0.05, o_d=0.28
)
disk1 = rs.DiskElement.from_geometry(
    n=4, material=rs.steel, width=0.07, i_d=0.05, o_d=0.28
)

stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
stfy = [0.8e7, 0.9e7, 1.0e7, 1.1e7]
freq = [400, 800, 1200, 1600]
bearing0 = rs.BearingElement(0, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
bearing1 = rs.BearingElement(6, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
oper_clearance_brg = [bearing0, bearing1]
rotor = rs.Rotor(shaft_elem, [disk0, disk1], oper_clearance_brg)

# coefficients for minimum clearance
stfx = [0.7e7, 0.8e7, 0.9e7, 1.0e7]
dampx = [2.0e3, 1.9e3, 1.8e3, 1.7e3]
freq = [400, 800, 1200, 1600]
bearing0 = rs.BearingElement(0, kxx=stfx, cxx=dampx, frequency=freq)
bearing1 = rs.BearingElement(6, kxx=stfx, cxx=dampx, frequency=freq)
min_clearance_brg = [bearing0, bearing1]

# coefficients for maximum clearance
stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
dampx = [2.8e3, 2.7e3, 2.6e3, 2.5e3]
freq = [400, 800, 1200, 1600]
bearing0 = rs.BearingElement(0, kxx=stfx, cxx=dampx, frequency=freq)
bearing1 = rs.BearingElement(6, kxx=stfx, cxx=dampx, frequency=freq)
max_clearance_brg = [bearing0, bearing1]

bearings = [min_clearance_brg, max_clearance_brg]
D = [0.35, 0.35]
H = [0.08, 0.08]
HP = [10000, 10000]
RHO_ratio = [1.11, 1.14]
RHOd = 30.45
RHOs = 37.65
oper_speed = 1000.0
config = rp.Config()

config.update_config(
    rotor_properties={
        "rotor_speeds": {
            "min_speed": 400,
            "max_speed": 1000,
            "oper_speed": 1000,
            "trip_speed": 1200,
            "unit": "rpm",
        }
    },
    bearings={
        "min_clearance": min_clearance_brg,
        "oper_clearance": min_clearance_brg,
        "max_clearance": max_clearance_brg,
    },
    run_campbell={
        "speed_range": np.linspace(100, 1200, 91),
        "frequency_units": "rpm",
    },
    run_unbalance_response={
        "probes": {
            "node": [0, 6],
            "orientation": [np.pi / 2, np.pi / 2],
            "unit": "rad",
        },
        "frequency_range": np.linspace(100, 1200, 101),
        "frequency_units": "rad/s",
        "plot_deflected_shape": {"speed": [815]},
    },
    plot_ucs={"stiffness_range": (6, 11)},
    stability_level1={
        "D": D,
        "H": H,
        "rated_power": HP,
        "rho_ratio": RHO_ratio,
        "rho_suction": RHOs,
        "rho_discharge": RHOd,
        "length_unit": "m",
        "power_unit": "hp",
        "density_unit": "kg/m**3",
    },
)
report = Report(
    rotor=rotor,
    config=config
)


results = report.run_report()

layout = report.generate_standard_layout_report(results)
html = layout.render_html_str()

with open("report.html", "w") as f:
    f.write(html)

``` 
The `report.html` file should look like [this](https://raw.githack.com/ross-rotordynamics/ross-report/master/examples/report.html).
