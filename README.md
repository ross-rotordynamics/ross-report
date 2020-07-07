# ROSS Report
This package aims to build a report for [ROSS](https://github.com/ross-rotordynamics/ross), a rotordynamics python package.

## Quick start to graphics submodule
In this first example we're going to see how to construct a simple `hello_world.html` with some text and an Plotly Figure.
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
    n=2, material=steel, width=0.07, i_d=0.05, o_d=0.28
)
disk1 = rs.DiskElement.from_geometry(
    n=4, material=steel, width=0.07, i_d=0.05, o_d=0.28
)

stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
stfy = [0.8e7, 0.9e7, 1.0e7, 1.1e7]
freq = [400, 800, 1200, 1600]
bearing0 = rs.BearingElement(0, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
bearing1 = rs.BearingElement(6, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)

rotor = rs.Rotor(shaft_elem, [disk0, disk1], [bearing0, bearing1])

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
report = Report(
        rotor=rotor,
        speed_range=(400, 1000),
        tripspeed=1200,
        bearing_stiffness_range=(5, 8),
        bearing_clearance_lists=bearings,
        speed_units="rad/s",
    )
    
D = [0.35, 0.35]
H = [0.08, 0.08]
HP = [10000, 10000]
RHO_ratio = [1.11, 1.14]
RHOd = 30.45
RHOs = 37.65
oper_speed = 1000.0
data = report.run(D, H, HP, oper_speed, RHO_ratio, RHOs, RHOd)
 
plot_rotor_1, ucs_fig_1, mode_fig_1 = report.assets_prep(data)["figs"]

text1 = """This is a report automatically generated using <a href= https://github.com/ross-rotordynamics/ross> ROSS</a>, a python package for rotordynamics analysis.
<br>Below there's a graphical representation of the rotor analyzed."""
text2 = """In this section the calculations carried out to evaluate the critical speed map and the rotor response to unbalance are described.
 The results of each calculation are shown at the end of this paragraph."""
text3 = """The undamped critical speed analysis is carried out according to API 617 7th edition para. 2.6.2.3. The rotor system as described in Appendix 1 is used. The bearings are represented by an equivalent spring constant between rotor and pedestals, which may then be considered as elastically mounted. Isotropic, linear bearing characteristics are assumed and no damping is considered present in the system. The stiffness range selected for the calculation is such to properly describe the behavior of the rotor and provide the required information to perform the next analysis steps. The actual stiffness range (achievable by adjusting bearing clearance) is much more limited and always inside the calculation range. The rotordynamic system is solved and the undamped lateral critical speeds are calculated as a function of support equivalent stiffness over the user defined stiffness range. The results are summarized in the critical speed maps as shown in the following pages. Superimposed on the same plot are the horizontal and vertical Bearing Clearance curves (Kxx and Kzz ) either for maximum and minimum Bearing Clearance. The intersections of the vertical Bearing Clearance and critical speed curves provide the undamped critical speed values and give, in a preliminary way, a rough estimation of the critical speed and Bearing Clearance range in operation. The 1st and 2nd mode shapes for maximum and minimum Bearing Clearance are also attached, with the only intent of mode shape identification. Therefore, the vibration amplitudes are normalized with respect to the maximum level.
"""
page1 = Page(
    content=[Text(text=text1),
             PlotlyFigure(figure=plot_rotor_1),
             Title(title="Critical Speed Analysis"),
             Text(text=text2),
             Title("Undamped Critical Speed Analysis"),
             Text(text=text3),
             ]
             )
page2 = Page(content=[PlotlyFigure(figure=ucs_fig_1), PlotlyFigure(figure=mode_fig_1),])
pages = [page1, page2]

html = Layout(pages=pages).render_html_str()

with open("report.html", "w") as f:
    f.write(html)

``` 
