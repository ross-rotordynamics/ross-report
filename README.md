# ROSS Report
![github actions](https://github.com/ross-rotordynamics/ross-report/workflows/Test/badge.svg)
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
fig.layout["title"]["text"] = "bars"
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
