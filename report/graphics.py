from pathlib import Path
import report
import base64 as b64
from plotly.graph_objs import Figure


class CSS:
    def __init__(self, path=Path(report.__file__).parent/"style.css"):
        self.path = Path(path)

    def __str__(self):
        with open(self.path) as css_file:
            css_code = css_file.read()
        return css_code

    def __repr__(self):
        return "CSS"


class Content:
    def __init__(self):
        pass

    def render_html_str(self):
        pass


class Text(Content):
    def __init__(self, text):
        super().__init__()
        self.text = text

    def render_html_str(self):
        html = (
            """<div class="row offset-2">\n
                    <div class="col-9  mt-3 mb-0 pb-0">\n
                        <p class="text-justify">\n
               """
            + self.text
            + "\n"
            + """
            </p>\n
        </div>\n
        </div>\n"""
        )
        return html


class Img(Content):
    """
    """

    def __init__(self, path_to_img):
        self.path_to_img = Path(path_to_img)

    def render_html_str(self):
        html = ""
        return html

    def __repr__(self):
        return f"{self.path_to_img}"


class Title(Content):
    """
    """

    def __init__(self, title):
        self.title = title

    def render_html_str(self):
        html = f"""
        <div class="row offset-2">
            <div class="col-9">
                <h1 class="h1">
                    {self.title}
                </h1>
            </div>
        </div>
"""

        return html


class Page:
    """
    """

    def __init__(
        self, content=None,
    ):
        for item in content:
            assert isinstance(item, Content) or isinstance(
                item, str
            ), "Every item of content must be either content or a string"
        self.content = content

    def render_html_str(self):
        html = ""
        for item in self.content:
            html += item.render_html_str()
        return html


class PlotlyFigure(Content):
    def __init__(self, figure):
        self.figure = figure

    def render_html_str(self):
        html = (
            """
        <div style="width: 750px;height: 350px" class="mx-auto">\n
        """
            + self.figure.to_html(full_html=False)
            + """\n
        </div>\n"""
        )
        return html


class Layout:
    """Report Layout
    """

    def __init__(self, css=CSS(), pages=None):
        assert (
            isinstance(pages, list) or isinstance(pages, Page) or pages is None
        ), "pages argument must be either a Page or a list of Pages."
        if isinstance(pages, list):
            for page in pages:
                assert isinstance(page, Page), "All items in page must be Pages."
        self.pages = pages

        assert isinstance(css, CSS), "css must be an CSS object."
        self.css = css

    def __repr__(self):
        rep = {}
        if isinstance(self.pages, list):
            for page in range(len(self.pages)):
                rep[self.pages[page]] = self.pages[page]
        elif isinstance(self.pages, Page):
            rep[self.pages] = self.pages.content
        return str(rep)

    def render_html_str(self):
        html = """
        <!DOCTYPE html>
        <html lang="en">
        <head>
            """+"""
            <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
            <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>
            """
        html = (
            html
            + "<style>"
            + str(self.css)
            + "</style>"
            + """
            <style>
            
                @media print { 
                    @page {
                        size: A4; /* DIN A4 standard, Europe */
                        margin: 1cm;
                        }
                    .offset-2{
                        margin-left: 3cm;
                        margin-right: 0;
                        margin-top: .5cm;
                    }
                    .mt-img{
                        margin-top: 5cm;
                    }
                    p{
                        width: 100ch;
                    }
                }
            </style>
            <meta charset="UTF-8">
            <title>ROSS Report</title>
            <link rel="icon" href="data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICAgICAgIHhtbG5zOm9zYj0iaHR0cDovL3d3dy5vcGVuc3dhdGNoYm9vay5vcmcvdXJpLzIwMDkvb3NiIgogICAgICAgIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIKICAgICAgICB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIgogICAgICAgIHhtbG5zOnJkZj0iaHR0cDovL3d3dy53My5vcmcvMTk5OS8wMi8yMi1yZGYtc3ludGF4LW5zIyIKICAgICAgICB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciCiAgICAgICAgeG1sbnM6c29kaXBvZGk9Imh0dHA6Ly9zb2RpcG9kaS5zb3VyY2Vmb3JnZS5uZXQvRFREL3NvZGlwb2RpLTAuZHRkIgogICAgICAgIHhtbG5zOmlua3NjYXBlPSJodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy9uYW1lc3BhY2VzL2lua3NjYXBlIgogICAgICAgIHdpZHRoPSIxNC4xMTA5OTFtbSIKICAgICAgICBoZWlnaHQ9IjEzLjM5MTAzNW1tIgogICAgICAgIHZpZXdCb3g9IjAgMCA0OS45OTk1NjkgNDcuNDQ4NTQ2IgogICAgICAgIGlkPSJzdmcyIgogICAgICAgIHZlcnNpb249IjEuMSIKICAgICAgICBpbmtzY2FwZTp2ZXJzaW9uPSIwLjkxIHIxMzcyNSIKICAgICAgICBzb2RpcG9kaTpkb2NuYW1lPSJyb3NzLWxvZ28uc3ZnIj4KICA8ZGVmcwogICAgIGlkPSJkZWZzNCI+CiAgICA8bGluZWFyR3JhZGllbnQKICAgICAgIGlkPSJsaW5lYXJHcmFkaWVudDQ0NjEiCiAgICAgICBvc2I6cGFpbnQ9InNvbGlkIj4KICAgICAgPHN0b3AKICAgICAgICAgc3R5bGU9InN0b3AtY29sb3I6I2Y4ZjhmODtzdG9wLW9wYWNpdHk6MTsiCiAgICAgICAgIG9mZnNldD0iMCIKICAgICAgICAgaWQ9InN0b3A0NDYzIiAvPgogICAgPC9saW5lYXJHcmFkaWVudD4KICAgIDxpbmtzY2FwZTpwZXJzcGVjdGl2ZQogICAgICAgc29kaXBvZGk6dHlwZT0iaW5rc2NhcGU6cGVyc3AzZCIKICAgICAgIGlua3NjYXBlOnZwX3g9IjAgOiA1MjYuMTgxMDQgOiAxIgogICAgICAgaW5rc2NhcGU6dnBfeT0iMCA6IDk5OS45OTk4NyA6IDAiCiAgICAgICBpbmtzY2FwZTp2cF96PSI3NDQuMDk0NDQgOiA1MjYuMTgxMDMgOiAxIgogICAgICAgaW5rc2NhcGU6cGVyc3AzZC1vcmlnaW49IjM3Mi4wNDcyIDogMzUwLjc4NzM1IDogMSIKICAgICAgIGlkPSJwZXJzcGVjdGl2ZTQxMzgiIC8+CiAgPC9kZWZzPgogIDxzb2RpcG9kaTpuYW1lZHZpZXcKICAgICBpZD0iYmFzZSIKICAgICBwYWdlY29sb3I9IiNmZmZmZmYiCiAgICAgYm9yZGVyY29sb3I9IiM2NjY2NjYiCiAgICAgYm9yZGVyb3BhY2l0eT0iMS4wIgogICAgIGlua3NjYXBlOnBhZ2VvcGFjaXR5PSIwLjAiCiAgICAgaW5rc2NhcGU6cGFnZXNoYWRvdz0iMiIKICAgICBpbmtzY2FwZTp6b29tPSIyLjgiCiAgICAgaW5rc2NhcGU6Y3g9Ii00Ljc5MDE2ODYiCiAgICAgaW5rc2NhcGU6Y3k9IjI2LjUzNDEyNCIKICAgICBpbmtzY2FwZTpkb2N1bWVudC11bml0cz0icHgiCiAgICAgaW5rc2NhcGU6Y3VycmVudC1sYXllcj0ibGF5ZXIxIgogICAgIHNob3dncmlkPSJ0cnVlIgogICAgIHNob3dib3JkZXI9ImZhbHNlIgogICAgIGZpdC1tYXJnaW4tdG9wPSIwIgogICAgIGZpdC1tYXJnaW4tbGVmdD0iMCIKICAgICBmaXQtbWFyZ2luLXJpZ2h0PSIwIgogICAgIGZpdC1tYXJnaW4tYm90dG9tPSIwIgogICAgIGlua3NjYXBlOndpbmRvdy13aWR0aD0iMTMyOSIKICAgICBpbmtzY2FwZTp3aW5kb3ctaGVpZ2h0PSI3NDQiCiAgICAgaW5rc2NhcGU6d2luZG93LXg9IjM3IgogICAgIGlua3NjYXBlOndpbmRvdy15PSIyNCIKICAgICBpbmtzY2FwZTp3aW5kb3ctbWF4aW1pemVkPSIxIiAvPgogIDxtZXRhZGF0YQogICAgIGlkPSJtZXRhZGF0YTciPgogICAgPHJkZjpSREY+CiAgICAgIDxjYzpXb3JrCiAgICAgICAgIHJkZjphYm91dD0iIj4KICAgICAgICA8ZGM6Zm9ybWF0PmltYWdlL3N2Zyt4bWw8L2RjOmZvcm1hdD4KICAgICAgICA8ZGM6dHlwZQogICAgICAgICAgIHJkZjpyZXNvdXJjZT0iaHR0cDovL3B1cmwub3JnL2RjL2RjbWl0eXBlL1N0aWxsSW1hZ2UiIC8+CiAgICAgICAgPGRjOnRpdGxlPjwvZGM6dGl0bGU+CiAgICAgIDwvY2M6V29yaz4KICAgIDwvcmRmOlJERj4KICA8L21ldGFkYXRhPgogIDxnCiAgICAgaW5rc2NhcGU6Z3JvdXBtb2RlPSJsYXllciIKICAgICBpZD0ibGF5ZXIyIgogICAgIGlua3NjYXBlOmxhYmVsPSJiYWNrZ3JvdW5kIgogICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0zLjU1MzM1MDgsNS44NTY3MDQ1KSIgLz4KICA8ZwogICAgIGlua3NjYXBlOmxhYmVsPSJMYXllciAxIgogICAgIGlua3NjYXBlOmdyb3VwbW9kZT0ibGF5ZXIiCiAgICAgaWQ9ImxheWVyMSIKICAgICB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMjM5LjY4MDcyLC03MDYuNTA1NDQpIj4KICAgIDxwYXRoCiAgICAgICBzdHlsZT0ib3BhY2l0eToxO2ZpbGw6IzU2NTY1NjtmaWxsLW9wYWNpdHk6MTtzdHJva2U6I2Q2MjcyODtzdHJva2Utd2lkdGg6MS4wMjA4NDQ4MjtzdHJva2UtbWl0ZXJsaW1pdDo0O3N0cm9rZS1kYXNoYXJyYXk6bm9uZTtzdHJva2UtZGFzaG9mZnNldDowO3N0cm9rZS1vcGFjaXR5OjAiCiAgICAgICBkPSJNIDI1IDQuMjM0Mzc1IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCA1LjUwOTc2NTYgMjMuNzI0NjA5IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCAyNSA0My4yMTQ4NDQgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDQ0LjQ5MDIzNCAyMy43MjQ2MDkgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDI1IDQuMjM0Mzc1IHogTSAyNy43NDAyMzQgMTMuOTEyMTA5IEEgMTMuMDcxNzQ0IDEzLjA3MTc0NCAwIDAgMSA0MC44MTI1IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDQwLjA1NjY0MSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMTQuNjY3OTY5IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDEzLjkxMjEwOSB6ICIKICAgICAgIGlkPSJwYXRoNDEzNiIKICAgICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKDIzOS42ODA3Miw3MDYuNTA1NDQpIiAvPgogIDwvZz4KPC9zdmc+Cg==">
        </head>
        <body>
    
        <div class="container-fluid">
            
            <div class="mt-4 offset-2 row">
                <div class="col-9">
                    <h1 class="h1">
                        ROSS Report
                    </h1>
                </div>

            <div class="col-3">
                <img id="ross-logo" src="data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICAgICAgIHhtbG5zOm9zYj0iaHR0cDovL3d3dy5vcGVuc3dhdGNoYm9vay5vcmcvdXJpLzIwMDkvb3NiIgogICAgICAgIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIKICAgICAgICB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIgogICAgICAgIHhtbG5zOnJkZj0iaHR0cDovL3d3dy53My5vcmcvMTk5OS8wMi8yMi1yZGYtc3ludGF4LW5zIyIKICAgICAgICB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciCiAgICAgICAgeG1sbnM6c29kaXBvZGk9Imh0dHA6Ly9zb2RpcG9kaS5zb3VyY2Vmb3JnZS5uZXQvRFREL3NvZGlwb2RpLTAuZHRkIgogICAgICAgIHhtbG5zOmlua3NjYXBlPSJodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy9uYW1lc3BhY2VzL2lua3NjYXBlIgogICAgICAgIHdpZHRoPSIxNC4xMTA5OTFtbSIKICAgICAgICBoZWlnaHQ9IjEzLjM5MTAzNW1tIgogICAgICAgIHZpZXdCb3g9IjAgMCA0OS45OTk1NjkgNDcuNDQ4NTQ2IgogICAgICAgIGlkPSJzdmcyIgogICAgICAgIHZlcnNpb249IjEuMSIKICAgICAgICBpbmtzY2FwZTp2ZXJzaW9uPSIwLjkxIHIxMzcyNSIKICAgICAgICBzb2RpcG9kaTpkb2NuYW1lPSJyb3NzLWxvZ28uc3ZnIj4KICA8ZGVmcwogICAgIGlkPSJkZWZzNCI+CiAgICA8bGluZWFyR3JhZGllbnQKICAgICAgIGlkPSJsaW5lYXJHcmFkaWVudDQ0NjEiCiAgICAgICBvc2I6cGFpbnQ9InNvbGlkIj4KICAgICAgPHN0b3AKICAgICAgICAgc3R5bGU9InN0b3AtY29sb3I6I2Y4ZjhmODtzdG9wLW9wYWNpdHk6MTsiCiAgICAgICAgIG9mZnNldD0iMCIKICAgICAgICAgaWQ9InN0b3A0NDYzIiAvPgogICAgPC9saW5lYXJHcmFkaWVudD4KICAgIDxpbmtzY2FwZTpwZXJzcGVjdGl2ZQogICAgICAgc29kaXBvZGk6dHlwZT0iaW5rc2NhcGU6cGVyc3AzZCIKICAgICAgIGlua3NjYXBlOnZwX3g9IjAgOiA1MjYuMTgxMDQgOiAxIgogICAgICAgaW5rc2NhcGU6dnBfeT0iMCA6IDk5OS45OTk4NyA6IDAiCiAgICAgICBpbmtzY2FwZTp2cF96PSI3NDQuMDk0NDQgOiA1MjYuMTgxMDMgOiAxIgogICAgICAgaW5rc2NhcGU6cGVyc3AzZC1vcmlnaW49IjM3Mi4wNDcyIDogMzUwLjc4NzM1IDogMSIKICAgICAgIGlkPSJwZXJzcGVjdGl2ZTQxMzgiIC8+CiAgPC9kZWZzPgogIDxzb2RpcG9kaTpuYW1lZHZpZXcKICAgICBpZD0iYmFzZSIKICAgICBwYWdlY29sb3I9IiNmZmZmZmYiCiAgICAgYm9yZGVyY29sb3I9IiM2NjY2NjYiCiAgICAgYm9yZGVyb3BhY2l0eT0iMS4wIgogICAgIGlua3NjYXBlOnBhZ2VvcGFjaXR5PSIwLjAiCiAgICAgaW5rc2NhcGU6cGFnZXNoYWRvdz0iMiIKICAgICBpbmtzY2FwZTp6b29tPSIyLjgiCiAgICAgaW5rc2NhcGU6Y3g9Ii00Ljc5MDE2ODYiCiAgICAgaW5rc2NhcGU6Y3k9IjI2LjUzNDEyNCIKICAgICBpbmtzY2FwZTpkb2N1bWVudC11bml0cz0icHgiCiAgICAgaW5rc2NhcGU6Y3VycmVudC1sYXllcj0ibGF5ZXIxIgogICAgIHNob3dncmlkPSJ0cnVlIgogICAgIHNob3dib3JkZXI9ImZhbHNlIgogICAgIGZpdC1tYXJnaW4tdG9wPSIwIgogICAgIGZpdC1tYXJnaW4tbGVmdD0iMCIKICAgICBmaXQtbWFyZ2luLXJpZ2h0PSIwIgogICAgIGZpdC1tYXJnaW4tYm90dG9tPSIwIgogICAgIGlua3NjYXBlOndpbmRvdy13aWR0aD0iMTMyOSIKICAgICBpbmtzY2FwZTp3aW5kb3ctaGVpZ2h0PSI3NDQiCiAgICAgaW5rc2NhcGU6d2luZG93LXg9IjM3IgogICAgIGlua3NjYXBlOndpbmRvdy15PSIyNCIKICAgICBpbmtzY2FwZTp3aW5kb3ctbWF4aW1pemVkPSIxIiAvPgogIDxtZXRhZGF0YQogICAgIGlkPSJtZXRhZGF0YTciPgogICAgPHJkZjpSREY+CiAgICAgIDxjYzpXb3JrCiAgICAgICAgIHJkZjphYm91dD0iIj4KICAgICAgICA8ZGM6Zm9ybWF0PmltYWdlL3N2Zyt4bWw8L2RjOmZvcm1hdD4KICAgICAgICA8ZGM6dHlwZQogICAgICAgICAgIHJkZjpyZXNvdXJjZT0iaHR0cDovL3B1cmwub3JnL2RjL2RjbWl0eXBlL1N0aWxsSW1hZ2UiIC8+CiAgICAgICAgPGRjOnRpdGxlPjwvZGM6dGl0bGU+CiAgICAgIDwvY2M6V29yaz4KICAgIDwvcmRmOlJERj4KICA8L21ldGFkYXRhPgogIDxnCiAgICAgaW5rc2NhcGU6Z3JvdXBtb2RlPSJsYXllciIKICAgICBpZD0ibGF5ZXIyIgogICAgIGlua3NjYXBlOmxhYmVsPSJiYWNrZ3JvdW5kIgogICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0zLjU1MzM1MDgsNS44NTY3MDQ1KSIgLz4KICA8ZwogICAgIGlua3NjYXBlOmxhYmVsPSJMYXllciAxIgogICAgIGlua3NjYXBlOmdyb3VwbW9kZT0ibGF5ZXIiCiAgICAgaWQ9ImxheWVyMSIKICAgICB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMjM5LjY4MDcyLC03MDYuNTA1NDQpIj4KICAgIDxwYXRoCiAgICAgICBzdHlsZT0ib3BhY2l0eToxO2ZpbGw6IzU2NTY1NjtmaWxsLW9wYWNpdHk6MTtzdHJva2U6I2Q2MjcyODtzdHJva2Utd2lkdGg6MS4wMjA4NDQ4MjtzdHJva2UtbWl0ZXJsaW1pdDo0O3N0cm9rZS1kYXNoYXJyYXk6bm9uZTtzdHJva2UtZGFzaG9mZnNldDowO3N0cm9rZS1vcGFjaXR5OjAiCiAgICAgICBkPSJNIDI1IDQuMjM0Mzc1IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCA1LjUwOTc2NTYgMjMuNzI0NjA5IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCAyNSA0My4yMTQ4NDQgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDQ0LjQ5MDIzNCAyMy43MjQ2MDkgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDI1IDQuMjM0Mzc1IHogTSAyNy43NDAyMzQgMTMuOTEyMTA5IEEgMTMuMDcxNzQ0IDEzLjA3MTc0NCAwIDAgMSA0MC44MTI1IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDQwLjA1NjY0MSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMTQuNjY3OTY5IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDEzLjkxMjEwOSB6ICIKICAgICAgIGlkPSJwYXRoNDEzNiIKICAgICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKDIzOS42ODA3Miw3MDYuNTA1NDQpIiAvPgogIDwvZz4KPC9zdmc+Cg=="/>
            </div>
    
        </div>
        """
        )
        if isinstance(self.pages, list):
            for page in range(len(self.pages)):
                html += self.pages[page].render_html_str()
                if page != len(self.pages) - 1:
                    html += """<<p class="page-break-before" style="page-break-before: always" /> \n>"""
        elif isinstance(self.pages, Page):
            html += self.pages.render_html_str()

        html += """
        
        </div>
        </body>
        </html>
        """
        return html
