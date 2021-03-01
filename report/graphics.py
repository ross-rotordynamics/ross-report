import base64 as b64
from pathlib import Path

from plotly.graph_objs import Figure

import report


class CSS:
    def __init__(self, path=Path(report.__file__).parent / "style.css"):
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
    def __init__(self, text, style=""):
        super().__init__()
        self.text = text
        self.style = style

    def render_html_str(self):
        html = f"""
                <div class="row offset-2">\n
                    <div class="col-9 mt-3 mb-0 pb-0">\n
                        <p style="{self.style}"class="text-justify">\n
                            {self.text} \n
                        </p>\n
                    </div>\n
                </div>\n
               """
        return html


class Img(Content):
    """"""

    def __init__(self, path_to_img):
        self.path_to_img = Path(path_to_img)

    def render_html_str(self):
        html = ""
        return html

    def __repr__(self):
        return f"{self.path_to_img}"


class Title(Content):
    """"""

    def __init__(self, title):
        self.title = title
        self.id = title

    def render_html_str(self):
        html = f"""
        <div class="row offset-2">
            <div class="col-9">
                <h1 class="h1"><a id="{self.id}"></a>
                    {self.title}
                </h1>
            </div>
        </div>
"""

        return html


class Page:
    """"""

    def __init__(
        self,
        content=None,
    ):

        for item in content:
            assert isinstance(item, Content) or isinstance(
                item, str
            ), "Every item of content must be either content or a string"
        self.content = content
        self._titles = []
        self._figures = []
        self._tables = []
        content = []
        for item in self.content:
            if isinstance(item, str):
                content.append(Text(item))
            else:
                content.append(item)
            if isinstance(item, Title):
                self._titles.append(item)
            if isinstance(item, PlotlyFigure) or isinstance(item, Img):
                self._figures.append(item)
            if isinstance(item, Table):
                self._tables.append(item)

    def render_html_str(self, figures_list, tables_list):
        html = ""
        for item in self.content:
            if isinstance(item, PlotlyFigure):
                figure_numb = len(figures_list)
                item = PlotlyFigure(
                    item.figure, item.width, legend=item.legend, id="Figure " + f"{figure_numb}"
                )
                figures_list.append(item)
            elif isinstance(item, Table):
                table_numb = len(tables_list)
                item = Table(item.table, item.width, legend=item.legend, id="Table " + f"{table_numb}")
                tables_list.append(item)
            html += item.render_html_str()

        return html, figures_list, tables_list


class PlotlyFigure(Content):
    def __init__(self, figure, width=900, id="", legend=""):
        self.figure = figure
        self.figure.update_layout(width=width)
        self.id = id
        self.width = width
        self.legend = legend

    def render_html_str(self, legend=""):

        html = f"""
        <div style="width: {self.figure.layout["width"]}px;height: {self.figure.layout["height"]}px" class="mx-auto" id="{self.id}">\n
            {self.figure.to_html(full_html=False)}
        </div>
        """

        fig_id = self.id.split(' ')
        fig_id[-1] = str(int(fig_id[-1])+1)
        fig_id = ' '.join(fig_id)

        html += f"""
        <div style="width: {self.figure.layout["width"]}px;height: auto;text-align: center"  class="mx-auto">\n
            <p class="text-center" style="font-size:11px;">{fig_id} {'- ' + self.legend if self.legend != "" else ""}</p>
        </div>
        """

        return html


class Table(Content):
    def __init__(self, pandas_data_frame, width, id="", legend=""):
        self.table = pandas_data_frame
        self.width = width
        self.id = id
        self.legend = legend

    def render_html_str(self, legend=""):
        html = ""
        table_id = self.id.split(' ')
        table_id[-1] = str(int(table_id[-1]) + 1)
        table_id = ' '.join(table_id)

        html += f"""
                        <div style="width: {self.width}px;" class="mx-auto" id="{self.id}">\n
                            <p class="text-center" style="font-size:11px;">{table_id} {'- ' + self.legend if self.legend != "" else ""}</p>
                        </div>
                        """

        html += f"""
        <div style="width: {self.width}px;" class="mx-auto" id="{self.id}">\n
            {self.table.to_html(classes="table table-striped table-hover table-responsive")}
        </div>
        """
        html = html.replace("&amp;#", "&#")

        return html


class Listing(Content):
    def __init__(self, items):
        self.items = items

    def __str__(self):
        return str(self.render_html_str())

    def render_html_str(self):
        html = """
        <div class="offset-2 mb-4 row">
                <div class="col-9">
                    <ui>
        """

        for item in self.items:
            html += f"""             <li style = "color: #555555; text-align: -webkit-match-parent;">{item}</li>\n"""

        html += """
                </ui>
            </div>
        </div>
        """
        return html


class Link(Content):
    def __init__(self, title, href, style="", internal=True):
        self.href = href
        self.title = title
        self.internal = internal
        self.style = style

    def render_html_str(self):
        return f"""<a style="{self.style}" class="text-justify" href = "{self._internal()}{self.href}">{self.title}</a>"""

    def _internal(self):
        if self.internal:
            return "#"
        else:
            return ""

    def __str__(self):
        return self.render_html_str()

    def __repr__(self):
        return self.render_html_str()


class Layout:
    """Report Layout"""

    def __init__(
        self,
        summary=True,
        tables_list_ref=True,
        figures_list_ref=True,
        css=CSS(),
        pages=None,
        main_title="ROSS Report",
    ):
        assert (
            isinstance(pages, list) or isinstance(pages, Page) or pages is None
        ), "pages argument must be either a Page or a list of Pages."
        if isinstance(pages, list):
            for page in pages:
                assert isinstance(page, Page), "All items in page must be Pages."
            self.pages = pages
        else:
            self.pages = [pages]

        assert isinstance(css, CSS), "css must be an CSS object."
        self.css = css

        assert isinstance(main_title, str), "main_title must be a str."
        self.main_title = main_title

        self.summary = summary

        self.figures_list = []
        self.figures_list_ref = figures_list_ref

        self.tables_list = []
        self.tables_list_ref = tables_list_ref

    def __repr__(self):
        rep = {}
        if isinstance(self.pages, list):
            for page in range(len(self.pages)):
                rep[self.pages[page]] = self.pages[page]
        elif isinstance(self.pages, Page):
            rep[self.pages] = self.pages.content
        return str(rep)

    def summary_renderer(self):
        layout_titles = []
        for page in self.pages:
            for title in page._titles:
                layout_titles.append(Link(title=title.title, href=title.title))
        if self.figures_list_ref:
            layout_titles.append(Link(title="Figures List", href="Figures List"))
        if self.tables_list_ref:
            layout_titles.append(Link(title="Tables List", href="Tables List"))
        summary = f"""
            <div class="mt-4 pb-4 offset-2 row">
                <div class="col-9">
                    <h1 class="h1">
                        Summary 
                    </h1>
                </div>
            </div>
        """
        summary += str(Listing(layout_titles))
        return summary

    def figures_list_renderer(self):
        content = [Title("Figures List")]
        image_titles = []
        for figure in range(len(self.figures_list)):
            image_titles.append(
                Link(title=f"Figure {figure + 1}" + (" - "+self.figures_list[figure].legend if self.figures_list[figure].legend != "" else ""), href=f"Figure {figure}")
            )
        content.append(Listing(image_titles))

        return content

    def tables_list_renderer(self):
        content = [Title("Tables List")]
        tables_titles = []
        for table in range(len(self.tables_list)):
            tables_titles.append(
                Link(title=f"Table {table + 1}"+ (" - "+self.tables_list[table].legend if self.tables_list[table].legend != "" else ""), href=f"Table {table}")
            )
        content.append(Listing(tables_titles))

        return content

    def render_pages(self, tables_list_ref=True, figures_list_ref=True):
        html = ""

        if self.summary is True:

            if isinstance(self.pages, list):
                figures_list = []
                tables_list = []

                for page in range(len(self.pages)):
                    rendered_page = self.pages[page].render_html_str(
                        figures_list=figures_list, tables_list=tables_list
                    )

                    html += rendered_page[0]
                    figures_list = rendered_page[1]
                    tables_list = rendered_page[2]

                    if page != len(self.pages) - 1:
                        html += """<<p class="page-break-before" style="page-break-before: always" /> \n>"""

                    if figures_list_ref and page == 0:
                        self.figures_list = figures_list
                        rendered_page = Page(
                            content=self.figures_list_renderer()
                        ).render_html_str(
                            figures_list=figures_list, tables_list=tables_list
                        )
                        html += rendered_page[0]

                    if tables_list_ref and page == 0:
                        self.tables_list = tables_list
                        rendered_page = Page(
                            content=self.tables_list_renderer()
                        ).render_html_str(
                            figures_list=figures_list, tables_list=tables_list
                        )
                        html += rendered_page[0]

            elif isinstance(self.pages, Page):
                html += self.pages.render_html_str()

        return html

    def render_html_str(self):
        rendered_pages = self.render_pages(
            figures_list_ref=self.figures_list_ref, tables_list_ref=self.tables_list_ref
        )
        summary = self.summary_renderer()
        html = (
            """
                <!DOCTYPE html>
                <html lang="en">
                <head>
                <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
                    """
            + """
            <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js" integrity="sha384-DfXdz2htPH0lsSSs5nCTpuj/zy4C+OGpamoFVy38MVBnE+IbbVYUew+OrCXaRkfj" crossorigin="anonymous"></script>
            <script src="https://cdn.jsdelivr.net/npm/popper.js@1.16.1/dist/umd/popper.min.js" integrity="sha384-9/reFTGAW83EW2RDu2S0VKaIzap3H66lZH81PoYlFhbGU+6BZp6G7niu735Sk7lN" crossorigin="anonymous"></script>
            <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/js/bootstrap.min.js" integrity="sha384-B4gt1jrGC7Jh4AgTPSdUtOBvfO8shuf57BaghqFfPlYxofvL8/KUEfYiJOMMV+rV" crossorigin="anonymous"></script>
            <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js"></script>
            <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js"></script>
            <script type="text/javascript" id="MathJax-script" async
              src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js">
            </script>
            """
        )
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
                    a,li{
                        margin-top: 0;
                        margin-bottom: 1rem;
                        color: #555555;
                    }
                }
            </style>
            <meta charset="UTF-8">
            <title>ROSS Report</title>
            <link rel="icon" href="data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICAgICAgIHhtbG5zOm9zYj0iaHR0cDovL3d3dy5vcGVuc3dhdGNoYm9vay5vcmcvdXJpLzIwMDkvb3NiIgogICAgICAgIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIKICAgICAgICB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIgogICAgICAgIHhtbG5zOnJkZj0iaHR0cDovL3d3dy53My5vcmcvMTk5OS8wMi8yMi1yZGYtc3ludGF4LW5zIyIKICAgICAgICB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciCiAgICAgICAgeG1sbnM6c29kaXBvZGk9Imh0dHA6Ly9zb2RpcG9kaS5zb3VyY2Vmb3JnZS5uZXQvRFREL3NvZGlwb2RpLTAuZHRkIgogICAgICAgIHhtbG5zOmlua3NjYXBlPSJodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy9uYW1lc3BhY2VzL2lua3NjYXBlIgogICAgICAgIHdpZHRoPSIxNC4xMTA5OTFtbSIKICAgICAgICBoZWlnaHQ9IjEzLjM5MTAzNW1tIgogICAgICAgIHZpZXdCb3g9IjAgMCA0OS45OTk1NjkgNDcuNDQ4NTQ2IgogICAgICAgIGlkPSJzdmcyIgogICAgICAgIHZlcnNpb249IjEuMSIKICAgICAgICBpbmtzY2FwZTp2ZXJzaW9uPSIwLjkxIHIxMzcyNSIKICAgICAgICBzb2RpcG9kaTpkb2NuYW1lPSJyb3NzLWxvZ28uc3ZnIj4KICA8ZGVmcwogICAgIGlkPSJkZWZzNCI+CiAgICA8bGluZWFyR3JhZGllbnQKICAgICAgIGlkPSJsaW5lYXJHcmFkaWVudDQ0NjEiCiAgICAgICBvc2I6cGFpbnQ9InNvbGlkIj4KICAgICAgPHN0b3AKICAgICAgICAgc3R5bGU9InN0b3AtY29sb3I6I2Y4ZjhmODtzdG9wLW9wYWNpdHk6MTsiCiAgICAgICAgIG9mZnNldD0iMCIKICAgICAgICAgaWQ9InN0b3A0NDYzIiAvPgogICAgPC9saW5lYXJHcmFkaWVudD4KICAgIDxpbmtzY2FwZTpwZXJzcGVjdGl2ZQogICAgICAgc29kaXBvZGk6dHlwZT0iaW5rc2NhcGU6cGVyc3AzZCIKICAgICAgIGlua3NjYXBlOnZwX3g9IjAgOiA1MjYuMTgxMDQgOiAxIgogICAgICAgaW5rc2NhcGU6dnBfeT0iMCA6IDk5OS45OTk4NyA6IDAiCiAgICAgICBpbmtzY2FwZTp2cF96PSI3NDQuMDk0NDQgOiA1MjYuMTgxMDMgOiAxIgogICAgICAgaW5rc2NhcGU6cGVyc3AzZC1vcmlnaW49IjM3Mi4wNDcyIDogMzUwLjc4NzM1IDogMSIKICAgICAgIGlkPSJwZXJzcGVjdGl2ZTQxMzgiIC8+CiAgPC9kZWZzPgogIDxzb2RpcG9kaTpuYW1lZHZpZXcKICAgICBpZD0iYmFzZSIKICAgICBwYWdlY29sb3I9IiNmZmZmZmYiCiAgICAgYm9yZGVyY29sb3I9IiM2NjY2NjYiCiAgICAgYm9yZGVyb3BhY2l0eT0iMS4wIgogICAgIGlua3NjYXBlOnBhZ2VvcGFjaXR5PSIwLjAiCiAgICAgaW5rc2NhcGU6cGFnZXNoYWRvdz0iMiIKICAgICBpbmtzY2FwZTp6b29tPSIyLjgiCiAgICAgaW5rc2NhcGU6Y3g9Ii00Ljc5MDE2ODYiCiAgICAgaW5rc2NhcGU6Y3k9IjI2LjUzNDEyNCIKICAgICBpbmtzY2FwZTpkb2N1bWVudC11bml0cz0icHgiCiAgICAgaW5rc2NhcGU6Y3VycmVudC1sYXllcj0ibGF5ZXIxIgogICAgIHNob3dncmlkPSJ0cnVlIgogICAgIHNob3dib3JkZXI9ImZhbHNlIgogICAgIGZpdC1tYXJnaW4tdG9wPSIwIgogICAgIGZpdC1tYXJnaW4tbGVmdD0iMCIKICAgICBmaXQtbWFyZ2luLXJpZ2h0PSIwIgogICAgIGZpdC1tYXJnaW4tYm90dG9tPSIwIgogICAgIGlua3NjYXBlOndpbmRvdy13aWR0aD0iMTMyOSIKICAgICBpbmtzY2FwZTp3aW5kb3ctaGVpZ2h0PSI3NDQiCiAgICAgaW5rc2NhcGU6d2luZG93LXg9IjM3IgogICAgIGlua3NjYXBlOndpbmRvdy15PSIyNCIKICAgICBpbmtzY2FwZTp3aW5kb3ctbWF4aW1pemVkPSIxIiAvPgogIDxtZXRhZGF0YQogICAgIGlkPSJtZXRhZGF0YTciPgogICAgPHJkZjpSREY+CiAgICAgIDxjYzpXb3JrCiAgICAgICAgIHJkZjphYm91dD0iIj4KICAgICAgICA8ZGM6Zm9ybWF0PmltYWdlL3N2Zyt4bWw8L2RjOmZvcm1hdD4KICAgICAgICA8ZGM6dHlwZQogICAgICAgICAgIHJkZjpyZXNvdXJjZT0iaHR0cDovL3B1cmwub3JnL2RjL2RjbWl0eXBlL1N0aWxsSW1hZ2UiIC8+CiAgICAgICAgPGRjOnRpdGxlPjwvZGM6dGl0bGU+CiAgICAgIDwvY2M6V29yaz4KICAgIDwvcmRmOlJERj4KICA8L21ldGFkYXRhPgogIDxnCiAgICAgaW5rc2NhcGU6Z3JvdXBtb2RlPSJsYXllciIKICAgICBpZD0ibGF5ZXIyIgogICAgIGlua3NjYXBlOmxhYmVsPSJiYWNrZ3JvdW5kIgogICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0zLjU1MzM1MDgsNS44NTY3MDQ1KSIgLz4KICA8ZwogICAgIGlua3NjYXBlOmxhYmVsPSJMYXllciAxIgogICAgIGlua3NjYXBlOmdyb3VwbW9kZT0ibGF5ZXIiCiAgICAgaWQ9ImxheWVyMSIKICAgICB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMjM5LjY4MDcyLC03MDYuNTA1NDQpIj4KICAgIDxwYXRoCiAgICAgICBzdHlsZT0ib3BhY2l0eToxO2ZpbGw6IzU2NTY1NjtmaWxsLW9wYWNpdHk6MTtzdHJva2U6I2Q2MjcyODtzdHJva2Utd2lkdGg6MS4wMjA4NDQ4MjtzdHJva2UtbWl0ZXJsaW1pdDo0O3N0cm9rZS1kYXNoYXJyYXk6bm9uZTtzdHJva2UtZGFzaG9mZnNldDowO3N0cm9rZS1vcGFjaXR5OjAiCiAgICAgICBkPSJNIDI1IDQuMjM0Mzc1IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCA1LjUwOTc2NTYgMjMuNzI0NjA5IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCAyNSA0My4yMTQ4NDQgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDQ0LjQ5MDIzNCAyMy43MjQ2MDkgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDI1IDQuMjM0Mzc1IHogTSAyNy43NDAyMzQgMTMuOTEyMTA5IEEgMTMuMDcxNzQ0IDEzLjA3MTc0NCAwIDAgMSA0MC44MTI1IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDQwLjA1NjY0MSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMTQuNjY3OTY5IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDEzLjkxMjEwOSB6ICIKICAgICAgIGlkPSJwYXRoNDEzNiIKICAgICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKDIzOS42ODA3Miw3MDYuNTA1NDQpIiAvPgogIDwvZz4KPC9zdmc+Cg==">
        </head>
        <body>
    """
            + f"""
        <div class="container-fluid">
            
            <div class="mt-4 offset-2 row">
                <div class="col-9">
                    <h1 style= "font-size: 2.5rem" class="h1">
                        {self.main_title}
                    </h1>
                </div>

            <div class="col-3">
                <a href="https://github.com/ross-rotordynamics/ross-report">
                    <img id="ross-logo" src="data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjwhLS0gQ3JlYXRlZCB3aXRoIElua3NjYXBlIChodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy8pIC0tPgoKPHN2ZwogICAgICAgIHhtbG5zOm9zYj0iaHR0cDovL3d3dy5vcGVuc3dhdGNoYm9vay5vcmcvdXJpLzIwMDkvb3NiIgogICAgICAgIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIKICAgICAgICB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIgogICAgICAgIHhtbG5zOnJkZj0iaHR0cDovL3d3dy53My5vcmcvMTk5OS8wMi8yMi1yZGYtc3ludGF4LW5zIyIKICAgICAgICB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciCiAgICAgICAgeG1sbnM6c29kaXBvZGk9Imh0dHA6Ly9zb2RpcG9kaS5zb3VyY2Vmb3JnZS5uZXQvRFREL3NvZGlwb2RpLTAuZHRkIgogICAgICAgIHhtbG5zOmlua3NjYXBlPSJodHRwOi8vd3d3Lmlua3NjYXBlLm9yZy9uYW1lc3BhY2VzL2lua3NjYXBlIgogICAgICAgIHdpZHRoPSIxNC4xMTA5OTFtbSIKICAgICAgICBoZWlnaHQ9IjEzLjM5MTAzNW1tIgogICAgICAgIHZpZXdCb3g9IjAgMCA0OS45OTk1NjkgNDcuNDQ4NTQ2IgogICAgICAgIGlkPSJzdmcyIgogICAgICAgIHZlcnNpb249IjEuMSIKICAgICAgICBpbmtzY2FwZTp2ZXJzaW9uPSIwLjkxIHIxMzcyNSIKICAgICAgICBzb2RpcG9kaTpkb2NuYW1lPSJyb3NzLWxvZ28uc3ZnIj4KICA8ZGVmcwogICAgIGlkPSJkZWZzNCI+CiAgICA8bGluZWFyR3JhZGllbnQKICAgICAgIGlkPSJsaW5lYXJHcmFkaWVudDQ0NjEiCiAgICAgICBvc2I6cGFpbnQ9InNvbGlkIj4KICAgICAgPHN0b3AKICAgICAgICAgc3R5bGU9InN0b3AtY29sb3I6I2Y4ZjhmODtzdG9wLW9wYWNpdHk6MTsiCiAgICAgICAgIG9mZnNldD0iMCIKICAgICAgICAgaWQ9InN0b3A0NDYzIiAvPgogICAgPC9saW5lYXJHcmFkaWVudD4KICAgIDxpbmtzY2FwZTpwZXJzcGVjdGl2ZQogICAgICAgc29kaXBvZGk6dHlwZT0iaW5rc2NhcGU6cGVyc3AzZCIKICAgICAgIGlua3NjYXBlOnZwX3g9IjAgOiA1MjYuMTgxMDQgOiAxIgogICAgICAgaW5rc2NhcGU6dnBfeT0iMCA6IDk5OS45OTk4NyA6IDAiCiAgICAgICBpbmtzY2FwZTp2cF96PSI3NDQuMDk0NDQgOiA1MjYuMTgxMDMgOiAxIgogICAgICAgaW5rc2NhcGU6cGVyc3AzZC1vcmlnaW49IjM3Mi4wNDcyIDogMzUwLjc4NzM1IDogMSIKICAgICAgIGlkPSJwZXJzcGVjdGl2ZTQxMzgiIC8+CiAgPC9kZWZzPgogIDxzb2RpcG9kaTpuYW1lZHZpZXcKICAgICBpZD0iYmFzZSIKICAgICBwYWdlY29sb3I9IiNmZmZmZmYiCiAgICAgYm9yZGVyY29sb3I9IiM2NjY2NjYiCiAgICAgYm9yZGVyb3BhY2l0eT0iMS4wIgogICAgIGlua3NjYXBlOnBhZ2VvcGFjaXR5PSIwLjAiCiAgICAgaW5rc2NhcGU6cGFnZXNoYWRvdz0iMiIKICAgICBpbmtzY2FwZTp6b29tPSIyLjgiCiAgICAgaW5rc2NhcGU6Y3g9Ii00Ljc5MDE2ODYiCiAgICAgaW5rc2NhcGU6Y3k9IjI2LjUzNDEyNCIKICAgICBpbmtzY2FwZTpkb2N1bWVudC11bml0cz0icHgiCiAgICAgaW5rc2NhcGU6Y3VycmVudC1sYXllcj0ibGF5ZXIxIgogICAgIHNob3dncmlkPSJ0cnVlIgogICAgIHNob3dib3JkZXI9ImZhbHNlIgogICAgIGZpdC1tYXJnaW4tdG9wPSIwIgogICAgIGZpdC1tYXJnaW4tbGVmdD0iMCIKICAgICBmaXQtbWFyZ2luLXJpZ2h0PSIwIgogICAgIGZpdC1tYXJnaW4tYm90dG9tPSIwIgogICAgIGlua3NjYXBlOndpbmRvdy13aWR0aD0iMTMyOSIKICAgICBpbmtzY2FwZTp3aW5kb3ctaGVpZ2h0PSI3NDQiCiAgICAgaW5rc2NhcGU6d2luZG93LXg9IjM3IgogICAgIGlua3NjYXBlOndpbmRvdy15PSIyNCIKICAgICBpbmtzY2FwZTp3aW5kb3ctbWF4aW1pemVkPSIxIiAvPgogIDxtZXRhZGF0YQogICAgIGlkPSJtZXRhZGF0YTciPgogICAgPHJkZjpSREY+CiAgICAgIDxjYzpXb3JrCiAgICAgICAgIHJkZjphYm91dD0iIj4KICAgICAgICA8ZGM6Zm9ybWF0PmltYWdlL3N2Zyt4bWw8L2RjOmZvcm1hdD4KICAgICAgICA8ZGM6dHlwZQogICAgICAgICAgIHJkZjpyZXNvdXJjZT0iaHR0cDovL3B1cmwub3JnL2RjL2RjbWl0eXBlL1N0aWxsSW1hZ2UiIC8+CiAgICAgICAgPGRjOnRpdGxlPjwvZGM6dGl0bGU+CiAgICAgIDwvY2M6V29yaz4KICAgIDwvcmRmOlJERj4KICA8L21ldGFkYXRhPgogIDxnCiAgICAgaW5rc2NhcGU6Z3JvdXBtb2RlPSJsYXllciIKICAgICBpZD0ibGF5ZXIyIgogICAgIGlua3NjYXBlOmxhYmVsPSJiYWNrZ3JvdW5kIgogICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKC0zLjU1MzM1MDgsNS44NTY3MDQ1KSIgLz4KICA8ZwogICAgIGlua3NjYXBlOmxhYmVsPSJMYXllciAxIgogICAgIGlua3NjYXBlOmdyb3VwbW9kZT0ibGF5ZXIiCiAgICAgaWQ9ImxheWVyMSIKICAgICB0cmFuc2Zvcm09InRyYW5zbGF0ZSgtMjM5LjY4MDcyLC03MDYuNTA1NDQpIj4KICAgIDxwYXRoCiAgICAgICBzdHlsZT0ib3BhY2l0eToxO2ZpbGw6IzU2NTY1NjtmaWxsLW9wYWNpdHk6MTtzdHJva2U6I2Q2MjcyODtzdHJva2Utd2lkdGg6MS4wMjA4NDQ4MjtzdHJva2UtbWl0ZXJsaW1pdDo0O3N0cm9rZS1kYXNoYXJyYXk6bm9uZTtzdHJva2UtZGFzaG9mZnNldDowO3N0cm9rZS1vcGFjaXR5OjAiCiAgICAgICBkPSJNIDI1IDQuMjM0Mzc1IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCA1LjUwOTc2NTYgMjMuNzI0NjA5IEEgMTkuNDg5NTgxIDE5LjQ4OTU4MSAwIDAgMCAyNSA0My4yMTQ4NDQgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDQ0LjQ5MDIzNCAyMy43MjQ2MDkgQSAxOS40ODk1ODEgMTkuNDg5NTgxIDAgMCAwIDI1IDQuMjM0Mzc1IHogTSAyNy43NDAyMzQgMTMuOTEyMTA5IEEgMTMuMDcxNzQ0IDEzLjA3MTc0NCAwIDAgMSA0MC44MTI1IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDQwLjA1NjY0MSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMTQuNjY3OTY5IDI2Ljk4NDM3NSBBIDEzLjA3MTc0NCAxMy4wNzE3NDQgMCAwIDEgMjcuNzQwMjM0IDEzLjkxMjEwOSB6ICIKICAgICAgIGlkPSJwYXRoNDEzNiIKICAgICAgIHRyYW5zZm9ybT0idHJhbnNsYXRlKDIzOS42ODA3Miw3MDYuNTA1NDQpIiAvPgogIDwvZz4KPC9zdmc+Cg=="/>
                </a>
            </div>
    
        </div>
        """
        )

        html += summary
        html += rendered_pages
        html += """
        
        </div>
        </body>
        </html>
        """
        return html
