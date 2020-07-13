# fmt: off
from collections.abc import Iterable
from copy import copy, deepcopy

import numpy as np
import pandas as pd
from plotly import express as px
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema

from ross.bearing_seal_element import BearingElement, SealElement
from ross.disk_element import DiskElement
from ross.materials import steel
from ross.plotly_theme import tableau_colors
from ross.rotor_assembly import Rotor
from ross.shaft_element import ShaftElement

# fmt: on

# set Plotly palette of colors
colors1 = px.colors.qualitative.Dark24

__all__ = ["Report", "report_example"]


class Report:
    """Report according to standard analysis.

    - Perform unbalance response
    - Perform Stability_level1 analysis
    - Apply Level 1 Screening Criteria
    - Perform Stability_level2 analysis

    Parameters
    ----------
    rotor : object
        A rotor built from rotor_assembly.
    config : object
        An instance of class Config() with the analyses configurations.

    Attributes
    ----------
    rotor_type: str
        Defines if the rotor is between bearings or overhung
    disk_nodes: list
        List of disk between bearings or overhung (depending on the
        rotor type)

    Returns
    -------
    A Report object

    Examples
    --------
    >>> import ross as rs
    >>> rotor = rs.rotor_example()
    >>> config = rs.Config()
    >>> report = rs.Report(rotor=rotor, config=config)
    >>> report.rotor_type
    'between_bearings'
    """

    def __init__(self, rotor, config):
        self.rotor = rotor

        # check if rotor is between bearings, single or double overhung
        # fmt: off
        if(
            all(i > min(rotor.df_bearings["n"]) for i in rotor.df_disks["n"]) and
            all(i < max(rotor.df_bearings["n"]) for i in rotor.df_disks["n"])
        ):
            rotor_type = "between_bearings"
            disk_nodes = [
                i for i in rotor.df_disks["n"] if(
                    i > min(rotor.df_bearings["n"]) and
                    i < max(rotor.df_bearings["n"])
                )
            ]
        elif(
            any(i < min(rotor.df_bearings["n"]) for i in rotor.df_disks["n"]) and
            all(i < max(rotor.df_bearings["n"]) for i in rotor.df_disks["n"])
        ):
            rotor_type = "single_overhung_l"
            disk_nodes = [
                i for i in rotor.df_disks["n"] if i < min(rotor.df_bearings["n"])
            ]
        elif(
            all(i > min(rotor.df_bearings["n"]) for i in rotor.df_disks["n"]) and
            any(i > max(rotor.df_bearings["n"]) for i in rotor.df_disks["n"])
        ):
            rotor_type = "single_overhung_r"
            disk_nodes = [
                i for i in rotor.df_disks["n"] if i > max(rotor.df_bearings["n"])
            ]
        elif(
            any(i < min(rotor.df_bearings["n"]) for i in rotor.df_disks["n"]) and
            any(i > max(rotor.df_bearings["n"]) for i in rotor.df_disks["n"])
        ):
            rotor_type = "double_overhung"
            disk_nodes = [
                i for i in rotor.df_disks["n"] if(
                        i < min(rotor.df_bearings["n"]) or
                        i > max(rotor.df_bearings["n"])
                )
            ]
        # fmt: on

        self.rotor_type = rotor_type
        self.disk_nodes = disk_nodes

        machine_options = ["compressor", "turbine", "axial_flow"]
        machine_type = config.rotor_properties.rotor_id.type
        if machine_type not in machine_options:
            raise ValueError(
                "rotor_id.type is set to {}. Please choose between {}.".format(
                    machine_type, machine_options
                )
            )

        if config.rotor_properties.rotor_id.tag is None:
            config.update_config(rotor_properties=dict(rotor_id=dict(tag=rotor.tag)))

        self.tag = config.rotor_properties.rotor_id.tag
        self.config = config

    def _rotor_instance(self, rotor, bearing_list):
        """Build an instance of an auxiliary rotor with different bearing clearances.

        Parameters
        ----------
        rotor : object
            A rotor built from rotor_assembly.
        bearing_list : list
            List with the bearing elements.

        Returns
        -------
        aux_rotor : Rotor.object
            Returns a rotor object copy with different bearing clearance.

        Example
        -------
        >>> import ross as rs
        >>> stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
        >>> damp = [2.8e3, 2.7e3, 2.6e3, 2.5e3]
        >>> freq = [400, 800, 1200, 1600]
        >>> bearing0 = rs.BearingElement(0, kxx=stfx, cxx=damp, frequency=freq)
        >>> bearing1 = rs.BearingElement(6, kxx=stfx, cxx=damp, frequency=freq)
        >>> bearings = [bearing0, bearing1]
        >>> rotor = rs.rotor_example()
        >>> report = rs.report_example()
        >>> aux_rotor = report._rotor_instance(rotor, bearings)
        """
        sh_elm = rotor.shaft_elements
        dk_elm = rotor.disk_elements
        pm_elm = rotor.point_mass_elements
        tag = rotor.tag

        aux_rotor = Rotor(sh_elm, dk_elm, bearing_list, pm_elm, tag=tag)

        return aux_rotor

    def run_report(self):
        """Run rotordynamics report.

        This method runs the rotordynamics analyses and prepare the results to
        generate the PDF report.

        Returns
        -------
        fig_ucs : list
            List with undamped critical speed map figures.
        fig_mode_shape : list
            List with mode shape figures.
        fig_unbalance : list
            List with unbalance response figures.
        df_unbalance : dataframe
            Dataframe for the unbalance response informations.
        fig_a_lvl1 : list
            List with "Applied Cross-Coupled Stiffness" (stability level 1) figures.
        fig_b_lvl1 : list
            List with "CSR vs. Mean Gas Density" (stability level 1) figures.
        df_lvl2 : dataframe
            Dataframe for the stability level 2 informations.
        summaries : pd.Dataframe
            Dataframes with a summary of stability level 1 and 2 analyses.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> # to run the report analysis, use:
        >>> # report.run_report()
        """
        fig_ucs = []
        fig_mode_shape = []
        fig_unbalance = []
        fig_a_lvl1 = []
        fig_b_lvl1 = []
        df_unbalance = []
        summaries = []

        rotor0 = self.rotor

        # static analysis
        static = rotor0.run_static()
        fig_static = [
            static.plot_free_body_diagram,
            static.plot_deformation,
            static.plot_shearing_force,
            static.plot_bending_moment,
        ]

        for k, bearings in self.config.bearings.__dict__.items():
            if not isinstance(bearings, Iterable):
                raise ValueError(
                    "{} option must be a list of bearing elements".format(k)
                )

            self.rotor = self._rotor_instance(rotor0, bearings)

            # undamped critical speed map
            fig_ucs.append(
                self._plot_ucs(
                    stiffness_range=self.config.plot_ucs.stiffness_range,
                    num=self.config.plot_ucs.num,
                    num_modes=self.config.plot_ucs.num_modes,
                    synchronous=self.config.plot_ucs.synchronous,
                )
            )

            for i, mode in enumerate([1, 3]):
                # mode shape figures
                fig_mode_shape.append(self._mode_shape(mode))

                # unbalance response figures and dataframe
                fig, _dict = self._unbalance_response(mode)
                fig_unbalance.append(fig)
                df = pd.DataFrame(_dict).astype(object)
                df_unbalance.append(df)

            # stability level 1 figures
            figs = self._stability_level_1()
            fig_a_lvl1.append(figs[0])
            fig_b_lvl1.append(figs[1])

            # stability level 2 dataframe
            df_lvl2 = self._stability_level_2()

            # API summary tables
            summaries.append(self._summary())

        df_unbalance = pd.concat(df_unbalance)

        self.rotor = rotor0

        return {
            "fig_ucs": fig_ucs,
            "fig_mode_shape": fig_mode_shape,
            "fig_unbalance": fig_unbalance,
            "df_unbalance": df_unbalance,
            "fig_a_lvl1": fig_a_lvl1,
            "fig_b_lvl1": fig_b_lvl1,
            "df_lvl2": df_lvl2,
            "summaries": summaries,
        }

    def assets_prep(self, dashboard_data):

        # Rotor figure
        plot_rotor = self.rotor.plot_rotor()
        plot_rotor.layout["height"] = 350
        plot_rotor.layout["width"] = 750

        plot_rotor.layout["xaxis"]["autorange"] = True
        plot_rotor.layout["yaxis"]["autorange"] = True
        plot_rotor.update_layout()

        # Undamped critical speed figure
        ucs_fig = make_subplots(
            cols=2,
            rows=1,
            x_title="Bearing Stiffness",
            y_title="Critical Speed",
            subplot_titles=["Undamped Critical Speed Map", "", "", ""],
        )

        ucs_row = 1
        ucs_col = 1
        for fig in dashboard_data["fig_ucs"]:
            for j in fig.data:
                ucs_fig.add_trace(j, row=ucs_row, col=ucs_col)
            if ucs_col == 1:
                ucs_col = 2

        ucs_fig.layout["height"] = 350
        ucs_fig.layout["width"] = 750
        ucs_fig.update_layout()

        # Undamped mode shape
        mode_fig = make_subplots(
            cols=2,
            rows=2,
            x_title="Rotor length",
            y_title="Non dimensional deformation",
            subplot_titles=["Undamped Mode Shape", "", "", ""],
        )
        mode_row = 1
        mode_col = 1

        for fig in dashboard_data["fig_mode_shape"]:
            for j in fig.data:
                mode_fig.add_trace(j, row=mode_row, col=mode_col)
            if mode_col == 1:
                mode_col = 2
            else:
                mode_row = 2
                mode_col = 1

        mode_fig.layout["height"] = 750
        mode_fig.layout["width"] = 750
        mode_fig.update_layout()

        figs = plot_rotor, ucs_fig, mode_fig
        html_figs = []
        for fig in figs:
            html_figs.append(fig.to_html(full_html=False))

        return {"figs": figs, "html_figs": html_figs}

    def _plot_ucs(self):
        """Plot undamped critical speed map.

        This method will plot the undamped critical speed map for a given range
        of stiffness values. If the range is not provided, the bearing
        stiffness at rated speed will be used to create a range.

        Returns
        -------
        fig : Plotly graph_objects.Figure()
            The figure object with the plot.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> fig = report._plot_ucs()
        """
        fig = self.rotor.plot_ucs(
            stiffness_range=self.config.plot_ucs.stiffness_range,
            num_modes=self.config.plot_ucs.num_modes,
            num=self.config.plot_ucs.num,
            synchronous=self.config.plot_ucs.synchronous,
        )

        return fig

    def _static_forces(self):
        """Calculate the bearing reaction forces.

        Returns
        -------
        Fb : list
            Bearing reaction forces.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> report._static_forces()
        array([44.09320349, 44.09320349])
        """
        # get reaction forces on bearings
        self.rotor.run_static()
        Fb = list(self.rotor.bearing_forces_nodal.values())
        Fb = np.array(Fb) / 9.8065

        return Fb

    def _unbalance_forces(self, mode):
        """Calculate the unbalance forces.

        The unbalance forces are calculated base on the rotor type:
            between_bearings :
                The unbalance forces derives from the reaction bearing forces.
            single_overung_l :
                The unbalance forces derives from the disk's masses on the
                shaft left end.
            single_overung_r :
                The unbalance forces derives from the disk's masses on the
                shaft right end.
            double_overung :
                The unbalance forces derives from the disk's masses on the
                shaft left and right ends.

        Parameters
        ----------
        mode : int
            n'th mode shape.

        Returns
        -------
        force : list
            Unbalancing forces.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> report._unbalance_forces(mode=0)
        array([0.04479869])
        """
        if self.config.rotor_properties.rotor_speeds.unit == "rpm":
            N = self.config.rotor_properties.rotor_speeds.max_speed
        elif self.config.rotor_properties.rotor_speeds.unit == "rad/s":
            N = self.config.rotor_properties.rotor_speeds.max_speed * 30 / np.pi

        if mode > 3:
            raise ValueError(
                "This module calculates only the response for the first "
                "two backward and forward modes. "
            )

        # get reaction forces on bearings
        if self.rotor_type == "between_bearings":
            Fb = self._static_forces()
            if mode == 0 or mode == 1:
                force = [max(6350e-6 * np.sum(Fb) / N, 254e-6 * np.sum(Fb))]

            if mode == 2 or mode == 3:
                force = [max(6350e-6 * f / N, 254e-6 * f) for f in Fb]

        # get disk masses
        elif self.rotor_type == "single_overhung_l":
            Wd = [
                disk.m
                for disk in self.rotor.disk_elements
                if disk.n < min(self.rotor.df_bearings["n"])
            ]
            Ws = [
                sh.m
                for sh in self.rotor.shaft_elements
                if sh.n_l < min(self.rotor.df_bearings["n"])
            ]
            W3 = np.sum(Wd + Ws)

            force = [6350e-6 * W3 / N]

        elif self.rotor_type == "single_overhung_r":
            Wd = [
                disk.m
                for disk in self.rotor.disk_elements
                if disk.n > max(self.rotor.df_bearings["n"])
            ]
            Ws = [
                sh.m
                for sh in self.rotor.shaft_elements
                if sh.n_r > max(self.rotor.df_bearings["n"])
            ]
            W3 = np.sum(Wd + Ws)

            force = [6350e-6 * W3 / N]

        elif self.rotor_type == "double_overhung":
            Wd_l = [
                disk.m
                for disk in self.rotor.disk_elements
                if disk.n < min(self.rotor.df_bearings["n"])
            ]
            Ws_l = [
                sh.m
                for sh in self.rotor.shaft_elements
                if sh.n_l < min(self.rotor.df_bearings["n"])
            ]
            Wd_r = [
                disk.m
                for disk in self.rotor.disk_elements
                if disk.n > max(self.rotor.df_bearings["n"])
            ]
            Ws_r = [
                sh.m
                for sh in self.rotor.shaft_elements
                if sh.n_r > max(self.rotor.df_bearings["n"])
            ]
            W3 = np.array([np.sum(Wd_l + Ws_l), np.sum(Wd_r + Ws_r)])

            force = 6350e-6 * W3 / N

        force = 2 * np.array(force)

        return force

    def _unbalance_response(self, mode):
        """Evaluate the unbalance response for the rotor.

        This analysis takes the critical speeds of interest, calculates the
        position and weight of the required unbalance and performs the analysis
        including:
         - Check if vibration at MCS is below the limit with the applied weight;
         - Check if the clearances are ok if the vibration deteriorate to the
         limit level;

        Parameters
        ----------
        mode : int
            n'th mode shape.
        samples : int
            Number of samples to generate de frequency range.

        Returns
        -------
        subplots : Plotly graph_objects.make_subplots()
            Plotly figure with Amplitude vs Frequency and Phase vs Frequency plots.
        unbalance_dict : dict
            A dictionary with information about simulation parameters to be
            displayed in the report. The dictionary contains:
                - Mode number;
                - Critical frequencies;
                - Amplification factors;
                - Separation margins (actual and required);
                - Unbalance stations;
                - Unbalance weights;
                - Unbalance phases;

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> fig, unbalance_dict = report._unbalance_response(mode=0)
        """
        maxspeed = self.config.rotor_properties.rotor_speeds.max_speed
        minspeed = self.config.rotor_properties.rotor_speeds.min_speed
        speed_factor = self.config.rotor_properties.rotor_speeds.speed_factor

        freq_range = self.config.run_unbalance_response.frequency_range
        modes = self.config.run_unbalance_response.modes
        cluster_points = self.config.run_unbalance_response.cluster_points
        num_modes = self.config.run_unbalance_response.num_modes
        num_points = self.config.run_unbalance_response.num_points
        rtol = self.config.run_unbalance_response.rtol

        if freq_range is None and not cluster_points:
            freq_range = np.linspace(0, speed_factor * maxspeed, 201)

        # returns de nodes where forces will be applied
        self._mode_shape(mode)
        node_min = self.node_min
        node_max = self.node_max
        nodes = [int(node) for sub_nodes in [node_min, node_max] for node in sub_nodes]

        _magnitude = self._unbalance_forces(mode)

        phase = []
        phase_angle = 0
        for node in nodes:
            phase.append(phase_angle)
            phase_angle += np.pi

        # fmt: off
        response = self.rotor.run_unbalance_response(
            nodes, _magnitude, phase, freq_range, modes,
            cluster_points, num_modes, num_points, rtol,
        )
        # fmt: on

        mag = response.magnitude

        plot = make_subplots(
            rows=2, cols=2, specs=[[{}, {"type": "polar", "rowspan": 2}], [{}, None]]
        )

        probe_nodes = self.config.run_unbalance_response.probes.node
        probe_orien = self.config.run_unbalance_response.probes.orientation
        unbalance_dict = {
            "probe {}".format(i + 1): None for i in range(len(probe_nodes))
        }

        k = 0
        for n, o in zip(probe_nodes, probe_orien):
            dof = self.rotor.number_dof * n + 1
            fig = response.plot(dof)

            for j, data in enumerate(fig["data"]):
                data.line.color = list(tableau_colors.keys())[k]
                data.marker.color = list(tableau_colors.keys())[k]
                data.name = "Probe Node {}".format(n)
                data.legendgroup = "Probe Node {}".format(n)
                data.showlegend = True if j == 0 else False
                plot.add_trace(data)

            plot.update_layout(fig.layout)

            _dict = {
                "Probe node": [n],
                "Probe orientation": [o],
                "Critical frequencies": [],
                "Amplification factor": [],
                "Separation margin - ACTUAL": [],
                "Separation margin - REQUIRED": [],
                "Unbalance station(s)": nodes,
                "Unbalance weight(s)": _magnitude.tolist(),
                "Unbalance phase(s)": phase,
            }

            magnitude = mag[dof]
            idx_max = argrelextrema(magnitude, np.greater)[0].tolist()
            wn = freq_range[idx_max]

            for i, peak in enumerate(magnitude[idx_max]):
                peak_n = 0.707 * peak
                peak_aux = np.linspace(peak_n, peak_n, len(freq_range))

                idx = np.argwhere(np.diff(np.sign(peak_aux - magnitude))).flatten()
                idx = np.sort(np.append(idx, idx_max[i]))

                # if speed range is not long enough to catch the magnitudes
                try:
                    idx_aux = [
                        list(idx).index(idx_max[i]) - 1,
                        list(idx).index(idx_max[i]) + 1,
                    ]
                    idx = idx[idx_aux]
                except IndexError:
                    idx = [list(idx).index(idx_max[i]) - 1, len(freq_range) - 1]

                # Amplification Factor (AF) - API684 - SP6.8.2.1
                AF = wn[i] / (freq_range[idx[1]] - freq_range[idx[0]])

                # Separation Margin (SM) - API684 - SP6.8.2.10
                if AF > 2.5 and wn[i] < minspeed:
                    SM = min([16, 17 * (1 - 1 / (AF - 1.5))]) / 100
                    SMspeed = wn[i] * (1 + SM)
                    SM_ref = (minspeed - wn[i]) / wn[i]

                elif AF > 2.5 and wn[i] > maxspeed:
                    SM = min([26, 10 + 17 * (1 - 1 / (AF - 1.5))]) / 100
                    SMspeed = wn[i] * (1 - SM)
                    SM_ref = (wn[i] - maxspeed) / maxspeed

                else:
                    SM = None
                    SM_ref = None
                    SMspeed = None

                _dict["Amplification factor"].append(AF)
                _dict["Separation margin - ACTUAL"].append(SM)
                _dict["Separation margin - REQUIRED"].append(SM_ref)
                _dict["Critical frequencies"].append(wn[i])

            unbalance_dict["probe {}".format(k + 1)] = _dict

            # amplitude limit (A1) - API684 - SP6.8.2.11
            A1 = 25.4 * np.sqrt(12000 / (30 * maxspeed / np.pi)) * 1e-6

            Amax = max(mag[dof])

            # Scale Factor (Scc) - API684 - SP6.8.2.11 / API617 - 4.8.2.11
            Scc = max(A1 / Amax, 0.5)
            Scc = min(Scc, 6.0)

            customdata = [minspeed, maxspeed]
            plot.add_trace(
                go.Scatter(
                    x=[minspeed, maxspeed, maxspeed, minspeed, minspeed],
                    y=[0, 0, max(mag[dof]), max(mag[dof]), 0],
                    customdata=customdata * 5,
                    mode="lines",
                    opacity=0.3,
                    fill="toself",
                    fillcolor=tableau_colors["green"],
                    line=dict(width=1.5, color=tableau_colors["green"]),
                    name="Operation Speed Range",
                    legendgroup="Operation Speed Range",
                    hoveron="points+fills",
                    showlegend=True if i == 0 else False,
                    hoverlabel=dict(bgcolor=tableau_colors["green"]),
                    hovertemplate=(
                        f"<b>min. speed: {customdata[0]:.1f}</b><br>"
                        + f"<b>max. speed: {customdata[1]:.1f}</b>"
                    ),
                ),
                row=1,
                col=1,
            )
            plot.add_trace(
                go.Scatter(
                    x=[minspeed, maxspeed],
                    y=[A1, A1],
                    mode="lines",
                    line=dict(width=2.0, color="black", dash="dashdot"),
                    name="Av1 - Mechanical test vibration limit",
                    showlegend=True if i == 0 else False,
                    hoverinfo="none",
                ),
                row=1,
                col=1,
            )
            plot.add_annotation(
                x=(minspeed + maxspeed) / 2,
                y=A1,
                axref="x",
                ayref="y",
                xshift=0,
                yshift=10,
                text="<b>Av1</b>",
                font=dict(size=18),
                showarrow=False,
                row=1,
                col=1,
            )

            k += 1

        return plot, unbalance_dict

    def _mode_shape(self, mode):
        """Evaluate the mode shapes for the rotor.

        This analysis presents the vibration mode for each critical speed.
        The importance is to locate the critical node, where the displacement
        is the greatest, then apply loads for unbalance response (stability
        level 1)

        Parameters
        ----------
        mode : int
            the n'th vibration mode

        Attributes
        ----------
        node_min : int
            Nodes where the maximum displacements occur
        node_max : int
            Nodes where the minimum displacements occur

        Returns
        -------
        fig : Plotly graph_objects.Figure()
            The figure object with the plot.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> fig = report._mode_shape(mode=0)
        >>> report.node_min
        array([], dtype=float64)
        >>> report.node_max
        array([3.])
        """
        nodes_pos = self.rotor.nodes_pos
        df_bearings = self.rotor.df_bearings
        df_disks = self.rotor.df_disks
        speed = self.config.rotor_properties.rotor_speeds.oper_speed

        modal = self.rotor.run_modal(speed=speed)
        xn, yn, zn, xc, yc, zc_pos, nn = modal.calc_mode_shape(mode=mode)

        # reduce 3D view to 2D view
        vn = np.zeros(len(zn))
        for i in range(len(zn)):
            theta = np.arctan(xn[i] / yn[i])
            vn[i] = xn[i] * np.sin(theta) + yn[i] * np.cos(theta)

        # remove repetitive values from zn and vn
        idx_remove = []
        for i in range(1, len(zn)):
            if zn[i] == zn[i - 1]:
                idx_remove.append(i)
        zn = np.delete(zn, idx_remove)
        vn = np.delete(vn, idx_remove)

        node_min = np.array([])
        node_max = np.array([])

        if self.rotor_type == "between_bearings":

            aux_idx_max = argrelextrema(vn, np.greater)[0].tolist()
            aux_idx_min = argrelextrema(vn, np.less)[0].tolist()

            # verification of rigid modes
            if len(aux_idx_max) == 0 and len(aux_idx_min) == 0:
                idx_max = np.argmax(vn)
                idx_min = np.argmin(vn)

                # corrects the index by the removed points
                for i in idx_remove:
                    if idx_min > i:
                        idx_min += 1
                    if idx_max > i:
                        idx_max += 1
                node_max = np.round(np.array([idx_max]) / nn)
                node_min = np.round(np.array([idx_min]) / nn)

            if len(aux_idx_min) != 0:
                idx_min = np.where(vn == min(vn[aux_idx_min]))[0].tolist()

                # corrects the index by the removed points
                for i in idx_remove:
                    if idx_min[0] > i:
                        idx_min[0] += 1
                node_min = np.round(np.array(idx_min) / nn)

            if len(aux_idx_max) != 0:
                idx_max = np.where(vn == max(vn[aux_idx_max]))[0].tolist()

                # corrects the index by the removed points
                for i in idx_remove:
                    if idx_max[0] > i:
                        idx_max[0] += 1
                node_max = np.round(np.array(idx_max) / nn)

        elif self.rotor_type == "double_overhung":
            node_max = [max(df_disks["n"])]
            node_min = [min(df_disks["n"])]

        elif self.rotor_type == "single_overhung_l":
            node_min = [min(df_disks["n"])]

        elif self.rotor_type == "single_overhung_r":
            node_max = [max(df_disks["n"])]

        nodes_pos = np.array(nodes_pos)
        rpm_speed = (30 / np.pi) * modal.wn[mode]

        self.node_min = node_min
        self.node_max = node_max

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=zn,
                y=vn,
                mode="lines",
                line=dict(width=4, color=colors1[3]),
                name="<b>Mode {}</b><br><b>Speed = {:.1f} RPM</b>".format(
                    mode, rpm_speed
                ),
                hovertemplate="Axial position: %{x:.2f}<br>Deformation: %{y:.2f}",
            )
        )
        fig.add_trace(
            go.Scatter(
                x=nodes_pos,
                y=np.zeros(len(nodes_pos)),
                mode="lines",
                line=dict(width=4, color=colors1[5], dash="dashdot"),
                name="centerline",
                hoverinfo="none",
                showlegend=False,
            )
        )
        fig.add_trace(
            go.Scatter(
                x=nodes_pos[df_bearings["n"]],
                y=np.zeros(len(df_bearings)),
                mode="markers",
                marker=dict(size=12, color=colors1[5]),
                name="bearing_node",
                showlegend=False,
                hovertemplate="Bearing Position: %{x:.2f}",
            )
        )

        pos0 = nodes_pos[min(df_bearings["n"])]
        pos1 = nodes_pos[max(df_bearings["n"])]
        fig.add_annotation(
            x=np.mean(nodes_pos[df_bearings["n"]]),
            y=0,
            axref="x",
            ayref="y",
            xshift=0,
            yshift=20,
            text="<b>Bearing Span = {:.2f}</b>".format(pos1 - pos0),
            font=dict(size=18),
            showarrow=False,
        )

        for node in nodes_pos[df_bearings["n"]]:
            fig.add_trace(
                go.Scatter(
                    x=[node, node],
                    y=[-2, 2],
                    mode="lines",
                    line=dict(width=2.5, color=colors1[5], dash="dash"),
                    name="Span",
                    legendgroup="Span",
                    hoverinfo="none",
                    showlegend=False,
                )
            )

        fig.update_xaxes(
            title_text="<b>Rotor lenght</b>",
            title_font=dict(family="Arial", size=20),
            tickfont=dict(size=16),
            gridcolor="lightgray",
            showline=True,
            linewidth=2.5,
            linecolor="black",
            mirror=True,
        )
        fig.update_yaxes(
            title_text="<b>Non dimensional deformation</b>",
            title_font=dict(family="Arial", size=20),
            tickfont=dict(size=16),
            range=[-2, 2],
            gridcolor="lightgray",
            showline=True,
            linewidth=2.5,
            linecolor="black",
            mirror=True,
        )
        fig.update_layout(
            width=1200,
            height=900,
            plot_bgcolor="white",
            hoverlabel_align="right",
            title=dict(
                text="<b>Undamped Mode Shape</b>".format(node), font=dict(size=20)
            ),
        )

        return fig              

    def stability_level_1(self, D, H, HP, oper_speed, RHO_ratio, RHOs, RHOd, unit="m"):
        """Stability analysis level 1.

        This analysis consider a anticipated cross coupling QA based on
        conditions at the normal operating point and the cross-coupling
        required to produce a zero log decrement, Q0.

        Components such as seals and impellers are not considered in this
        analysis.

        Parameters
        ----------
        D: list
            Impeller diameter, m (in.),
            Blade pitch diameter, m (in.),
        H: list
            Minimum diffuser width per impeller, m (in.),
            Effective blade height, m (in.),
        HP: list
            Rated power per stage/impeller, W (HP),
        oper_speed: float
            Operating speed, rpm,
        RHO_ratio: list
            Density ratio between the discharge gas density and the suction
            gas density per impeller (RHO_discharge / RHO_suction),
            kg/m3 (lbm/in.3),
        RHOs: float
            Suction gas density in the first stage, kg/m3 (lbm/in.3).
        RHOd: float
            Discharge gas density in the last stage, kg/m3 (lbm/in.3),
        unit: str, optional
            Adopted unit system. Options are "m" (meter) and "in" (inch)
            Default is "m"

        Attributes
        ----------
        condition: bool
            False: Stability Level 1 satisfies the analysis;
            True: Stability Level 2 is required.

        Return
        ------
        fig1 : Plotly graph_objects.Figure()
            Applied Cross-Coupled Stiffness vs. Log Decrement plot.
        fig2 : Plotly graph_objects.Figure()
            CSR vs. Mean Gas Density plot.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> fig1, fig2 = report.stability_level_1(D=[0.35, 0.35],
        ...                                       H=[0.08, 0.08],
        ...                                       HP=[10000, 10000],
        ...                                       RHO_ratio=[1.11, 1.14],
        ...                                       RHOd=30.45,
        ...                                       RHOs=37.65,
        ...                                       oper_speed=1000.0)
        >>> report.Qa
        23022.32142857143
        """
        steps = 11
        if unit == "m":
            C = 9.55
        elif unit == "in":
            C = 63.0
        else:
            raise TypeError("choose between meters (m) or inches (in)")

        if len(D) != len(H):
            raise Exception("length of D must be the same of H")

        Qa = 0.0
        cross_coupled_array = np.array([])
        # Qa - Anticipated cross-coupling for compressors - API 684 - SP6.8.5.6
        if self.machine_type == "compressor":
            Bc = 3.0
            Dc, Hc = D, H
            for i, disk in enumerate(self.rotor.disk_elements):
                if disk.n in self.disk_nodes:
                    qi = HP[i] * Bc * C * RHO_ratio[i] / (Dc[i] * Hc[i] * oper_speed)
                    Qi = np.linspace(0, 10 * qi, steps)
                    cross_coupled_array = np.append(cross_coupled_array, Qi)
                    Qa += qi

        # Qa - Anticipated cross-coupling for turbines - API 684 - SP6.8.5.6
        if self.machine_type == "turbine" or self.machine_type == "axial_flow":
            Bt = 1.5
            Dt, Ht = D, H
            for i, disk in enumerate(self.rotor.disk_elements):
                if disk.n in self.disk_nodes:
                    qi = (HP[i] * Bt * C) / (Dt[i] * Ht[i] * oper_speed)
                    Qi = np.linspace(0, 10 * qi, steps)
                    cross_coupled_array = np.append(cross_coupled_array, Qi)
                    Qa += qi

        # Defining cross-coupling range to 10*Qa - API 684 - SP6.8.5.8
        Qi = np.linspace(0, 10 * Qa, steps)
        cross_coupled_array = np.append(cross_coupled_array, Qi)
        cross_coupled_array = cross_coupled_array.reshape(
            [len(self.disk_nodes) + 1, steps]
        ).T

        log_dec = np.zeros(len(cross_coupled_array))

        # remove disks and seals from the rotor model
        bearing_list = [
            copy(b)
            for b in self.rotor.bearing_elements
            if not isinstance(b, SealElement)
        ]

        # Applying cross-coupling on rotor mid-span
        if self.rotor_type == "between_bearings":
            for i, Q in enumerate(cross_coupled_array[:, -1]):
                bearings = [copy(b) for b in bearing_list]

                # cross-coupling introduced at the rotor mid-span
                n = np.round(np.mean(self.rotor.nodes))
                cross_coupling = BearingElement(n=int(n), kxx=0, cxx=0, kxy=Q, kyx=-Q)
                bearings.append(cross_coupling)

                aux_rotor = Rotor(
                    shaft_elements=self.rotor.shaft_elements,
                    disk_elements=[],
                    bearing_elements=bearings,
                    rated_w=self.rotor.rated_w,
                )
                modal = aux_rotor.run_modal(speed=oper_speed * np.pi / 30)
                non_backward = modal.whirl_direction() != "Backward"
                log_dec[i] = modal.log_dec[non_backward][0]

        # Applying cross-coupling for each disk - API 684 - SP6.8.5.9
        else:
            for i, Q in enumerate(cross_coupled_array[:, :-1]):
                bearings = [copy(b) for b in bearing_list]
                # cross-coupling introduced at overhung disks
                for n, q in zip(self.disk_nodes, Q):
                    cross_coupling = BearingElement(n=n, kxx=0, cxx=0, kxy=q, kyx=-q)
                    bearings.append(cross_coupling)

                aux_rotor = Rotor(
                    shaft_elements=self.rotor.shaft_elements,
                    disk_elements=[],
                    bearing_elements=bearings,
                    rated_w=self.rotor.rated_w,
                )
                modal = aux_rotor.run_modal(speed=oper_speed * np.pi / 30)
                non_backward = modal.whirl_direction() != "Backward"
                log_dec[i] = modal.log_dec[non_backward][0]

        # verifies if log dec is greater than zero to begin extrapolation
        cross_coupled_Qa = cross_coupled_array[:, -1]
        if log_dec[-1] >= 0:
            g = interp1d(
                cross_coupled_Qa, log_dec, fill_value="extrapolate", kind="linear"
            )
            stiff = cross_coupled_Qa[-1] * (1 + 1 / (len(cross_coupled_Qa)))
            while g(stiff) > 0:
                log_dec = np.append(log_dec, g(stiff))
                cross_coupled_Qa = np.append(cross_coupled_Qa, stiff)
                stiff += cross_coupled_Qa[-1] / (len(cross_coupled_Qa))
            Q0 = cross_coupled_Qa[-1]

        else:
            idx = min(range(len(log_dec)), key=lambda i: abs(log_dec[i]))
            Q0 = cross_coupled_Qa[idx]

        # Find value for log_dec corresponding to Qa
        log_dec_a = log_dec[np.where(cross_coupled_Qa == Qa)][0]

        # CSR - Critical Speed Ratio
        crit_speed = self.rotor.run_modal(speed=self.maxspeed).wn[0]
        CSR = self.maxspeed / crit_speed

        # RHO_mean - Average gas density
        RHO_mean = (RHOd + RHOs) / 2
        RHO = np.linspace(0, RHO_mean * 5, 501)

        # CSR_boundary - function to define the CSR boundaries
        CSR_boundary = np.piecewise(
            RHO,
            [RHO <= 16.53, RHO > 16.53, RHO == 60, RHO > 60],
            [2.5, lambda RHO: (-0.0115 * RHO + 2.69), 2.0, 0.0],
        )

        # Plotting area

        fig1 = go.Figure()

        fig1.add_trace(
            go.Scatter(
                x=cross_coupled_Qa,
                y=log_dec,
                mode="lines",
                showlegend=False,
                hoverinfo="none",
            )
        )
        fig1.add_trace(
            go.Scatter(
                x=[Qa],
                y=[log_dec_a],
                mode="markers",
                name="<b>Qa: Anticipated cross-coupling</b>",
                hoverinfo="none",
            )
        )
        fig1.add_annotation(
            x=Qa,
            y=log_dec_a,
            axref="x",
            ayref="y",
            xshift=15,
            yshift=15,
            text="<b>Qa</b>",
            showarrow=False,
        )
        fig1.update_xaxes(
            title_text="<b>Applied Cross-Coupled Stiffness, Q (N/m)</b>",
            rangemode="nonnegative",
        )
        fig1.update_yaxes(title_text="<b>Log Dec</b>", rangemode="nonnegative")
        fig1.update_layout(
            title=dict(
                text=(
                    "<b>Applied Cross-Coupled Stiffness vs. Log Decrement</b><br>"
                    + "<b>(API 684 - SP 6.8.5.10)</b>"
                )
            )
        )

        fig2 = go.Figure()
        fig2.add_annotation(
            x=RHO_mean,
            y=CSR,
            axref="x",
            ayref="y",
            xshift=40,
            yshift=0,
            text="<b>{}</b>".format(self.tag),
            showarrow=False,
        )

        for text, x, y in zip(["Region A", "Region B"], [30, 60], [1.20, 2.75]):
            fig2.add_annotation(
                x=x,
                y=y,
                axref="x",
                ayref="y",
                xshift=0,
                yshift=0,
                text=f"<b>{text}</b>",
                opacity=0.4,
                showarrow=False,
            )

        fig2.add_trace(
            go.Scatter(
                x=RHO,
                y=CSR_boundary,
                mode="lines",
                showlegend=False,
                hoverinfo="none",
                xaxis="x",
            )
        )
        fig2.add_trace(
            go.Scatter(
                x=0.062428 * RHO,
                y=CSR_boundary,
                mode="lines",
                showlegend=False,
                hoverinfo="none",
                xaxis="x2",
            )
        )
        fig2.add_trace(
            go.Scatter(
                x=[RHO_mean],
                y=[CSR],
                mode="markers",
                name="<b>CSR: Critical Speed Ratio</b>",
                hoverinfo="none",
                xaxis="x",
            )
        )

        fig2.update_xaxes(mirror=True)
        fig2.update_yaxes(
            title_text="<b>Maximum Critical Speed Ratio</b>",
            rangemode="nonnegative",
            domain=[0.1, 1],
        )
        fig2.update_layout(
            xaxis=dict(
                title_text="<b>kg/m³</b>",
                rangemode="nonnegative",
                overlaying="x2",
                anchor="y",
            ),
            xaxis2=dict(
                title_text="<b>lb/ft³</b>",
                rangemode="nonnegative",
                anchor="free",
                side="bottom",
                position=0,
            ),
            title=dict(
                text=(
                    "<b>CSR vs. Mean Gas Density</b><br>"
                    + "<b>(API 684 - SP 6.8.5.10)</b>"
                )
            ),
        )

        # Level 1 screening criteria - API 684 - SP6.8.5.10
        idx = min(range(len(RHO)), key=lambda i: abs(RHO[i] - RHO_mean))

        if self.machine_type == "compressor":
            if Q0 / Qa < 2.0:
                condition = True

            if log_dec_a < 0.1:
                condition = True

            if 2.0 < Q0 / Qa < 10.0 and CSR > CSR_boundary[idx]:
                condition = True

            else:
                condition = False

        if self.machine_type == "turbine" or self.machine_type == "axial flow":
            if log_dec_a < 0.1:
                condition = True

            else:
                condition = False

        # updating attributes
        self.Q0 = Q0
        self.Qa = Qa
        self.log_dec_a = log_dec_a
        self.CSR = CSR
        self.Qratio = Q0 / Qa
        self.crit_speed = crit_speed
        self.MCS = self.maxspeed
        self.RHO_gas = RHO_mean
        self.condition = condition

        return fig1, fig2

    def stability_level_2(self):
        """Stability analysis level 2.

        For the level 2 stability analysis additional sources that contribute
        to the rotor stability shall be considered such as:
        a)  labyrinth seals;
        b)  damper seals;
        c)  impeller/blade flow aerodynamic effects;
        d)  internal friction.

        Returns
        -------
        df_logdec: pd.DataFrame
            A dataframe relating the logarithmic decrement for each case analyzed.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> dataframe = report.stability_level_2()
        """
        # Build a list of seals
        seal_list = [
            copy(b) for b in self.rotor.bearing_elements if isinstance(b, SealElement)
        ]

        bearing_list = [
            copy(b)
            for b in self.rotor.bearing_elements
            if not isinstance(b, SealElement)
        ]

        log_dec_seal = []
        log_dec_disk = []
        log_dec_full = []
        data_seal = {}
        data_disk = {}
        data_rotor = {}

        # Evaluate log dec for each component - Disks
        if len(self.rotor.disk_elements):
            for disk in self.rotor.disk_elements:
                aux_rotor = Rotor(
                    shaft_elements=self.rotor.shaft_elements,
                    disk_elements=[disk],
                    bearing_elements=bearing_list,
                    rated_w=self.maxspeed,
                )
                modal = aux_rotor.run_modal(speed=self.maxspeed)
                non_backward = modal.whirl_direction() != "Backward"
                log_dec_disk.append(modal.log_dec[non_backward][0])

            # Evaluate log dec for group bearings + disks
            disk_tags = [
                "Shaft + Bearings + " + disk.tag for disk in self.rotor.disk_elements
            ]

            # Evaluate log dec for group bearings + all disks
            if len(self.rotor.disk_elements) > 1:
                all_disks_tag = " + ".join(
                    [disk.tag for disk in self.rotor.disk_elements]
                )
                disk_tags.append("Shaft + Bearings + " + all_disks_tag)

                aux_rotor = Rotor(
                    shaft_elements=self.rotor.shaft_elements,
                    disk_elements=self.rotor.disk_elements,
                    bearing_elements=bearing_list,
                    rated_w=self.maxspeed,
                )
                modal = aux_rotor.run_modal(speed=self.maxspeed)
                non_backward = modal.whirl_direction() != "Backward"
                log_dec_disk.append(modal.log_dec[non_backward][0])

            data_disk = {"tags": disk_tags, "log_dec": log_dec_disk}

        # Evaluate log dec for each component - Seals
        if len(seal_list):
            for seal in seal_list:
                bearings_seal = deepcopy(bearing_list)
                bearings_seal.append(seal)

                aux_rotor = Rotor(
                    shaft_elements=self.rotor.shaft_elements,
                    disk_elements=[],
                    bearing_elements=bearings_seal,
                    rated_w=self.maxspeed,
                )
                modal = aux_rotor.run_modal(speed=self.maxspeed)
                non_backward = modal.whirl_direction() != "Backward"
                log_dec_seal.append(modal.log_dec[non_backward][0])

            seal_tags = ["Shaft + Bearings + " + seal.tag for seal in seal_list]

            if len(seal_list) > 1:
                # Evaluate log dec for group bearings + seals
                all_seals_tag = " + ".join([seal.tag for seal in seal_list])
                seal_tags.append("Shaft + Bearings + " + all_seals_tag)

                aux_rotor = Rotor(
                    shaft_elements=self.rotor.shaft_elements,
                    disk_elements=[],
                    bearing_elements=self.rotor.bearing_elements,
                    rated_w=self.maxspeed,
                )
                modal = aux_rotor.run_modal(speed=self.maxspeed)
                non_backward = modal.whirl_direction() != "Backward"
                log_dec_seal.append(modal.log_dec[non_backward][0])

            data_seal = {"tags": seal_tags, "log_dec": log_dec_seal}

        # Evaluate log dec for all components
        modal = self.rotor.run_modal(speed=self.maxspeed)
        non_backward = modal.whirl_direction() != "Backward"
        log_dec_full.append(modal.log_dec[non_backward][0])
        rotor_tags = [self.tag]

        data_rotor = {"tags": rotor_tags, "log_dec": log_dec_full}

        df_logdec_disk = pd.DataFrame(data_disk)
        df_logdec_seal = pd.DataFrame(data_seal)
        df_logdec_full = pd.DataFrame(data_rotor)
        df_logdec = pd.concat([df_logdec_disk, df_logdec_seal, df_logdec_full])
        df_logdec = df_logdec.reset_index(drop=True)

        self.df_logdec_disk = df_logdec_disk
        self.df_logdec_seal = df_logdec_seal
        self.df_logdec_full = df_logdec_full
        self.df_logdec = df_logdec

        return df_logdec

    def summary(self):
        """Return datarfreames for Report summary.

        This method will create dataframes with relevant info about the report.

        Returns
        -------
        df_stab_lvl1 : pd.DataFrame
            Dataframe with stability level 1 results
        df_stab_lvl2 : pd.DataFrame
            Dataframe with stability level 2 results

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> stability1 = report.stability_level_1(D=[0.35, 0.35],
        ...                          H=[0.08, 0.08],
        ...                          HP=[10000, 10000],
        ...                          RHO_ratio=[1.11, 1.14],
        ...                          RHOd=30.45,
        ...                          RHOs=37.65,
        ...                          oper_speed=1000.0)
        >>> stability2 = report.stability_level_2()
        >>> df_lvl1, df_lvl2 = report.summary()
        """
        stab_lvl1_data = dict(
            tags=[self.tag],
            machine_type=[self.machine_type],
            Q0=[self.Q0],
            Qa=[self.Qa],
            log_dec_a=[self.log_dec_a],
            Qratio=[self.Qratio],
            crti_speed=[self.crit_speed],
            MCS=[self.MCS],
            CSR=[self.CSR],
            RHO_gas=[self.RHO_gas],
        )
        stab_lvl2_data = dict(
            tags=self.df_logdec["tags"], logdec=self.df_logdec["log_dec"]
        )

        df_stab_lvl1 = pd.DataFrame(stab_lvl1_data)
        df_stab_lvl2 = pd.DataFrame(stab_lvl2_data)

        return df_stab_lvl1, df_stab_lvl2

    def plot_summary(self):
        """Plot the report .

        This method will create tables to be presented in the report.

        Returns
        -------
        fig : Plotly graph_objects.make_subplots()
            The figure object with the tables.

        Example
        -------
        >>> import ross as rs
        >>> report = rs.report_example()
        >>> stability1 = report.stability_level_1(D=[0.35, 0.35],
        ...                          H=[0.08, 0.08],
        ...                          HP=[10000, 10000],
        ...                          RHO_ratio=[1.11, 1.14],
        ...                          RHOd=30.45,
        ...                          RHOs=37.65,
        ...                          oper_speed=1000.0)
        >>> stability2 = report.stability_level_2()
        >>> table = report.plot_summary()
        """
        stab_lvl1_data, stab_lvl2_data = self.summary()
        for var in stab_lvl1_data.columns[2:]:
            stab_lvl1_data[str(var)] = np.round(stab_lvl1_data[str(var)], 3)

        stab_lvl2_data["logdec"] = np.round(stab_lvl2_data["logdec"], 4)

        stab_lvl1_titles = [
            "<b>Rotor Tag</b>",
            "<b>Machine Type</b>",
            "<b>Q_0</b>",
            "<b>Q_A</b>",
            "<b>log dec (δ)</b>",
            "<b>Q_0 / Q_A</b>",
            "<b>1st Critical Spped</b>",
            "<b>MCS</b>",
            "<b>CSR</b>",
            "<b>Gas Density</b>",
        ]
        stab_lvl2_titles = ["<b>Components</b>", "<b>Log. Dec.</b>"]

        fig = make_subplots(
            rows=2,
            cols=1,
            specs=[[{"type": "table"}], [{"type": "table"}]],
            subplot_titles=["<b>Stability Level 1</b>", "<b>Stability Level 2</b>"],
        )

        colors = ["#ffffff", "#c4d9ed"]
        cell_colors = [colors[i % 2] for i in range(len(stab_lvl1_data["tags"]))]
        fig.add_trace(
            go.Table(
                header=dict(
                    values=stab_lvl1_titles,
                    font=dict(family="Verdana", size=14, color="white"),
                    line=dict(color="#1e4162", width=1.5),
                    fill=dict(color="#1e4162"),
                    align="center",
                ),
                cells=dict(
                    values=[stab_lvl1_data[str(var)] for var in stab_lvl1_data.columns],
                    font=dict(family="Verdana", size=14, color="#12263b"),
                    line=dict(color="#c4d9ed", width=1.5),
                    fill=dict(color=[cell_colors * len(stab_lvl1_data["tags"])]),
                    align="center",
                    height=25,
                ),
            ),
            row=1,
            col=1,
        )

        cell_colors = [colors[i % 2] for i in range(len(stab_lvl2_data["tags"]))]
        fig.add_trace(
            go.Table(
                header=dict(
                    values=stab_lvl2_titles,
                    font=dict(family="Verdana", size=14, color="white"),
                    line=dict(color="#1e4162", width=1.5),
                    fill=dict(color="#1e4162"),
                    align="center",
                ),
                cells=dict(
                    values=[stab_lvl2_data[str(var)] for var in stab_lvl2_data.columns],
                    font=dict(family="Verdana", size=14, color="#12263b"),
                    line=dict(color="#c4d9ed", width=1.5),
                    fill=dict(color=[cell_colors * len(stab_lvl2_data["tags"])]),
                    align="center",
                    height=25,
                ),
            ),
            row=2,
            col=1,
        )

        return fig


def report_example():
    """Build a report example.

    This function returns an instance of a simple report from a rotor
    example. The purpose of this is to make available a simple model
    so that doctest can be written using this.

    Returns
    -------
    An instance of a report object.

    Examples
    --------
    >>> import ross as rs
    >>> report = rs.report_example()
    >>> report.rotor_type
    'between_bearings'
    """
    i_d = 0
    o_d = 0.05
    n = 6
    L = [0.25 for _ in range(n)]

    shaft_elem = [
        ShaftElement(
            l,
            i_d,
            o_d,
            material=steel,
            shear_effects=True,
            rotary_inertia=True,
            gyroscopic=True,
        )
        for l in L
    ]

    disk0 = DiskElement.from_geometry(
        n=2, material=steel, width=0.07, i_d=0.05, o_d=0.28
    )
    disk1 = DiskElement.from_geometry(
        n=4, material=steel, width=0.07, i_d=0.05, o_d=0.28
    )

    stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
    stfy = [0.8e7, 0.9e7, 1.0e7, 1.1e7]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(0, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
    bearing1 = BearingElement(6, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)

    rotor = Rotor(shaft_elem, [disk0, disk1], [bearing0, bearing1])

    # coefficients for minimum clearance
    stfx = [0.7e7, 0.8e7, 0.9e7, 1.0e7]
    dampx = [2.0e3, 1.9e3, 1.8e3, 1.7e3]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(0, kxx=stfx, cxx=dampx, frequency=freq)
    bearing1 = BearingElement(6, kxx=stfx, cxx=dampx, frequency=freq)
    min_clearance_brg = [bearing0, bearing1]

    # coefficients for maximum clearance
    stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
    dampx = [2.8e3, 2.7e3, 2.6e3, 2.5e3]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(0, kxx=stfx, cxx=dampx, frequency=freq)
    bearing1 = BearingElement(6, kxx=stfx, cxx=dampx, frequency=freq)
    max_clearance_brg = [bearing0, bearing1]

    bearings = [min_clearance_brg, max_clearance_brg]
    return Report(
        rotor=rotor,
        speed_range=(400, 1000),
        tripspeed=1200,
        bearing_stiffness_range=(5, 8),
        bearing_clearance_lists=bearings,
        speed_units="rad/s",
    )
