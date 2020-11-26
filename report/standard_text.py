"""Standard text for ROSS Report."""

__all__ = ["ReportTemplateText"]


class ReportTemplateText:
    """Create text instances to ROSS REPORT.

    Each argument holds a dictionary related to a single Report section.
    The arguments are collected and set to the HTML file according to the content
    organization.

    Parameters
    ----------
    report : rp.Report object
        An instance from Report class. The texts collects some arguments from the
        report object to customize the informations for each rotor model and its
        respective analyses.
    """

    def __init__(self, report):
        self.introduction = dict(
            intro=(
                """
                ROSS - Rotordynamics Open-Source Software is a library written in Python
                for rotordynamic analyses. It's developed by Petrobrás and Federal
                University of Rio de Janeiro.
                It allows the construction of rotor models and their
                numerical simulation. Shaft elements, as a default, are modeled with the
                Timoshenko beam theory, which considers shear and rotary inertia
                effects, and discretized by means of the Finite Element Method. Disks
                are assumed to be rigid bodies, thus their strain energy is not taken
                into account. And bearings/seals are included as linear
                stiffness/damping coefficients.

                ROSS carries out several different analyses which include:
                1. Static analysis;
                2. Modal analysis - Natural frequencies and mode shapes determination;
                3. Damped and undemped critical speed analysis;
                4. Unbalance response analysis;
                5. Time domain linear analysis;
                6. Stability analysis of the rotor.
                """
            )
        )

        self.static_analysis = dict(
            intro=(
                """
                The static analysis calculates the shaft deformation, shearing forces
                and bending moments for the rotor, given its self weight (shaft and
                other couplings).
                Figure XXXX shows the free-body diagram representation, where Fd stands
                for Disk wieght and Fb stands for bearing reaction forces.
                Figure XXXX shows the shaft static deformation (gyroscopic effect not
                included - speed = 0 RPM).
                Figures XXXX and XXXX show the diagrams for shearing force and
                bending moment respectively.
                """
            )
        )

        self.undamped_critical_speed_map = dict(
            intro=(
                """
                The Critical Speed Map determines approximated values of the natural
                frequencies as function of the bearings stiffness coefficients.
                The intersections of these bearings curves with the natural frequencies
                curves define the undamped critical speeds. The horizontal dashed lines
                represent the rotor operation speeds.
                """
            ),
            min_clearance=(
                """
                Figure XXXX shows the Undamped Critical Speed Map under the minimum
                clearance of the bearings.
                """
            ),
            rated_clearance=(
                """
                Figure XXXX shows the Undamped Critical Speed Map under the rated
                clearance of the bearings.
                """
            ),
            max_clearance=(
                """
                Figure XXXX shows the Undamped Critical Speed Map under the maximum
                clearance of the bearings.
                """
            ),
        )

        self.damped_critical_speed_map = dict(
            intro=(
                """
                The Damped Critical Speed Map, also called Campbell Diagram determines
                approximated values of the damped natural frequencies as function of
                the rotor speed. The intersections of each harmonic curve with the
                natural frequencies scatter curves define the critical speeds.
                Furthermore, the damping level of each mode, measured by the logarithm
                decrement, is presented with a color scale.
                """
            ),
            min_clearance=(
                """
                Figure XXXX shows the Campbell diagram under the minimum clearance of
                the bearings.
                """
            ),
            rated_clearance=(
                """
                Figure XXXX shows the Campbell diagram under the rated clearance of the
                bearings.
                """
            ),
            max_clearance=(
                """
                Figure XXXX shows the Campbell diagram under the maximum clearance of
                the bearings.
                """
            ),
        )

        self.mode_shapes = dict(
            intro=(
                """The mode shapes are calculate through the rotor modal analysis. The
                results present the 2d shapes with respective natural frequencies, whirl
                direction and the log dec. The modal analysis is performed to the rotor
                operation speed.
                """
            ),
            min_clearance=(
                f"""
                Figure XXXX shows first two mode shape of {report.tag} rotor, for the
                minimum clearance.
                """
            ),
            rated_clearance=(
                f"""
                Figure XXXX shows first two mode shape of {report.tag} rotor, for the
                rated clearance.
                """
            ),
            max_clearance=(
                f"""
                Figure XXXX shows first two mode shape of {report.tag} rotor, for the
                maximum clearance.
                """
            ),
        )

        self.unbalance_response = dict(
            intro=(
                """The Unbalance Response Analysis represents the rotor synchronous
                excitation due to rotor unbalance. The results present a diagram with
                plots of amplitude and phase versus frequency and a polar plot of
                amplitude versus phase The setup of unbalance positions, weights and
                phases is defined according to API 684 SP6.8.2.7 and SP6.8.2.8, which is
                based on the machine type, bearings and couplings positions and the
                mode shapes configurations. The Amplification Factors,
                Separation Margins and Scale Factors are defined by
                API684 - SP6.8.2.1, SP6.8.2.10 and SP6.8.2.11 respectively.
                """
            ),
            min_clearance=(
                """
                The unbalance response diagram is shown in Figure XXXX and Table XXXX
                show a brief results summary under the minimum clearance of the
                bearings. The amplitude in all curves are calculated for the nodes and
                orientations for each probe selected in the analysis.
                """
            ),
            rated_clearance=(
                """
                The unbalance response diagram is shown in Figure XXXX and Table XXXX
                show a brief results summary under the rated clearance of the bearings.
                The amplitude in all curves are calculated for the nodes and
                orientations for each probe selected in the analysis.
                """
            ),
            max_clearance=(
                """
                The unbalance response diagram is shown in Figure XXXX and Table XXXX
                show a brief results summary under the maximum clearance of the
                bearings. The amplitude in all curves are calculated for the nodes and
                orientations for each probe selected in the analysis.
                """
            ),
        )

        self.deflected_shape = dict(
            intro=(
                """The deflected shape analysis presets results for the 3d shape rotor
                deformation due the applied imbalance for a given speed, the 2d shape of
                the absolute value for the major axis and the bending moment diagram.
                """
            ),
            min_clearance=(
                f"""
                The plots of deflected shapes for speed {report.config.run_unbalance_response.plot_deflected_shape.speed}
                with minimum clearance are shown below.
                """
            ),
            rated_clearance=(
                f"""
                The plots of deflected shapes for speed {report.config.run_unbalance_response.plot_deflected_shape.speed}
                with rated clearance are shown below.
                """
            ),
            max_clearance=(
                f"""
                The plots of deflected shapes for speed {report.config.run_unbalance_response.plot_deflected_shape.speed}
                with maximum clearance are shown below.
                """
            ),
        )

        self.level_1_analysis = dict(
            intro=(
                """
                The Stability Level 1 Analysis determines the natural frequencies and
                the corresponding logarithmic decrements (log decs) of the damped
                rotor/support system using a complex value analysis.
                This analysis is performed with a varying amount of cross coupling
                introduced at the rotor mid-span for between bearing rotors or at the
                center of gravity of the stage or impeller for single overhung rotors.
                For double overhung rotors, the cross coupling shall be placed at each
                stage or impeller concurrently and shall reflect the ratio of the
                anticipated cross coupling (qa, calculated for each impeller or stage).

                The anticipated cross coupling, QA, present in the rotor is defined by
                the following procedures:

                For centrifugal compressors:

                    INSERT FORMULA

                HP is the rated power per impeller, Nm/s (HP);
                Bc is 3;
                C is 9.55 (63);
                ρd is the discharge gas density per impeller, kg/m3 (lbm/ft3);
                ρs is the suction gas density per impeller, kg/m3 (lbm/ft3);
                Dc is the impeller diameter, mm (in.);
                Hc is the minimum of diffuser or impeller discharge width per impeller, mm (in.);
                Nr is the normal operating speed for calculation of aerodynamic excitation (rpm);
                qa is the cross coupling for each individual impeller, kN/mm (klbf/in).

                For axial flow rotors:

                    INSERT FORMULA

                Bt is 1.5;
                Dt is the blade pitch diameter, mm (in.);
                Ht is the effective blade height, mm (in.).
                """
            ),
            min_clearance=(
                """
                The result of level I stability （plot of Applied Cross-coupled
                stiffness vs logarithmic decrement） is shown as Figure XXXX. The plot
                shows the relationship of the logarithmic decrement and Cross-coupled
                stiffness of rotor, Qa is anticipated cross coupling stiffness, Q0 is
                the amount of the applied cross coupling required to produce a zero
                logarithmic decrement (where the curve crosses the abscissa).
                Figure XXXX is a screening criteria relating the Critical Speed Ratio
                (CSR) and the average gas density. If the screening point is located
                on Region B, further stabiliy analysis is required.
                """
            ),
            rated_clearance=(
                """
                The result of level I stability （plot of Applied Cross-coupled
                stiffness vs logarithmic decrement） is shown as Figure XXXX. The plot
                shows the relationship of the logarithmic decrement and Cross-coupled
                stiffness of rotor, Qa is anticipated cross coupling stiffness, Q0 is
                the amount of the applied cross coupling required to produce a zero
                logarithmic decrement (where the curve crosses the abscissa).
                Figure XXXX is a screening criteria relating the Critical Speed Ratio
                (CSR) and the average gas density. If the screening point is located
                on Region B, further stabiliy analysis is required.
                """
            ),
            max_clearance=(
                """
                The result of level I stability （plot of Applied Cross-coupled
                stiffness vs logarithmic decrement） is shown as Figure XXXX. The plot
                shows the relationship of the logarithmic decrement and Cross-coupled
                stiffness of rotor, Qa is anticipated cross coupling stiffness, Q0 is
                the amount of the applied cross coupling required to produce a zero
                logarithmic decrement (where the curve crosses the abscissa).
                Figure XXXX is a screening criteria relating the Critical Speed Ratio
                (CSR) and the average gas density. If the screening point is located
                on Region B, further stabiliy analysis is required.
                """
            ),
        )

        self.level_2_analysis = dict(
            intro=(
                """
                The Stability Level 2 analysis shall include the dynamic characteristics
                of all components that, somehow, affects the stability behaviour of the
                rotor machine. These dynamic effects replace the anticipated cross
                coupling, QA, calculated in the Stability Level 1 Analysis.
                """
            ),
            min_clearance=(
                """
                The Table XXXX below shows the log decrement for several rotor
                configurations, considering each component individually and the full
                rotor model. Each row present a configuration and the respective log
                decrement calculated for the maximum continuous speed.
                """
            ),
            rated_clearance=(
                """
                The Table XXXX below shows the log decrement for several rotor
                configurations, considering each component individually and the full
                rotor model. Each row present a configuration and the respective log
                decrement calculated for the maximum continuous speed.
                """
            ),
            max_clearance=(
                """
                The Table XXXX below shows the log decrement for several rotor
                configurations, considering each component individually and the full
                rotor model. Each row present a configuration and the respective log
                decrement calculated for the maximum continuous speed.
                """
            ),
        )

        self.conclusion = dict(
            Pass=(
                f"""
                According to the analysis, the rotor {report.tag} complies with the
                API 617-2014 Standard, the lateral vibration and stability analysis of
                the rotor are acceptable.
                """
            ),
            Not_Pass=(
                f"""
                According to the analysis, the rotor {report.tag} do not complies with
                the API 617-2014 Standard, the lateral vibration and / or stability
                analysis of the rotor are not acceptable. Further verifications are
                required
                """
            ),
        )

    def __getitem__(self, option):
        """Return the value for a given option from the dictionary.

        Parameters
        ----------
        option : str
            A dictionary key corresponding to the text section.

        Raises
        ------
        KeyError
            Raises an error if the parameter doesn't belong to the dictionary.

        Returns
        -------
        Return the value for the given key.
        """
        if option not in self.__dict__.keys():
            raise KeyError("Option '{}' not found.".format(option))

        return self.__dict__[option]
