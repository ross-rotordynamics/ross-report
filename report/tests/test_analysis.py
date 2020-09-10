import numpy as np
import pytest
from numpy.testing import assert_allclose

from report.analysis import Report
from report.config import Config
from ross.bearing_seal_element import BearingElement
from ross.disk_element import DiskElement
from ross.materials import steel
from ross.rotor_assembly import Rotor
from ross.shaft_element import ShaftElement


@pytest.fixture
def report0():
    # rotor type: between bearings
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
    oper_clerance_brg = [bearing0, bearing1]

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

    config = Config()
    config.update_config(
        rotor_properties={
            "rotor_speeds": {
                "min_speed": 400,
                "max_speed": 1000,
                "oper_speed": 800,
                "trip_speed": 1200,
                "unit": "rad/s",
            },
            "rotor_id": {"type": "compressor"},
        },
        bearings={
            "oper_clearance": oper_clerance_brg,
            "min_clearance": min_clearance_brg,
            "max_clearance": max_clearance_brg,
        },
        run_campbell={"speed_range": np.linspace(0, 1500, 51)},
        run_unbalance_response={
            "probes": {
                "node": [1, 4],
                "orientation": [np.pi / 2, np.pi / 2],
                "unit": "rad",
            },
            "frequency_range": np.linspace(0, 1500, 101),
            "plot_deflected_shape": {"speed": [615]},
        },
        plot_ucs={"stiffness_range": (5, 8)},
        stability_level1={
            "D": [0.35, 0.35],
            "H": [0.08, 0.08],
            "rated_power": [6000, 8000],
            "rho_ratio": [1.11, 1.14],
            "rho_suction": 30.45,
            "rho_discharge": 37.65,
            "length_unit": "m",
            "power_unit": "hp",
            "density_unit": "kg/m**3",
        },
    )

    return Report(rotor, config)


@pytest.fixture
def report1():
    # rotor type: single overhung
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
        n=0, material=steel, width=0.07, i_d=0.05, o_d=0.28
    )

    stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
    stfy = [0.8e7, 0.9e7, 1.0e7, 1.1e7]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(2, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
    bearing1 = BearingElement(6, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
    oper_clerance_brg = [bearing0, bearing1]

    rotor = Rotor(shaft_elem, [disk0], [bearing0, bearing1])

    # coefficients for minimum clearance
    stfx = [0.7e7, 0.8e7, 0.9e7, 1.0e7]
    dampx = [2.0e3, 1.9e3, 1.8e3, 1.7e3]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(2, kxx=stfx, cxx=dampx, frequency=freq)
    bearing1 = BearingElement(6, kxx=stfx, cxx=dampx, frequency=freq)
    min_clearance_brg = [bearing0, bearing1]

    # coefficients for maximum clearance
    stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
    dampx = [2.8e3, 2.7e3, 2.6e3, 2.5e3]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(2, kxx=stfx, cxx=dampx, frequency=freq)
    bearing1 = BearingElement(6, kxx=stfx, cxx=dampx, frequency=freq)
    max_clearance_brg = [bearing0, bearing1]

    config = Config()
    config.update_config(
        rotor_properties={
            "rotor_speeds": {
                "min_speed": 400,
                "max_speed": 1000,
                "oper_speed": 800,
                "trip_speed": 1200,
                "unit": "rad/s",
            },
            "rotor_id": {"type": "turbine"},
        },
        bearings={
            "oper_clearance": oper_clerance_brg,
            "min_clearance": min_clearance_brg,
            "max_clearance": max_clearance_brg,
        },
        run_campbell={"speed_range": np.linspace(0, 1500, 51)},
        run_unbalance_response={
            "probes": {
                "node": [1, 4],
                "orientation": [np.pi / 2, np.pi / 2],
                "unit": "rad",
            },
            "frequency_range": np.linspace(0, 1500, 101),
            "plot_deflected_shape": {"speed": [615]},
        },
        plot_ucs={"stiffness_range": (5, 8)},
        stability_level1={
            "D": [0.35, 0.35],
            "H": [0.08, 0.08],
            "rated_power": [10000, 10000],
            "rho_ratio": [1.11, 1.14],
            "rho_suction": 30.45,
            "rho_discharge": 37.65,
            "length_unit": "m",
            "power_unit": "hp",
            "density_unit": "kg/m**3",
        },
    )

    return Report(rotor, config)


@pytest.fixture
def report2():
    # rotor type: double overhung
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
        n=0, material=steel, width=0.07, i_d=0.05, o_d=0.28
    )
    disk1 = DiskElement.from_geometry(
        n=6, material=steel, width=0.07, i_d=0.05, o_d=0.28
    )

    stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
    stfy = [0.8e7, 0.9e7, 1.0e7, 1.1e7]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(2, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
    bearing1 = BearingElement(4, kxx=stfx, kyy=stfy, cxx=2e3, frequency=freq)
    oper_clerance_brg = [bearing0, bearing1]

    rotor = Rotor(shaft_elem, [disk0, disk1], [bearing0, bearing1])

    # coefficients for minimum clearance
    stfx = [0.7e7, 0.8e7, 0.9e7, 1.0e7]
    dampx = [2.0e3, 1.9e3, 1.8e3, 1.7e3]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(2, kxx=stfx, cxx=dampx, frequency=freq)
    bearing1 = BearingElement(4, kxx=stfx, cxx=dampx, frequency=freq)
    min_clearance_brg = [bearing0, bearing1]

    # coefficients for maximum clearance
    stfx = [0.4e7, 0.5e7, 0.6e7, 0.7e7]
    dampx = [2.8e3, 2.7e3, 2.6e3, 2.5e3]
    freq = [400, 800, 1200, 1600]
    bearing0 = BearingElement(2, kxx=stfx, cxx=dampx, frequency=freq)
    bearing1 = BearingElement(4, kxx=stfx, cxx=dampx, frequency=freq)
    max_clearance_brg = [bearing0, bearing1]

    config = Config()
    config.update_config(
        rotor_properties={
            "rotor_speeds": {
                "min_speed": 3820.0,
                "max_speed": 9550.0,
                "oper_speed": 8000.0,
                "trip_speed": 10500.0,
                "unit": "rpm",
            },
            "rotor_id": {"type": "axial_flow"},
        },
        bearings={
            "oper_clearance": oper_clerance_brg,
            "min_clearance": min_clearance_brg,
            "max_clearance": max_clearance_brg,
        },
        run_campbell={"speed_range": np.linspace(0, 1500, 51)},
        run_unbalance_response={
            "probes": {
                "node": [1, 4],
                "orientation": [np.pi / 2, np.pi / 2],
                "unit": "rad",
            },
            "frequency_range": np.linspace(0, 1500, 101),
            "plot_deflected_shape": {"speed": [615]},
        },
        plot_ucs={"stiffness_range": (5, 8)},
        stability_level1={
            "D": [0.35, 0.35],
            "H": [0.08, 0.08],
            "rated_power": [6000, 8000],
            "rho_ratio": [1.11, 1.14],
            "rho_suction": 30.45,
            "rho_discharge": 37.65,
            "length_unit": "m",
            "power_unit": "hp",
            "density_unit": "kg/m**3",
        },
    )

    return Report(rotor, config)


def test_initial_attributes(report0, report1, report2):
    assert report0.rotor_type == "between_bearings"
    assert report0.disk_nodes == [2, 4]
    assert report0.tag == "Rotor 0"
    assert report1.rotor_type == "single_overhung_l"
    assert report1.disk_nodes == [0]
    assert report1.tag == "Rotor 0"
    assert report2.rotor_type == "double_overhung"
    assert report2.disk_nodes == [0, 6]
    assert report2.tag == "Rotor 0"


def test_report_static_forces(report0, report1, report2):
    F_0 = report0._static_forces()
    F_1 = report1._static_forces()
    F_2 = report2._static_forces()
    assert_allclose(F_0[0], 44.09320348748279, atol=1e-6)
    assert_allclose(F_0[1], 44.09320348748279, atol=1e-6)
    assert_allclose(F_1[0], 66.13980523122419, atol=1e-6)
    assert_allclose(F_1[1], -10.5440269209198, atol=1e-6)
    assert_allclose(F_2[0], 44.09320348748279, atol=1e-6)
    assert_allclose(F_2[1], 44.09320348748279, atol=1e-6)


def test_unbalance_forces(report0, report1, report2):
    Uforce_00 = report0._unbalance_forces(mode=0)
    assert_allclose(Uforce_00, [0.04479869], atol=1e-6)

    Uforce_02 = report0._unbalance_forces(mode=2)
    assert_allclose(Uforce_02, [0.02239935, 0.02239935], atol=1e-6)

    Uforce_10 = report1._unbalance_forces(mode=0)
    assert_allclose(Uforce_10, [2.0450646e-05], atol=1e-8)

    Uforce_12 = report1._unbalance_forces(mode=2)
    assert_allclose(Uforce_12, [2.0450646e-05], atol=1e-8)

    Uforce_20 = report2._unbalance_forces(mode=0)
    assert_allclose(Uforce_20, np.array([2.0450646e-05, 2.0450646e-05]), atol=1e-6)

    Uforce_22 = report2._unbalance_forces(mode=2)
    assert_allclose(Uforce_22, np.array([2.0450646e-05, 2.0450646e-05]), atol=1e-6)


def test_report_mode_shape(report0, report1, report2):
    n1, n2 = report0._mode_shape(mode=1)
    nodes = [int(node) for sub_nodes in [n1, n2] for node in sub_nodes]
    assert nodes == [3]

    n1, n2 = report1._mode_shape(mode=1)
    nodes = [int(node) for sub_nodes in [n1, n2] for node in sub_nodes]
    assert nodes == [0]

    n1, n2 = report1._mode_shape(mode=3)
    nodes = [int(node) for sub_nodes in [n1, n2] for node in sub_nodes]
    assert nodes == [0]

    n1, n2 = report2._mode_shape(mode=1)
    nodes = [int(node) for sub_nodes in [n1, n2] for node in sub_nodes]
    assert nodes == [0, 6]

    n1, n2 = report2._mode_shape(mode=3)
    nodes = [int(node) for sub_nodes in [n1, n2] for node in sub_nodes]
    assert nodes == [0, 6]


def test_stability_level1(report0, report1, report2):
    _ = report0._stability_level_1()

    assert_allclose(report0.Q0, 19021.9367934373, atol=1e-4)
    assert_allclose(report0.Qa, 2113.54853260414442, atol=1e-4)
    assert_allclose(report0.log_dec_a, 0.014520423503078583, atol=1e-4)
    assert_allclose(report0.CSR, 8.83742092096016, atol=1e-4)
    assert_allclose(report0.Qratio, 9.0, atol=1e-4)
    assert_allclose(report0.crit_speed, 113.15518508670297, atol=1e-4)
    assert_allclose(report0.MCS, 1000.0, atol=1e-4)
    assert_allclose(report0.rho_gas, 34.05, atol=1e-4)
    assert report0.condition == True

    _ = report1._stability_level_1()

    assert_allclose(report1.Q0, 25570.065206064668, atol=1e-4)
    assert_allclose(report1.Qa, 669.69218396836, atol=1e-4)
    assert_allclose(report1.log_dec_a, 0.031775752903274904, atol=1e-6)
    assert_allclose(report1.CSR, 9.927312211419821, atol=1e-4)
    assert_allclose(report1.Qratio, 38.1818181818182, atol=1e-4)
    assert_allclose(report1.crit_speed, 100.73220008630899, atol=1e-4)
    assert_allclose(report1.MCS, 1000.0, atol=1e-4)
    assert_allclose(report1.rho_gas, 34.05, atol=1e-4)
    assert report1.condition == True

    _ = report2._stability_level_1()

    assert_allclose(report2.Q0, 2434436.079545, atol=1e-4)
    assert_allclose(report2.Qa, 895.3125000000002, atol=1e-4)
    assert_allclose(report2.log_dec_a, 0.19275430850265318, atol=1e-4)
    assert_allclose(report2.CSR, 11.718265284342298, atol=1e-4)
    assert_allclose(report2.Qratio, 2719.0909090908012, atol=1e-4)
    assert_allclose(report2.crit_speed, 85.34314910322335, atol=1e-4)
    assert_allclose(report2.MCS, 9550.0, atol=1e-4)
    assert_allclose(report2.rho_gas, 34.05, atol=1e-4)
    assert report2.condition == False


def test_stability_level2(report0, report1, report2):
    df0 = report0._stability_level_2()
    df1 = report1._stability_level_2()

    assert_allclose(
        df0["log_dec"].tolist(),
        [
            0.011399321894782566,
            0.011399321895234833,
            0.010612141367754908,
            0.010612141367232914,
        ],
        atol=1e-6,
    )
    assert_allclose(
        df1["log_dec"].tolist(),
        [0.02449532360217691, 0.024495323586707394],
        atol=1e-6,
    )
