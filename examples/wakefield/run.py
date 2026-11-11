"""Plasma wake behind a laser pulse or a particle bunch in a static box.

    python -m pypic.launch -n 8 examples/wakefield/run.py --config examples/wakefield/configs/laser_linear.yaml
    python examples/wakefield/analyze.py output/wakefield/laser_linear

``--set drive=false`` runs the same box without the driver (null test).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plan import params_from_dict, plan_wakefield  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("wakefield", plan_wakefield, params_from_dict, doc=__doc__)
