"""Relativistic unmagnetised pair shock: jump conditions and Fermi acceleration.

python -m pypic.launch -n 8 examples/relativistic_shock/run.py --config examples/relativistic_shock/configs/pair_2d.yaml
python examples/relativistic_shock/analyze.py output/relativistic_shock/pair_2d --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_shock  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("relativistic_shock", plan_shock, params_from_dict, doc=__doc__)
