"""Bow shock of a magnetised planet in the solar wind (3-D): absorbing sphere with a dipole in an open box.

python -m pypic.launch -n 8 examples/bow_shock/run.py --config examples/bow_shock/configs/gate.yaml
python examples/bow_shock/analyze.py output/bow_shock/gate --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_bow_shock  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("bow_shock", plan_bow_shock, params_from_dict, doc=__doc__)
