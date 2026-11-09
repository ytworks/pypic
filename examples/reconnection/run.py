"""Relativistic pair magnetic reconnection (double Harris sheet): particle acceleration.

python -m pypic.launch -n 8 examples/reconnection/run.py --config examples/reconnection/configs/pair_sigma10.yaml
python examples/reconnection/analyze.py output/reconnection/pair_sigma10 --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_reconnection  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("reconnection", plan_reconnection, params_from_dict, doc=__doc__)
