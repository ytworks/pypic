"""Langmuir oscillations: Landau damping, Dawson's exact cold oscillation, wave breaking.

python -m pypic.launch -n 8 examples/landau_damping/run.py --config examples/landau_damping/configs/landau.yaml
python examples/landau_damping/analyze.py output/landau_damping/landau --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_langmuir  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("landau_damping", plan_langmuir, params_from_dict, doc=__doc__)
