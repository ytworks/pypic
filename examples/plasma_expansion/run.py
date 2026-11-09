"""Plasma expansion into vacuum (isothermal, Mora 2003): ion acceleration by the sheath field.

python -m pypic.launch -n 8 examples/plasma_expansion/run.py --config examples/plasma_expansion/configs/isothermal.yaml
python examples/plasma_expansion/analyze.py output/plasma_expansion/isothermal --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_expansion  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("plasma_expansion", plan_expansion, params_from_dict, doc=__doc__)
