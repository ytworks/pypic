"""An Alfvénic wave carried by a flow through an inflow face, optionally into a piston shock.

python -m pypic.launch -n 2 examples/alfven_inflow/run.py --config examples/alfven_inflow/configs/gate.yaml
python examples/alfven_inflow/analyze.py output/alfven_inflow/gate
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_alfven_inflow  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("alfven_inflow", plan_alfven_inflow, params_from_dict, doc=__doc__)
