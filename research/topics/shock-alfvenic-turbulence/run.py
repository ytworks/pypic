"""Run one configuration of the topic (planner: plan.py).

    python -m pypic.launch -n 8 research/topics/shock-alfvenic-turbulence/run.py \
        --config research/topics/shock-alfvenic-turbulence/configs/<name>.yaml --out <dir>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_topic  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("shock_alfvenic_turbulence", plan_topic, params_from_dict, doc=__doc__)
