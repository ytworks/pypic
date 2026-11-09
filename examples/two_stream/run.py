"""3-D two-stream instability experiment.

Single GPU:   python examples/two_stream/run.py --config examples/two_stream/configs/ei.yaml
Multi GPU:    python -m pypic.launch -n 8 examples/two_stream/run.py --config .../ep.yaml
Plan only:    python examples/two_stream/run.py --config ... --plan-only
Overrides:    --set u_drift=0.3 --set solver=implicit --set overrides.ppc=32 --set output.codec=deflate
Output:       output/two_stream/<config>/ (run.json, scalars.h5, reduced/, fields/, tracks/)
"""
from pypic.setup.planner import params_from_dict, plan_two_stream
from pypic.setup.runner import main_for

if __name__ == "__main__":
    main_for("two_stream", plan_two_stream, params_from_dict, doc=__doc__)
