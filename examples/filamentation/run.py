"""Current-filamentation (Weibel) instability of counter-streaming beams.

python -m pypic.launch -n 8 examples/filamentation/run.py --config examples/filamentation/configs/pair.yaml
python examples/filamentation/analyze.py output/filamentation/pair --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_filamentation  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("filamentation", plan_filamentation, params_from_dict, doc=__doc__)
