"""Blast wave around a magnetised body (3-D): a debris shell expands into a magnetised plasma at rest.

python -m pypic.launch -n 8 examples/blast_wave/run.py --config examples/blast_wave/configs/gate.yaml
python examples/blast_wave/analyze.py output/blast_wave/gate
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_blast_wave  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("blast_wave", plan_blast_wave, params_from_dict, doc=__doc__)
