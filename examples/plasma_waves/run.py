"""Thermal plasma waves: omega-k spectra of a (magnetised) Maxwellian plasma vs linear theory.

python examples/plasma_waves/run.py --config examples/plasma_waves/configs/unmagnetized.yaml
python -m pypic.launch -n 8 examples/plasma_waves/run.py --config .../perpendicular.yaml
python examples/plasma_waves/analyze.py output/plasma_waves/unmagnetized --movie
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_waves  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("plasma_waves", plan_waves, params_from_dict, doc=__doc__)
