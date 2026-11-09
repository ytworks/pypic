"""Test-particle orbits with exact relativistic solutions.

python examples/particle_orbits/run.py --config examples/particle_orbits/configs/plane_wave.yaml
python examples/particle_orbits/run.py --config ... --set dt_scale=0.5 --out output/particle_orbits/plane_wave_half
python examples/particle_orbits/analyze.py output/particle_orbits/plane_wave \
    --refined output/particle_orbits/plane_wave_half
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plan import params_from_dict, plan_orbits  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

if __name__ == "__main__":
    main_for("particle_orbits", plan_orbits, params_from_dict, doc=__doc__)
