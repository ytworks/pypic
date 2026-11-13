"""Runner of the dipole runs (dipole-design.md): examples/blast_wave with its body and a dipole that dominates
the whole box.

python -m pypic.launch -n 8 research/topics/spherical-shock-injection/run_dipole.py --config <config> --out <dir>

The planner of the example sizes the box by the stopping radius R_B of the uniform field B0. Here B0 is kept
negligible against the dipole (the dipole radius, where the two are equal, lies three box half-widths out), so R_B,
the Alfven Mach number and the beta that the planner prints refer to B0 and mean nothing; the local values are
in dipole-design.md. This runner lowers the planner's limit on the half-width (in R_B) and changes nothing else:
`half_width` and `run` of the config give the box and the length directly.
"""
import os
import sys

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "examples", "blast_wave")
sys.path.insert(0, os.path.abspath(EXAMPLE))

import plan  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

HALF_WIDTH_MIN = 0.001  # in R_B of the negligible B0

if __name__ == "__main__":
    plan.HALF_WIDTH_MIN = HALF_WIDTH_MIN
    main_for("blast_wave", plan.plan_blast_wave, plan.params_from_dict, doc=__doc__)
