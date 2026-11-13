"""Runner of the strong-shock runs (strong-shock-run.md): examples/blast_wave in a box smaller than 1.3 R_B.

python -m pypic.launch -n 8 research/topics/spherical-shock-injection/run_piston.py --config <config> --out <dir>

The planner of the example sizes the box by the stopping radius R_B = R_m M_A^(2/3) and refuses a half-width
below 1.3 R_B: room for the front until the debris has stopped. A strong shock needs the opposite: a debris
cloud so heavy that it does not slow down inside the box (R_m of the order of the half-width), and a run that
ends before the front reaches the faces. This runner lowers that one limit and changes nothing else; the
config gives `half_width` and `run` (both in units of R_B) so that the run ends inside the box. The faces are
open and the ambient plasma keeps entering, as in the example.
"""
import os
import sys

EXAMPLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "examples", "blast_wave")
sys.path.insert(0, os.path.abspath(EXAMPLE))

import plan  # noqa: E402

from pypic.setup.runner import main_for  # noqa: E402

HALF_WIDTH_MIN = 0.3  # in R_B; the example's 1.3 holds the front of a cloud that stops, this one does not stop

if __name__ == "__main__":
    plan.HALF_WIDTH_MIN = HALF_WIDTH_MIN
    main_for("blast_wave", plan.plan_blast_wave, plan.params_from_dict, doc=__doc__)
