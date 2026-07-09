"""Reasoning-guided structural repair.

The reasoning step reads the structured residual summary and infers a mechanism:
a knockdown that raises a target the model lowers implies a sign flip or an
intermediate repressor; a residual that is non-monotone in two regulators implies
replacing an additive term with a product-term gate. It emits a structural edit.
Re-fit, re-diagnose, re-read, and iterate to a plateau or a budget. Implemented in
Step 3, where it escalates to a logic gate on a non-monotone module.
"""
raise NotImplementedError
