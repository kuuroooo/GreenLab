Guideline applicability (fixed and indentation-correct)

Applied in separate files:
- G1: Factor repeated expressions in backprop into variables. File: two_hidden_layers_neural_network_G1.py
- G3: Early termination & cached loop-end in training. File: two_hidden_layers_neural_network_G3.py
- G6: Remove redundant buffers; reuse computed outputs. File: two_hidden_layers_neural_network_G6.py
- G7: Vectorized batch prediction. File: two_hidden_layers_neural_network_G7.py

Combined:
- ALL version: two_hidden_layers_neural_network_ALL.py (G1+G3+G6+G7)

Not applied (not meaningful for this single-thread NumPy script):
- G4, G9, G12, G14

All files use fixed example input in example().
