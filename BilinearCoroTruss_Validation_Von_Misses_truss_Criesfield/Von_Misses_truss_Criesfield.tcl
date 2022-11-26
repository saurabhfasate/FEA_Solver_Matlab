# ===
# Von_Misses_truss_Criesfield
# ===

# ==============================================
# Define Nodes and Elements
# ==============================================
node 1  0.0  0.0 0
node 2 -10.0 0.0 0.
node 3  10.0 0.0 0.
node 4   0.0 0.0 0.5

element CoroTruss 1 2 4 1 5e6 10.0124 1000 0
element CoroTruss 2 3 4 1 5e6 10.0124 1000 0
element Spring 3 1 4 2000

# ==============================================
# Define Fixity
# ==============================================
fix 1 1 1 1
fix 2 1 1 1
fix 3 1 1 1
fix 4 0 1 0

# mass $nodeTag (ndf $massValues)
# mass 2 0.2533 0. 0.
