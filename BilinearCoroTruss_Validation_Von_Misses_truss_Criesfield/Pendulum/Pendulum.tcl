# ===
# Von_Misses_truss_Criesfield
# ===

# ==============================================
# Define Nodes and Elements
# ==============================================
node 1  0.0 0.0 0
node 2  2.0 0.0 0.

element BilinearCoroTruss 1 1 2 1 2.0e6 2 1000 0

# ==============================================
# Define Fixity
# ==============================================
fix 1 1 1 1
fix 2 0 1 0

# mass $nodeTag (ndf $massValues)
mass 2 1.0 0. 1.
