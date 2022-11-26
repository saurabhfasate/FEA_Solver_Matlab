model BasicBuilder -ndm 3 -ndf 3

# ===
# Von_Misses_truss_Criesfield
# ===

# ==============================================
# Define Nodes and Elements
# ==============================================
node 1  0.0  0.0 -9.5
node 2 -10.0 0.0 0.0
node 3  10.0 0.0 0.0
node 4   0.0 0.0 0.5


# uniaxialMaterial Elastic $matTag $E 
uniaxialMaterial Elastic 1 5.0e6
uniaxialMaterial Elastic 2 10000

# element corotTruss $eleTag $iNode $jNode $A $matTag <-rho $rho> <-cMass $cFlag> <-doRayleigh $rFlag>
element corotTruss 1 2 4 1 1
element corotTruss 2 3 4 1 1
# element truss      3 1 4 1 2
# ==============================================
# Define Fixity
# ==============================================
fix 1 1 1 1
fix 2 1 1 1
fix 3 1 1 1
fix 4 1 1 0

# Set a parameter for the axial load
set P 2000;                

# # Create a Plain load pattern with a Linear TimeSeries
pattern Plain 1 "Linear" {
        # Create nodal loads at nodes 2
        #    nd       FX   FY  FZ 
        load 4   0.0 0.0 [expr -$P]  
}
recorder Node -file "dispopensees.txt" -time -node 4 -dof 1 3 disp
recorder Node -file "reactopensees.txt" -time -node 1 2 3 -dof 1 3 reaction

# initialize in case we need to do an initial stiffness iteration
# initialize
# ------------------------------
# End of model generation
# ------------------------------
# ------------------------------
# Start of analysis generation
# ------------------------------
system FullGeneral
constraints Plain
numberer Plain
# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-6 and a max number of iterations of 10
test NormDispIncr 1.0e-6 50 0
# Create the solution algorithm, a Newton-Raphson algorithm
algorithm Newton
# Create the integration scheme, the LoadControl scheme using steps of 0.1 
integrator LoadControl 0.025
# Create the analysis object
analysis Static
# # ------------------------------
# # End of analysis generation
# # ------------------------------
# # ------------------------------
# # Finally perform the analysis
# # ------------------------------
# # perform the gravity load analysis, requires 10 steps to reach the load level
analyze 40
puts "Gravity load analysis completed";

# print node  1 2 3 4 5
# print ele 1 2 3 4
# print -node 2