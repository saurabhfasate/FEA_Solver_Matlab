# ===
# SDOF_AKChopra Sinusoidal Force
# ===

# ==============================================
# Define Nodes and Elements
# ==============================================
node 1 0.0 0.0 0.
node 2 1.0 0.0 0.

element BilinearSpring 1 1 2 10 0.75 0
# ==============================================
# Define Fixity
# ==============================================
fix 1 1 1 1
fix 2 0 1 1

# # mass $nodeTag (ndf $massValues)
# set M 0.2533
mass 2 0.2533 0. 0.

# # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# # Time History ANALYSIS 
# # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# # Damping
# # ----------------
# set T1 1.0

# set w1 [expr 2*3.14/$T1]

# set z 0.05
# set a0 [expr 2*$z*$w1]
# set a1 [expr 2*$z/$w1]
# puts "$a0 $a1"
# rayleigh $a0 $a1 0. 0.

# timeSeries Path 1 -time {0.0 0.1 0.2 0.3 0.4 0.5 0.6} -values {0.0 5.0 8.66 10.0 8.66 5.0 0.0}
# pattern Plain 1 1 {
# load 2 1.0 0. ;
# }

# set dt 0.1
# set nsteps [expr 10]

# recorder Node -file dydisp.txt -time -node 2 -dof 1 2 disp

# # Run a transient analysis with no loading (to record free vibration)
# # -----------------------------------------------------------------------
# system FullGeneral; # System of equations solver
# constraints Plain; # Constraint handler
# numberer Plain; # DOF numberer
# algorithm Newton
# test NormDispIncr 1.0e-9 100 0
# # test EnergyIncr 1.0e-12 100 0; # Convergence test tolerance maxIter displayCode
# integrator Newmark 0.5 0.25; 
# analysis Transient;# Create the analysis object
# set ok [analyze $nsteps $dt]

# if {$ok == 0} {
   # puts "Time History analysis completed SUCCESSFULLY";
# } else {
   # puts "Time History analysis FAILED";    
# }
