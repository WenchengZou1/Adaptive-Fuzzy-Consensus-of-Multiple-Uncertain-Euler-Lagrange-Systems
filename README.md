Adaptive Fuzzy Consensus of Multiple Uncertain Euler-Lagrange Systems With Sampled-Data Output Interactions

Abstract—In this paper, the consensus problem for a class
of multi-agent systems is investigated, where each agent is
described by an Euler-Lagrange system and only sampled-data
output interaction among agents is allowed. In addition, the
Euler-Lagrange system considered possesses a higher degree of
uncertainty; specifically, the regression matrix is also unknown.
The intricate heterogeneity, inherent uncertainties, and stringent
constraints of information interactions compel us to develop a
new protocol. The protocol is synthesized by organically integrating 
theories and techniques such as graph theory, sampled-data
cooperative control, fixed-time control, and fuzzy approximation.
Each agent is equipped with a first-order difference trajectory
generator that updates exclusively at sampling instants using the
output information received from its neighbors. The virtual trajectory 
serves as the tracking target for the corresponding agent
systems output. A fixed-time fuzzy adaptive tracker is proposed
for each agent to track the generated trajectory. The problem of
analyzing the protocol’s effectiveness is deconstructed into two
coupled subproblems: a synchronization problem for a perturbed
first-order differentiator and a practical fixed-time tracking
problem for a single Euler-Lagrange system. This problem is
adequately addressed primarily based on the Lyapunov function
method. Finally, the effectiveness of the proposed protocol is
validated through numerical simulation and ROS experiment.
