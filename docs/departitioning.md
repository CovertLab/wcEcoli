* implement gillespie
* enforce time step independence
** go through processes and identify which need to be reworked
*** complexation - implement with gillespie
*** two component + equilibrium - combine and replace with gillespie
*** metabolism - apply John's prototype to meet demand (PI controller)
* remove partitioning
** implement updates as messages
** implement update filters to process updates
** if invalid state is reached, throw updates away and recalculate with half time step
** else if valid state is reached, apply updates and continue