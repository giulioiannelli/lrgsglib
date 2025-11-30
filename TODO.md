# TODO List

### General 
- [ ] Adjust the `Chronometer` for each dynamical and computational objects.
- [ ] Adjust the logging and in general the verbosity of outputs, from notebooks to programs
- [ ] 

### Implement Different Dynamics for Ising Model
- [ ] Add alternative update rules beyond Glauber dynamics
  - [ ] Metropolis-Hastings dynamics
  - [ ] Heat-bath (Gibbs) dynamics
  - [ ] Kawasaki dynamics (spin exchange)
  - [ ] Swendsen-Wang cluster updates
- [ ] Refactor `IsingDynamics` class to support multiple dynamics types
- [ ] Add parameter to select dynamics type in initialization
- [ ] Update C/C++ backend implementations if needed
- [ ] Add tests for each dynamics type
- [ ] Document performance characteristics of each method
  
### Restyle Serializer
- [ ] Use the L3D_Recon_srun as template.

### Documentation
- [ ] Add more examples to docstrings
- [ ] Create tutorial notebooks for common workflows
- [ ] Document all dynamics options once implemented

### Testing
- [ ] Performance benchmarking suite

---

*Last updated: 2025-11-06*
