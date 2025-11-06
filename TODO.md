# TODO List

### 1. Implement Different Dynamics for Ising Model
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

### 2. Documentation
- [ ] Add more examples to docstrings
- [ ] Create tutorial notebooks for common workflows
- [ ] Document all dynamics options once implemented

### 3. Testing
- [ ] Performance benchmarking suite

---

*Last updated: 2025-11-06*
