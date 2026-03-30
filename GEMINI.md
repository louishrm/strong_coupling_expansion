# Context: Strong Coupling Expansion (Hubbard Model)

## 1. Domain Scope
- **Physics:** Perturbative expansion of the pure Hubbard model on a 2D square lattice in powers of hopping $t$.
- **Method:** Strong Coupling / High-Temperature Expansion.
- **Key Constraint:** **NO Wick's Theorem.** Diagrams are products of atomic cumulants (connected Green's functions), not Gaussian contractions.
- **Objective:** Compute free energy coefficients.

## 2. Tech Stack
- **Language:** C++ (Standard C++17 or higher recommended).
- **Core Library:** **TRIQS** (Toolbox for Research on Interacting Quantum Systems).
- **Integration:** Monte Carlo methods for integral evaluation.
- **Data Structures:** - Diagrams represented as **Graphs**.
  - Cumulants are stored and retrieved via **Memoization**.

## 3. Implementation Rules
- **Memoization:** Ensure the memoization key generation (hashing of graphs/diagrams) is collision-free and extremely fast. 
- **Monte Carlo:**
  - Ensure the random number generator is statistically sound (e.g., `std::mt19937` or TRIQS generators).
  - Hot loops must be zero-overhead.

## 4. Optimization Priorities (Speed > Memory)
- **Memory:** Minimize heap allocations inside the Monte Carlo loop. Use pre-allocated buffers or stack memory.
- **Passing:** Pass heavy objects (graphs, large arrays) by `const &`.
- **Parallelism:** If suggesting loops, consider OpenMP (`#pragma omp`) thread safety, specifically regarding the memoization cache.
- **Inlining:** Encourage inlining for small graph traversal helper functions.

## 5. Interaction Guidelines
- **Analysis:** Before suggesting a fix, analyze the existing code style. Do not rewrite working logic into a different paradigm unless explicitly asked.
- **Math:** Use LaTeX for physics explanations.
- **Verification:** When modifying code, explicitly verify that the changes pass the tests.
- **Changes** When suggesting changes, first make sure that you explain in good detail what change you are making, how you are doing it and why you are doing it. 

## 6. Code Style
- **Classes and Structs** names are always in PascalCase, header files only include the declaration of members and functions. Implementation is always done in .cpp including constructor. Members referenced using this->. 

## 7. TODO: Diagram Logic Refactor (Diagram2 & Vertex)

### 7.1 Vertex Hierarchy & Management
- [ ] **Vertex Cache System**: 
  - At order $n$, maintain a vector of `VertexType` (global cache) for all possible vertex orders from $1$ to $n/2$.
  - Ensure `Diagram2` uses these shared `VertexType` pointers to create `VertexInstance` objects.
  - Implement the "sticky note" logic: `VertexInstance` checks its local cache first, then delegating to `VertexType`'s global cache.

### 7.2 Staggered-Dimer Lattice Embedding ($N_{sites}=2$)
- [ ] **Triangular Lattice Walk**:
  - Replace `Graph::compute_free_multiplicity` with a lattice walk that embeds the graph on a staggered triangular lattice (formed by tiling dimers on the square lattice).
  - Enforce connectivity constraints: A hopping line must connect site $i$ of the outgoing dimer to site $j$ of the incoming dimer where $i \neq j$.
- [ ] **Configuration Extraction**:
  - Collect all valid sequences of fermion creation/annihilation operators (spin/site indices) resulting from the embeddings.
  - For $N_{sites}=1$, simplify this to the sequence of fermion operators for spin species.

### 7.3 Symmetry & Orbit Reduction
- [ ] **Dual-Layer Symmetry Engine**:
  - **Internal (Physics)**: Implement the SU(2) spin-flip and the Dimer Site-Swap ($1 \leftrightarrow 2$) symmetry of the cumulant $C(1,2) = C(2,1)$.
- [ ] **Weight Management**:
  - Ensure the `OrbitalConfiguration.weight` correctly compounds both the spatial multiplicity and the spin/site multiplicity.
  - Divide by the Graph's internal automorphism count to avoid overcounting topological symmetries.

### 7.4 MCMC Performance Optimization
- [ ] **Incremental Re-evaluation (Dirty Tracking)**:
  - Optimize the MCMC loop to propose only one time-change ($d\tau$) at a time.
  - Identify the at most 2 vertices affected by the time change and mark their `VertexInstance` as `dirty`.
  - In `Diagram2::evaluate`, re-evaluate only the `dirty` vertices using the hierarchical cache (Local -> Global -> Solver).

