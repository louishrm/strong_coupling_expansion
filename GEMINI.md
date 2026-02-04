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

## 7. Todo
@./TODO.md
