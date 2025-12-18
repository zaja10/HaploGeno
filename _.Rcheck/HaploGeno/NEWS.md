# HaploGeno 0.3.0

- **New Features**:
  - Added `filter_monomorphic()` to automatically identify and exclude zero-variance markers, improving analysis stability.
  - Implemented `save_project()` and `load_haplo_project()` for robust persistence of analysis objects and their associated file-backed matrices.
  - Added a demo dataset (`load_demo_data()`) to facilitate immediate testing and visualization.
- **Improvements**:
  - Enhanced `print()` method to provide a detailed summary of the `HaploObject` state.
  - Hardened `load_pheno()` to enforce dimension consistency between genotypes and phenotypes.
  - Refined imputation logic to ensure integer-rounded values.

# HaploGeno 0.2.0

- **Performance**:
  - Optimized haplotype encoding using C++ (`encode_block_fast`), significantly reducing processing time for large datasets.
  - Integrated Conjugate Gradient solver for KRR to handle large sample sizes efficiently.
- **Methodology**:
  - Implemented LD-based block definition (`define_blocks_ld`) with customizable $r^2$ thresholds.
  - Added support for local GEBV calculation and significance testing.

# HaploGeno 0.1.0

- **Initial Release**:
  - Basic `HaploObject` structure with `bigstatsr` backend.
  - Support for fixed-window haplotype blocking.
  - Initial implementation of Kernel Ridge Regression for global prediction.
