# Vendored third-party libraries.
#
# Each is compiled as an OBJECT library so its objects can be baked
# directly into sc_expansion_c. This avoids the need to install / export
# the dependency separately (which would otherwise be required because
# sc_expansion_c is itself a static library that's exported).

# ----- bliss 0.77 -----
# Canonical labeling and automorphism group computation.
# Upstream: https://users.aalto.fi/~tjunttil/bliss/  (LGPL-2.1)
set(BLISS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bliss/src)
add_library(bliss_obj OBJECT
  ${BLISS_DIR}/abstractgraph.cc
  ${BLISS_DIR}/bliss_C.cc
  ${BLISS_DIR}/defs.cc
  ${BLISS_DIR}/digraph.cc
  ${BLISS_DIR}/graph.cc
  ${BLISS_DIR}/orbit.cc
  ${BLISS_DIR}/partition.cc
  ${BLISS_DIR}/uintseqhash.cc
  ${BLISS_DIR}/utils.cc
)
set_property(TARGET bliss_obj PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(bliss_obj PUBLIC $<BUILD_INTERFACE:${BLISS_DIR}>)
# bliss is third-party: silence its warnings, keep it optimized
target_compile_options(bliss_obj PRIVATE -w -O3 -DNDEBUG)
