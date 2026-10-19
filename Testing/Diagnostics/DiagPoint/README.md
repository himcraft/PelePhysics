# Fixed-point diagnostic tests

These tests compile the production `DiagPoint` and the PelePhysics
Diagnostics base against real AMReX, without chemistry or the full flow solver.
Use a compiler supported by the checked-out AMReX (currently GCC >= 11 or
Clang >= 14).

```sh
cmake -S Testing/Diagnostics/DiagPoint -B tmp_build_dir/diagpoint-2d \
  -DAMReX_SPACEDIM=2 -DAMReX_MPI=ON -DAMReX_OMP=OFF
cmake --build tmp_build_dir/diagpoint-2d -j 8
ctest --test-dir tmp_build_dir/diagpoint-2d --output-on-failure
```

Repeat with a separate build directory and `AMReX_SPACEDIM=3`. Set
`AMReX_MPI=OFF` to run the serial variants. MPI tests launch two ranks and
deliberately assign the sampled fine cell to the non-I/O rank. Each scenario
has a separate working directory; only its own `probe.dat` fixture is removed
at test startup.

Fixtures verify coarse/fine selection, field order, physical boundaries,
nonzero physical origins and domain indices, regridding, restart append and
duplicate suppression, malformed/incompatible histories, invalid inputs,
and rejection of fully covered EB cells through the actual `vfrac` field.
It also verifies sampling without a volume-fraction field; EB rejection uses
the configurable `volume_fraction_field`. The standalone fixture does not
construct an EB geometry or exercise a GPU. Production initialization and
checkpoint-restart tests belong to the host solver (PeleC provides one in
`Tests/DiagPoint`).
