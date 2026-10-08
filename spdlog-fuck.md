# Report: The spdlog CMake Target Collision (`spdlog` vs `spdlog::spdlog`)

## Executive summary

`cmake --build build` for `multiplierless-cpp` failed at the **configure** step (not the compile
step) with a CMake ALIAS-target collision:

```
CMake Error at build/_deps/spdlog-src/CMakeLists.txt:187 (add_library):
  add_library cannot create ALIAS target "spdlog::spdlog" because another
  target with the same name already exists.
```

The failure is caused by **EllAlgo v1.6.9** (`luk036/ellalgo-cpp`), which unconditionally
`FetchContent`-declares spdlog and hardcodes the plain `spdlog` target name in its own `CMakeLists.txt`.
When the consumer resolves spdlog from a **system/installed package** (only the namespaced
`spdlog::spdlog` target exists), EllAlgo re-fetches spdlog from source and the two targets collide.

A series of increasingly elaborate workarounds were applied in this repo before the true root cause
was identified. The final state is minimal: spdlog is declared with a plain `CPMAddPackage` exactly
like every other dependency, plus a one-line guard that keeps spdlog source-built (which is what
EllAlgo assumes). The real fix belongs upstream in EllAlgo and is tracked as
[ellalgo-cpp#7](https://github.com/luk036/ellalgo-cpp/issues/7).

## Symptoms

With a system-installed spdlog present (Termux: `/data/data/com.termux/files/usr/lib/cmake/spdlog`)
and the build cache setting `CPM_USE_LOCAL_PACKAGES=ON`, configure failed with a cascade of errors,
all rooted in spdlog's own `CMakeLists.txt` being processed a second time:

- `add_library cannot create target "spdlog" ...`
- `add_library cannot create ALIAS target "spdlog::spdlog" because another target with the same name already exists.`
- `add_library cannot create ALIAS target "spdlog::spdlog_header_only" ...`

In `projgeom-cpp` the same underlying situation surfaced differently: its own CMakeLists called
`get_target_property(SPDLOG_COMPILE_DEFS spdlog ...)` and failed with
`non-existent target "spdlog"` because only `spdlog::spdlog` (imported) existed.

## Background: two CMake representations of spdlog

| spdlog obtained via                | Targets that exist                           |
| ---------------------------------- | -------------------------------------------- |
| built from source (CPM download)   | `spdlog` **and** `spdlog::spdlog` (alias)    |
| system/installed package (find_package) | only `spdlog::spdlog` (imported)          |

- `spdlog::spdlog` is the namespaced, imported target exported by an installed spdlog package. The
  `::` signals an external/imported target and avoids name collisions. Consumers are *supposed* to
  link this.
- `spdlog` is the plain target that spdlog's own `CMakeLists.txt` creates via `add_library(spdlog
  ...)` when it is **built from source**, plus `add_library(spdlog::spdlog ALIAS spdlog)`.

Everything downstream here referenced the **plain** name, which only exists in the source-built case.

## Root cause

`ellalgo-cpp` v1.6.9 `CMakeLists.txt` lines 50-64:

```cmake
# spdlog (bundles fmt) - used by the library's logger
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.17.0
)
set(SPDLOG_INSTALL OFF)
FetchContent_MakeAvailable(spdlog)

get_target_property(SPDLOG_INCLUDE_DIRS spdlog INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(SPDLOG_COMPILE_DEFS spdlog INTERFACE_COMPILE_DEFINITIONS)
```

Two defects:

1. **Unconditional `FetchContent_MakeAvailable(spdlog)`.** There is no
   `find_package(spdlog CONFIG QUIET)` or `if(TARGET spdlog)` guard. If spdlog is already provided by
   the consumer, EllAlgo re-fetches and re-adds it; spdlog's own `add_library(spdlog::spdlog ALIAS
   spdlog)` then collides with the existing imported/alias target.
2. **Hardcoded plain `spdlog` target name.** An installed spdlog only exports `spdlog::spdlog`, so
   even if the re-fetch were suppressed, the two `get_target_property(... spdlog ...)` calls fail.

**Why `CPM_USE_LOCAL_PACKAGES=ON` triggers it:** CPM's documented opt-in makes `CPMAddPackage`
resolve packages via `find_package` first. For spdlog it then returns an imported `spdlog::spdlog`
and — crucially — does **not** register the FetchContent "populated" override that CPM installs when
it *builds* a package from source. A raw `FetchContent` consumer like EllAlgo therefore sees spdlog as
unpopulated and re-fetches it, causing the collision.

## Timeline of attempted fixes

1. **Toggle `CPM_USE_LOCAL_PACKAGES OFF` around spdlog only** (`set(...)` / `unset(...)`).
   Worked, but mutates a global option around a single package; judged "ugly and confusing".

2. **Bridge target + FetchContent populated-marking.** Create a plain `spdlog` INTERFACE target that
   forwards to the imported `spdlog::spdlog` (copying include dirs and compile definitions), and mark
   the spdlog FetchContent as already-populated so EllAlgo's `FetchContent_MakeAvailable(spdlog)`
   becomes a no-op.
   - First implemented by calling CPM's private `cpm_override_fetchcontent()`.
   - Then replaced with direct `define_property`/`set_property` of FetchContent's internal global
     properties `_FetchContent_spdlog_{sourceDir,binaryDir,populated}` to avoid depending on a
     CPM-internal function. Two non-obvious facts surfaced here:
     - `FetchContent_SetPopulated()` (CMake >= 3.24) is documented as **dependency-provider-only**
       and is explicitly unsupported outside that context, so it cannot be used.
     - `FetchContent_GetProperties()` reads the populated flag with `get_property(... DEFINED)`,
       which only returns true if the property was registered via `define_property` first — a bare
       `set_property` silently does not register it.

3. **Made the blocks byte-identical across `multiplierless-cpp` and `projgeom-cpp`.** Verified
   identical via `diff`; however, the populated-marking is *load-bearing* only in multiplierless
   (EllAlgo uses raw FetchContent) and is *inert* in projgeom (Fractions re-declares spdlog via
   `CPMAddPackage`, which CPM short-circuits through its own already-added registry before any
   FetchContent logic runs). Removing the marking in projgeom was verified to change nothing.

4. **"Keep only `spdlog::spdlog`".** 
   - projgeom: fully namespaced-only — removed the entire bridge; its own `get_target_property`
     now reads `spdlog::spdlog`.
   - multiplierless: our own references now use `spdlog::spdlog`, but the bridge *had* to remain
     because EllAlgo hardcodes the plain `spdlog` target.

5. **Realization: spdlog is not the problem — EllAlgo is.** All the workarounds were compensating
   for EllAlgo's defect. The complexity collapsed once the real root cause was named.

## Final solution

In both `multiplierless-cpp` and `projgeom-cpp`, all spdlog special-casing was deleted. spdlog is
now a plain `CPMAddPackage`, identical in shape to EllAlgo/Csd/Ginger/Fractions:

```cmake
# spdlog is used by the compiled logger source. SPDLOG_FMT_EXTERNAL is always set so spdlog uses the
# external fmt and avoids duplicate symbols.
CPMAddPackage(
  NAME spdlog
  GIT_TAG v1.17.0
  GITHUB_REPOSITORY gabime/spdlog
  OPTIONS "SPDLOG_INSTALL YES" "SPDLOG_FMT_EXTERNAL YES"
)
```

plus a single defensive line placed after the CPM include (`CMakeLists.txt:57-60`):

```cmake
# Always build dependencies from source. EllAlgo declares spdlog itself via FetchContent and assumes
# the source-built spdlog targets; the system-package opt-in (CPM_USE_LOCAL_PACKAGES) would swap in an
# imported spdlog::spdlog and break that assumption, so spdlog needs no special handling here.
set(CPM_USE_LOCAL_PACKAGES OFF)
```

With spdlog built from source, the plain `spdlog` target and the `spdlog::spdlog` alias both exist,
CPM marks the fetch populated, and EllAlgo's `FetchContent_MakeAvailable(spdlog)` is a no-op.

Our own compile-definition propagation reads the namespaced target
(`CMakeLists.txt:191`):

```cmake
get_target_property(SPDLOG_COMPILE_DEFS spdlog::spdlog INTERFACE_COMPILE_DEFINITIONS)
```

### Residual EllAlgo workaround

`CMakeLists.txt:99-104` still works around a related EllAlgo quirk: EllAlgo compiles spdlog headers
in its library but does not link fmt, and with `SPDLOG_FMT_EXTERNAL` those headers include fmt, so
fmt's include dir must be added to EllAlgo manually:

```cmake
if(TARGET EllAlgo)
  get_target_property(FMT_INCLUDE_DIRS fmt::fmt INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(EllAlgo PRIVATE ${FMT_INCLUDE_DIRS})
endif()
```

## Verification

- `multiplierless-cpp`: configure + full build + `./build/MultiplierlessTests`
  — **62 test cases, 199 assertions — 100% pass (0 failed, 0 skipped)**.
- `projgeom-cpp`: configure + full build + `./build/ProjGeomTests`
  — **96 test cases, 129 assertions — 100% pass (0 failed, 0 skipped)**.
- spdlog builds from source (`-- CPM: Adding package spdlog@1.17.0`, `-- Build spdlog: 1.17.0`); the
  system spdlog is no longer pulled in.

## Issues created

- **Upstream bug (the actual fix belongs here):**
  [ellalgo-cpp#7](https://github.com/luk036/ellalgo-cpp/issues/7) — "CMakeLists: unconditional
  FetchContent of spdlog collides with consumer-provided spdlog". Suggested patch: guard the fetch
  with `find_package(spdlog CONFIG QUIET)` (falling back to FetchContent only if absent) and use the
  namespaced `spdlog::spdlog` target.
- **Follow-up in this repo:**
  [multiplierless-cpp#5](https://github.com/luk036/multiplierless-cpp/issues/5) — "Follow-up: remove
  EllAlgo spdlog workarounds once ellalgo-cpp#7 is fixed". Once upstream is fixed:
  1. delete `set(CPM_USE_LOCAL_PACKAGES OFF)` (restores the README's system/vcpkg/Conan opt-in),
  2. re-evaluate the fmt include-dir workaround (lines 99-104),
  3. verify both `-DCPM_USE_LOCAL_PACKAGES=ON` and default configures build and pass tests.

## Lessons learned

- **Diagnose before patching.** The first error trace pointed straight at `build/_deps/spdlog-src`
  being re-added by EllAlgo's `FetchContent` — that was the root cause. Several rounds of
  increasingly clever workarounds in this repo would have been unnecessary had the upstream defect
  been reported first.
- **Third-party targets come in two flavors.** A plain target (created by source builds) and a
  namespaced imported target (from installed packages) are not interchangeable in code that
  hardcodes a plain name.
- **CPM's local-package path skips the FetchContent override.** Opting into `CPM_USE_LOCAL_PACKAGES`
  is incompatible with dependencies that assume source-built targets.
- **`FetchContent_SetPopulated` is provider-only**, and FetchContent's internal `DEFINED` property
  semantics require `define_property` — details that make hand-rolling the "populated" override
  fragile and a strong signal the problem belongs upstream.
