// Debug-only: spectral_fact_root tests deferred to Release suite.
// SIGABRT in ginger's Aberth solver (checked-iterator / assertion
// in MSVC Debug CRT).  The root-based function is verified by
// Release builds; the FFT-based spectral_fact tests run in Debug.
