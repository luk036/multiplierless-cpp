set_version("1.0.6")
set_license("MIT")
set_policy("check.target_package_licenses", false)
set_languages("c++20")

add_rules("mode.debug", "mode.release", "mode.coverage")
add_requires("nlohmann_json", "doctest", "fftw", "fmt", "spdlog", "cxxopts")

if is_mode("release") then
    set_optimize("fastest")
end

if is_mode("coverage") then
    add_cxflags("-ftest-coverage", "-fprofile-arcs", {force = true})
end

if is_plat("linux") then
    set_warnings("all", "error")
elseif is_plat("windows") then
    add_cxflags("/EHsc /utf-8 /W4 /WX /wd4702", {force = true})
    add_syslinks("ws2_32")
end

-- sibling projects: EllAlgo (Arr), Ginger (root finding), Csd (multiplier synthesis)
local ellalgo_dir = path.join(os.projectdir(), "../ellalgo-cpp")
local csd_dir = path.join(os.projectdir(), "../csd-cpp")
local ginger_dir = path.join(os.projectdir(), "../ginger-cpp")

-- RapidCheck from CMake build (optional)
local mode_dir = is_mode("release") and "Release" or "Debug"
local rc_dir = path.join(os.projectdir(), "build", "_deps", "rapidcheck-src")
local rc_lib_dir = path.join(os.projectdir(), "build", "_deps", "rapidcheck-build", mode_dir)
local rc_lib = is_plat("windows")
    and path.join(rc_lib_dir, "rapidcheck.lib")
    or path.join(rc_lib_dir, "librapidcheck.a")

target("Ginger")
    set_kind("static")
    -- aberth.hpp publicly includes <lds/lds.hpp>, so the Lds include dir must be public
    add_includedirs(path.join(ginger_dir, "include"), {public = true})
    add_includedirs(path.join(os.projectdir(), "../lds-cpp/include"), {public = true})
    add_files(path.join(ginger_dir, "source/aberth.cpp"),
              path.join(ginger_dir, "source/autocorr.cpp"),
              path.join(ginger_dir, "source/rootfinding.cpp"))
    set_group("Dependencies")

target("Csd")
    set_kind("static")
    set_languages("c++17")
    add_includedirs(path.join(csd_dir, "include"), {public = true})
    add_files(path.join(csd_dir, "source/csd.cpp"),
              path.join(csd_dir, "source/csd_multiplier.cpp"),
              path.join(csd_dir, "source/lcsre.cpp"))
    set_group("Dependencies")

target("EllAlgo")
    set_kind("static")
    add_includedirs(path.join(ellalgo_dir, "include"), {public = true})
    add_files(path.join(ellalgo_dir, "source/*.cpp"))
    add_packages("fmt", "spdlog")
    set_group("Dependencies")

target("Multiplierless")
    set_kind("static")
    add_includedirs("include", {public = true})
    add_includedirs(path.join(ellalgo_dir, "include"), {public = true})
    add_files("source/*.cpp")
    add_deps("Ginger", "Csd", "EllAlgo")
    add_packages("fftw", "fmt", "spdlog")

target("MultiplierlessTests")
    add_deps("Multiplierless", "EllAlgo")
    add_includedirs("include")
    add_files("test/source/*.cpp")
    add_packages("fftw", "fmt", "spdlog", "doctest")
    add_tests("default")

    if os.isdir(rc_dir) and os.isfile(rc_lib) then
        add_includedirs(path.join(rc_dir, "include"))
        add_linkdirs(rc_lib_dir)
        add_links("rapidcheck")
        add_defines("RAPIDCHECK_H")
    end

target("MultiplierlessStandalone")
    add_deps("Multiplierless", "EllAlgo")
    add_includedirs("include")
    add_files("standalone/source/main.cpp")
    add_packages("cxxopts", "fmt", "spdlog", "fftw")

target("LoggerExample")
    add_deps("Multiplierless")
    add_includedirs("include")
    add_files("standalone/source/logger_example.cpp")
    add_packages("fmt", "spdlog", "fftw")

target("FirDesign")
    add_deps("Multiplierless", "EllAlgo")
    add_includedirs("include")
    add_files("standalone/source/fir_design.cpp")
    add_packages("nlohmann_json", "fmt", "spdlog", "fftw")
