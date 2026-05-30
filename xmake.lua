set_languages("c++20")
set_policy("build.optimization.lto", true)

add_rules("mode.debug", "mode.release", "mode.coverage")

if is_mode("release") then
	set_optimize("fastest")
end

if is_mode("coverage") then
	add_cxflags("-ftest-coverage", "-fprofile-arcs", {force = true})
end

if is_plat("linux") then
	set_warnings("all", "error")
	-- add_cxflags("-Wconversion", {force = true})
elseif is_plat("windows") then
	add_cxflags("/EHsc /utf-8 /W4 /WX /wd4459 /wd4819 /wd4996 /wd4267 /wd4244", {force = true})
end

-- Local dependencies from CMake build/_deps (offline-friendly)
local deps = path.join(os.projectdir(), "build/_deps")
local mode_dir = is_mode("release") and "Release" or "Debug"

-- EllAlgo-cpp include path (core Arr lives here now)
local ellalgo_inc = path.join(os.projectdir(), "../ellalgo-cpp/include")

-- FFTW from conda environment (adjust for your system)
local fftw_prefix = "D:/scoop/apps/miniconda3/current/envs/cppflow/Library"

-- Header-only: doctest
local doctest_dir = path.join(deps, "doctest-src")
local doctest_h = path.join(doctest_dir, "doctest", "doctest.h")

-- Header-only: cxxopts
local cxxopts_dir = path.join(deps, "cxxopts-src")
local cxxopts_h = path.join(cxxopts_dir, "include", "cxxopts.hpp")

-- Compiled: fmt
local fmt_dir = path.join(deps, "fmt-src")
local fmt_lib_dir = path.join(deps, "fmt-build", mode_dir)

-- Compiled: spdlog
local spdlog_dir = path.join(deps, "spdlog-src")
local spdlog_lib_dir = path.join(deps, "spdlog-build", mode_dir)

-- EllAlgo source
local ellalgo_dir = path.join(deps, "ellalgo-src")

-- RapidCheck from CMake build (optional)
local rc_dir = path.join(deps, "rapidcheck-src")
local rc_lib_dir = path.join(deps, "rapidcheck-build", mode_dir)
local rc_lib = is_plat("windows")
	and path.join(rc_lib_dir, "rapidcheck.lib")
	or  path.join(rc_lib_dir, "librapidcheck.a")

target("EllAlgo")
	set_kind("static")
	add_includedirs(path.join(ellalgo_dir, "include"), {public = true})
	add_includedirs(path.join(fmt_dir, "include"), {public = true})
	add_includedirs(path.join(spdlog_dir, "include"), {public = true})
	add_files(path.join(ellalgo_dir, "source/*.cpp"))
	set_group("Dependencies")

target("Multiplierless")
	set_kind("static")
	add_includedirs("include", {public = true})
	add_includedirs(ellalgo_inc, {public = true})
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"), {public = true})
	add_includedirs(path.join(fftw_prefix, "include"), {public = true})
	add_includedirs(path.join(fmt_dir, "include"), {public = true})
	add_includedirs(path.join(spdlog_dir, "include"), {public = true})
	add_files("source/*.cpp")
	add_deps("EllAlgo")
	add_linkdirs(path.join(fftw_prefix, "lib"))
	add_linkdirs(fmt_lib_dir)
	add_linkdirs(spdlog_lib_dir)
	add_links("fftw3", "fftw3f", "fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end

target("MultiplierlessTests")
	set_kind("binary")
	add_deps("Multiplierless", "EllAlgo")
	add_includedirs("include")
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"))
	add_includedirs(path.join(fftw_prefix, "include"))
	add_includedirs(path.join(fmt_dir, "include"))
	add_includedirs(path.join(spdlog_dir, "include"))
	if os.isfile(doctest_h) then
		add_includedirs(path.join(doctest_dir))
	end
	add_files("test/source/*.cpp")
	add_linkdirs(path.join(fftw_prefix, "lib"))
	add_linkdirs(fmt_lib_dir)
	add_linkdirs(spdlog_lib_dir)
	add_links("fftw3", "fftw3f", "fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end
	add_tests("default")

	-- RapidCheck from CMake build/_deps
	if os.isdir(rc_dir) and os.isfile(rc_lib) then
		add_includedirs(path.join(rc_dir, "include"))
		add_linkdirs(rc_lib_dir)
		add_links("rapidcheck")
		add_defines("RAPIDCHECK_H")
	end

target("MultiplierlessStandalone")
	set_kind("binary")
	add_deps("Multiplierless", "EllAlgo")
	add_includedirs("include")
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"))
	if os.isfile(cxxopts_h) then
		add_includedirs(path.join(cxxopts_dir, "include"))
	end
	add_files("standalone/source/main.cpp")
	add_linkdirs(fmt_lib_dir)
	add_linkdirs(spdlog_lib_dir)
	add_links("fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end

target("LoggerExample")
	set_kind("binary")
	add_deps("Multiplierless")
	add_includedirs("include")
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"))
	add_files("standalone/source/logger_example.cpp")
	add_linkdirs(fmt_lib_dir)
	add_linkdirs(spdlog_lib_dir)
	add_links("fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end
