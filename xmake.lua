set_languages("c++20")
-- set_policy("build.optimization.lto", true)

add_rules("mode.debug", "mode.release", "mode.coverage")

add_requires("nlohmann_json")
add_requires("doctest", { alias = "doctest" })
add_requires("fmt", { alias = "fmt" })
add_requires("spdlog", { alias = "spdlog" })
add_requires("cxxopts", { alias = "cxxopts" })

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
	add_cxflags("/EHsc /utf-8 /W4 /WX /wd5285 /wd4996 /wd4267 /wd4244", {force = true})
end

-- EllAlgo-cpp include path (core Arr lives here now)
local ellalgo_inc = path.join(os.projectdir(), "../ellalgo-cpp/include")

-- FFTW from conda environment (adjust for your system)
local fftw_prefix = "D:/scoop/apps/miniconda3/current/envs/cppflow/Library"

-- EllAlgo source (sibling project)
local ellalgo_dir = path.join(os.projectdir(), "../ellalgo-cpp")

-- RapidCheck from CMake build (optional)
local mode_dir = is_mode("release") and "Release" or "Debug"
local rc_dir = path.join(os.projectdir(), "build", "_deps", "rapidcheck-src")
local rc_lib_dir = path.join(os.projectdir(), "build", "_deps", "rapidcheck-build", mode_dir)
local rc_lib = is_plat("windows")
	and path.join(rc_lib_dir, "rapidcheck.lib")
	or  path.join(rc_lib_dir, "librapidcheck.a")

local csd_dir = path.join(os.projectdir(), "../csd-cpp")
local ginger_dir = path.join(os.projectdir(), "../ginger-cpp")
local fftw_prefix = "D:/scoop/apps/miniconda3/current/envs/cppflow/Library"

target("Ginger")
	set_kind("static")
	set_languages("c++17")
	add_includedirs(path.join(ginger_dir, "include"), {public = true})
	add_files(path.join(ginger_dir, "source/aberth.cpp"),
	          path.join(ginger_dir, "source/autocorr.cpp"),
	          path.join(ginger_dir, "source/rootfinding.cpp"))
	set_group("Dependencies")

target("Csd")
	set_kind("static")
	set_languages("c++14")
	add_includedirs(path.join(csd_dir, "include"), {public = true})
	add_files(path.join(csd_dir, "source/csd.cpp"), path.join(csd_dir, "source/csd_multiplier.cpp"),
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
	add_includedirs(ellalgo_inc, {public = true})
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"), {public = true})
	add_includedirs(path.join(fftw_prefix, "include"), {public = true})
	add_files("source/*.cpp")
	add_deps("Ginger", "Csd", "EllAlgo")
	add_linkdirs(path.join(fftw_prefix, "lib"))
	add_links("fftw3", "fftw3f")
	add_packages("fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end

target("MultiplierlessTests")
	set_kind("binary")
	add_deps("Multiplierless", "EllAlgo")
	add_includedirs("include")
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"))
	add_files("test/source/*.cpp")
	add_links("fftw3", "fftw3f")
	add_packages("fmt", "spdlog", "doctest")
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
	add_files("standalone/source/main.cpp")
	add_packages("cxxopts", "fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end

target("LoggerExample")
	set_kind("binary")
	add_deps("Multiplierless")
	add_includedirs("include")
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"))
	add_files("standalone/source/logger_example.cpp")
	add_packages("fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end

target("FirDesign")
	set_kind("binary")
	add_deps("Multiplierless", "EllAlgo")
	add_includedirs("include")
	add_includedirs(path.join(os.projectdir(), "build/PackageProjectInclude"))
	add_files("standalone/source/fir_design.cpp")
	add_packages("nlohmann_json")
	add_packages("fmt", "spdlog")
	if is_plat("windows") then
		add_syslinks("ws2_32")
	end
