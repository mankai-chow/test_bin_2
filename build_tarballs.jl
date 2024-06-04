# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder, Pkg

name = "FuzzifiED"
version = v"0.5.10"

# Collection of sources required to complete build
sources = [
    DirectorySource("./src/")
]

# Bash recipe for building across all platforms
script = raw"""
cd $WORKSPACE/srcdir
if [[ ${nbits} == 32 ]]; then
    gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./cfs.f90
    gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./bs.f90
    gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./op.f90
    gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./diag.f90
    gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./diag_re.f90
    gfortran -fPIC -shared -larpack -fopenmp -L ${libdir} -o ${libdir}/libfuzzified.$dlext ./*.o
else 
    gfortran -fPIC -fdefault-integer-8 -larpack -fopenmp -L ${libdir} -c ./cfs.f90
    gfortran -fPIC -fdefault-integer-8 -larpack -fopenmp -L ${libdir} -c ./bs.f90
    gfortran -fPIC -fdefault-integer-8 -larpack -fopenmp -L ${libdir} -c ./op.f90
    gfortran -fPIC -fdefault-integer-8 -larpack -fopenmp -L ${libdir} -c ./diag.f90
    gfortran -fPIC -fdefault-integer-8 -larpack -fopenmp -L ${libdir} -c ./diag_re.f90
    gfortran -fPIC -fdefault-integer-8 -shared -larpack -fopenmp -L ${libdir} -o ${libdir}/libfuzzified.$dlext ./*.o
fi
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = supported_platforms()
platforms = expand_gfortran_versions(platforms)

# The products that we will ensure are always built
products = [
    LibraryProduct("libfuzzified", :LibpathFuzzifiED)
]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae"); platforms=filter(!Sys.isbsd, platforms)),
    Dependency(PackageSpec(name="LLVMOpenMP_jll", uuid="1d63c593-3942-5779-bab2-d838dc0a180e"); platforms=filter(Sys.isbsd, platforms)),
    Dependency(PackageSpec(name="Arpack_jll", uuid="68821587-b530-5797-8361-c406ea357684"))
]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies; julia_compat="1.6")
