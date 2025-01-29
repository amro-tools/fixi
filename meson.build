project('fixi', 'cpp', version:'0.1', default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

compiler = meson.get_compiler('cpp')

cpp_args = []
cpp_args += compiler.get_supported_arguments([
  '-Wno-unused-local-typedefs',  # Ignore unused local typedefs warnings
  '-Wno-array-bounds',           # Suppress out-of-bounds array access warnings
  '-ffast-math',                 # Enable faster, non-IEEE math calculations
  '-fno-finite-math-only',       # Allow Inf and NaN
  # These are based on recommendations from
  # https://youtu.be/_enXuIxuNV4?si=LvtMqPwJ6jYbDY66
  '-fno-semantic-interposition', # Assume no interposition for module functions
  '-fno-plt',                    # Avoid PLT for function calls within shared libs
  '-Bsymbolic',                  # Resolve symbols to local definitions
  '-fopenmp'
])

if cpp_args.length() > 0
  message('Adding compiler flags', cpp_args)
endif

inc = include_directories([['include', 'thirdparty']])

deps = []

fixi_dep = declare_dependency(include_directories : inc)
omp_dep = dependency('openmp')

deps += fixi_dep
deps += omp_dep


executable('fixi_cli', 'src/main.cpp', include_directories : inc, cpp_args : cpp_args, dependencies : deps)

tests = [
    ['Test_Dimer', 'tests/test_dimer.cpp'],
    ['Test_Water', 'tests/test_water.cpp'],
    ['Test_Buckets', 'tests/test_buckets.cpp'],
]

test_inc = []
test_inc += inc

Catch2 = dependency('Catch2', method : 'cmake', modules : ['Catch2::Catch2WithMain', 'Catch2::Catch2'])

foreach t : tests
  exe = executable(t.get(0), t.get(1),
    dependencies : [Catch2, deps],
    include_directories : test_inc,
    cpp_args : cpp_args
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach

# Python module

py = import('python').find_installation('python', modules: ['numpy'])

python_deps = deps + [
  py.dependency(),
  dependency('pybind11'),
]

# for the bindings
py.extension_module(
  'fixicpp',
  sources : [
    'python_bindings/bindings.cpp'
  ],
  # include_directories: ,
  dependencies: python_deps,
  cpp_args : cpp_args,
  install: true,
  subdir: 'pyfixi/'
)

# pyflowy main package
py.install_sources([
    'pyfixi/__init__.py',
    'pyfixi/constraints.py'
  ],
  pure: false, # install next to compiled extension
  subdir: 'pyfixi/'
)

