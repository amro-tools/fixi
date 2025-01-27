project('fixi', 'cpp', version:'0.1', default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

compiler = meson.get_compiler('cpp')

inc = include_directories([['include', 'thirdparty']])

fixi_dep = declare_dependency(include_directories : inc)

executable('fixi_cli', 'src/main.cpp', include_directories : inc)

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
    dependencies : [Catch2],
    include_directories : test_inc
  )
  test(t.get(0), exe, workdir : meson.project_source_root())
endforeach