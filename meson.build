project('fixi', 'cpp', version:'0.1', default_options : ['warning_level=3', 'cpp_std=c++20', 'optimization=3'])

compiler = meson.get_compiler('cpp')

inc = include_directories([['include', 'thirdparty']])

fixi_dep = declare_dependency(include_directories : inc)

executable('fixi_cli', 'src/main.cpp', include_directories : inc)