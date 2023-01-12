py_binary(
    name = "init_data",
    srcs = ["init_data.py"],
    deps = [":init_data_lib"],
)

py_library(
    name = "init_data_lib",
    srcs = ["init_data.py"],
    deps = [
        ":org_hs_eg_db_lookup",
        ":tcga_cnv_parser",
    ],
)

py_library(
    name = "org_hs_eg_db_lookup",
    srcs = ["input_data_adapters/org_hs_eg_db_lookup.py"],
)

py_test(
    name = "org_hs_eg_db_lookup_tests",
    srcs = ["input_data_adapters/org_hs_eg_db_lookup_tests.py"],
    data = [
    ],
    deps = [":org_hs_eg_db_lookup"],
)

py_library(
    name = "tcga_cnv_parser",
    srcs = ["input_data_adapters/tcga_cnv_parser.py"],
)

py_test(
    name = "tcga_cnv_parser_tests",
    srcs = ["input_data_adapters/tcga_cnv_parser_tests.py"],
    deps = [":init_data"],
)
