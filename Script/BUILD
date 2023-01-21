py_binary(
    name = "init_data",
    srcs = [
        "__init__.py",
        "init_data.py",
    ],
    deps = [":input_data_adapters"],
)

py_library(
    name = "init_data_lib",
    srcs = [
        "__init__.py",
        "init_data.py",
    ],
    deps = [":input_data_adapters"],
)

py_library(
    name = "input_data_adapters",
    srcs = [
        "input_data_adapters/__init__.py",
        "input_data_adapters/org_hs_eg_db_lookup.py",
        "input_data_adapters/tcga_cnv_parser.py",
    ],
)

py_test(
    name = "org_hs_eg_db_lookup_tests",
    size = "small",
    srcs = ["input_data_adapters/org_hs_eg_db_lookup_tests.py"],
    data = [
    ],
    deps = [":input_data_adapters"],
)

py_test(
    name = "tcga_cnv_parser_tests",
    size = "small",
    srcs = ["input_data_adapters/tcga_cnv_parser_tests.py"],
    data = [":test_data"],
    deps = [":input_data_adapters"],
)

filegroup(
    name = "test_data",
    srcs = [
        "input_data_adapters/test_data/phenotype-01.csv",
        "input_data_adapters/test_data/threshold-01.csv",
    ],
)
