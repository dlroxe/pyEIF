py_binary(
    name = "init_data",
    srcs = ["init_data.py"],
    deps = [":init_data_lib"],
)

py_library(
    name = "init_data_lib",
    srcs = ["init_data.py"],
    deps = [
        ":hs_data_lookup",
        ":org_hs_eg_db_lookup",
    ],
)

py_test(
    name = "init_data_tests",
    srcs = ["init_data_tests.py"],
    deps = [":init_data"],
)

py_library(
    name = "hs_data_lookup",
    srcs = ["input_data_adapters/hs_data_lookup.py"],
)

py_test(
    name = "hs_data_lookup_tests",
    srcs = ["input_data_adapters/hs_data_lookup_tests.py"],
    data = [
        "input_data_adapters/test_data/hs_data/first_20_Hs.data",
        "input_data_adapters/test_data/hs_data/first_3_short_Hs.data",
        "input_data_adapters/test_data/hs_data/short_Hs.data",
    ],
    deps = [":hs_data_lookup"],
)

py_library(
    name = "org_hs_eg_db_lookup",
    srcs = ["input_data_adapters/org_hs_eg_db_lookup.py"],
)
