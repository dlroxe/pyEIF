py_binary(
    name = "init_data",
    srcs = ["init_data.py"],
    deps = [":init_data_lib"],
)

py_library(
    name = "init_data_lib",
    srcs = ["init_data.py"],
    deps = [
        ":entrez_lookup",
        ":org_hs_eg_db_lookup",
    ],
)

py_test(
    name = "init_data_tests",
    srcs = ["init_data_tests.py"],
    deps = [":init_data"],
)

py_library(
    name = "entrez_lookup",
    srcs = ["entrez_lookup.py"],
)

py_test(
    name = "entrez_lookup_tests",
    srcs = ["entrez_lookup_tests.py"],
    data = [
        "test_data/entrez_data/first_20_Hs.data",
        "test_data/entrez_data/first_3_short_Hs.data",
        "test_data/entrez_data/short_Hs.data",
    ],
    deps = [":entrez_lookup"],
)

py_library(
    name = "org_hs_eg_db_lookup",
    srcs = ["org_hs_eg_db_lookup.py"],
)
