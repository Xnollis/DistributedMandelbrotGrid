aux_source_directory(. ALL_SRCS_LIST)
add_subdirectory(mpn)
add_library(mpir ${ALL_SRCS_LIST})
target_link_libraries(mpir
mpir_mpn #在子目录中另外生成的一个子Lib
)