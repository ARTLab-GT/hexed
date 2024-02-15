file(GLOB DOC_FILES CONFIGURE_DEPENDS *.dox)
configure_file(config.in doc_config)
file(MAKE_DIRECTORY html)
file(COPY header.png DESTINATION html)
file(COPY header_background.png DESTINATION html)
file(COPY ref_coords.svg DESTINATION html)
file(COPY phys_coords.svg DESTINATION html)
add_custom_command(
    OUTPUT doxygen_output.txt
    COMMAND python3 ${CMAKE_SOURCE_DIR}/script/vis_benchmark.py ${CMAKE_SOURCE_DIR}/benchmark.txt False
    COMMAND doxygen doc_config > doxygen_output.txt
    DEPENDS
        doc_config
        hexed
        ${CMAKE_SOURCE_DIR}/python/hexed_utils.py
        style.css
        ${DOC_FILES}
        ${CMAKE_SOURCE_DIR}/script/vis_benchmark.py
        ${CMAKE_SOURCE_DIR}/benchmark.txt
)
add_custom_target(
    documentation ALL
    DEPENDS doxygen_output.txt
)
