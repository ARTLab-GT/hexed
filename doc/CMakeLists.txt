file(GLOB DOC_FILES CONFIGURE_DEPENDS *.dox)
configure_file(config.in doc_config)
file(COPY header.png DESTINATION html)
file(COPY header_background.png DESTINATION html)
add_custom_command(
    OUTPUT doxygen_output.txt
    COMMAND doxygen doc_config > doxygen_output.txt
    DEPENDS
        doc_config
        hexed
        ${CMAKE_SOURCE_DIR}/python/hexed_utils.py
        header.png
        style.css
        ${DOC_FILES}
)
add_custom_target(
    documentation ALL
    DEPENDS doxygen_output.txt
)
