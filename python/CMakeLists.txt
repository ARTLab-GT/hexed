# configure python files
get_target_property(CPPYY_INCLUDE hexed INCLUDE_DIRECTORIES)
if (use_occt)
    get_target_property(OCCT_INCLUDE_DIRECTORIES TKBRep INTERFACE_INCLUDE_DIRECTORIES)
    set(CPPYY_INCLUDE ${CPPYY_INCLUDE};${OCCT_INCLUDE_DIRECTORIES})
endif()
configure_file(hexed_py_config.py.in hexed_py_config.py)
configure_file(hexed_python.py hexed.py)
