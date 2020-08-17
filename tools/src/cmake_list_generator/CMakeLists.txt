# Add the cmake_list_generator
add_custom_target(cmake_gen)
add_custom_command(TARGET cmake_gen
	 	   POST_BUILD
		   COMMAND ${CMAKE_COMMAND} -E copy
		           "${IRL_TOOLS_SRC_DIR}/cmake_list_generator/cmake_list_generator.py"
			   "${IRL_TOOLS_DIR}/bin")
