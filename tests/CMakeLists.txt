file(
	GLOB_RECURSE
	source_files_headers
	inc/*
)

file(
	GLOB_RECURSE
	source_files_source
	src/*
)

include_directories(
	inc/
	$(CMAKE_SOURCE_DIR)
)

# En-t�tes externes de la biblioth�que.
source_group(
	"Headers"
	FILES
	$(source_files_headers)
)

source_group(
	"Sources"
	FILES
	$(source_files_source)
)

add_executable(Tests ${source_files_headers} ${source_files_source})
set_target_properties(Tests PROPERTIES LINKER_LANGUAGE CXX)
# target_link_libraries(Tests Formal_Maths)