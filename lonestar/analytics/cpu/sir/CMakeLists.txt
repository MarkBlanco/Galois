add_executable(sir sir.cpp)
add_dependencies(apps sir)
target_link_libraries(sir PRIVATE Galois::shmem lonestar)

#add_test_scale(small1 bfs "${BASEINPUT}/reference/structured/rome99.gr")
#add_test_scale(small2 bfs "${BASEINPUT}/scalefree/rmat10.gr")
#add_test_scale(web bfs "${BASEINPUT}/random/r4-2e26.gr")
