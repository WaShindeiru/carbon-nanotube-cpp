cmake_minimum_required(VERSION 3.14)
project(md_fec CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

add_executable(md_fec prog_full.cpp)
target_compile_options(md_fec PRIVATE -Ofast)
target_link_libraries(md_fec PRIVATE nlohmann_json::nlohmann_json)
