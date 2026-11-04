# IPC Toolkit (https://github.com/ipc-sim/ipc-toolkit)
# License: MIT

if(TARGET ipc::toolkit)
    return()
endif()

message(STATUS "Third-party: creating target 'ipc::toolkit'")

include(CPM)
CPMAddPackage(
	NAME "ipc-toolkit"
	GIT_REPOSITORY https://github.com/fsichetti/ipc-toolkit.git
	GIT_TAG "3c5caa9a187002b331b5963500bc2bf5c471f03c"
	OPTIONS "IPC_TOOLKIT_WITH_CUDA OFF"
)
