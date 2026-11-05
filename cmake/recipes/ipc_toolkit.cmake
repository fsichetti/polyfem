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
	GIT_TAG "a21e439db035616899f2b98c1799ea2356d18893"
	OPTIONS "IPC_TOOLKIT_WITH_CUDA OFF"
)
