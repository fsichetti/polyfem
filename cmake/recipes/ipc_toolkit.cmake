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
	GIT_TAG "af46b2ceeb5b4d0b73f370ca3ee831c3ce9fb8a7"
	OPTIONS "IPC_TOOLKIT_WITH_CUDA OFF"
)
