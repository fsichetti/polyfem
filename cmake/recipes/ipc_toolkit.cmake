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
	GIT_TAG "b7094f5c3b5451261e75217b7db11e7175b9f752"
	OPTIONS "IPC_TOOLKIT_WITH_CUDA OFF"
)
