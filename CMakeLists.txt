cmake_minimum_required(VERSION 2.8)
project(LCSExtractor)

find_package(VTK)
if (VTK_FOUND)
 	# include (${VTK_USE_FILE})
 	message("VTK found")
else()
 	message(ERROR "VTK not found")
endif()

include_directories(${VTK_INCLUDE_DIRS})

find_package(CUDA)
if (CUDA_FOUND)
	include(FindCUDA)
	set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -gencode arch=compute_35,code=sm_35)
	message("CUDA found")
else()
	message("CUDA not found")
endif()
#include(${CUDA_USE_FILE})

#if (VTK_FOUND)
#	include(${VTK_USE_FILE})
#else (VTK_FOUND)
#	message(FATAL_ERROR, "Cannot build the executable without VTK. Please set the VTK_DIR")
#endif (VTK_FOUND)

#cuda_compile(lcsBigBlockInitializationForPositions lcsBigBlockInitializationForPositions.cu)
#cuda_compile(lcsBigBlockInitializationForVelocities lcsBigBlockInitializationForVelocities.cu)
#cuda_compile(lcsBlockedTracingOfRK4 lcsBlockedTracingOfRK4.cu)
#cuda_compile(lcsCollectActiveParticlesForNewInterval lcsCollectActiveParticlesForNewInterval.cu)
#cuda_compile(lcsCollectActiveParticlesForNewRun lcsCollectActiveParticlesForNewRun.cu)
#cuda_compile(lcsExclusiveScanForInt lcsExclusiveScanForInt.cu)
#cuda_compile(lcsGetGroupsForBlocks lcsGetGroupsForBlocks.cu)
#cuda_compile(lcsGetStartOffsetInParticles lcsGetStartOffsetInParticles.cu)
#cuda_compile(lcsInitialCellLocation lcsInitialCellLocation.cu)
#cuda_compile(lcsRedistributeParticles lcsRedistributeParticles.cu)
#cuda_compile(lcsTetrahedronBlockIntersection lcsTetrahedronBlockIntersection.cu)

cuda_add_executable(extractor
	lcsGeometry.cpp
        utility.cpp
        algebra.cpp
        main.cpp
)


target_link_libraries(extractor ${VTK_LIBRARIES})
