PROJECT(fox Fortran)

SET(CMAKE_VERBOSE_MAKEFILE FALSE)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#======================#
# set another compiler #
#======================#
# to set another compiler than the default, 
# do the usual from within the build directory:
# $ FC=g95 
# $ rm -rf * && cmake ../
# $ make 

#=====================#
# verbose make output #
#=====================#
# to see out all commands used during building,
# do as follows:
# $ rm -rf * && cmake ../
# $ make VERBOSE=1

#======================#
# speed up compilation #
#======================#
# to speed up compilation by using multiple cores/processors,
# do the following:
# $ rm -rf * && cmake ../
# $ make -j

#==============================#
# get the machine architecture #
#==============================#

SET(OS ${CMAKE_SYSTEM})
SET(MACHINE_ARCH ${CMAKE_SYSTEM_PROCESSOR})

MESSAGE(STATUS "determining operating system and architecture:")
MESSAGE("   -> your operating system is : ${OS}")
MESSAGE("   -> your architecture is     : ${MACHINE_ARCH}")


#==========================#
# check if m4 is available #
#==========================#

MESSAGE(STATUS "Searching for m4 scripting language")
FIND_PROGRAM(M4_PROGRAM m4)
IF(${M4_PROGRAM} MATCHES "NOTFOUND")
  MESSAGE("   -> WARNING : could not find m4, cannot generate fortran files from m4 sources")
ELSE(${M4_PROGRAM} MATCHES "NOTFOUND")
  MESSAGE("   -> ${M4_PROGRAM}")
ENDIF(${M4_PROGRAM} MATCHES "NOTFOUND")

#============================#
# set some other definitions #
#============================#

# file name extensions
SET(m4_ext .m4)
SET(f90pp_ext .F90)
SET(f90_ext .f90)

###===================================###
###===================================###
### determine some preprocessor flags ###
###===================================###
###===================================###

#===========================================================#
# determine end-of-line character based on operating system #
#===========================================================#

INCLUDE(./cmake/DetermineEOL.cmake)

#============================================#
# determine if and how FLUSH intrinsic works #
#============================================#

INCLUDE(./cmake/CheckFlushIntrinsic.cmake)

#=====================================================#
# check for 'associated in restricted expression' bug #
#=====================================================#

INCLUDE(./cmake/CheckAssociatedBug.cmake)

#============================================#
# determine if and how ABORT intrinsic works #
#============================================#

INCLUDE(./cmake/CheckAbortIntrinsic.cmake)

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FPPFLAGS}")


ADD_SUBDIRECTORY(fsys)
ADD_SUBDIRECTORY(utils)
ADD_SUBDIRECTORY(common)
ADD_SUBDIRECTORY(wxml)
ADD_SUBDIRECTORY(wcml)
ADD_SUBDIRECTORY(wkml)
ADD_SUBDIRECTORY(sax)
ADD_SUBDIRECTORY(dom)
ADD_SUBDIRECTORY(examples)

