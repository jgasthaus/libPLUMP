AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}
AM_CXXFLAGS = $(BOOST_CPPFLAGS)

lib_LTLIBRARIES = libplump.la
libplump_la_SOURCES = libplump/context_tree.cc \
					  libplump/hpyp_model.cc \
					  libplump/hpyp_parameters.cc \
					  libplump/hpyp_restaurants.cc \
					  libplump/pyp_sample.cc \
					  libplump/random.cc \
					  libplump/serialization.cc \
					  libplump/stirling.cc \
					  libplump/switching_restaurant.cc

library_includedir=$(includedir)/libplump
library_include_HEADERS = libplump/config.h \
                          libplump/context_tree.h \
                          libplump/hpyp_model.h \
                          libplump/hpyp_parameters.h \
                          libplump/hpyp_restaurants.h \
                          libplump/hpyp_restaurant_interface.h \
                          libplump/libplump.h \
                          libplump/mini_map.h \
                          libplump/node_manager.h \
                          libplump/node_manager_interface.h \
                          libplump/pool.h \
                          libplump/pyp_sample.h \
                          libplump/random.h \
                          libplump/serialization.h \
                          libplump/stirling.h \
                          libplump/subseq.h \
                          libplump/switching_restaurant.h \
                          libplump/utils.h

bin_PROGRAMS = score_file gradient_test
score_file_SOURCES = utils/score_file.cc
score_file_LDADD = libplump.la
score_file_LDFLAGS = $(BOOST_LDFLAGS) \
                     $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) \
					 $(BOOST_SERIALIZATION_LIB) $(BOOST_IOSTREAMS_LIB)
gradient_test_SOURCES = utils/gradient_test.cc
gradient_test_LDADD = libplump.la
gradient_test_LDFLAGS = $(BOOST_LDFLAGS) \
                     $(BOOST_SYSTEM_LIB) $(BOOST_PROGRAM_OPTIONS_LIB) $(BOOST_FILESYSTEM_LIB) \
					 $(BOOST_SERIALIZATION_LIB) $(BOOST_IOSTREAMS_LIB)
