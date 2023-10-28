#---Remove built elements
CLEAN_FORT = rm -f *.mod
CLEAN_LIBS = rm -f *.a
CLEAN_OBJS = rm -rf *.o *.gcda *.gcno __pycache__
CLEAN_EXE = rm -rf *.dSYM *.pyc .pytest_cache
CLEAN_RESULTS = rm -f oft_in oft_in.xml *.results mesh.*.h5 vector_dump.*.h5 scalar_dump.*.h5 *.xmf *.rst *.hist *.err dump.dat
