GCC = gcc -Wall -Werror -Wextra -std=c11 -D_GNU_SOURCE
TST_LIBS = -lcheck
OS := $(shell uname -s)


ifeq ($(OS),Linux)
		CC+=-D LINUX_OS
		TST_LIBS+=-lsubunit -lrt -lm -lpthread
else
	ifeq ($(OS),Darwin)
		CC+=-D MAC_OS
		TST_LIBS+=-lm -lpthread
	endif
endif


all: clean s21_math.a test gcov_report

clean:
	rm -rf *.a && rm -rf *.o *.exe
	rm -rf *.info && rm -rf *.gcda && rm -rf *.gcno &&  rm -rf *.gcov *.html *.css *.cfg
	rm -rf report/ && rm -rf *.out

s21_math.a: clean
	$(GCC) -c s21_math.c -o ./s21_math.o
	ar rcs s21_math.a s21_math.o
	ranlib s21_math.a


test: s21_math.a
	$(GCC) tests.c -L. s21_math.a $(TST_LIBS) -o ./tests.exe
	$(GCC) --coverage tests.c s21_math.c s21_math.a -L. s21_math.a $(TST_LIBS) -o ./tests.exe
	./tests.exe

gcov_report: s21_math.a test
	./tests.exe
	gcov -a -f s21_math.gcda
	gcovr --html-details -o report.html

check:
	cppcheck *.h *.c
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extension=c *.c *.h
