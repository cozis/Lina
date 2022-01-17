gcc tests/test.c  src/lina.c -o test  -Wall -Wextra -g -Isrc/
gcc tests/test_loader.c src/lina.c -o test_loader -Wall -Wextra -g -Isrc/

gcc tests/my_test.c src/lina.c -o my_test -Wall -Wextra -g -Isrc/
