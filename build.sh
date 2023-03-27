gcc tests/test.c src/lina.c src/qr.c -o test -Wall -Wextra -g -Isrc/ -lm
gcc tests/test_loader.c src/lina.c src/qr.c -o test_loader -Wall -Wextra -g -Isrc/ -lm