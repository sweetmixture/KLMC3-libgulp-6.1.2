#include <stdio.h>

int main() {
    char filename[] = "/work/ta159/ta159/ta159wkjee/Software/klmc3_compiler_tests/KLMC3-libgulp-6.1.2-gnu/Src/_build_libgulp_mpierr/rm/abc/example.txt";  // Replace with your file's path

    // Attempt to delete the file
    if (remove(filename) == 0) {
        printf("File '%s' deleted successfully.\n", filename);
    } else {
        perror("Error deleting the file");
    }

    return 0;
}

