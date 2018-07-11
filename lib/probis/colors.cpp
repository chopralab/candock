#include "colors.h"

bool Colors::is_visited(int c) {
    //  // ce zelimo upostevati vse aminokisline proteina
    //  return true;

    for (int i = 0; i < size; i++)
        if (color[i] == c) return true;
    return false;
}

void Colors::copy(Colors c) {
    for (int i = 0; i < c.size; i++) color[i] = c.color[i];
    size = c.size;
}
