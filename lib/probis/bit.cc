#include "bit.h"

bool Bit::get(int i, int j) {
  return e[i/WORD][j]&1<<(i%WORD);
}

void Bit::set_one(int i, int j) {
  e[i/WORD][j] = e[i/WORD][j]|1<<(i%WORD);
} 

void Bit::set_zero(int i, int j) {
  e[i/WORD][j] = e[i/WORD][j]&(~(1<<(i%WORD)));
} 

void Bit::set_zero_fast(int i, int j) {
  e[i][j] = 0;
} 
