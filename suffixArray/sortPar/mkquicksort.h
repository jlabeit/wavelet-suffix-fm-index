#include <algorithm>
#include <sdsl/int_vector.hpp>

template<class T, class C, int alphabet_size>
inline C pick_pivot(T* first, T* last, C* text, T depth) {
	// Take care not to go past the end of the text	
	C a = text[*first + depth];
	C b = text[*(last-1) + depth];
	C c = text[*((first + last)/ 2)];
	if (a > b && a < c || a < b && a > c)
	       return a;
	if (b > a && b < c || b < a && b > c)
		return b;
	return c;	
}

template<class T, class C, int alphabet_size>
void mkquicksortB(T* first, T* last, C* text, T depth,  bit_vector& end_flags) {
		if (last - first < THRESHOLD) {
			return;
		}	
		C pivot = pick_pivot(first, last, text, depth);
		T* i = first;
		T* j = last-1;
		while (i < j) {

		}
}
