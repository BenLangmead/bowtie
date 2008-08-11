#ifndef LH3_ALGO_HH
#define LH3_ALGO_HH

// algorithm, C++ version. adapted from algo.h of C version.

typedef struct
{
	size_t left,right;
} ALGO_QSortStack;

template<class ALGO_TYPE>
void algo_sort(size_t n, ALGO_TYPE a[])
{
	size_t s, t, i, j, k;
	ALGO_QSortStack *top, *stack;
	ALGO_TYPE rp, swap_tmp;

	if (n == 0) return;
/*	stack = (ALGO_QSortStack*)malloc(sizeof(ALGO_QSortStack) * (size_t)((sizeof(size_t)*log(n)/ALGO_LOG2)+2)); */
	stack = (ALGO_QSortStack*)malloc(sizeof(ALGO_QSortStack) * ((sizeof(size_t)*64)+2));

	top = stack; s = 0; t = n-1;
	while (1) {
		if (s < t) {
			i = s; j = t; k = (i+j)>>1; rp = a[k];
			swap_tmp = a[k]; a[k] = a[t]; a[t] = swap_tmp;
			do {
				do { ++i; } while (a[i] < rp);
				do { --j; } while (j && rp < a[j]);
				swap_tmp = a[i]; a[i] = a[j]; a[j] = swap_tmp;
			} while (i < j);
			swap_tmp = a[i]; a[i] = a[j]; a[j] = swap_tmp;
			swap_tmp = a[i]; a[i] = a[t]; a[t] = swap_tmp;
			if (i-s > t-i) {
				if (i-s > 9) { top->left = s; top->right = i-1; ++top; }
				if (t-i > 9) s = i+1;
				else s = t;
			} else {
				if (t-i > 9) { top->left = i+1; top->right = t; ++top; }
				if (i-s > 9) t = i-1;
				else t = s;
			}
		} else {
			if (top == stack) {
				free(stack);
				for (i = 1; i < n; ++i)
					for (j = i; j > 0 && a[j] < a[j-1]; --j) {
						swap_tmp = a[j]; a[j] = a[j-1]; a[j-1] = swap_tmp;
					}
				return;
			} else { --top; s = top->left; t = top->right; }
		}
	}
}

template<class TYPE>
inline void algo_swap(TYPE &a, TYPE &b)
{
	TYPE c;
	c = a; a = b; b = c;
}

template<class TYPE>
TYPE algo_ksmall(size_t n, TYPE array[], size_t k)
/* Return the kth smallest value in array array[0..n-1], The input array will be rearranged
 * to have this value in array[k-1], with all smaller elements moved to arr[0..k-2] (in
 * arbitrary order) and all larger elements in arr[k..n] (also in arbitrary order) */
{
	TYPE *arr, a;
	size_t i, ir, j, l, mid;

	if (k == 0) k = 1;
	arr = array - 1;
	l = 1;
	ir = n;
	for (;;) {
		if (ir <= l + 1) { /* Active partition contains 1 or 2 elements */
			if (ir == l + 1 && arr[ir] < arr[l]) /* Case of 2 elements */
				algo_swap(arr[l], arr[ir]);
			return arr[k];
		} else {
			mid = (l + ir) >> 1;
			algo_swap(arr[mid], arr[l+1]);
			if (arr[ir] < arr[l]) algo_swap(arr[l], arr[ir]);
			if (arr[ir] < arr[l+1]) algo_swap(arr[l+1], arr[ir]);
			if (arr[l+1] < arr[l]) algo_swap(arr[l], arr[l+1]);
			i = l + 1; /* initialize pointers for partitioning */
			j = ir;
			a = arr[l+1]; /* partition element */
			for (;;) { /* beginning of innermost loop */
				do ++i; while (arr[i] < a); /* scan up to find element > a */
				do --j; while (a < arr[j]); /* scan down to find element < a */
				if (j < i) break; /* Pointers crossed. Partitioning complete. */
				algo_swap(arr[i], arr[j]);
			}
			arr[l+1] = arr[j];	/* insert partitioning element */
			arr[j] = a;
			if (j >= k) ir = j - 1; /* Keep active the partition that contains the kth element */
			if (j <= k) l = i;
		}
	}
}

template<class TYPE>
void algo_shuffle(int n, TYPE *array)
{
	int i;
	TYPE tmp;
	for (i = n - 1; i > 0; --i) {
		int rand_ind = (int)(drand48() * i);
		tmp = array[i]; array[i] = array[rand_ind]; array[rand_ind] = tmp;
	}
}

// Heap related functions (binary heap)

template<class ALGO_TYPE>
inline void algo_heap_adjust(ALGO_TYPE l[], int i, int n)
{
	int k;
	ALGO_TYPE tmp = l[i];
	for (;;) {
		k = (i << 1) + 1; // the left child
		if (k >= n) { l[i] = tmp; return; }
		if (k < n - 1 && l[k] < l[k+1]) ++k;
		if (tmp < l[k]) { l[i] = l[k]; i = k; }
		else { l[i] = tmp; return; }
	}
}
template<class ALGO_TYPE>
void algo_heap_make(ALGO_TYPE l[], int lsize)
{
	for (int i = (lsize >> 1) - 1; i >= 0; --i)
		algo_heap_adjust(l, i, lsize);
}
template<class ALGO_TYPE>
void algo_heap_sort(ALGO_TYPE l[], int lsize)
{
	ALGO_TYPE swap_tmp;
	for (int i = lsize - 1; i > 0; --i) {
		ALGO_TYPE swap_tmp = l[0];
		l[0] = l[i]; l[i] = swap_tmp;
		algo_heap_adjust(l, 0, i);
	}
}

#endif /* LH3_ALGO_HH */
