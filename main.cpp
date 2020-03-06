#include <iostream>
#include <algorithm>

using namespace std;

template <class T>
void printArray(T a[], int n) {
    cout << "[ ";
    for (int i = 0; i < n; i++) {
        cout << a[i] << ((i < n-1)?", ":" ]\n\n");
    }
}

// Functions to generate random populated numeric arrays
void createRandomDoubleArray(double *a, int n) {
    for (int i = 0; i <n; i++) {
        a[i] = drand48();
    }
}

void createRandomFloatArray(float *a, int n) {
    for (int i = 0; i <n; i++) {
        a[i] = static_cast<float> (arc4random()) / static_cast<float> (RAND_MAX);
    }
}

void createRandomIntArray(int *a, int n) {
    for (int i = 0; i <n; i++) {
        a[i] = arc4random();
    }
}

void createRandomPositiveIntArray(int *a, int n, int upper_bound) {
    for (int i = 0; i < n; i++) {
        a[i] = arc4random_uniform(upper_bound);
    }
}

template <class T>
T getMax(T *a, int n) {
    T max = a[0];
    for (int i = 1; i < n; i++) {
        if (a[i] > max) max = a[i];
    }
    return max;
}

template <class T>
T getMin(T *a, int n) {
    T min = a[0];
    for (int i = 1; i < n; i++) {
        if (a[i] < min) min = a[i];
    }
    return min;
}

template <class T>
void exchangeSort(T *a, int n) {
    for (int i = 0; i < n -1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (a[i] > a[j]) {
                swap(a[i], a[j]);
            }
        }
    }
}

template <class T>
void selectionSort(T *a, int n) {
    for (int i = 0; i < n - 1; i++) {
        int jMin = i;
        for (int j = i + 1; j < n; j++) {
            if (a[j] < a[jMin]) {
                jMin = j;
            }
        }
        if (jMin != i) {
            swap(a[i], a[jMin]);
        }
    }
}

template <class T>
void insertionSort(T *a, int n) {
    int i;
    T temp;
    for (int k = 1; k < n; k++) {
        int order = 0;
        temp = a[k];
        i = k - 1;
        while ((i >= 0) && !order) {
            if(a[i] > temp) {
                a[i + 1] = a[i];
            } else {
                order = 1;
            }
        }
        a[i + 1] = temp;
    }
}

template <class T>
void bubbleSort(T *a, int n) {
    int unsorted = 1;
    int lap, j;
    for(lap = 0; lap < n - 1 && unsorted; lap++) {
        unsorted = 0;
        for(j = 0; j < n - lap - 1; j++) {
            if (a[j] > a[j + 1]) {
                swap(a[j], a[j + 1]);
                unsorted = 1;
            }
        }

    }
}

template <class T>
void merge(T *a, int first, int mid, int last) {
    int i, j, k;
    int n1 = mid - first + 1;
    int n2 = last - mid;

    T left[n1], right[n2];

    for(i = 0; i < n1; i++) {
        left[i] = a[first + i];
    }

    for(j = 0; j < n2; j++) {
        right[j] = a[mid + 1 + j];
    }

    i = 0;
    j = 0;
    k = first;
    while (i < n1 && j < n2) {
        if (left[i] <= right[j]) {
            a[k] = left[i];
            i++;
        } else {
            a[k] = right[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        a[k] = left[i];
        i++;
        k++;
    }

    while (j < n2) {
        a[k] = right[j];
        j++;
        k++;
    }
}

template <class T>
void mergeSort(T *a, int first, int last) {
    if (first < last) {
        int mid = (first + last) / 2;
        mergeSort(a, first, mid);
        mergeSort(a, mid +1, last);
        merge(a, first, mid, last);
    }
}

template <class T>
void shellSort(T *a, int n) {
    // We need to track *leap* (integer) and we need indices to the cells we are going to compare.
    int leap;

    // Our leap is the half of the array length (n/2) and while it is greater than 0, we check the
    // array cells.  Every time we process the array, we divide leap by 2 and begin again.
    leap = n / 2;
    while (leap > 0) {
        for (int i = leap; i < n; i ++) {
            int j = i - leap;
            while (j >= 0) {
                int k = j + leap;
                if (a[j] <= a[k]) {
                    j--;
                } else {
                    swap(a[j], a[k]);
                    j -= leap;
                }
            }
            leap /= 2;
        }
    }
}

template <class T>
void quickSort(T *a, int first, int last) {
    // We need a T variable to hold the pivot element and we also need indices to the cells we are going
    // to compare.
    int i, j, middle;
    T pivot;

    // We obtain the index to the middle element by adding first and last, and then dividing by two. We
    // can now get the value of the pivot element (the one at the middle index).
    middle = (first + last) / 2;
    pivot = a[middle];

    // Then we duplicate the values of our first and last indices.
    i = first;
    j = last;

    // Next, while our i index (the first index duplicate) is less than or equals to our j index (the
    // last index duplicate):
    //      We increase our i index while the element on that position is less than the pivot.
    //      We decrease our j index while the element on that position is greater than the pivot.
    //      If index i is less than or equals to index j, we swap the elements at those positions
    //      and then we increment i by 1 and decrement j by 1.
    do {
        while (a[i] < pivot) i++;
        while (a[j] > pivot) j--;

        if (i <= j) {
            swap(a[i], a[j]);
            i++;
            j--;
        }
    } while (i <= j);
    // Once we have rearrenged our arrays elements, if index first is less than index j we recursively
    // call quickSort sending the array, first, and j as last.
    if (first < j) {
        quickSort(a, first, j);
    }

    // Then, if i is less than last, we recursively call quickSort sending the array, i as first, and last.
    if (i < last) {
        quickSort(a, i, last);
    }
}

// This function helps us "build" the Max Heap (a complete binary tree) from the T *a array of size n.
// The parameter node points to the "root" of our tree.
template <class T>
void heapify(T *a, int n, int node) {
    // We need to identify the pointer to the node with the greatest value.  We initially set this
    // pointer to the node parameter.
    int largest = node;
    // We also need pointers to the left and right children of the node.  They are computed as
    // two times node plus one (left) and two times node plus two (right).
    int left = node * 2 + 1;
    int right = node * 2 + 2;

    // We need to check if any of the root's children is greater than the root.
    // If our computed left is less than the size of our array (n) and the element in the left
    // position is larger than the element in the largest position, we set the largest pointer to left.
    if (left < n && a[left] > a[largest]) {
        largest = left;
    }

    // Then, if our computed right is less than the size of our array and the element is the right
    // position is larger than the element in the largest position, we set the largest pointer to right.
    if (right < n && a[right] > a[largest]) {
        largest = right;
    }
    // If the pointer to the largest element is different to the pointer to the root element (node), then
    // we swap their values and we recursively call again to heapify sending the array, its size, and the
    // pointer to the largest element as or new root.
    if (largest != node) {
        swap(a[largest], a[node]);
        heapify(a, n, largest);
    }
}

template <class T>
void heapSort(T *a, int n) {
    // We run the heapify algorithm sending as root node pointers from (n/2)-1 to 0; that is, from the
    // middle of the array to the first element.  This will help us "build" the binary tree for the
    // first time.
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(a, n, i);
    }
    // Once we have the tree, we have the tree, we are sure that we have the largest element in the
    // first cell of our array (index 0).
    // Now, we go in a loop from the last position of our array (n-1) to 0 and we:
    //      Swap the elements in the position [0] with the last unused position in the aray.
    //      Heapify the rest of the array having our root in cell 0.
    for (int i = n - 1; i >= 0; i--) {
        swap(a[0], a[n-1]);
        heapify(a, i, 0);
    }
}

template <class T>
void bucketSort(T *a, int n) {
    // We need a vector of class T, of n elements, that will hold the buckets we need to sort our array.
    vector<T> bucket[n];

    // We store each element of our array in their correspondent bucket (n*a[i]).
    for (int i = 0; i < n; i++) {
        int j = n * a[i];
        bucket[j].push_back(a[i]);
    }

    // Then, we sort each bucket.  For simplicity, we will use the sort function of the algorithm library,
    // sending each buckets begin() and end() properties.
    for (int i = 0; i < n; i++) {
        sort(bucket[i].begin(), bucket[i].end());
    }
    // We retrieve each element from each bucket and store them in order in the array.
    int index = 0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < bucket[i].size(); j++) {
            a[index++] = bucket[i][j];
        }
    }
}

template <class T>
void countingSort(T *a, int n) {
    // We need two arrays, one with size maximum value of array a plus 1, and the other the size of the
    // original array to hold the arranged data.
    T max = getMax(a, n);
    int count[max + 1];
    T output[n];

    // The count array is used to know how many elements there are of any given index (count[a[i]]), so
    // we initialise its cells to 0 and then we process the a array to count a[i] occurrences
    // (i from 0 to n-1).
    for (int i = 0; i <= max; i++) count[i] = 0;
    for (int i = 0; i < n; i++) {
        count[a[i]]++;
    }
    // Then we process the count array from 1 to max, accumulating the index count, as in
    // count[i] += count[i -1]
    for (int i = 1; i < n; i++) {
        count[i] += count[i -1];
    }
    // Once we have the count array, we enter a loop for an index having values from n-1 to 0, modifying
    // the output array such that output[count[a[i]]-1] = a[i], and decreasing count[a[i]] by 1.
    for (int i = 0; i < n; i++) {
        output[count[a[i]] - 1] = a[i];
        count[a[i]]--;
    }
    // Finally, we copy element by element the output array into the a array.
    for (int i = 0; i < n; i++) {
        a[i] = output[i];
    }
}

template <class T>
void countingSortForRadix(T *a, int n, int place) {
    // This special version of counting sort always use a ten-elements array (for digits 0..9)
    // We still need an output array of the same size of the a array (n), and a count array of
    // ten elements as stated above.
    const int max = 10;
    int count[max];
    T output[n];

    // We initialize our count array cells to 0.
    for (int i = 0; i < max; i++) count[i] = 0;
    // Then, for every item in our array a, we check for the digit in its "place" place and
    // count it accordingly (count[(a[i]/place) % 10]++)
    for (int i = 0; i < n; i++) {
        count[(a[i]/place) % 10]++;
    }
    // Next we process the count array from 1 to max, accumulating the index count, as in
    // count[i] += count[i - 1]
    for (int i = 1; i < max; i++) {
        count[i] += count[i - 1];
    }
    // Now we enter a loop for an index having values from n-1 to 0, so we fill the output array
    // accordingly to output[count[(a[i] / place) % 10]-1] = a[i], and then decreasing
    // count[(a[i] / place) % 10] by one
    for (int i = 0; i < n; i++) {
        output[count[(a[i]/place) % 10] - 1] = a[i];
        count[(a[i]/place) % 10]--;
    }
    // Finally, we copy element by element the output array into the a array.
    for (int i = 0; i < n; i++) {
        a[i] = output[i];
    }
}

template <class T>
void radixSort(T *a, int n) {
    // We need to know the max value within the a array, so we can know when to stop
    // processing such array.
    int max = getMax(a, n);

    // We enter a loop to process digit place by digit place, in a LSD radix sort.  Our loop
    // begins at place 1, will keep working while max/place is greater than 0, and place will
    // increase ten-fold each loop.  Within our loop, we will call radix sort sending the
    // a array, the n size and the current place (1, 10, 100, etc.)
    for (int place = 1; max / place > 0; place *= 10) {
        countingSortForRadix(a, n, place);
    }
}

int main() {
    const int size = 20;
    int myArray[size];

    createRandomPositiveIntArray(myArray, size, 100);
    printArray(myArray, size);
    countingSort(myArray, size);
    printArray(myArray, size);
    return 0;
}
