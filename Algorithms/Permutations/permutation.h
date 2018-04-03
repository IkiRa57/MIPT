#include <iostream>

template<typename T>
void reverse(T *a, T *b) {
    while (b - a > 0) {
        T temp = *a;
        *a = *b;
        *b = temp;
        ++a;
        --b;
    }
}

template<typename T>
void swap(T &a, T &b) {
    T temp = a;
    a = b;
    b = temp;
}

class Permutation {
private:
    int length;
    int *arr;
public:
    Permutation(int n) {
        length = n;
        arr = new int[length];
        for (int i = 0; i < n; ++i) {
            arr[i] = i;
        }
    }

    Permutation(int n, int *b) {
        length = n;
        arr = new int[length];
        for (int i = 0; i < n; ++i) {
            arr[i] = b[i];
        }
    }

    Permutation(const Permutation &x) {
        length = x.length;
        arr = new int[length];
        for (int i = 0; i < length; ++i) {
            arr[i] = x.arr[i];
        }
    }

    ~Permutation() {
        delete[] arr;
    }

    Permutation &operator=(const Permutation &a) {

        if (this == &a) {
            return *this;
        }
        Permutation temp(a);
        swap(length, temp.length);
        swap(arr, temp.arr);
        return *this;
    }

    const Permutation operator*(const Permutation &a) const {
        Permutation temp(length);
        for (int i = 0; i < length; ++i) {
            temp.arr[i] = arr[a.arr[i]];
        }
        return temp;
    }

    Permutation &operator*=(const Permutation &a) {
        return (*this = *this * a);
    }

    Permutation &operator++() {
        for (int i = length - 2; i >= 0; --i) {
            if (arr[i] < arr[i + 1]) {
                for (int x = length - 1; x >= i + 1; --x) {
                    if (arr[x] > arr[i]) {
                        int temp = arr[i];
                        arr[i] = arr[x];
                        arr[x] = temp;
                        break;
                    }
                }
                int j = i + 1, k = length - 1;
                reverse(arr + j, arr + k);
                break;
            }
        }
        return *this;
    }

    Permutation &operator--() {
        for (int i = length - 2; i >= 0; --i) {
            if (arr[i] > arr[i + 1]) {
                for (int x = length - 1; x >= i + 1; --x) {
                    if (arr[i] > arr[x]) {
                        int temp = arr[i];
                        arr[i] = arr[x];
                        arr[x] = temp;
                        break;
                    }
                }
                int j = i + 1, k = length - 1;
                reverse(arr + j, arr + k);
                break;
            }
        }
        return *this;
    }

    Permutation operator--(int) {
        Permutation temp = *this;
        --(*this);
        return temp;
    }

    Permutation operator++(int) {
        Permutation temp = *this;
        ++(*this);
        return temp;
    }

    const Permutation next() const {
        Permutation temp(length, arr);
        return ++temp;
    }

    const Permutation previous() const {
        Permutation temp(length, arr);
        return --temp;
    }

    bool operator==(const Permutation &m) const {
        if (this->length != m.length) {
            return false;
        }
        for (int i = 0; i < length; ++i) {
            if (arr[i] != m.arr[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Permutation &m) const {
        return !(*this == m);
    }

    bool operator<(const Permutation &m) const {
        for (int i = 0; i < length; ++i) {
            if (arr[i] < m.arr[i]) {
                return true;
            }
            if (length < m.length)
                return true;
            if (length > m.length) {
                return false;
            }
        }
        return false;
    }

    bool operator<=(const Permutation &m) const {
        return (*this < m || *this == m);
    }

    bool operator>=(const Permutation &m) const {
        return !(*this < m);
    }

    bool operator>(const Permutation &m) const {
        return (m < *this);
    }

    int operator[](int i) const {
        return arr[i];
    }

    Permutation inverse() const {
        Permutation temp(length);
        for (int i = 0; i < length; ++i) {
            temp.arr[this->arr[i]] = i;
        }
        return temp;
    }

    template<typename T>
    void operator()(T *b) {
        T *c = new T[length];
        for (int i = 0; i < length; ++i) {
            c[i] = b[i];
        }
        for (int i = 0; i < length; ++i) {
            b[this->arr[i]] = c[i];
        }
        delete[] c;
    }

    friend std::ostream &operator<<(std::ostream &os, const Permutation &p);
};

std::ostream &operator<<(std::ostream &os, const Permutation &p) {
    for (int i = 0; i < p.length; ++i) {
        os << p.arr[i] << ' ';
    }
    return os;
}
