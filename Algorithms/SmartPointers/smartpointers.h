#include <bits/stdc++.h>

using namespace std;

template<class T>
class UniquePtr {
private:
    T *ptr;
public:
    typedef T value_type;
    typedef T *pointer;
    typedef T &reference;

    UniquePtr() : ptr(nullptr) {}

    UniquePtr(pointer p) {
        ptr = p;
    }

    UniquePtr(UniquePtr &&other) : ptr(other.ptr) {
        other.ptr = nullptr;
    }

    UniquePtr(const UniquePtr &other) = delete;

    UniquePtr &operator=(UniquePtr &&other) {
        delete ptr;
        ptr = other.ptr;
        other.ptr = nullptr;
        return *this;
    }

    UniquePtr &operator=(const UniquePtr &other) = delete;

    ~UniquePtr() {
        delete ptr;
    }

    pointer get() const {
        return ptr;
    }

    reference operator*() const {
        return *ptr;
    }

    pointer operator->() const {
        return ptr;
    }

    void swap(UniquePtr &other) {
        std::swap(ptr, other.ptr);
    }

    pointer release() {
        pointer res(ptr);
        ptr = nullptr;
        return res;
    }

    void reset(pointer nptr = pointer()) {
        if (ptr != nptr) {
            pointer res = ptr;
            ptr = nptr;
            delete res;
        }
    }
};

struct SharedCount {
    size_t use_count;
    size_t weak_count;

    SharedCount() : use_count(0), weak_count(0) {}
};

template<typename T>
class WeakPtr;

template<typename T>
class SharedPtr {
private:
    T *ptr;
    SharedCount *s_count;

    void clean() {
        if (!use_count()) {
            delete ptr;
        }
        if (s_count && !use_count() && !s_count->weak_count) {
            delete s_count;
        }
        ptr = nullptr;
        s_count = nullptr;
    }

    void increment_count() {
        if (!s_count) {
            return;
        }
        ++s_count->use_count;
    }

    void decrement_count() {
        if (!s_count) {
            return;
        }
        if (s_count->use_count) {
            --s_count->use_count;
        }
    }

public:
    typedef T value_type;
    typedef T *pointer;
    typedef T &reference;

    SharedPtr() : ptr(nullptr), s_count(nullptr) {

    }

    SharedPtr(pointer p) : ptr(p) {
        s_count = new SharedCount;
        increment_count();
    }

    SharedPtr(const SharedPtr &other) : ptr(other.ptr), s_count(other.s_count) {
        increment_count();

    }

    SharedPtr(const WeakPtr<T> &other) : ptr(other.ptr), s_count(other.s_count) {
        if (!s_count) {
            s_count = new SharedCount();
        }
        increment_count();
    }

    SharedPtr(SharedPtr &&other) : ptr(other.ptr), s_count(other.s_count) {
        other.ptr = nullptr;
        other.s_count = nullptr;
    }

    SharedPtr &operator=(const SharedPtr &other) {
        decrement_count();
        clean();
        ptr = other.ptr;
        s_count = other.s_count;
        increment_count();
        return *this;

    }

    SharedPtr &operator=(SharedPtr &&other) {
        decrement_count();
        clean();
        ptr = other.ptr;
        s_count = other.s_count;
        other.ptr = nullptr;
        other.s_count = nullptr;
        return *this;
    }

    ~SharedPtr() {
        decrement_count();
        clean();
    }

    int use_count() const {
        return s_count ? s_count->use_count : 0;
    }

    pointer get() const {
        return ptr;
    }

    reference operator*() const {
        return *ptr;
    }

    pointer operator->() const {
        return ptr;
    }

    void swap(SharedPtr &other) {
        std::swap(s_count, other.s_count);
        std::swap(ptr, other.ptr);
    }

    void reset(pointer nptr = pointer()) {
        if (ptr != nptr) {
            decrement_count();
            clean();
            ptr = nptr;
            s_count = new SharedCount();
            increment_count();
        }
    }

    template<typename>
    friend
    class WeakPtr;
};

template<typename T>
class WeakPtr {
private:
    T *ptr;
    SharedCount *s_count;

    void clean() {
        if (s_count && !use_count() && !s_count->weak_count) {
            delete s_count;
        }
        ptr = nullptr;
        s_count = nullptr;
    }

    void increment_count() {
        if (s_count) {
            ++s_count->weak_count;
        }
    }

    void decrement_count() {
        if (s_count && s_count->weak_count) {
            --s_count->weak_count;
        }
    }

public:
    WeakPtr() : ptr(nullptr), s_count(nullptr) {
    }

    WeakPtr(const WeakPtr &other) : ptr(other.ptr), s_count(other.s_count) {
        increment_count();
    }

    WeakPtr(WeakPtr &&other) : ptr(other.ptr), s_count(other.s_count) {
        other.ptr = nullptr;
        other.s_count = nullptr;
    }

    WeakPtr(const SharedPtr<T> &other) : ptr(other.ptr), s_count(other.s_count) {
        increment_count();
    }

    WeakPtr &operator=(const WeakPtr &other) {
        decrement_count();
        ptr = other.ptr;
        s_count = other.s_count;
        increment_count();
        return *this;
    }

    WeakPtr &operator=(WeakPtr &&other) {
        decrement_count();
        clean();
        ptr = other.ptr;
        s_count = other.s_count;
        other.ptr = nullptr;
        other.s_count = nullptr;
        return *this;
    }

    WeakPtr &operator=(const SharedPtr<T> &other) {
        decrement_count();
        clean();
        ptr = other.ptr;
        s_count = other.s_count;
        increment_count();
        return *this;
    }

    ~WeakPtr() {
        decrement_count();
        clean();
    }

    int use_count() const {
        return s_count ? s_count->use_count : 0;
    }

    bool expired() const {
        return !use_count();
    }

    SharedPtr<T> lock() const {
        if (expired()) {
            return SharedPtr<T>();
        }
        return SharedPtr<T>(*this);
    }

    void swap(WeakPtr &other) {
        std::swap(s_count, other.s_count);
        std::swap(ptr, other.ptr);
    }

    void reset() {
        decrement_count();
        clean();
    }

    template<typename>
    friend
    class SharedPtr;
};

