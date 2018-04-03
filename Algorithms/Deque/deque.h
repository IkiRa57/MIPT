#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <exception>
#include <deque>

template<class T>
class Deque {
private:
    struct block {
        T *chunk;

        block() {
            chunk = new T[512];
        }

        block(const block &another) {
            chunk = new T[512];
            std::copy(another.chunk, another.chunk + 512, chunk);
        }

        ~block() {
            delete[] chunk;
        }
    };

    block **arr;
    size_t deque_size;
    size_t capacity;
    int blockNum;
    int head;
    int headBlock;
    int tail;
    int tailBlock;

    void Resize() {
        bool resize_type = (deque_size == capacity);
        int new_blockNum = resize_type ? 2 * blockNum : blockNum / 2;
        block **Newarr = new block *[new_blockNum];
        int put_place = tail;
        int put_Block = 0;
        Newarr[0] = new block();
        for (unsigned i = 0; i < deque_size; ++i) {
            Newarr[put_Block]->chunk[put_place] = arr[(tailBlock + (i % 512 + tail >= 512) + i / 512) %
                                                      blockNum]->chunk[(i + tail) % 512];
            if (put_place == 511) {
                put_Block = (put_Block + 1) % new_blockNum;
                Newarr[put_Block] = new block();
            }
            put_place = (put_place + 1) % 512;
        }
        delete[] arr;
        for (int i = put_Block + 1; i < new_blockNum; ++i) {
            Newarr[i] = new block();
        }
        arr = Newarr;
        tailBlock = 0;
        headBlock = put_Block;
        capacity *= (new_blockNum) / blockNum;
        blockNum = new_blockNum;
    }

    template<bool is_const_iterator = true>
    class OwnIterator : public std::iterator<std::random_access_iterator_tag, T, long long,
            typename std::conditional<is_const_iterator, const T *, T *>::type,
            typename std::conditional<is_const_iterator, const T &, T &>::type> {
    private:
        const Deque<T> &container;
        int currBlock;
        int offset;

        void get_iterator(int s) {
            int k = ((s % 512 + container.tail) >= 512);
            currBlock = (container.tailBlock + k + s / 512) % container.blockNum;
            offset = (container.tail + s) % 512;
        }

    public:
        int get_position() const {
            int k = !((offset % 512) >= container.tail);

            int block_count = ((currBlock) + container.blockNum - container.tailBlock - k) % container.blockNum;
            int offset_count = (offset + 512 - container.tail) % 512;
            return 512 * block_count + offset_count;
        }

        OwnIterator() {}

        OwnIterator(const Deque<T> &some_deque, int place) : container(some_deque) {
            get_iterator(place);
        }


        OwnIterator(const OwnIterator<false> &another)
                : container(another.container),
                  currBlock(another.currBlock),
                  offset(another.offset) {}

        OwnIterator &operator=(const OwnIterator &another) {
            if (container != another.container) {
                throw (std::logic_error("Invalid assigment of iterators"));
            }
            currBlock = another.currBlock;
            offset = another.offset;
            return *this;
        }

        typedef typename std::conditional<is_const_iterator, const T &, T &>::type reference_type;
        typedef typename std::conditional<is_const_iterator, const T *, T *>::type pointer_type;

        bool operator==(const OwnIterator &another) const {
            return (&container == &another.container &&
                    currBlock == another.currBlock &&
                    offset == another.offset);
        }

        bool operator!=(const OwnIterator &another) const {
            return !(*this == another);
        }

        bool operator<(const OwnIterator &another) const {
            return get_position() < another.get_position();
        }

        bool operator>(const OwnIterator &another) const {
            return get_position() > another.get_position();
        }

        bool operator<=(const OwnIterator &another) const {
            return !(*this > another);
        }

        bool operator>=(const OwnIterator &another) const {
            return !(*this < another);
        }

        reference_type operator*() const {
            return container.arr[currBlock]->chunk[offset];
        }

        pointer_type operator->() const {
            return container.arr[currBlock]->chunk + offset;
        }

        reference_type operator[](int n) const {
            return *(*this + n);
        }

        long long operator-(const OwnIterator &another) const {

            return get_position() - another.get_position();
        }

        OwnIterator &operator++() {
            if (offset == 511) {
                currBlock = (currBlock + 1) % container.blockNum;
            }
            offset = (offset + 1) % 512;
            return *this;
        }

        OwnIterator &operator++(int n) {
            OwnIterator copied(*this);
            ++this;
            return copied;
        }

        OwnIterator &operator--() {
            if (offset == 0) {
                currBlock = (currBlock + container.blockNum - 1) % container.blockNum;
            }
            offset = (offset + 511) % 512;
            return *this;
        }

        OwnIterator &operator--(int n) {
            OwnIterator copied(*this);
            --this;
            return copied;
        }

        OwnIterator operator+=(int n) {
            get_iterator(get_position() + n);
            return *this;
        }

        OwnIterator operator-=(int n) {
            get_iterator(get_position() - n);
            return *this;
        }

        OwnIterator operator+(int n) const {
            OwnIterator temp(*this);
            temp += n;
            return temp;
        }

        OwnIterator operator-(int n) const {
            OwnIterator temp(*this);
            temp -= n;
            return temp;
        }

        friend class OwnIterator<true>;
    };


public:
    typedef OwnIterator<true> const_iterator;
    typedef std::reverse_iterator <const_iterator> const_reverse_iterator;
    typedef OwnIterator<false> iterator;
    typedef std::reverse_iterator <iterator> reverse_iterator;

    Deque() : deque_size(0), capacity(512), blockNum(1), head(0), headBlock(0), tail(0), tailBlock(0) {
        arr = new block *[1];
        arr[0] = new block();
    }

    Deque(const Deque &another) : deque_size(another.deque_size),
                                  capacity(another.capacity),
                                  blockNum(another.blockNum),
                                  head(another.head),
                                  headBlock(another.headBlock),
                                  tail(another.tail),
                                  tailBlock(another.tailBlock) {

        arr = new block *[another.blockNum];
        for (int i = 0; i < another.blockNum; ++i) {
            arr[i] = new block(*another.arr[i]);
        }
    }

    Deque operator=(const Deque &another) {
        if (this != &another) {
            Deque<T> copied(another);
            capacity = copied.capacity;
            deque_size = copied.deque_size;
            blockNum = copied.blockNum;
            head = copied.head;
            headBlock = copied.headBlock;
            tail = copied.tail;
            tailBlock = copied.tailBlock;
            delete[] arr;
            arr = new block *[copied.blockNum];
            for (int i = 0; i < copied.blockNum; ++i) {
                if (copied.arr[i]) {
                    arr[i] = new block(*copied.arr[i]);
                }
            }
        }
        return *this;
    }

    ~Deque() {
        delete[] arr;
    }

    bool operator==(const Deque &another) const {
        if (this == &another) {
            return true;
        } else {
            if (deque_size != another.deque_size) {
                return false;
            }
            return std::equal(begin(), end(), another.begin());
        }
    }

    bool operator!=(const Deque &another) const {
        return !(*this == another);
    }

    void push_back(T n) {
        arr[headBlock]->chunk[head] = n;
        if (head == 511) {
            headBlock = (headBlock + 1) % blockNum;
        }
        head = (head + 1) % 512;
        ++deque_size;
        if (deque_size == capacity) {
            Resize();
        }
    }

    void pop_back() {
        if (deque_size <= capacity / 4 && capacity >= 512 * 4) {
            Resize();
        }
        if (head == 0) {
            headBlock = (headBlock - 1 + blockNum) % blockNum;
        }
        head = (head + 511) % 512;
        --deque_size;
    }

    void push_front(T n) {
        if (deque_size == 0) {
            arr[tailBlock]->chunk[tail] = n;
            if (head == tail) {
                if (head == 511) {
                    headBlock = (headBlock + 1) % blockNum;
                }
                head = (head + 1) % 512;
            }
            ++deque_size;
        } else {
            if (tail == 0) {
                tailBlock = (tailBlock - 1 + blockNum) % blockNum;
            }
            tail = (tail + 511) % 512;
            arr[tailBlock]->chunk[tail] = n;
            ++deque_size;
            if (deque_size == capacity) {
                Resize();
            }
        }
    }

    void pop_front() {
        if (deque_size <= capacity / 4 && capacity >= 512 * 4) {
            Resize();
        }
        if (tail == 511) {
            tailBlock = (tailBlock + 1) % blockNum;
        }
        tail = (tail + 1) % 512;
        --deque_size;
    }

    T &operator[](const int s) {
        int k = (s % 512 + tail >= 512);
        int sBlock = ((tailBlock + k + s / 512) % blockNum);

        return arr[sBlock]->chunk[(s + tail) % 512];
    }

    T operator[](const int s) const {
        int k = (s % 512 + tail >= 512);
        int sBlock = ((tailBlock + k + s / 512) % blockNum);
        return arr[sBlock]->chunk[(s + tail) % 512];
    }

    bool empty() const {
        return !deque_size;
    }

    size_t size() const {
        return deque_size;
    }

    T &back() {
        return operator[](deque_size - 1);
    }

    T back() const {
        T Back = operator[](deque_size - 1);
        return Back;
    }

    T &front() {
        return operator[](0);
    }

    T front() const {
        T Front = operator[](0);
        return Front;
    }

    iterator begin() {
        return iterator(*this, 0);
    }

    const iterator begin() const {
        return iterator(*this, 0);
    }

    iterator end() {
        return iterator(*this, deque_size);
    }

    const iterator end() const {
        return iterator(*this, deque_size);
    }

    reverse_iterator rbegin() {
        return std::reverse_iterator<iterator>(end());
    }

    const reverse_iterator rbegin() const {
        return std::reverse_iterator<iterator>(end());
    }

    reverse_iterator rend() {
        return std::reverse_iterator<iterator>(begin());
    }

    const reverse_iterator rend() const {
        return std::reverse_iterator<iterator>(begin());
    }

    const_iterator cbegin() {
        return const_iterator(*this, 0);
    }

    const const_iterator cbegin() const {
        return const_iterator(*this, 0);
    }

    const_iterator cend() {
        return const_iterator(*this, deque_size);
    }

    const const_iterator cend() const {
        return const_iterator(*this, deque_size);
    }

    const_reverse_iterator crbegin() {
        return std::reverse_iterator<const_iterator>(cend());
    }

    const const_reverse_iterator crbegin() const {
        return std::reverse_iterator<const_iterator>(cend());
    }

    const_reverse_iterator crend() {
        return std::reverse_iterator<const_iterator>(cbegin());
    }

    const const_reverse_iterator crend() const {
        return std::reverse_iterator<const_iterator>(cbegin());
    }


    void printDeque() {

        int i = 0;
        int currBlock = tailBlock;
        int k = tail;
        while (i < deque_size) {
            if (k != 511) {
                std::cout << arr[currBlock]->chunk[k] << ' ';
                ++i;
            } else {
                std::cout << arr[currBlock]->chunk[k] << ' ';
                currBlock = (currBlock + 1) % blockNum;

                ++i;
            }
            k = (k + 1) % 512;
        }
    }

    friend std::ostream &operator<<(std::ostream &out, const Deque &deque) {
        for (unsigned i = 0; i < deque.size(); ++i) {
            out << deque[i] << ' ';
        }
        return out;
    }
};
