#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <exception>
#include <iomanip>
#include <cassert>
#include <cmath>

#define Finite Field

class Rational;

class BigInteger {
private:
    std::vector<unsigned> number;
    bool sign;
    static const int BASE = 1000000000;
    static const int BASE_DEGREE = 9;
public:
    BigInteger() : sign(true) {}

    BigInteger(int n) {
        sign = (n >= 0);
        n = std::abs(n);
        do {
            number.push_back(n % BASE);
            n = n / BASE;
        } while (n);
    }

    BigInteger(const BigInteger &x) {
        sign = x.sign;
        number = x.number;
    }

    void reverse() {
        unsigned int i, j = number.size() - 1;
        for (i = 0; i < j; ++i) {
            unsigned temp = number[i];
            number[i] = number[j];
            number[j] = temp;
            --j;
        }
    }

    std::string toString() const {
        std::string s;
        if (!sign) {
            s += '-';
        }
        s += std::to_string(number[number.size() - 1]);
        for (int i = number.size() - 2; i >= 0; --i) {
            std::string a = std::to_string(number[i]);
            for (unsigned j = 0; j < 9 - a.size(); ++j) {
                s += '0';
            }
            s += a;
        }
        return s;
    }

    void check() {
        while (!number.empty() && !number[number.size() - 1]) {
            number.pop_back();
        }
        if (number.empty()) {
            number.push_back(0);
            sign = 1;
        }
    }

    explicit operator bool() const {
        return (number.size() == 1 && number[0] == 0) ? 0 : 1;
    }

    BigInteger &operator=(const BigInteger &a) {
        sign = a.sign;
        number = a.number;
        check();
        return *this;
    }

    //Bool functions 
    bool operator==(const BigInteger &another) const {

        if (this == &another) {
            return true;
        }
        if (number.size() != another.number.size() || sign != another.sign) {
            return false;
        } else {
            for (int i = another.number.size() - 1; i >= 0; --i) {
                if (another.number[i] != number[i]) {
                    return false;
                }
            }
            return true;
        }
    }

    bool operator!=(const BigInteger &another) const {
        return !(*this == another);
    }

    bool operator<(const BigInteger &another) const {
        if (sign != another.sign) {
            return sign < another.sign;
        } else {
            if (number.size() > another.number.size()) {
                return !sign;
            } else if (number.size() < another.number.size()) {
                return sign;
            } else {
                for (int i = another.number.size() - 1; i >= 0; --i) {
                    if (another.number[i] > number[i]) {
                        return another.sign;
                    }
                    if (another.number[i] < number[i]) {
                        return !another.sign;
                    }
                }
            }
        }
        return false;
    }

    bool operator>=(const BigInteger &another) const {
        return !(*this < another);
    }

    bool operator>(const BigInteger &another) const {
        return (another < *this);
    }

    bool operator<=(const BigInteger &another) const {
        return !(*this > another);
    }

    //Arithmetic operators

    const BigInteger abs() const {
        return (sign) ? *this : -*this;
    }

    const BigInteger gcd(const BigInteger &a) const {
        BigInteger temp1 = abs(), temp2 = a.abs(), temp;
        while (temp1) {
            temp2 %= temp1;
            temp = temp2;
            temp2 = temp1;
            temp1 = temp;
        }
        return temp2;
    }

    const BigInteger operator-() const {
        BigInteger temp;
        temp.sign = !sign;
        temp.number = number;
        return temp;
    }

    BigInteger &operator+=(const BigInteger &another) {
        if (sign == another.sign) {
            bool carry = 0;
            if (number.size() >= another.number.size()) {
                for (unsigned int i = 0; i < number.size(); ++i) {
                    int k = (i >= another.number.size()) ? 0 : another.number[i];
                    int x = number[i] + k;
                    number[i] = ((x + carry) % BASE);
                    carry = (x + carry) / BASE;
                    if (i > another.number.size() && carry == 0) {
                        break;
                    }
                }
                number.push_back(carry);
            } else {
                for (unsigned int i = 0; i < another.number.size(); ++i) {
                    int k = (i >= number.size()) ? 0 : number[i];
                    int x = another.number[i] + k;
                    if (i < number.size()) {
                        number[i] = ((x + carry) % BASE);
                        carry = (x + carry) / BASE;
                    } else {
                        number.push_back((x + carry) % BASE);
                        carry = (x + carry) / BASE;
                    }
                    if (i > another.number.size() && carry == 0) {
                        break;
                    }
                }
                number.push_back(carry);
            }
        } else if (abs() >= another.abs()) {
            bool carry = 0;
            for (unsigned int i = 0; i < number.size(); ++i) {
                int k = (i >= another.number.size()) ? 0 : another.number[i];
                int x = number[i] - carry - k;
                carry = 0;
                if (x < 0) {
                    x = x + BASE;
                    carry = 1;
                }
                number[i] = (x % BASE);
            }
        } else {
            sign = !sign;
            bool carry = 0;
            for (unsigned int i = 0; i < another.number.size(); ++i) {
                int k = (number.size() <= i) ? 0 : number[i];
                int x = BASE + another.number[i] - carry - k;
                if (i < number.size()) {
                    number[i] = (x % BASE);
                } else {
                    number.push_back(x % BASE);
                }
                carry = !(x / BASE);
            }
        }
        check();
        return *this;
    }

    BigInteger &operator-=(const BigInteger &another) {
        return (*this += (-another));
    }

    BigInteger &operator++() {
        return *this += 1;
    }

    BigInteger &operator--() {
        return *this -= 1;
    }

    BigInteger operator++(int) {
        BigInteger c = *this;
        *this += 1;
        return c;
    }

    BigInteger operator--(int) {
        BigInteger c = *this;
        *this -= 1;
        return c;
    }

    BigInteger &operator*=(const BigInteger &another) {
        BigInteger temp;
        for (unsigned int i = 0; i < another.number.size(); ++i) {
            BigInteger temp2;
            for (unsigned int j = 0; j < i; ++j) {
                temp2.number.push_back(0);
            }
            int carry = 0;
            for (unsigned int l = 0; l < number.size() || carry; ++l) {
                int k = (l >= number.size()) ? 0 : number[l];
                long long x = carry + (long long) another.number[i] * k;
                temp2.number.push_back(x % BASE);
                carry = x / BASE;
            }
            temp += temp2;
        }
        temp.sign = (another.sign == sign);
        temp.check();
        *this = temp;
        return *this;
    }

    BigInteger operator/=(const BigInteger &another) {
        if (another.abs() > abs()) {
            *this = 0;
            return *this;
        }
        BigInteger temp;
        unsigned int k = another.number.size();
        BigInteger temp1;
        BigInteger degrees2[30];
        int powers2[30];
        degrees2[0] = another.abs();
        powers2[0] = 1;
        for (int i = 1; i < 30; ++i) {
            degrees2[i] = degrees2[i - 1];
            degrees2[i] += degrees2[i - 1];
            powers2[i] = powers2[i - 1];
            powers2[i] += powers2[i - 1];
        }
        for (unsigned int i = 0; i < k; ++i) {
            temp1.number.push_back(number[number.size() - k + i]);
        }
        if (temp1 < another.abs()) {
            temp1.reverse();
            temp1.number.push_back(number[number.size() - 1 - k]);
            k++;
            temp1.reverse();
        }
        while (number.size() >= k) {
            int t = 29;
            int i = 0;
            while (t >= 0) {
                if (temp1 >= degrees2[t]) {
                    temp1 -= degrees2[t];
                    i += powers2[t];
                }
                --t;
            }
            temp.number.push_back(i);
            if (number.size() >= k) {
                temp1.check();
                temp1.reverse();
                temp1.number.push_back(number[number.size() - 1 - k]);
                temp1.check();
                temp1.reverse();
            }
            temp1.check();
            ++k;
        }
        temp.sign = (sign == another.sign);
        temp.reverse();
        *this = temp;
        return *this;
    }

    BigInteger &operator%=(const BigInteger &a) {
        BigInteger temp = *this;
        temp /= a;
        temp *= a;
        *this -= temp;
        return *this;
    }

    void addzero(const int n) {
        reverse();
        for (int i = 1; i <= n; ++i) {
            number.push_back(0);
        }
        reverse();
    }

    friend std::istream &operator>>(std::istream &is, BigInteger &n);

    friend class Rational;
};

std::istream &operator>>(std::istream &is, BigInteger &n) {
    std::string s;
    is >> s;
    n.sign = (s[0] != '-');
    n.number.clear();
    int k = 0;
    int d = 0;
    for (int i = s.size() - 1; i >= 0; --i) {
        if (s[i] != '-') {
            d *= 10;
            d += s[i] - '0';
            ++k;
            if (!(k % 9) || !i || s[i - 1] == '-') {
                int f = 0;
                while (k) {
                    f *= 10;
                    f += d % 10;
                    d /= 10;
                    --k;
                }
                n.number.push_back(f);
            }
        }
    }
    n.check();
    return is;
}

std::ostream &operator<<(std::ostream &os, const BigInteger &n) {
    os << n.toString();
    return os;
}

const BigInteger operator+(const BigInteger &b, const BigInteger &a) {
    BigInteger temp = a;
    temp += b;
    return temp;
}

const BigInteger operator-(const BigInteger &b, const BigInteger &a) {
    BigInteger temp = b;
    temp -= a;
    return temp;
}

const BigInteger operator*(const BigInteger &b, const BigInteger &a) {
    BigInteger temp = b;
    temp *= a;
    return (temp);
}

const BigInteger operator/(const BigInteger &b, const BigInteger &another) {
    BigInteger temp = b;
    temp /= another;

    return temp;
}

const BigInteger operator%(const BigInteger &b, const BigInteger &a) {
    BigInteger temp = b;
    temp %= a;
    return temp;
}

class Rational {
private:
    BigInteger up;
    BigInteger down;
public:
    Rational() {}

    Rational(int a, int b) {
        if (a == 0) {
            up = 0;
            down = 1;
        } else if ((a > 0) != (b > 0)) {
            up = -a;
            down = b;
        } else {
            up = a;
            down = b;
        }
        normalize();
    }

    Rational(const BigInteger &a, const BigInteger &b) {
        up = (a.sign == b.sign) ? a.abs() : -a.abs();
        down = b.abs();
        normalize();
    }

    Rational(const BigInteger &a) {
        up = a;
        down = 1;
    }

    Rational(int a) {
        up = a;
        down = 1;
    }

    void normalize() {
        if (!up) {
            down = 1;
        } else {
            if (down < 0) {
                down = -down;
                up = -up;
            }
            BigInteger a = up.gcd(down);
            BigInteger temp = a;
            up /= a;
            down /= a;
        }
    }

    std::string toString() const {
        std::string s;
        s += up.toString();
        if (down - 1) {
            s += '/';
            s += down.toString();
        }
        return s;
    }

    explicit operator bool() const {
        return !(up == 0);
    }

    Rational &operator=(const Rational &a) {
        up = a.up;
        down = a.down;
        return *this;
    }

    bool operator==(const Rational &another) const {
        if (this == &another) {
            return true;
        }
        return (up * another.down == down * another.up);
    }

    bool operator!=(const Rational &another) const {
        return !(*this == another);
    }

    bool operator<(const Rational &another) const {
        return ((up * another.down) - (down * another.up) < 0);
    }

    bool operator>=(const Rational &another) const {
        return !(*this < another);
    }

    bool operator>(const Rational &another) const {
        return (another < *this);
    }

    bool operator<=(const Rational &another) const {
        return !(*this > another);
    }

    //Arithmetic operations
    const Rational operator-() const {
        return Rational(-up, down);
    }

    Rational &operator+=(const Rational &another) {
        up *= another.down;
        up += another.up * down;
        down *= another.down;
        normalize();
        return *this;
    }

    Rational &operator-=(const Rational &another) {
        *this += -another;
        return *this;
    }

    Rational &operator++() {
        *this += 1;
        return *this;
    }

    Rational &operator--() {
        return *this -= 1;
    }

    Rational operator++(int) {
        Rational temp = *this;
        *this -= 1;
        return temp;
    }

    Rational operator--(int) {
        Rational temp = *this;
        *this -= 1;
        return temp;
    }

    Rational &operator*=(const Rational &another) {
        up *= another.up;
        down *= another.down;
        normalize();
        return *this;
    }

    Rational &operator/=(const Rational &another) {
        up *= another.down;
        down *= another.up;
        normalize();
        return *this;
    }

    std::string asDecimal(size_t precision = 0) const {
        BigInteger a, b;
        std::string s;
        if (up < 0) {
            s += "-";
        }
        s += (up.abs() / down).toString();
        if (precision == 0) {
            return s;
        }
        s += '.';
        a = up.abs() % down;
        b = down;
        for (unsigned int k = 0; k < precision; ++k) {
            a *= 10;
            int i = 0;
            while (a > b) {
                a -= b;
                i++;
            }
            s += i + '0';
            a %= b;
        }
        return s;
    }

    double toDouble() const {
        double t = 0;
        BigInteger a, b;
        a = up / down;
        for (int i = a.number.size() - 1; i >= 0; --i) {
            t *= up.BASE;
            t += a.number[i];
        }
        a = up % down;
        a *= up.BASE;
        a /= down;
        double k = 0;
        for (int i = a.number.size() - 1; i >= 0; --i) {
            k *= up.BASE;
            k += a.number[i];
        }
        k /= up.BASE;
        t += k;
        if (up < 0) {
            t = -t;
        }
        return t;
    }

    explicit operator double() {
        return toDouble();
    }

    friend std::istream &operator>>(std::istream &is, Rational &n) {
        is >> n.up;
        n.down = 1;
        n.normalize();
        return is;
    }
};

std::ostream &operator<<(std::ostream &os, const Rational &n) {
    os << n.toString();
    return os;
}

const Rational operator+(const Rational &b, const Rational &a) {
    Rational temp(b);
    temp += a;
    return temp;
}

const Rational operator-(const Rational &b, const Rational &a) {
    return (b + (-a));
}

const Rational operator*(const Rational &b, const Rational &a) {
    Rational temp(b);
    temp *= a;
    return temp;
}

const Rational operator/(const Rational &b, const Rational &a) {
    Rational temp(b);
    temp /= a;
    return temp;
}


enum is_prime {
    undefined, prime, composite
};

template<int N>
class Finite {
private:
    static is_prime NumberType;
    long long value;

    static bool isPrime(const int number) {
        if (number == 1) {
            return false;
        }
        for (int i = 2; i <= ceil(sqrt(number)) && i < number; ++i) {
            if ((number % i) == 0) {
                return false;
            }
        }
        return true;
    }

    static int mod(long long n) {
        int x = (n % N + N) % N;
        return x;
    }

    static int BinPow(int a, int n) {
        if (n == 0) {
            return 1;
        }
        if (n % 2 == 1) {
            long long b = BinPow(a, n - 1);
            return mod(b * a);
        } else {
            long long b = BinPow(a, n / 2);
            return mod(b * b);
        }
    }

    static int Inverse(int a, int n) {
        return BinPow(a, n - 2);
    }

public:

    Finite(int number = 0) : value(mod(number)) {
        if (NumberType == undefined) {
            NumberType = isPrime(N) ? prime : composite;
        }
    }

    Finite(const Finite &another) : value(mod(another.value)) {
        if (NumberType == undefined) {
            NumberType = isPrime(N) ? prime : composite;
        }
    }

    explicit operator bool() const {
        return !value;
    }

    bool operator==(const Finite &another) const {
        return (value == another.value);
    }

    bool operator!=(const Finite &another) const {
        return !(value == another.value);
    }

    Finite &operator=(const Finite &another) {
        value = another.value;
        return *this;
    }

    Finite operator+(const Finite &another) const {
        Finite temp(value + another.value);
        return temp;
    }

    Finite operator-(const Finite &another) const {
        Finite temp(value - another.value);
        return temp;
    }

    Finite operator-() const {
        Finite temp(-value);
        return temp;
    }

    Finite &operator+=(const Finite &another) {
        value = mod(value + another.value);
        return *this;
    }

    Finite &operator++() {
        return *this += 1;
    }

    Finite &operator--() {
        return *this -= 1;
    }

    Finite operator++(int t) {
        Finite temp = *this;
        ++(*this);
        return temp;
    }

    Finite operator--(int t) {
        Finite temp = *this;
        --(*this);
        return temp;
    }

    Finite operator-=(const Finite &another) {
        value = mod(value - another.value);
        return *this;
    }

    Finite operator*(const Finite &another) const {
        Finite temp(mod(value * another.value));
        if (value * another.value % N != temp.value) {
        }
        return temp;
    }

    Finite &operator*=(const Finite &another) {
        value = mod(value * another.value);
        return *this;
    }

    Finite operator/(const Finite &another) const {
        if (NumberType != prime) {
            throw std::logic_error("Division in ring");
        }
        Finite temp(mod(value * BinPow(another.value, N - 2)));
        if ((temp.value * another.value) % N != value) {
        }
        return temp;
    }

    Finite &operator/=(const Finite &another) {
        if (NumberType != prime) {
            throw std::logic_error("Division in ring");
        }
        value = mod(value * Inverse(another.value, N));
        return *this;
    }

    friend std::istream &operator>>(std::istream &in, Finite &number) {
        in >> number.value;
        return in;
    }

    friend std::ostream &operator<<(std::ostream &out, const Finite &number) {
        out << number.value;
        return out;
    }
};

template<int N>
is_prime Finite<N>::NumberType = undefined;

constexpr int nextPow(int n) {
    return !n ? 1 : !(n & (n - 1)) ? n : (1 << (32 - __builtin_clz(n)));
}

template<unsigned int N, unsigned int M, typename Field = Rational>
class Matrix {
private:
    static constexpr int Pow = nextPow(N > M ? N : M);
    Field **array;
public:
    Matrix() {
        array = new Field *[N]();
        for (unsigned int i = 0; i < N; ++i) {
            array[i] = new Field[M]();
        }
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                array[i][j] = (Field) 0;
            }
        }
        if (N == M) {
            for (unsigned i = 0; i < N; ++i) {
                array[i][i] = (Field) 1;
            }
        }
    }

    Matrix(const Matrix &another) {
        array = new Field *[N]();
        for (unsigned int i = 0; i < N; ++i) {
            array[i] = new Field[M]();
        }
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] = another.array[i][j];
            }
        }
    }

    Matrix(std::vector <std::vector<Field>> &vect) {
        array = new Field *[N]();
        for (unsigned int i = 0; i < N; ++i) {
            array[i] = new Field[M]();
        }
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] = (Field) vect[i][j];
            }
        }
    }

    Matrix(std::vector <std::vector<int>> &vect) {
        array = new Field *[N]();
        for (unsigned int i = 0; i < N; ++i) {
            array[i] = new Field[M]();
        }
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] = (Field) vect[i][j];
            }
        }
    }


    Field *&operator[](int i) {
        return array[i];
    }

    const Field *operator[](int i) const {
        return array[i];
    }

    Matrix &operator=(const Matrix<N, M, Field> &another) {
        if (this == &another) {
            return *this;
        }
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] = another.array[i][j];
            }
        }
        return *this;
    }

    template<unsigned int K, unsigned int L>
    bool operator==(const Matrix<K, L, Field> &another) const {
        if (K != N || L != M) {
            return false;
        }
        if (this == &another) {
            return true;
        }
        if (std::min(rank(), another.rank()) == 0 && std::max(rank(), another.rank()) == 4)
            return true;

        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                if (array[i][j] != another.array[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    template<unsigned int K, unsigned int L>
    bool operator!=(const Matrix<K, L, Field> &another) const {
        return !(*this == another);
    }

    Matrix &operator+=(const Matrix<N, M, Field> &another) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] += another.array[i][j];
            }
        }
        return *this;
    }

    Matrix &operator-=(const Matrix<N, M, Field> &another) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] -= another.array[i][j];
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix<N, M, Field> &another) const {
        Matrix temp(*this);
        temp += another;
        return temp;
    }

    Matrix operator-(const Matrix<N, M, Field> &another) const {
        Matrix temp(*this);
        temp -= another;
        return temp;
    }

    Matrix &operator*=(const Field &scalar) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                array[i][j] *= scalar;
            }
        }
        return *this;
    }

    Matrix &operator/=(const Field &scalar) {
        return *this *= (1 / scalar);
    }

    Matrix operator*(const Field &scalar) const {
        Matrix temp(*this);
        temp *= scalar;
        return temp;
    }

    Matrix operator/(const Field &scalar) const {
        Matrix temp(*this);
        temp /= scalar;
        return temp;
    }

    Matrix<M, N, Field> transposed() const {
        Matrix<M, N, Field> trans;
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                trans[j][i] = array[i][j];
            }
        }
        return trans;
    }

    void addString(int i, int j, const Field coof) {
        for (unsigned int k = 0; k < M; ++k) {
            array[i][k] += coof * array[j][k];
        }
    }

    void addColumn(int i, int j, const Field coof) {
        for (unsigned int k = 0; k < N; ++k) {
            array[k][i] += coof * array[k][j];
        }
    }

    void multiplyString(int i, const Field coof) {
        for (unsigned int k = 0; k < M; ++k) {
            array[i][k] *= coof;
        }
    }

    void divideString(int i, const Field coof) {
        for (unsigned int k = 0; k < M; ++k) {
            array[i][k] /= coof;
        }
    }

    Matrix GaussTrans() const {
        Matrix temp(*this);
        unsigned int stair = 0, height = 0;
        while (stair < M && height < N) {

            bool halfemptyColumn = false;
            if (temp.array[height][stair] == 0) {
                halfemptyColumn = true;
                for (unsigned int i = height + 1; i < N; ++i) {

                    if (temp.array[i][stair] != 0) {
                        temp.addString(height, i, 1);
                        halfemptyColumn = false;
                        break;
                    }
                }
            }
            if (halfemptyColumn) {
                ++stair;
                continue;
            } else {
                for (unsigned int i = height + 1; i < N; ++i) {
                    temp.addString(i, height, -temp.array[i][stair] / temp.array[height][stair]);
                }
                ++height;
                ++stair;
            }
        }
        return temp;
    }

    Field det() const {
        static_assert(N == M, "Determinant of non-square matrix");
        Matrix gauss = GaussTrans();
        Field deter = gauss[0][0];
        for (unsigned int i = 1; i < N; ++i) {
            deter *= gauss[i][i];
        }
        return deter;
    }

    Field trace() const {
        static_assert(N == M, "Trace of non-square matrix");
        Field sum = 0;
        for (unsigned int i = 0; i < N; ++i) {
            sum += array[i][i];
        }
        return sum;
    }

    int rank() const {
        Matrix gauss = GaussTrans();
        for (unsigned int i = 0; i < N; ++i) {
            bool emptyString = true;
            for (unsigned int j = 0; j < M; ++j) {
                if (gauss[i][j] != 0) {
                    emptyString = false;
                    break;
                }
            }
            if (emptyString) {
                return i;
            }
        }
        return N;
    }

    void invert() {
        static_assert(N == M, "Can't invert non-square matrix");
        Matrix<N, 2 * N, Field> temp;
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < N; ++j) {
                temp[i][j] = array[i][j];
                temp[i][N + j] = (i == j) ? 1 : 0;
            }
        }
        temp = temp.GaussTrans();
        for (unsigned int i = N - 1; i != 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                temp.addString(j, i, -temp[j][i] / temp[i][i]);
            }
        }
        for (unsigned int i = 0; i < N; ++i) {
            temp.divideString(i, temp[i][i]);
        }
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < N; ++j) {
                array[i][j] = temp[i][j + N];
            }
        }
    }

    Matrix inverted() const {
        Matrix temp = *this;
        temp.invert();
        return temp;
    }

    std::vector <Field> getRow(unsigned int k) const {
        std::vector <Field> a;
        for (unsigned i = 0; i < M; ++i) {
            a.push_back(array[k][i]);
        }
        return a;
    }

    std::vector <Field> getColumn(unsigned int k) const {
        std::vector <Field> a;
        for (unsigned i = 0; i < N; ++i) {
            a.push_back(array[i][k]);
        }
        return a;
    }

    template<unsigned int K>
    Matrix<N, K, Field> operator*(const Matrix<M, K, Field> &other) const {
        Matrix<N, K, Field> solution;
        if (N <= 64 && M <= 64 && K <= 64) {
            Matrix<N, K, Field> answer;
            for (unsigned i = 0; i < N; ++i) {
                for (unsigned j = 0; j < M; ++j) {
                    Field sum = 0;
                    for (unsigned int k = 0; k < M; ++k) {
                        sum += array[i][k] * other[k][j];
                    }
                    answer[i][j] = sum;
                }
            }
            return answer;
        }

        Matrix<(Pow > other.Pow) ? Pow : other.Pow, (Pow > other.Pow) ? Pow : other.Pow, Field> A, B, C;

        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                A[i][j] = array[i][j];
            }
        }

        for (unsigned int i = 0; i < M; ++i) {
            for (unsigned int j = 0; j < K; ++j) {
                B[i][j] = other[i][j];
            }
        }

        Matrix<((Pow > other.Pow) ? Pow : other.Pow) / 2, ((Pow > other.Pow) ? Pow : other.Pow) / 2, Field>
                A11, A12, A21, A22, B11, B12, B21, B22, C11, C12, C21, C22, P1, P2, P3, P4, P5, P6, P7;
        const int k = ((Pow > other.Pow) ? Pow : other.Pow) / 2;
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                A11[i][j] = A[i][j];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                A12[i][j] = A[i][j + k];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                A21[i][j] = A[i + k][j];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                A22[i][j] = A[i + k][j + k];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                B11[i][j] = B[i][j];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                B12[i][j] = B[i][j + k];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                B21[i][j] = B[i + k][j];
            }
        }
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                B22[i][j] = B[i + k][j + k];
            }
        }

        P1 = (A11 + A22) * (B11 + B22);
        P2 = (A21 + A22) * B11;
        P3 = A11 * (B12 - B22);
        P4 = A22 * (B21 - B11);
        P5 = (A11 + A12) * B22;
        P6 = (A21 - A11) * (B11 + B12);
        P7 = (A12 - A22) * (B21 + B22);
        C11 = P1 + P4 - P5 + P7;
        C12 = P3 + P5;
        C21 = P2 + P4;
        C22 = P1 - P2 + P3 + P6;
        for (unsigned int i = 0; i < k; ++i) {
            for (unsigned int j = 0; j < k; ++j) {
                C[i][j] = C11[i][j];
                C[i][j + k] = C12[i][j];
                C[i + k][j] = C21[i][j];
                C[i + k][j + k] = C22[i][j];
            }
        }
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < K; ++j) {
                solution[i][j] = C[i][j];
            }
        }
        return solution;
    }

    template<unsigned int K>
    Matrix<N, M, Field> operator*=(const Matrix<M, K, Field> &other) {
        assert(N == M && M == K);
        return *this = *this * other;
    }

    friend std::istream &operator>>(std::istream &in, Matrix &matr) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                in >> matr[i][j];
            }
        }
        return in;
    }

    friend std::ostream &operator<<(std::ostream &out, const Matrix &matr) {
        for (unsigned int i = 0; i < N; ++i) {
            for (unsigned int j = 0; j < M; ++j) {
                out << std::setprecision(2) << matr.array[i][j] << ' ';
            }
            out << std::endl;
        }
        return out;
    }

    template<unsigned N1, unsigned M1, typename Field1>
    friend
    class Matrix;
};

template<unsigned int N1, typename Field1 = Rational>
using SquareMatrix = Matrix<N1, N1, Field1>;

template<unsigned N1, unsigned M1, typename Field = Rational>
Matrix<N1, M1, Field> operator*(const Field &scalar, const Matrix<N1, M1, Field> &a) {
    return a * scalar;
}
