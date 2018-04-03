#include <vector>
#include <string>
#include <iostream>

class Rational;

class BigInteger {
private:
    std::vector<unsigned short> number;
    bool sign;
    static const int BASE = 10;
public:
    BigInteger() : sign(true) {}

    BigInteger(int n) {
        sign = (n >= 0);
        n = std::abs(n);
        do {
            number.push_back(n % BASE);
            n = n / 10;
        } while (n);
    }

    BigInteger(const BigInteger &x) {
        sign = x.sign;
        number = x.number;
    }

    void reverse() {
        unsigned int i, j = number.size() - 1;
        for (i = 0; i < j; ++i) {
            unsigned short temp = number[i];
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
        for (int i = number.size() - 1; i >= 0; --i) {
            s += number[i] + '0';
        }
        return s;
    }

    void check() {
        while (!number.empty() && number[number.size() - 1] == 0) {
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
        if (number.size() != another.number.size() || sign != another.sign) {
            return false;
        }
        for (int i = another.number.size() - 1; i >= 0; --i) {
            if (another.number[i] != number[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const BigInteger &another) const {
        return !(*this == another);
    }

    bool operator<(const BigInteger &another) const {
        if (sign != another.sign) {
            return sign < another.sign;
        }
        if (number.size() > another.number.size()) {
            return !sign;
        }
        if (number.size() < another.number.size()) {
            return sign;
        }
        for (int i = another.number.size() - 1; i >= 0; --i) {
            if (another.number[i] > number[i]) {
                return another.sign;
            }
            if (another.number[i] < number[i]) {
                return !another.sign;
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
        BigInteger temp1 = abs(), temp2 = a.abs();
        BigInteger temp;
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
                int x = carry + another.number[i] * k;
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

    BigInteger &operator/=(const BigInteger &another) {

        if (another.abs() > abs()) {
            *this = 0;
            return *this;
        } else {
            BigInteger temp;
            unsigned int k = another.number.size();
            BigInteger temp1;
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
                int i = 0;
                while (temp1 >= another.abs()) {
                    temp1 -= another.abs();
                    ++i;
                }
                temp.number.push_back(i);
                if (number.size() >= k) {
                    temp1.check();
                    temp1.reverse();
                    temp1.number.push_back(number[number.size() - 1 - k]);
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

std::ostream &operator<<(std::ostream &os, const BigInteger &n) {
    os << n.toString();
    return os;
}

std::istream &operator>>(std::istream &is, BigInteger &n) {
    std::string s;
    is >> s;
    n.sign = (s[0] == '-') ? false : true;
    n.number.clear();
    for (int i = s.size() - 1; i >= 0; --i) {
        if (s[i] != '-') {
            n.number.push_back(s[i] - '0');
        }
    }
    n.check();
    return is;
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
            BigInteger temp = up.gcd(down);
            up /= temp;
            down /= temp;
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

    Rational &operator=(const Rational &a) {
        up = a.up;
        down = a.down;
        return *this;
    }

    bool operator==(const Rational &another) const {
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

    double to_double() const {
        double t = 0;
        BigInteger a, b;
        a = up / down;
        for (int i = a.number.size() - 1; i >= 0; --i) {
            t *= 10;
            t += a.number[i];
        }
        a = up % down;
        a.addzero(10);
        a /= down;
        double k = 0;
        for (int i = a.number.size() - 1; i >= 0; --i) {
            k *= 10;
            k += a.number[i];
        }
        for (int i = 0; i < 10; ++i) {
            k /= 10;
        }
        t += k;
        if (up < 0) {
            t = -t;
        }
        return t;
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



