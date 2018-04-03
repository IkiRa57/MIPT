
#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <typeinfo>
#include <cassert>

class Line;

const double EPS__ = 1e-6;

struct Point {
    double x;
    double y;

    Point() {}

    Point(double a, double b) {
        x = a;
        y = b;
    }

    Point &operator=(const Point &another) {
        x = another.x;
        y = another.y;
        return *this;
    }

    const Point operator+(const Point &another) const {
        Point temp;
        temp.x = (fabs(x + another.x) < EPS__) ? 0 : x + another.x;
        temp.y = (fabs(y + another.y) < EPS__) ? 0 : y + another.y;
        return temp;
    }

    const Point operator-() const {
        Point temp;
        temp.x = -(x);
        temp.y = -(y);
        return temp;
    }

    const Point operator-(const Point &another) const {
        return (*this + (-another));
    }

    const Point operator*(const double &a) const {
        Point temp;
        temp.x = a * x;
        temp.y = a * y;
        return temp;
    }

    const Point operator/(const double &a) const {
        Point temp;
        temp.x = x / a;
        temp.y = y / a;
        return temp;
    }

    bool operator==(const Point &another) const {
        return ((fabs(x - another.x) < EPS__) && (fabs(y - another.y) < EPS__));
    }

    bool operator!=(const Point &another) const {
        return !(*this == another);
    }

    const Point rotate(const Point center, double angle) const {
        Point temp1, temp;
        temp1 = *this - center;
        temp.x = center.x + temp1.x * (cos(angle)) - temp1.y * (sin(angle));
        if (fabs(temp.x) < EPS__) temp.x = 0;
        temp.y = center.y + temp1.x * (sin(angle)) + temp1.y * (cos(angle));
        if (fabs(temp.y) < EPS__) temp.y = 0;
        return temp;
    }

    const Point scale(const Point center, double coefficient) const {
        Point temp;
        temp.x = center.x + (x - center.x) * coefficient;
        if (fabs(temp.x) < EPS__) temp.x = 0;
        temp.y = center.y + (y - center.y) * coefficient;
        if (fabs(temp.y) < EPS__) temp.y = 0;
        return temp;
    }

    const Point reflex(const Point center) {
        Point temp;
        temp.x = 2 * center.x - x;
        temp.y = 2 * center.y - y;
        return temp;
    }

    const double distance(const Point &a) const {
        return sqrt((x - a.x) * (x - a.x) + (y - a.y) * (y - a.y));
    }

    const double length() const {
        return sqrt(double(x * x + y * y));
    }

    const double length2() const {
        return (x * x + y * y);
    }

    const Point reflex(Line axis);
};

const Point operator*(const double &a, const Point &b) {
    Point temp;
    temp.x = a * b.x;
    temp.y = a * b.y;
    return temp;
}

std::ostream &operator<<(std::ostream &os, const Point &n) {
    double a = n.x;
    double b = n.y;
    os << '(' << a << ", " << b << ')';
    return os;
}

class Line {
private:
    double A;
    double B;
    double C;
public:
    Line() {}

    Line(Point a, Point b) {
        if (a.x * b.y == a.y * b.x) {
            C = 0;
            A = a.y - b.y;
            B = -a.x + b.x;
            if (A < 0) {
                A = -A;
                B = -B;
                C = -C;
            }
            double r = sqrt(A * A + B * B);
            A /= r;
            B /= r;
        } else {
            C = -1;
            A = (b.y - a.y) / (a.x * b.y - b.x * a.y);
            B = (a.x - b.x) / (a.x * b.y - b.x * a.y);
            if (A < 0) {
                A = -A;
                B = -B;
                C = -C;
            }
            double r = sqrt(A * A + B * B);
            A /= r;
            B /= r;
            C /= r;

        }
    }

    Line(Point a, double coof) {
        Point temp;
        temp.x = a.x + 1;
        temp.y = a.y + coof;
        *this = Line(a, temp);
    }

    Line(double k, double b) {
        Point temp(0, b);
        *this = Line(temp, k);
    }

    bool isHigher(const Point &point) const {
        return ((A != 0) && -(B * point.y + C) / A >= point.x);
    }

    friend std::ostream &operator<<(std::ostream &os, const Line &n);

    bool operator==(const Line &another) const {
        return (A == another.A && B == another.B && C == another.C);
    }

    bool operator!=(const Line &another) const {
        return !(*this == another);
    }

    Point normal() const {
        Point temp(A / sqrt(A * A + B * B), B / sqrt(A * A + B * B));
        return temp;
    }

    friend const Point Point::reflex(Line axis);
};

const Point Point::reflex(Line axis) {
    Point temp;
    double xq = (axis.A * x + axis.B * y + axis.C) / (axis.A * axis.A + axis.B * axis.B);
    temp.x = x - 2 * xq * axis.A;
    temp.y = y - 2 * xq * axis.B;
    return temp;
}

std::ostream &operator<<(std::ostream &os, const Line &n) {
    double a = n.A;
    double b = n.B;
    double c = n.C;
    os << a << "*x + " << b << "*y + " << c << " = 0";
    return os;
}

class Shape {
public:
    virtual ~Shape() {}

    virtual bool operator==(const Shape &another) const = 0;

    virtual bool operator!=(const Shape &another) const = 0;

    virtual double perimeter() const = 0;

    virtual double area() const = 0;

    virtual bool isCongruentTo(const Shape &another) const = 0;

    virtual bool isSimilarTo(const Shape &another) const = 0;

    virtual bool containsPoint(const Point &p) const = 0;

    virtual void rotate(Point center, double angle) = 0;

    virtual void reflex(Point center) = 0;

    virtual void reflex(Line axis) = 0;

    virtual void scale(Point center, double coefficient) = 0;
};

class Vector {
private:
    double x;
    double y;
public:
    Vector() {}

    Vector(const Point &first, const Point &second) {
        x = second.x - first.x;
        y = second.y - first.y;
    }

    Vector(const Point &A) {
        x = A.x;
        y = A.y;
    }

    const double det(const Vector &A) const {
        return (x * A.y - y * A.x);
    }

    const double operator*(const Vector &A) const {
        return (x * A.x + y * A.y);
    }

    const double length() const {
        return sqrt(x * x + y * y);
    }

    const double sqrlength() const {
        return (x * x + y * y);
    }

    const double angle(const Vector &A) const {
        return (A * (*this) / (A.length() * length()));
    }
};

class Ellipse : public Shape {
protected:
    Point F1;
    Point F2;
    double sum;
    double a;
    double b;
public:
    Ellipse() {}

    Ellipse(Point &f1, Point &f2, double Sum) {
        F1 = f1;
        F2 = f2;
        sum = Sum;
        a = Sum / 2;
        Vector xq(f1, f2);
        b = sqrt(a * a - xq.sqrlength() / 4);
    }

    double area() const {
        return M_PI * a * b;
    }

    double perimeter() const {
        return M_PI * (3 * (a + b) - sqrt((3 * a + b) * (3 * b + a)));
    }

    bool operator==(const Shape &another) const {
        const Ellipse *other = dynamic_cast<const Ellipse *>(&another);
        if (!other) {
            return false;
        }
        return (focuses() == other->focuses() && sum == other->sum);
    }

    bool operator!=(const Shape &another) const {
        return !(*this == another);
    }

    bool isCongruentTo(const Shape &another) const {
        const Ellipse *other = dynamic_cast<const Ellipse *>(&another);
        if (!other) {
            return false;
        }
        return (a == other->a && b == other->b);
    }

    bool isSimilarTo(const Shape &another) const {
        const Ellipse *other = dynamic_cast<const Ellipse *>(&another);

        if (!other) {
            return false;
        }
        return (eccentricity() == other->eccentricity());
    }

    bool containsPoint(const Point &point) const {
        return (F1.distance(point) + F2.distance(point) < sum);
    }

    std::pair <Point, Point> focuses() const {
        return (std::make_pair(F1, F2));
    }

    const Point center() const {
        return (F1 + F2) / 2;
    }

    double eccentricity() const {
        return (sqrt(1 - b * b / (a * a)));
    }

    std::pair <Line, Line> directrixes() {
        Point A = F1.scale(center(), a / eccentricity());
        Line D1(A, (F1.y - F2.y) / (F2.x - F1.x));
        Point B = F2.scale(center(), a / eccentricity());
        Line D2(B, (F1.y - F2.y) / (F2.x - F1.x));
        return (std::make_pair(D1, D2));
    }

    void rotate(Point center, double angle) {
        F1 = F1.rotate(center, angle);
        F2 = F2.rotate(center, angle);
    }

    void reflex(Point center) {
        F1 = F1.reflex(center);
        F2 = F2.reflex(center);
    }

    void reflex(Line axis) {
        F1 = F1.reflex(axis);
        F2 = F2.reflex(axis);
    }

    void scale(Point center, double coefficient) {
        F1 = F1.scale(center, coefficient);
        F2 = F2.scale(center, coefficient);
        sum *= fabs(coefficient);
        a *= fabs(coefficient);
        b *= fabs(coefficient);
    }

    friend std::ostream &operator<<(std::ostream &os, const Ellipse &n);
};

std::ostream &operator<<(std::ostream &os, const Ellipse &n) {
    os << n.focuses().first << ' ' << n.focuses().second << ' ' << n.sum << std::endl;
    return os;
}

class Circle : public Ellipse {
public:
    Circle() {}

    Circle(Point &center, double radius) : Ellipse(center, center, 2 * radius) {}

    double radius() const {
        return (sum / 2.);
    }

    double perimeter() const {
        return 2 * M_PI * radius();
    }

    double area() const {
        return M_PI * radius() * radius();
    }

};


std::vector<int> z_function(std::vector <Point> array) {
    int n = array.size();
    std::vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        if (i <= r) {
            z[i] = std::min(r - i + 1, z[i - l]);
        }
        while (i + z[i] < n && array[z[i]] == array[i + z[i]]) {
            ++z[i];
        }
        if (i + z[i] - 1 > r) {
            l = i;
            r = i + z[i] - 1;
        }
    }
    return z;
}


class Polygon : public Shape {
protected:
    std::vector <Point> vert;

    void init(const Point &p) {
        vert.push_back(p);
    }

    template<class... Args>
    void init(const Point &p, Args... args) {
        vert.push_back(p);
        init(args...);
    }

public:
    Polygon() {}

    Polygon(const std::vector <Point> &vert2) : vert(vert2.begin(), vert2.end()) {}

    template<class... Args>
    Polygon(Args... args) {
        init(args...);
    }

    std::vector <Point> getVertices() {
        return vert;
    }

    const int verticesCount() const {
        return vert.size();
    }

    double area() const {
        double ar = 0;
        for (size_t i = 0; i < vert.size(); ++i) {
            Vector S(vert[i]);
            Vector S2(vert[(i + 1) % vert.size()]);
            ar += S.det(S2);
        }
        return fabs(ar / 2);
    }

    double perimeter() const {
        double P = 0;
        for (size_t i = 0; i < vert.size(); ++i) {
            Vector S(vert[i], vert[(i + 1) % vert.size()]);
            P += S.length();
        }
        return P;
    }

    bool isConvex() const {
        int k = 1;
        int l = 0;

        for (size_t i = 0; i < vert.size(); ++i) {
            Vector A(vert[i], vert[(i + 1) % vert.size()]);
            Vector B(vert[(i + 1) % vert.size()], vert[(i + 2) % vert.size()]);
            if (A.det(B) > 0) {
                k = 0;
            }
            if (A.det(B) < 0) {
                l = 1;
            }
        }
        return !(k < l);
    }

    bool operator==(const Shape &another) const {
        const Polygon *another2 = dynamic_cast<const Polygon *>(&another);
        if (!another2) {
            return false;
        }
        std::vector <Point> a = another2->vert;
        std::vector <Point> b;
        for (size_t i = 0; i < another2->vert.size(); ++i) {
            b.push_back(another2->vert[another2->vert.size() - 1 - i]);
        }
        if (another2->vert.size() != vert.size()) {
            return false;
        } else {
            for (size_t i = 0; i < vert.size(); ++i) {
                a.push_back(vert[i]);
            }
            for (size_t i = 0; i < vert.size(); ++i) {
                a.push_back(vert[i]);
            }
            int n = a.size();
            std::vector<int> z(n);
            z = z_function(a);
            for (size_t i = 0; i < 2 * vert.size(); ++i) {
                if (z[vert.size() + i] >= (int) vert.size()) {
                    return true;
                }
            }
            for (size_t i = 0; i < vert.size(); ++i) {
                b.push_back(vert[i]);
            }
            for (size_t i = 0; i < vert.size(); ++i) {
                b.push_back(vert[i]);
            }
            std::vector<int> y;
            y = z_function(b);
            for (size_t i = 0; i < 2 * vert.size(); ++i) {
                if (y[vert.size() + i] >= (int) vert.size()) {
                    return true;
                }
            }
        }
        return false;
    }

    bool operator!=(const Shape &another) const {
        return !(*this == another);
    }

    bool isCongruentTo(const Shape &another) const {
        const Polygon *another2 = dynamic_cast<const Polygon *>(&another);
        if (!another2) {
            return false;
        }
        if (vert.size() != another2->vert.size()) {
            return false;
        }
        int n = vert.size();

        std::vector <Point> another3;
        for (int i = 0; i < n; ++i) {
            another3.push_back(another2->vert[n - 1 - i]);
        }
        for (int order = 0; order < 2; ++order) {
            std::vector <Point> temp = (order) ? another3 : another2->vert;
            for (int i = 0; i < n; ++i) {
                std::vector<double> angles;
                std::vector<double> another_angles;
                std::vector<double> ratio;
                bool k = true;
                for (size_t j = 1; j < temp.size(); ++j) {
                    Vector a(vert[j - 1], vert[(j)]);
                    Vector b(temp[(i + j - 1) % n], temp[(i + j) % n]);
                    ratio.push_back(a.sqrlength() / b.sqrlength());
                    Polygon A(vert[(j - 1)], vert[j], vert[(j + 1) % n]);
                    Polygon B(temp[(i + j - 1) % n], temp[(i + j) % n], temp[(i + j + 1) % n]);
                    angles.push_back(fabs(A.area() / a.sqrlength()) < EPS__ ? 0 : (A.area() / a.sqrlength()));
                    another_angles.push_back(fabs(B.area() / b.sqrlength()) < EPS__ ? 0 : (B.area() / b.sqrlength()));
                }
                for (size_t j = 0; j < angles.size(); ++j) {
                    if (fabs(ratio[j] - 1) > EPS__ || (angles[j] - another_angles[j]) > EPS__) {
                        k = false;
                    }
                }
                if (k) {
                    return true;
                }
            }
        }
        return false;
    }

    bool isSimilarTo(const Shape &another) const {
        std::vector <Point> other_vert(vert.size());
        double coefficient = another.perimeter() / perimeter();
        for (size_t i = 0; i < vert.size(); ++i) {
            other_vert[i] = vert[i] * coefficient;
        }
        Polygon other(other_vert);
        return other.isCongruentTo(another);
    }

    bool containsPoint(const Point &point) const {
        int k = 0;
        for (size_t i = 0; i < vert.size(); ++i) {
            k += (std::max(vert[i].x, vert[(i + 1) % vert.size()].x) <= point.x
                  && std::max(vert[i].y, vert[(i + 1) % vert.size()].y) >= point.y
                  && std::min(vert[i].y, vert[(i + 1) % vert.size()].y) < point.y);
        }
        return (k % 2);
    }

    void reflex(Point center) {
        std::vector <Point> vect;
        for (size_t i = 0; i < vert.size(); ++i) {
            vect.push_back(2 * center - vert[i]);
        }
        vert = vect;
    }

    void rotate(Point center, double angle) {
        std::vector <Point> vect;
        for (size_t i = 0; i < vert.size(); ++i) {
            vect.push_back(vert[i].rotate(center, angle));
        }
        vert = vect;
    }

    void reflex(Line axis) {
        std::vector <Point> vect;
        for (size_t i = 0; i < vert.size(); ++i) {
            vect.push_back(vert[i].reflex(axis));
        }
        vert = vect;
    }

    void scale(Point center, double coefficient) {
        std::vector <Point> vect;
        for (size_t i = 0; i < vert.size(); ++i) {
            vect.push_back(vert[i].scale(center, coefficient));
        }
        vert = vect;
    }

    friend std::ostream &operator<<(std::ostream &os, const Polygon &n);
};

std::ostream &operator<<(std::ostream &os, const Polygon &n) {
    for (size_t i = 0; i < n.vert.size(); ++i) {
        os << n.vert[i] << ' ';
    }
    return os;
}

class Rectangle : public Polygon {
public:
    Rectangle() {}

    template<typename Type>
    Rectangle(const Point &a, const Point &b, Type coef2) {
        double coef = static_cast<double>(coef2);
        vert.resize(4);
        coef = (coef >= 1) ? coef : 1 / coef;
        double proportion = (coef * coef) / (1 + coef * coef);
        Point n = Line(a, b).normal();
        Point d = b - a;
        vert[0] = a;
        vert[1] = a + proportion * d + (d.length() * coef / (1 + coef * coef)) * n;
        vert[2] = b;
        vert[3] = b - proportion * d - (d.length() * coef / (1 + coef * coef)) * n;
    }

    template<class... Args>
    Rectangle(Args... args): Polygon(args...) {}

    Rectangle(const std::vector <Point> &vert) : Polygon(vert) {
        assert(vert.size() == 4);
    }

    Point center() {
        return ((vert[0] + vert[2]) / 2);
    }

    std::pair <Line, Line> diagonals() const {
        return std::make_pair(Line(vert[0], vert[2]), Line(vert[1], vert[3]));
    };
};

class Square : public Rectangle {
public:
    Square() {}

    Square(const Point &a, const Point &b) : Rectangle(a, b, 1) {}

    const Circle circumscribedCircle() const {
        Point center = (vert[0] + vert[2]) / 2;
        double R = (vert[2] - vert[0]).length() / 2;
        return Circle(center, R);
    }

    const Circle inscribedCircle() const {
        Point center = (vert[0] + vert[2]) / 2;
        double r = (vert[1] - vert[0]).length() / 2;
        return Circle(center, r);
    }
};

class Triangle : public Polygon {
public:
    Triangle() {}

    template<class... Args>
    Triangle(Args... args): Polygon(args...) {}

    Triangle(const std::vector <Point> &vert) : Polygon(vert) {
        assert(vert.size() == 3);
    }

    const Point centroid() const {
        return (vert[0] / 3 + vert[1] / 3 + vert[2] / 3);
    }

    const Circle circumscribedCircle() const {
        const Point A = vert[0], B = vert[1], C = vert[2];
        double a2 = (B - C).length2(), b2 = (C - A).length2(), c2 = (A - B).length2();
        double xa = a2 * (b2 + c2 - a2), xb = b2 * (a2 + c2 - b2), xc = c2 * (a2 + b2 - c2);
        double X = xa + xb + xc;
        Point center = A * xa / X + B * xb / X + C * xc / X;
        double R = (A - center).length();
        return Circle(center, R);
    }

    const Circle inscribedCircle() const {
        Point A = vert[0], B = vert[1], C = vert[2];
        double a = (B - C).length(), b = (C - A).length(), c = (A - B).length();
        double X = a + b + c;
        Point center = A * a / X + B * b / X + C * c / X;
        double r = 2 * area() / perimeter();
        return Circle(center, r);
    }

    const Point orthocenter() const {
        Point A = vert[0], B = vert[1], C = vert[2];
        double a2 = (B - C).length2(), b2 = (C - A).length2(), c2 = (A - B).length2();
        double xa = (a2 + b2 - c2) * (a2 - b2 + c2), xb = (a2 + b2 - c2) * (-a2 + b2 + c2), xc =
                (a2 - b2 + c2) * (-a2 + b2 + c2);
        double X = xa + xb + xc;
        Point ort = A * xa / X + B * xb / X + C * xc / X;
        return ort;
    }

    const Line EulerLine() const {
        return Line(centroid(), ninePointsCircle().center());
    }

    const Circle ninePointsCircle() const {
        Point center = (orthocenter() + circumscribedCircle().center()) / 2.;
        double R = circumscribedCircle().radius() / 2.;
        return Circle(center, R);
    }
};
