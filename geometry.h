#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <algorithm>

using std::vector;
using std::pair;
using std::cout;
using std::endl;

class Line;
struct Point {
    double x = 0;
    double y = 0;
    Point() {}
    Point(double x, double y);
    bool operator==(const Point& other) const;
    Point& operator+=(const Point& other);
    Point& operator*=(double k);
    Point& operator-=(const Point& other);
    bool operator<(const Point& other) const;
    void rotate(const Point& center, double angle);
    void reflect(const Point& center);
    void reflect(const Line& axis);
    void scale(const Point& center, double k);
    Point get_perpendicular() const;
    double get_length() const;
    Point rotary_homothety(double angle, double k) const;
    static double oriented_angle_between_vectors(const Point& vec1, const Point& vec2);
    static double angle_between_vectors(const Point& vec1, const Point& vec2);
    static double scalar_product(const Point& p1, const Point& p2);
    static double skew_product(const Point& p1, const Point& p2);
    static bool is_on_one_line(const Point& p1, const Point& p2, const Point& p3);
    static bool is_point_on_ray(const Point& beginning_of_ray, const Point& point_on_ray, const Point& point);
    static bool is_point_on_segment(const Point& first, const Point& second, const Point& point);
    static bool is_vectors_parallel(const Point& vec1, const Point& vec2);
};
class Line {
    double A = 0;
    double B = 0;
    double C = 0;
public:
    Line() {}
    Line(double A, double B, double C);
    Line(Point p1, Point p2);
    Line(double k, double b);
    Line(Point p, double k);
    Point direction() const;
    Point get_first_point() const;
    Point get_second_point() const;
    double value_at_point(const Point& point) const;
    double sqAB() const;
    bool operator==(const Line& other) const;
    bool is_point_on_line(const Point& point) const;
    static Point intersect_lines(const Line& l1, const Line& l2);
    static Point projection_on_line(const Line& line, const Point& point);
    static Line middle_perpendicular(const Point& p1, const Point& p2);
    static Line bisector_of_angle(const Point& p1, const Point& p2, const Point& p3);
    static bool are_lines_parallel(const Line& line1, const Line& line2);
    static Line line_by_angle(const Point& point, double angle);
    double distance_point_line(const Point& point) const;
};
class Polygon;
class Ellipse;
Point operator-(const Point& p1, const Point& p2);
Point operator*(const Point& p1, double k);
Point operator+(const Point& p1, const Point& p2);
bool operator==(const Polygon& first, const Polygon& second);
bool operator==(const Ellipse& first, const Ellipse& second);

namespace Constants {
    const double EPS = 1e-10;
    const double PI = atan2(0, -1);
    const double DEGS_EQ_2PI = 360;
    const double DEGS_EQ_PI = 180;
    const double HALF = 0.5;
    const double RECTANGLE_NUMBER_OF_VERTEX = 4;
}

namespace GeometryFunctions {
    bool equal(double a, double b);
    double radians_to_degs(double angle_in_radians);
    double degs_to_radians(double angle_in_radians);
}

Point::Point(double x, double y): x(x), y(y) {}

bool Point::operator==(const Point& other) const {
    return GeometryFunctions::equal(x, other.x) && GeometryFunctions::equal(y, other.y);
}
Point& Point::operator+=(const Point& other) {
    x += other.x;
    y += other.y;
    return (*this);
}
Point& Point::operator*=(double k) {
    x *= k;
    y *= k;
    return (*this);
}
Point& Point::operator-=(const Point& other) {
    return (*this) += other * -1;
}
bool Point::operator<(const Point& other) const {
    return std::make_pair(x, y) < std::make_pair(other.x, other.y);
}
void Point::rotate(const Point& center, double angle) {
    angle = GeometryFunctions::degs_to_radians(angle);
    (*this) -= center;
    double rotated_x = cos(angle) * x - sin(angle) * y;
    double rotated_y = sin(angle) * x + cos(angle) * y;
    x = rotated_x, y = rotated_y;
    (*this) += center;
}
void Point::reflect(const Point& center) {
    (*this) = center * 2 - (*this);
}
void Point::reflect(const Line& axis) {
    reflect(Line::projection_on_line(axis, (*this)));
} 
void Point::scale(const Point& center, double k) {
    (*this) = ((*this) - center) * k + center;
}
Point Point::get_perpendicular() const {
    return Point(-y, x);    
}
double Point::get_length() const {
    return hypot(x, y);
}
Point Point::rotary_homothety(double angle, double k) const {
    Point copy_point = *this;
    copy_point.rotate(Point(0, 0), angle);
    return copy_point *= k;
}
double Point::oriented_angle_between_vectors(const Point& vec1, const Point& vec2) {
    return GeometryFunctions::radians_to_degs(atan2(Point::skew_product(vec1, vec2), Point::scalar_product(vec1, vec2)));    
}
double Point::angle_between_vectors(const Point& vec1, const Point& vec2) {
    return fabs(oriented_angle_between_vectors(vec1, vec2));    
}
double Point::scalar_product(const Point& p1, const Point& p2) {
    return p1.x * p2.x +  p1.y * p2.y;
}
double Point::skew_product(const Point& p1, const Point& p2) {
    return p1.x * p2.y -  p1.y * p2.x;
}
bool Point::is_on_one_line(const Point& p1, const Point& p2, const Point& p3) {
    return GeometryFunctions::equal(skew_product(p1 - p2, p1 - p3), 0); 
}
bool Point::is_point_on_ray(const Point& beginning_of_ray, const Point& point_on_ray, const Point& point) {
    return GeometryFunctions::equal(skew_product(point - beginning_of_ray, point_on_ray - beginning_of_ray), 0) && (GeometryFunctions::equal(scalar_product(point - beginning_of_ray, point_on_ray - beginning_of_ray), 0) || scalar_product(point - beginning_of_ray, point_on_ray - beginning_of_ray) > 0);
}
bool Point::is_point_on_segment(const Point& first, const Point& second, const Point& point) {
    return is_point_on_ray(first, second, point) && is_point_on_ray(second, first, point);
}
bool Point::is_vectors_parallel(const Point& vec1, const Point& vec2) {
    return GeometryFunctions::equal(skew_product(vec1, vec2), 0);
}

Point operator+(const Point& p1, const Point& p2) {
    return Point(p1.x + p2.x, p1.y + p2.y);
}

Point operator-(const Point& p1, const Point& p2) {
    return Point(p1.x - p2.x, p1.y - p2.y);
}

Point operator*(const Point& p1, double k) {
    return Point(p1.x * k, p1.y * k);
}

Point operator/(const Point& point, double coefficient) {
    return Point(point.x / coefficient, point.y / coefficient);
}

Line::Line(double A, double B, double C): A(A), B(B), C(C) {}
Line::Line(Point p1, Point p2): A(p2.y - p1.y), B(p1.x - p2.x), C(p1.y * p2.x - p2.y * p1.x) {}
Line::Line(double k, double b): Line(Point(0, b), Point(1, k + b)) {}
Line::Line(Point p, double k): Line(p, Point(p.x + 1, p.y + k)) {}
Point Line::direction() const {
    return Point(-B, A);
}
Point Line::get_first_point() const {
    if (!GeometryFunctions::equal(A, 0)) return Point(-C / A, 0);
    return Point(0, -C / B);
}
Point Line::get_second_point() const {
    return get_first_point() + direction();
}
double Line::value_at_point(const Point& point) const {
    return A * point.x + B * point.y + C;
}
double Line::sqAB() const {
    return sqrt(A * A + B * B);
}
bool Line::operator==(const Line& other) const {
    return Point::is_on_one_line(get_first_point(), get_second_point(), other.get_first_point()) && Point::is_on_one_line(get_first_point(), get_second_point(), other.get_second_point());
}
bool Line::is_point_on_line(const Point& point) const {
    return GeometryFunctions::equal(value_at_point(point), 0);
}
Point Line::intersect_lines(const Line& l1, const Line& l2) {
    double det = l1.A * l2.B - l1.B * l2.A;
    double det_x = (-l1.C) * l2.B - l1.B * (-l2.C);
    double det_y = -(-l1.C) * l2.A + l1.A * (-l2.C);
    return Point(det_x / det, det_y / det);
}
Point Line::projection_on_line(const Line& line, const Point& point) {
    return intersect_lines(line, Line(point, point + (line.get_first_point() - line.get_second_point()).get_perpendicular()));
}
Line Line::middle_perpendicular(const Point& p1, const Point& p2) {
    return Line((p1 + p2) * Constants::HALF, (p1 + p2) * Constants::HALF + (p1 - p2).get_perpendicular());
}
Line Line::bisector_of_angle(const Point& p1, const Point& p2, const Point& p3) {
    double angle = Point::angle_between_vectors(p1 - p2, p3 - p2) / 2;
    angle *= (Point::skew_product(p1 - p2, p3 - p2) < 0) ? (-1) : 1;
    return Line(p2, p2 + (p1 - p2).rotary_homothety(angle, 1));
}
bool Line::are_lines_parallel(const Line& line1, const Line& line2) {
    return Point::is_vectors_parallel(line1.direction(), line2.direction());
}
Line Line::line_by_angle(const Point& point, double angle) {
    Point dir(1, 0);
    dir.rotate(Point(0, 0), angle);
    return Line(point, point + dir);
}

double Line::distance_point_line(const Point& point) const {
    return fabs(value_at_point(point) / sqAB());
} 

class Shape {
public:
    virtual double perimeter() const = 0; 
    virtual double area() const = 0;
    virtual bool isCongruentTo(const Shape& other) const = 0;
    virtual bool operator==(const Shape& other) const = 0;
    virtual bool isSimilarTo(const Shape& other) const = 0;
    virtual bool containsPoint(const Point& point) const = 0;
    virtual void rotate(const Point& point, double angle) = 0;
    virtual void reflect(const Point& point) = 0;
    virtual void reflect(const Line& line) = 0;
    virtual void scale(const Point& point, double coefficient) = 0;

    virtual ~Shape() = default;
};


class Polygon: public Shape {
    friend bool operator==(const Polygon& first, const Polygon& second);
    void push_back_points(const Point& p) {
        points.push_back(p);
    }

    template<class... Points> 
    void push_back_points(const Point& p, Points... args) {
        points.push_back(p);
        push_back_points(args...);
    }
protected:
    std::vector<Point> points;
public:
    Polygon() {}
    Polygon(const vector<Point>& points): points(points) {}
    Polygon(const std::initializer_list<Point>& init_points): points(init_points) {}

    template<class... Points>
    Polygon(Points... args) {
        push_back_points(args...);
    }

    size_t verticesCount() const {
        return points.size();
    }
    const vector<Point>& getVertices() const {
        return points;
    }

    bool isConvex() const {
        if (verticesCount() <= 2) return false; 
        bool is_turn_left = 0;
        for (size_t pos = 0; pos < verticesCount(); pos++) {
            size_t pos1 = (pos + 1) % verticesCount();
            size_t pos2 = (pos + 2) % verticesCount();
            Point first_edge = points[pos1] - points[pos], second_edge = points[pos2] - points[pos1];
            if (pos != 0) {
                if ((Point::skew_product(first_edge, second_edge) < 0) == is_turn_left) {
                    return false;
                }
            }
            is_turn_left = !(Point::skew_product(first_edge, second_edge) < 0);
        }
        return true;
    }
    double perimeter() const override {
        double ans = 0;
        for (size_t i = 0; i < verticesCount(); i++) {
            ans += (points[i] - points[(i + 1) % verticesCount()]).get_length();
        }
        return ans;
    }
    double area() const override {
        double ans = 0;
        for (size_t i = 0; i < verticesCount(); i++) {
            ans += Point::skew_product(points[i], points[(i + 1) % verticesCount()]);
        }
        return fabs(ans) / 2;
    }
    static bool is_arrays_equal(const vector<Point>& first, const vector<Point>& second) {
        for (size_t i = 0; i < first.size(); i++) {
            size_t v1 = 0, v2 = i;
            bool is_match = true;
            for (size_t j = 0; j < first.size(); j++, v1++, (++v2) %= first.size()) {
                if (first[v1] != second[v2]) {
                    is_match = false;
                    break;
                }
            }
            if (is_match) {
                return true;
            }
        }
        return false;
    }
    static bool is_arrays_congruent_similar(const vector<Point>& first, const vector<Point>& second, bool is_congruent) {
        for (size_t i = 0; i < first.size(); i++) {
            size_t v1 = 0, v2 = i;
            bool is_match = true;
            bool is_ratio_negative = 0;
            for (size_t j = 0; j < first.size(); j++, v1++, (++v2) %= first.size()) {
                Point first_edge1 = first[v1] - first[(v1 + 1) % first.size()];
                Point first_edge2 = first[(v1 + 1) % first.size()] - first[(v1 + 2) % first.size()];
                Point second_edge1 = second[v2] - second[(v2 + 1) % first.size()];
                Point second_edge2 = second[(v2 + 1) % first.size()] - second[(v2 + 2) % first.size()];
                double skew_product_fe1_fe2 = Point::skew_product(first_edge1, first_edge2);
                double skew_product_se1_se2 = Point::skew_product(second_edge1, second_edge2);
                double fe1_len = first_edge1.get_length();
                double fe2_len = first_edge2.get_length();
                double se1_len = second_edge1.get_length();
                double se2_len = second_edge2.get_length();
                double angle_between_fe1_fe2 = Point::angle_between_vectors(first_edge1, first_edge2);
                double angle_between_se1_se2 = Point::angle_between_vectors(second_edge1, second_edge2);      
                double oriented_angle_between_fe1_fe2 = Point::oriented_angle_between_vectors(first_edge1, first_edge2);
                double oriented_angle_between_se1_se2 = Point::oriented_angle_between_vectors(second_edge1, second_edge2);               
                if (is_congruent) {
                    if (!GeometryFunctions::equal(se1_len, fe1_len)) {
                        is_match = false;
                        break;
                    }
                    if (!GeometryFunctions::equal(fabs(skew_product_fe1_fe2), fabs(skew_product_se1_se2)) || (j != 0 && (skew_product_fe1_fe2 * skew_product_se1_se2 < 0) != is_ratio_negative)) {
                        is_match = false;
                        break;
                    }
                    is_ratio_negative = (skew_product_fe1_fe2 * skew_product_se1_se2 < 0);
                } else {
                    if (!GeometryFunctions::equal(fe1_len * se2_len, fe2_len * se1_len)) {
                        is_match = false;
                        break;
                    }
                    if (!GeometryFunctions::equal(angle_between_fe1_fe2, angle_between_se1_se2) || (j != 0 && (oriented_angle_between_fe1_fe2 * oriented_angle_between_se1_se2 < 0) != is_ratio_negative)) {
                        is_match = false;
                        break;
                    }
                    is_ratio_negative = oriented_angle_between_fe1_fe2 * oriented_angle_between_se1_se2 < 0;
                }
            }
            if (is_match) {
                return true;
            }
        }
        return false;
    }
    bool isCongruentTo(const Shape& other) const override {
        const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other);
        if (!other_polygon) return false;
        if (other_polygon->verticesCount() != verticesCount()) {
            return false;
        }
        if (is_arrays_congruent_similar(points, other_polygon->points, true)) return true;
        vector<Point> copy_points = points;
        reverse(copy_points.begin(), copy_points.end());
        return is_arrays_congruent_similar(copy_points, other_polygon->points, true);
    }
    bool operator==(const Shape& other) const override {
        const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other);
        if (!other_polygon) return false;
        return (*this) == (*other_polygon);
    }
    bool isSimilarTo(const Shape& other) const override {
        const Polygon* other_polygon = dynamic_cast<const Polygon*>(&other);
        if (!other_polygon) return false;
        if (other_polygon->verticesCount() != verticesCount()) {
            return false;
        }
        if (is_arrays_congruent_similar(points, other_polygon->points, false)) return true;
        vector<Point> copy_points = points;
        std::reverse(copy_points.begin(), copy_points.end());
        return is_arrays_congruent_similar(copy_points, other_polygon->points, false);
    }

    Point get_edge(size_t vertex_number) const {
        return points[(vertex_number + 1) % verticesCount()] - points[vertex_number];
    }

    bool is_point_on_edge(const Point& point) const {
        for (size_t i = 0; i < verticesCount(); i++) {
            if (Point::is_point_on_segment(points[i], points[(i + 1) % verticesCount()], point)) {
                return true;
            }
        }
        return false;
    }

    bool containsPoint(const Point& point) const override {
        if (is_point_on_edge(point)) return true;
        Line suitable_line;
        for (size_t number_of_parts = 2; ; number_of_parts++) {
            bool found_suitable = true;
            for (size_t i = 0; i < number_of_parts; i++) {
                found_suitable = true;
                suitable_line = Line::line_by_angle(point, static_cast<double>(i) * Constants::DEGS_EQ_2PI / static_cast<double>(number_of_parts));
                for (size_t j = 0; j < verticesCount(); j++) {
                    if (suitable_line.is_point_on_line(points[j])) {
                        found_suitable = false;
                        break;
                    }
                }
                if (found_suitable) break;
            }
            if (found_suitable) break;
        }
        size_t number_of_intersections = 0;
        for (size_t i = 0; i < verticesCount(); i++) {
            Line edge_line = Line(points[i], points[(i + 1) % points.size()]);
            if (!Line::are_lines_parallel(edge_line, suitable_line) && Point::is_point_on_segment(points[i], points[(i + 1) % points.size()], Line::intersect_lines(edge_line, suitable_line)) && Point::is_point_on_ray(point, point + suitable_line.direction(), Line::intersect_lines(edge_line, suitable_line))) {
                number_of_intersections++;
            }
        }
        return number_of_intersections % 2 == 1;
    }

    void rotate(const Point& point, double angle) override {
        for (size_t i = 0; i < verticesCount(); i++) {
            points[i].rotate(point, angle);
        }
    }
    void reflect(const Point& point) override {
        for (size_t i = 0; i < verticesCount(); i++) {
            points[i].reflect(point);
        }
    }
    void reflect(const Line& line) override {
        for (size_t i = 0; i < verticesCount(); i++) {
            points[i].reflect(line);
        }
    }
    void scale(const Point& point, double coefficient) override {
        for (size_t i = 0; i < verticesCount(); i++) {
            points[i].scale(point, coefficient);
        }
    }
    void print() const {
        for (auto i : points) {
            cout << i.x << " " << i.y << endl;
        }
    }
};

bool operator==(const Polygon& first, const Polygon& second) {
    if (first.verticesCount() != second.verticesCount()) {
        return false;
    }
    if (Polygon::is_arrays_equal(first.points, second.points)) return true;
    vector<Point> copy_first_points = first.points;
    reverse(copy_first_points.begin(), copy_first_points.end());
    if (Polygon::is_arrays_equal(copy_first_points, second.points)) return true;
    return false;
}

class Ellipse: public Shape {
    friend bool operator==(const Ellipse& first, const Ellipse& second);
protected:
    double large_semi_axis = 0;
    Point focus1, focus2;
    double get_focal_length() const {
        return (focus1 - focus2).get_length() / 2;
    }
    double get_small_semi_axis() const {
        return sqrt(large_semi_axis * large_semi_axis - get_focal_length() * get_focal_length()); 
    }
public:
    Ellipse() {}
    Ellipse(const Point& focus1, const Point& focus2, double summ_of_lens): large_semi_axis(summ_of_lens / 2), focus1(focus1), focus2(focus2) {}

    pair<Point, Point> focuses() const {
        return std::make_pair(focus1, focus2);
    }
    Point center() const {
        return (focus1 + focus2) * Constants::HALF;
    }
    double eccentricity() const {
        return get_focal_length() / large_semi_axis;
    }
    pair<Line, Line> directrices() {
        Point cf1 = focus1 - center();
        Point cf2 = focus2 - center();
        cf1 *= large_semi_axis / (eccentricity() * get_focal_length());
        cf2 *= large_semi_axis / (eccentricity() * get_focal_length());
        Line dir1 = Line(center() + cf1, center() + cf1 + cf1.get_perpendicular());
        Line dir2 = Line(center() + cf2, center() + cf2 + cf2.get_perpendicular());
        return std::make_pair(dir1, dir2);
    }
    double perimeter() const override {
        return Constants::RECTANGLE_NUMBER_OF_VERTEX * large_semi_axis * std::comp_ellint_2(eccentricity()) ;
    }
    double area() const override {
        return Constants::PI * large_semi_axis * get_small_semi_axis();
    }
    bool isCongruentTo(const Shape& other) const override {
        const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other);
        if (!other_ellipse) return false;
        return GeometryFunctions::equal(large_semi_axis, other_ellipse->large_semi_axis) && GeometryFunctions::equal(get_small_semi_axis(), other_ellipse->get_small_semi_axis());
    }
    bool operator==(const Shape& other) const override {
        const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other);
        if (!other_ellipse) return false;
        return (*this) == (*other_ellipse);
    }
    bool isSimilarTo(const Shape& other) const override {
        const Ellipse* other_ellipse = dynamic_cast<const Ellipse*>(&other);
        if (!other_ellipse) return false;
        return GeometryFunctions::equal(large_semi_axis * other_ellipse->get_small_semi_axis(), get_small_semi_axis() * other_ellipse->large_semi_axis); 
    }
    bool containsPoint(const Point& point) const override {
        return ((focus1 - point).get_length() + (focus2 - point).get_length()) < 2 * large_semi_axis || GeometryFunctions::equal(((focus1 - point).get_length() + (focus2 - point).get_length()), 2 * large_semi_axis);
    }
    void rotate(const Point& point, double angle) override {
        focus1.rotate(point, angle);
        focus2.rotate(point, angle);
    }
    void reflect(const Point& point) override {
        focus1.reflect(point);
        focus2.reflect(point);
    }
    void reflect(const Line& line) override {
        focus1.reflect(line);
        focus2.reflect(line);
    }
    void scale(const Point& point, double coefficient) override {
        focus1.scale(point, coefficient);
        focus2.scale(point, coefficient);
        large_semi_axis *= coefficient;
    }
};

bool operator==(const Polygon&, const Ellipse&) {
    return false;
}

class Circle: public Ellipse {
public:
    Circle(const Point& center, double radius): Ellipse(center, center, 2 * radius) {}
    double radius() const {
        return large_semi_axis;
    }
    double perimeter() const override {
        return 2 * Constants::PI * radius();
    }
    double area() const override {
        return Constants::PI * radius() * radius();
    }
    bool containsPoint(const Point& point) const override {
        return (center() - point).get_length() < radius() || GeometryFunctions::equal((center() - point).get_length(), radius());
    }
};

bool operator==(const Ellipse& first, const Ellipse& second) {
    return GeometryFunctions::equal(first.large_semi_axis, second.large_semi_axis) && ((first.focus1 == second.focus1 && first.focus2 == second.focus2) || (first.focus2 == second.focus1 && first.focus1 == second.focus2));
}

class Rectangle: public Polygon {
public:
    Rectangle() {}
    Rectangle(const Point& p1, const Point& p2, double k) {
        k = (k < 1) ? (1 / k) : k;
        Point p3 = p1 + (p2 - p1).rotary_homothety(GeometryFunctions::radians_to_degs(-atan2(1, k)), k / sqrt(1 + k * k));
        Point p4 = p1 + (p2 - p1).rotary_homothety(GeometryFunctions::radians_to_degs(atan2(k, 1)), 1 / sqrt(1 + k * k));
        points = {p1, p3, p2, p4};
    }
    double diagonal() const {
        return (points[2] - points[0]).get_length();
    }
    double long_side() const {
        return std::max((points[1] - points[0]).get_length(), (points[2] - points[1]).get_length());
    }
    double short_side() const {
        return std::min((points[1] - points[0]).get_length(), (points[2] - points[1]).get_length());
    }
    Point center() const {
        return (points[0] + points[2]) * Constants::HALF;
    }
    pair<Line, Line> diagonals() const {
//NOLINTNEXTLINE
        return std::make_pair(Line(points[0], points[1]), Line(points[1], points[3]));
    }
    double perimeter() const override {
        return 2 * short_side() + 2 * long_side();
    }
    double area() const override {
        return short_side() * long_side();
    }
};

class Square: public Rectangle {
public:
    Square() {}
    Square(const Point& p1, const Point& p2): Rectangle(p1, p2, 1) {}
    double side() const {
        return short_side();
    }
    Circle circumscribedCircle() const {
        return Circle(center(), diagonal() / 2);
    }
    Circle inscribedCircle() const {
        return Circle(center(), side() / 2); 
    }
    double perimeter() const override {
        return Constants::RECTANGLE_NUMBER_OF_VERTEX * side();
    }
    double area() const override {
        return side() * side();
    }
};

class Triangle: public Polygon {
public:
    Triangle(const Point& p1, const Point& p2, const Point& p3): Polygon(p1, p2, p3) {}
    Circle circumscribedCircle() const {
        Line mp1 = Line::middle_perpendicular(points[0], points[1]);
        Line mp2 = Line::middle_perpendicular(points[0], points[2]);
        Point center = Line::intersect_lines(mp1, mp2);
        return Circle(center, (points[0] - center).get_length());
    }
    Circle inscribedCircle() const {
        Line bis1 = Line::bisector_of_angle(points[0], points[1], points[2]);
        Line bis2 = Line::bisector_of_angle(points[0], points[2], points[1]);
        Point center = Line::intersect_lines(bis1, bis2);
        return Circle(center, Line(points[0], points[1]).distance_point_line(center));
    }
    Point intersection_middle_perpendiculars() const {
        Line mp1 = Line::middle_perpendicular(points[0], points[1]);
        Line mp2 = Line::middle_perpendicular(points[0], points[2]);
        return Line::intersect_lines(mp1, mp2);
    }
    Point centroid() const {
//NOLINTNEXTLINE
        return ((points[2] - points[0]) + (points[1] - points[0])) / 3.0 + points[0];
    }
    Point orthocenter() const {
        return Line::intersect_lines(Line(points[0], points[0] + (points[2] - points[1]).get_perpendicular()), Line(points[1], points[1] + (points[2] - points[0]).get_perpendicular()));
    }
    Line EulerLine() const {
        return Line(intersection_middle_perpendiculars(), orthocenter());
    }
    Circle ninePointsCircle() const {
        return Triangle((points[0] + points[1]) * Constants::HALF, (points[1] + points[2]) * Constants::HALF, (points[0] + points[2]) * Constants::HALF).circumscribedCircle();
    }
    double perimeter() const override {
        return (points[0] - points[1]).get_length() + (points[1] - points[2]).get_length() + (points[2] - points[0]).get_length();
    }
    double area() const override {
        return Line(points[0], points[1]).distance_point_line(points[2]) * (points[1] - points[0]).get_length() / 2; 
    }
};

namespace GeometryFunctions {
    bool equal(double a, double b) {
        return fabs(a - b) < Constants::EPS;
    }
    double radians_to_degs(double angle_in_radians) {
        return (angle_in_radians / Constants::PI) * Constants::DEGS_EQ_PI;
    }
    double degs_to_radians(double angle_in_radians) {
        return (angle_in_radians / Constants::DEGS_EQ_PI) * Constants::PI;
    }
}