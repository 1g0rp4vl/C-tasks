#pragma once
#include <iostream>
#include <compare>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <string>
#include <random>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <cassert>
#include <string>

using std::vector;
using std::complex;
using std::cout;
using std::endl;
using std::cin;
using std::string;

class BigInteger;

class BigInteger {
public:
    constexpr static const double PI = 3.141592653589793238462643383279;
    static const size_t baseSim = 100000000;
    static const size_t baseMul = 100;
    static const int standardBase = 10;
    constexpr static const double half = 0.5;
private:
    size_t base = baseSim;
    bool is_negative = 0;
    vector<int> digits;
    int sign(int a) const {
        if (a < 0) return -1;
        if (a > 0) return 1;
        return 0;
    }
    void change_sign() {
        is_negative ^= true;
        is_negative = is_negative && !is_null();
    }
    bool is_null() const {
        return digits.size() == 1 && digits[0] == 0;
    }
    size_t get_base() const {
        return base;
    }
    static size_t exp(size_t number, size_t degree) {
        if (degree == 0) return 1;
        return number * exp(number, degree - 1);
    }
    static size_t reverse_bits(size_t number, size_t cnt_bits) {
        size_t new_number = 0;
        for (size_t i = 0; i < cnt_bits; i++, (new_number <<= 1)) {
            new_number += ((number & (1 << i)) > 0 ? 1 : 0);
        }
        new_number >>= 1;
        return new_number;
    }
    static size_t log2_ll(size_t n) {
        size_t cur_cnt = 0, cur_num = 1;
        while (cur_num < n) {
            cur_cnt++;
            cur_num *= 2;
        }
        return cur_cnt;
    } 
    static size_t log10(size_t number) {
        size_t cur_pow = 0, cur_num = 1;
        while (cur_num < number) {
            cur_pow++;
            cur_num *= standardBase;
        }
        return cur_pow;
    }
    static void fft(vector<complex<double>>& polynom, complex<double> fft_num) {
        size_t n = polynom.size();
        if (n == 1) return;
        for (size_t i = 0; i < n; i++) {
            if (reverse_bits(i, log2_ll(n)) < i) {
                std::swap(polynom[reverse_bits(i, log2_ll(n))], polynom[i]);
            }
        }
        for (size_t cur_len = 2; cur_len <= n; cur_len *= 2) {
            complex<double> current_fft_num = fft_num;
            for (size_t i = cur_len; i < n; i *= 2) {
                current_fft_num *= current_fft_num;
            }
            for (size_t i = 0; i < n; i += cur_len) {
                complex<double> fft_num_deg = complex<double>(1, 0);
                for (size_t j = i; j < i + cur_len / 2; j++) {
                    complex<double> u = polynom[j];
                    complex<double> v = polynom[j + cur_len / 2] * fft_num_deg;
                    polynom[j] = u + v;
                    polynom[j + cur_len / 2] = u - v;
                    fft_num_deg *= current_fft_num;
                }
            }
        }
    }
    static void invfft(vector<complex<double>>& polynom, complex<double> fft_num) {
        fft(polynom, conj(fft_num));
        for (size_t i = 0; i < polynom.size(); i++) {
            polynom[i] /= static_cast<double>(polynom.size());
        }
        return;
    }
    static void delete_first_nulls(vector<int>& arr) {
        while (arr.size() > 1 && arr.back() == 0) arr.pop_back();
    }
    static void change_base_vectors(size_t old_base, size_t new_base, vector<int>& digits) {
        size_t log10_new_base = log10(new_base);
        vector<uint8_t> to_base10;
        size_t log10_Base = log10(old_base);
        for (size_t i = 0; i < digits.size(); i++) {
            size_t counter = 0;
            while (digits[i] > 0) {
                counter++;
                to_base10.push_back(static_cast<uint8_t>(digits[i] % standardBase));
                digits[i] /= standardBase;
            }
            while (counter < log10_Base) {
                to_base10.push_back(0);
                counter++;
            }
        }
        digits.resize((to_base10.size() + log10_new_base - 1) / log10_new_base);
        for (size_t i = 0; i < to_base10.size(); i += log10_new_base) {
            digits[i / log10_new_base] = 0;
            for (size_t j = std::min(i + log10_new_base, to_base10.size()); (j--) > i;) {
                digits[i / log10_new_base] *= standardBase;
                digits[i / log10_new_base] += to_base10[j];
            }
        }
        delete_first_nulls(digits);
    }
    static vector<int> multiply_polynoms(vector<int> first_p, vector<int> second_p) {
        change_base_vectors(baseSim, baseMul, first_p);
        change_base_vectors(baseSim, baseMul, second_p);
        size_t cur_num = first_p.size() + second_p.size();
        size_t n = 1;
        while (n < cur_num) {
            n *= 2;
        }
        vector<complex<double>> compl_first(n), compl_second(n);
        std::copy(first_p.begin(), first_p.end(), compl_first.begin());
        std::copy(second_p.begin(), second_p.end(), compl_second.begin());
        complex<double> fft_num(cos(2 * PI / static_cast<double>(n)), sin(2 * PI / static_cast<double>(n)));
        fft(compl_first, fft_num);
        fft(compl_second, fft_num);
        for (size_t i = 0; i < n; i++) {
            compl_first[i] *= compl_second[i];
        }
        invfft(compl_first, fft_num);
        vector<int> result(n);   
        for (size_t i = 0; i < n; i++) {
            result[i] = static_cast<int>(floor(real(compl_first[i]) + half));
        }
        return result;
    }

public:
    void swap(BigInteger& other) {
        std::swap(base, other.base);
        digits.swap(other.digits);
        std::swap(is_negative, other.is_negative);
    }
    size_t size() const {
        return digits.size();
    }
    BigInteger(): is_negative(0), digits{0} {}
    BigInteger(const BigInteger& other) = default;
    BigInteger& operator=(const BigInteger& other) = default;
    BigInteger(int number): is_negative(number < 0) {
        if (number == 0) {
            digits.push_back(0);
            return;
        }
        number *= sign(number);
        while (number > 0) {
            digits.push_back(number % static_cast<int>(base));
            number /= static_cast<int>(base);
        }
    }
    BigInteger(const string& str): base(standardBase), is_negative((str.back() == '-') ? true : false) {
        for (size_t i = 0; i < ((str.back() == '-') ? str.size() - 1 : str.size()); i++) {
            digits.push_back(str[i] - '0');
        }
        change_base(baseSim);
    }
    std::strong_ordering operator<=> (const BigInteger& other) const {
        auto sign_ordering = other.is_negative <=> is_negative;
        if (sign_ordering != std::strong_ordering::equal) return sign_ordering;
        auto size_ordering = digits.size() <=> other.digits.size();
        if (size_ordering != std::strong_ordering::equal) return is_negative ? (0 <=> size_ordering) : size_ordering;
        auto res_ordering = std::strong_ordering::equal;
        for (size_t i = digits.size(); (i--) > 0;) {
            if (digits[i] != other.digits[i]) {
                res_ordering = digits[i] <=> other.digits[i];
                break;
            }
        }
        return is_negative ? (0 <=> res_ordering) : res_ordering;
    }

    bool operator== (const BigInteger& other) const = default;

    BigInteger& operator+=(const BigInteger& other) {
        digits.resize(std::max(size(), other.size()));
        int sign_of_mul_signs = ((is_negative == other.is_negative) ? 1 : (-1));
        for (size_t i = 0; i < size(); i++) {
            digits[i] += ((i < other.size()) ? other.digits[i] : 0) * sign_of_mul_signs;
            if (digits[i] == -1ll && i + 1 == size()) {
                break;
            }
            long long to_next = 0;
            if (digits[i] >= static_cast<int>(base)) {
                to_next++;
                digits[i] -= static_cast<int>(base);
            } else if (digits[i] < 0) {
                to_next--;
                digits[i] += static_cast<int>(base);
            }
            if (to_next == 0 && i >= other.size()) break;
            if (to_next != 0) {
                if (i + 1 == size()) {
                    digits.push_back(0);
                }
                digits[i + 1] += static_cast<int>(to_next);
            }
        }
        bool result_negative = 0;
        if (digits.back() == -1) {
            result_negative = 1;
            size_t first_not_null = 0;
            for (size_t i = 0; i < size(); i++) {
                if (digits[i] != 0) {
                    first_not_null = i;
                    break;
                }
            }
            if (first_not_null + 1 != size()) {
                digits[first_not_null] = static_cast<int>(base) - digits[first_not_null];
                for (size_t i = first_not_null + 1; i + 1 < size(); i++) {
                    digits[i] = static_cast<int>(base) - digits[i] - 1;
                }
                digits[size() - 1] = 0;
            } else {
                digits[first_not_null] = 1; 
            }
        }
        delete_first_nulls(digits);
        is_negative = (is_negative != result_negative) && !is_null();
        return *this;
    }

    BigInteger& operator-=(const BigInteger& other) {
        change_sign();
        *this += other;
        change_sign();
        return *this;
    }

    BigInteger& operator++() {
        return *this += 1;
    }

    BigInteger operator++(int) {
        BigInteger copy = *this;
        *this += 1;
        return copy;
    }

    BigInteger& operator--() {
        return *this -= 1;
    }

    BigInteger operator--(int) {
        BigInteger copy = *this;
        *this -= 1;
        return copy;
    }

    BigInteger operator-() const {
        BigInteger copy = *this;
        copy.change_sign();
        return copy;
    }

    BigInteger& operator*=(const BigInteger& other) {
        digits = multiply_polynoms(digits, other.digits);
        base = baseMul;
        for (size_t i = 0; i < digits.size(); i++) {
            if (digits[i] >= static_cast<int>(base)) {
                if (i + 1 == digits.size()) {
                    digits.push_back(0);
                }
                digits[i + 1] += digits[i] / static_cast<int>(base);
            }
            digits[i] %= static_cast<int>(base);
        }
        delete_first_nulls(digits);
        is_negative = (is_negative != other.is_negative) && !is_null();
        change_base(baseSim);
        return *this;
    } 

    BigInteger& operator/= (const BigInteger& other) {
        if (size() < other.size() || is_null()) {
            is_negative = 0;
            digits = {0};
            return *this;
        }
        if (this == &other) {
            is_negative = 0;
            digits = {1};
            return *this;
        }
        vector<int> ans;
        BigInteger copy_this = *this;
        copy_this.is_negative = 0;
        BigInteger copy_other = other;
        for (size_t i = size() - other.size() + 1; (i--) > 0;) {
            int left_border = 0, right_border = static_cast<int>(base);
            while (left_border < right_border - 1) {
                copy_other = other;
                int middle = (left_border + right_border) / 2; 
                copy_other.mul_figure(middle);
                copy_other.mul_base_deg(i);
                copy_other.is_negative = 0;
                if (copy_this >= copy_other) {
                    left_border = middle;
                } else {
                    right_border = middle;
                }
            }
            ans.push_back(left_border);
            copy_other = other;
            copy_other.mul_base_deg(i);
            copy_other.is_negative = 0;
            copy_other.mul_figure(left_border);
            copy_this -= copy_other;
        }
        std::reverse(ans.begin(), ans.end());
        digits = ans;
        delete_first_nulls(digits);
        is_negative = (is_negative != other.is_negative) && !is_null();
        return *this;
    }

    void change_base(size_t new_base) {
        change_base_vectors(base, new_base, digits);
        base = new_base;
    }

    BigInteger& operator%= (const BigInteger& other) {
        BigInteger copy_this = *this;
        copy_this /= other;
        copy_this *= other;
        *this -= copy_this;
        return *this;
    }

    void mul_figure(int figure) {
        unsigned long long carry = 0;
        for (size_t i = 0; i < size(); i++) {
            carry = carry + static_cast<unsigned long long>(digits[i]) * static_cast<unsigned long long>(figure);
            digits[i] = static_cast<int>(carry % base);
            carry /= base;
            if (carry != 0 && i + 1 == size()) {
                digits.push_back(0);
            }
        }
    }

    void mul_base_deg(size_t degree) {
        std::reverse(digits.begin(), digits.end());
        for (size_t i = 0; i < degree; i++) {
            digits.push_back(0);
        }
        std::reverse(digits.begin(), digits.end());
    }

    explicit operator bool() const {
        return !is_null();
    }

    string toString() const {
        string result = "";
        if (is_negative) {
            result += '-';
        }
        size_t log10_base = log10(base);
        for (size_t i = size(); (i--) > 0;) {
            string figure = std::to_string(digits[i]);
            if (i + 1 != size()) {
                std::reverse(figure.begin(), figure.end());
                while (figure.size() < log10_base) {
                    figure += '0';
                }
                std::reverse(figure.begin(), figure.end());
            }
            result += figure;
        }
        return result;
    } 
};

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger copy_first = first;
    return copy_first += second;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    BigInteger copy_first = first;
    return copy_first -= second;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger copy_first = first;
    return copy_first *= second;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger copy_first = first;
    return copy_first /= second;
}

BigInteger operator%(BigInteger first, const BigInteger& second) {
    BigInteger copy_first = first;
    return copy_first %= second;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

std::istream& operator>> (std::istream& in, BigInteger& number) {
    char c = '0';
    bool is_beginning = 1;
    string str = "";
    while (in.get(c)) {
        if (isspace(c)) {
            if (!is_beginning) {
                break;
            }
        } else {
            is_beginning = 0;
            str += c;
        }
    }
    reverse(str.begin(), str.end());
    number = BigInteger(str);
    return in;
}

BigInteger gcd(const BigInteger& first, const BigInteger& second) {
    BigInteger copy_first = first, copy_second = second;
    while (copy_first != 0) {
        copy_second %= copy_first;
        copy_first.swap(copy_second);
    }
    return copy_second;
}

BigInteger lcm(const BigInteger& first, const BigInteger& second) {
    return (first * second) / gcd(first, second);
}

class Rational {
    constexpr static const double standardBaseInv = 0.1;
    static const size_t enoughNumberOfDigits = 20;
    BigInteger numerator, denominator;
    bool is_null() {
        return numerator == 0;
    }
    void fix_fraction() {
        if (denominator < 0) {
            numerator *= -1;
            denominator *= -1;
        }
        bool is_numerator_negative = (numerator < 0);
        numerator *= is_numerator_negative ? (-1) : 1;
        BigInteger grecomdiv = gcd(numerator, denominator);
        numerator /= grecomdiv;
        denominator /= grecomdiv;
        numerator *= is_numerator_negative ? (-1) : 1;
    }

public:
    Rational(): numerator(0), denominator(1) {}
    Rational(const BigInteger& other): numerator(other), denominator(1) {}
    Rational(int other): numerator(other), denominator(1) {}
    Rational(const BigInteger& numerator, const BigInteger& denominator): numerator(numerator), denominator(denominator) {
        fix_fraction();
    }   

    Rational& operator+= (const Rational& other) {
        numerator = numerator * other.denominator + other.numerator * denominator;
        denominator = denominator * other.denominator;
        fix_fraction();
        return *this;
    }

    Rational& operator-= (const Rational& other) {
        numerator = numerator * other.denominator - other.numerator * denominator;
        denominator = denominator * other.denominator;
        fix_fraction();
        return *this;
    }

    Rational& operator*= (const Rational& other) {
        numerator = numerator * other.numerator;
        denominator = denominator * other.denominator;
        fix_fraction();
        return *this;
    }

    Rational& operator/= (const Rational& other) {
        *this = Rational(numerator * other.denominator, denominator * other.numerator);
        return *this;
    }

    Rational operator-() const {
        return Rational(-numerator, denominator);
    }

    std::strong_ordering operator<=>(const Rational& other) const {
        if (numerator == 0 && other.numerator == 0) {
            return std::strong_ordering::equal;
        } else if (numerator < 0 && other.numerator >= 0) {
            return std::strong_ordering::less;
        } else if (numerator >= 0 && other.numerator < 0) {
            return std::strong_ordering::greater;
        }
        Rational substraction = *this;
        substraction -= other;
        if (substraction.numerator < 0) return std::strong_ordering::less;
        if (substraction.numerator > 0) return std::strong_ordering::greater;
        return std::strong_ordering::equal; 
    }

    bool operator==(const Rational& other) const = default;

    string toString() const {
        string result = numerator.toString();
        if (denominator != 1) {
            result += '/';
            result += denominator.toString();
        }
        return result;
    }

    string asDecimal(uint32_t precision = 0) const {
        BigInteger copy_numerator = numerator;
        copy_numerator.change_base(BigInteger::standardBase);
        copy_numerator.mul_base_deg(precision);
        copy_numerator.change_base(BigInteger::baseSim);
        copy_numerator /= denominator;
        if (precision == 0) {
            return copy_numerator.toString();
        }
        string result = "";
        if (copy_numerator < 0) {
            result += '-';
            copy_numerator *= -1;
        }
        string copy_numerator_str = copy_numerator.toString();
        if (copy_numerator_str.size() <= precision) {
            result += "0.";
            for (size_t i = 0; i < precision - copy_numerator_str.size(); i++) {
                result += '0';
            }
            result += copy_numerator_str;
        } else {
            result += copy_numerator_str.substr(0, copy_numerator_str.size() - precision);
            result += '.';
            result += copy_numerator_str.substr(copy_numerator_str.size() - precision, precision);
        }
        return result;
    }

    explicit operator double() const {
        string str = asDecimal(enoughNumberOfDigits);
        char* str_chararray = new char[str.size()];
        std::copy(str.begin(), str.end(), str_chararray);
        char* end = nullptr;
        double result = strtod(str_chararray, &end);
        delete[] str_chararray;
        return result;
    }
};

Rational operator+(const Rational& first, const Rational& second) {
    Rational copy_first = first;
    return copy_first += second;
}

Rational operator-(const Rational& first, const Rational& second) {
    Rational copy_first = first;
    return copy_first -= second;
}

Rational operator*(const Rational& first, const Rational& second) {
    Rational copy_first = first;
    return copy_first *= second;
}

Rational operator/(const Rational& first, const Rational& second) {
    Rational copy_first = first;
    return copy_first /= second;
}

BigInteger operator""_bi(const char* arr_of_chars, size_t) {
    string str = arr_of_chars;
    reverse(str.begin(), str.end());
    return BigInteger(str);
}

BigInteger operator""_bi(unsigned long long number) {
    string str = std::to_string(number);
    reverse(str.begin(), str.end());
    return BigInteger(str);
}


//aboba