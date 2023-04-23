//исправил, то, что говорили, не понял насчет оператора +
#include <iostream>
#include <algorithm>
#include <cstring>

using std::copy;

class String {
    char* array_ = nullptr;
    size_t size_ = 0;
    size_t capacity_ = 0;

public:
    void set_end(size_t position) {
        array_[position] = '\0';
    }

    String(): array_(new char[1]) {
        set_end(0);
    }

    String(const String& another): array_(new char[another.capacity_ + 1]), size_(another.size_), capacity_(another.capacity_) {
        memmove(array_, another.array_, size_ + 1);
    } 

    String(const char* c_style_string): array_(new char[strlen(c_style_string) + 1]), size_(strlen(c_style_string)), capacity_(strlen(c_style_string)) {
        memmove(array_, c_style_string, size_);
        set_end(size_);
    }

    String(size_t n, char c): array_(new char[n + 1]), size_(n), capacity_(n) {
        memset(array_, c, n);
        set_end(size_);
    }

    const char& operator[](size_t index) const {
        return array_[index];
    }

    char& operator[](size_t index) {
        return array_[index];
    }

    size_t size() const {
        return size_;
    }

    size_t capacity() const {
        return capacity_;
    }

    size_t length() const {
        return size_;
    }

    void reallocate(size_t new_capacity) {
        char* new_array = new char[new_capacity + 1];
        memmove(new_array, array_, size_ + 1);
        capacity_ = new_capacity;
        delete[] array_;
        array_ = new_array;
    }


    void push_back(char c) {
        if (capacity_ == size_) {
            reallocate(2 * (size_ + 1));
        }
        array_[size_++] = c;
        set_end(size_);
    }

    void pop_back() {
        set_end(--size_);
    }

    const char& front() const {
        return array_[0];
    }

    const char& back() const {
        return array_[size_ - 1];
    }

    char& front() {
        return array_[0];
    }

    char& back() {
        return array_[size_ - 1];
    }

    const char* data() const {
        return array_;
    } 

    char* data() {
        return array_;
    }

    String& operator=(const String& another) {
        if (size_ > another.size_) {
            memmove(array_, another.array_, another.size_ + 1);
            size_ = another.size_;
        } else {
            String copy = another;
            swap(copy);
            //delete[] copy.array_;
        }
        return *this;
    }

    void swap(String& another) {
        std::swap(array_, another.array_);
        std::swap(size_, another.size_);
        std::swap(capacity_, another.capacity_);
    }

    String& operator+=(char c) {
        push_back(c);
        return *this;
    }

    String& operator+=(const String& another) {
        if (capacity_ - size_ < another.size_) {
            reallocate(size_ + another.size_);
        }
        memmove(array_ + size_, another.array_, another.size_ + 1);
        size_ += another.size_;
        return *this;
    }

    size_t search(const String& substring, bool reversed) const {
        for (size_t i = reversed ? (size_ - substring.size_ + 1) : 0; reversed ? ((i--) > 0) : (i + substring.size_ <= size_); reversed ? i : ++i) {
            bool is_substring = true;
            for (size_t j = i; j < i + substring.size_; j++) {
                if (array_[j] != substring[j - i]) {
                    is_substring = false;
                    break;
                }
            }
            if (is_substring) return i;
        }
        return length();
    }

    size_t find(const String& substring) const {
        return search(substring, false);
    }

    size_t rfind(const String& substring) const {
        return search(substring, true);
    }

    String substr(size_t start, size_t count) const {
        String substring = String(count, '\0');
        memmove(substring.array_, array_ + start, count);
        return substring;
    }

    bool empty() const {
        return size_ == 0;
    }

    void clear() {
        set_end(size_ = 0);
    }

    void shrink_to_fit() {
        if (capacity_ > size_) {
            reallocate(size_);
        }
    }

    ~String() {
        delete[] array_;
    }
};

bool operator==(const String& first, const String& second) {
    if (first.size() != second.size()) {
        return false;
    }
    for (size_t i = 0; i < first.size(); i++) {
        if (first[i] != second[i]) {
            return false;
        }
    }
    return true;
}

bool operator!=(const String& first, const String& second) {
    return !(first == second);
}

bool operator<(const String& first, const String& second) {
    for (size_t i = 0; i < std::min(first.size(), second.size()); i++) {
        if (first[i] < second[i]) {
            return true;
        } else if (first[i] > second[i]) {
            return false;
        }
    }
    return first.size() < second.size();
}

bool operator>(const String& first, const String& second) {
    return second < first;
}

bool operator>= (const String& first, const String& second) {
    return !(first < second);
}

bool operator<= (const String& first, const String& second) {
    return !(first > second);
}

String operator+(String first, const String& second) {
    first += second;
    return first;
}

String operator+(String string, char c) {
    string += c;
    return string;
}

String operator+(char c, const String& string) {
    return String(1, c) + string;
}

std::istream& operator>>(std::istream& in, String& string) {
    string.clear();
    char c;
    bool is_begining = 1;
    while (in.get(c)) {
        //std::cout << c << " " << is_begining << std::endl;
        if (isspace(c)) {
            if (!is_begining) {
                break;
            }
        } else {
            is_begining = false;
            string += c;
        }
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const String& string) {
    for (size_t i = 0; i < string.size(); i++) {
        out << string[i];
    }
    return out;
}
