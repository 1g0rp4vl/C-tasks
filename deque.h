    #include <iostream>

    template<typename T>
    class Deque {
        static const ssize_t size_of_block = 32; 

        struct Position {
            ssize_t ind1 = 0;
            ssize_t ind2 = 0;
            Position() {}
            Position(ssize_t ind1, ssize_t ind2): ind1(ind1), ind2(ind2) {}
            Position(ssize_t pos): ind1((pos < 0) ? ((pos - size_of_block + 1) / size_of_block) : (pos / size_of_block)), ind2(pos - size_of_block * ind1) {}
            Position& operator+=(ssize_t add_number) {
                ssize_t new_pos = ind1 * size_of_block + ind2 + add_number;
                ind1 = (new_pos < 0) ? ((new_pos - size_of_block + 1) / size_of_block) : (new_pos / size_of_block);
                ind2 = new_pos - size_of_block * ind1;
                return *this;
            }
            Position& operator-=(ssize_t add_number) {
                return (*this) += -add_number;
            }
            Position& operator++() {
                return (*this) += 1;
            }
            Position& operator--() {
                return (*this) -= 1;
            }
            Position operator++(int) {
                auto copy = *this;
                (*this) += 1;
                return copy;
            }
            Position operator--(int) {
                auto copy = *this;
                (*this) -= 1;
                return copy;
            }
            ssize_t full_position() const {
                return ind1 * size_of_block + ind2;
            }
            auto operator<=>(const Position& other) const {
                return std::make_pair(ind1, ind2) <=> std::make_pair(other.ind1, other.ind2);
            }
            bool operator==(const Position& other) const {
                return std::make_pair(ind1, ind2) == std::make_pair(other.ind1, other.ind2);
            }
            Position operator+(ssize_t number) const {
                Position copy = *this;
                copy += number;
                return copy;
            }
            ssize_t operator-(const Position& other) const {
                ssize_t full_pos_this = ind1 * size_of_block + ind2;
                ssize_t full_pos_other = other.ind1 * size_of_block + other.ind2;
                return full_pos_this - full_pos_other;
            }
        };

        template <typename Value>
        class base_iterator {
            template <typename Value2>
            friend class Deque<T>::base_iterator;
            T** block = nullptr;
            Position pos = Position(0, 0);
            Value* elem = nullptr;
        public:
            using difference_type = ptrdiff_t;
            using value_type = Value;
            using reference = value_type&;
            using iterator_category = std::random_access_iterator_tag;
            using pointer = Value*;
            base_iterator() = default;
            base_iterator(T** begin_block, Position _pos) {
                pos = _pos;
                block = begin_block ? (begin_block + pos.ind1) : nullptr;
                elem = block ? (*(block) + pos.ind2) : nullptr;
            }
            reference operator*() const {
                return *elem;
            }
            value_type* operator->() const {
                return elem;
            }
            base_iterator& operator++() {
                if (pos.ind2 == size_of_block - 1) {
                    elem = block ? *(++block) : nullptr;
                } else {
                    ++elem;
                }
                ++pos;
                return *this;
            }
            base_iterator operator++(int) {
                auto copy = *this;
                ++*this;
                return copy;   
            }
            base_iterator& operator--() {
                if (pos.ind2 != 0) {
                    --elem;
                } else if (pos.ind1 != 0) {
                    elem = *(--block) + size_of_block - 1;
                } else {
                    block ? --block : nullptr;
                    elem = nullptr;
                }
                --pos;
                return *this;
            }
            base_iterator operator--(int) {
                auto copy = *this;
                --*this;
                return copy;   
            }
            base_iterator& operator+=(difference_type diff) {
                Position new_pos = pos + diff;
                block += block ? (new_pos.ind1 - pos.ind1) : 0;
                if (new_pos.ind1 < 0) {
                    elem = nullptr;
                } else {
                    elem = block ? ((*block) + new_pos.ind2) : nullptr;
                }
                pos = new_pos;
                return *this;
            }
            base_iterator& operator-=(difference_type diff) {
                return *this += (-diff);
            }
            base_iterator operator+(difference_type diff) const {
                auto copy = *this;
                copy += diff;
                return copy;
            }
            base_iterator operator-(difference_type diff) const {
                return (*this) + (-diff);
            }
            auto operator<=>(const base_iterator& other) const {
                return pos <=> other.pos;
            }
            bool operator==(const base_iterator& other) const {
                return pos == other.pos;
            }
            value_type* get_ptr() const {
                return elem;
            }
            Position get_position() const {
                return pos;
            }
            operator base_iterator<const Value>() const {
                base_iterator<const Value> res;
                res.block = block;
                res.pos = pos;
                res.elem = elem;
                return res;
            }
            ptrdiff_t operator-(const base_iterator& other) const {
                return pos - other.pos;
            }
        };
    public:
        using iterator = base_iterator<T>;
        using const_iterator = base_iterator<const T>;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;
        iterator begin() {
            return iterator(arr_of_blocks, begin_pos);
        }
        const_iterator begin() const {
            return const_iterator(arr_of_blocks, begin_pos);
        }
        const_iterator cbegin() const {
            return const_iterator(arr_of_blocks, begin_pos);
        }
        iterator end() {
            return iterator(arr_of_blocks, begin_pos + sz);
        }
        const_iterator end() const {
            return const_iterator(arr_of_blocks, begin_pos + sz);
        }
        const_iterator cend() const {
            return const_iterator(arr_of_blocks, begin_pos + sz);
        }
        reverse_iterator rbegin() {
            return reverse_iterator(end());
        }
        const_reverse_iterator rbegin() const {
            return const_reverse_iterator(end());
        }
        const_reverse_iterator crbegin() const {
            return const_reverse_iterator(cend());
        }
        reverse_iterator rend() {
            return reverse_iterator(begin());
        }
        const_reverse_iterator rend() const {
            return const_reverse_iterator(begin());
        }
        const_reverse_iterator crend() const {
            return const_reverse_iterator(cbegin());
        }
        size_t sz_arr_of_blocks = 0;
        size_t begin_block = 0;
        size_t end_block = 0;
        size_t sz = 0;
        T** arr_of_blocks = nullptr;
        Position begin_pos = Position(0, 0);
        enum class type_of_contruction {
            default_c,
            copy_one_elem,
            copy_all_elems
        };
        template <type_of_contruction type>
        void constructor(T** arr, Position beg_pos, Position end_pos, size_t beg_bl, size_t end_bl, const T* ptr = nullptr, const_iterator it_other = iterator(nullptr, Position(0, 0))) {
            if (end_bl - beg_bl == 0) return;
            size_t i = beg_bl;
            iterator it;
            try {
                for (; i < end_bl; i++) {
                    arr[i] = reinterpret_cast<T*>(new char[sizeof(T) * size_of_block]);
                }
                for (it = iterator(arr, beg_pos); it != iterator(arr, end_pos); ++it) {
                    if constexpr (type == type_of_contruction::default_c) {
                        new (it.get_ptr()) T();
                    } else if constexpr (type == type_of_contruction::copy_one_elem) {
                    // std::cout << it_beg.pos.ind1 << " " << it_beg.pos.ind2 << std::endl;
                        new (it.get_ptr()) T(*ptr);
                    } else {
                        new (it.get_ptr()) T(*(it_other++));
                    }
                }
            } catch(...) {
                clear_memory(arr, iterator(arr, beg_pos), it, beg_bl, i);
                throw;
            }
        }
        void clear_memory(T** arr, iterator it_beg, iterator it_end, size_t beg_bl, size_t ind) {
            for (; it_beg != it_end; ++it_beg) {
                (it_beg.get_ptr())->~T();
            }
            for (size_t j = beg_bl; j < ind; j++) {
                delete[] reinterpret_cast<char*>(arr[j]);
            }
            delete[] arr;
        }
        size_t capacity_beg() const {
            return begin_pos - Position(0, 0);
        }
        size_t capacity_end() const {
            if (sz_arr_of_blocks == 0) return 0;
            return Position(sz_arr_of_blocks - 1, size_of_block - 1) - (begin_pos + sz);
        }
    public:
        Deque() {}
        Deque(int sz): sz_arr_of_blocks(sz == 0 ? 0 : (sz + size_of_block) / size_of_block), begin_block(0), end_block(begin_block + sz_arr_of_blocks), sz(sz), arr_of_blocks(sz == 0 ? nullptr : (new T*[sz_arr_of_blocks])) {
            constructor<type_of_contruction::default_c>(arr_of_blocks, begin_pos, begin_pos + sz, begin_block, end_block);
        }
        Deque(int sz, const T& el): sz_arr_of_blocks(sz == 0 ? 0 : (sz + size_of_block) / size_of_block), begin_block(0), end_block(begin_block + sz_arr_of_blocks), sz(sz), arr_of_blocks(sz == 0 ? nullptr : (new T*[sz_arr_of_blocks])) { 
            //std::cout << (arr_of_blocks == nullptr) << std::endl;
            constructor<type_of_contruction::copy_one_elem>(arr_of_blocks, begin_pos, begin_pos + sz, begin_block, end_block, &el);
        }
        Deque(const Deque& other): sz_arr_of_blocks(other.sz_arr_of_blocks), begin_block(other.begin_block), end_block(other.end_block), sz(other.sz), arr_of_blocks(sz_arr_of_blocks == 0 ? nullptr : new T*[sz_arr_of_blocks]), begin_pos(other.begin_pos) {
            if (sz_arr_of_blocks == 0) return;
            constructor<type_of_contruction::copy_all_elems>(arr_of_blocks, begin_pos, begin_pos + sz, begin_block, end_block, nullptr, other.begin());
        }
        Deque<T>& operator=(const Deque& other) {
            T** new_arr_of_blocks = (other.sz_arr_of_blocks == 0) ? nullptr : new T*[other.sz_arr_of_blocks];
            constructor<type_of_contruction::copy_all_elems>(new_arr_of_blocks, other.begin_pos, other.begin_pos + other.sz, other.begin_block, other.end_block, nullptr, other.begin());
            clear_memory(arr_of_blocks, begin(), end(), begin_block, end_block);
            std::swap(arr_of_blocks, new_arr_of_blocks);
            sz = other.sz;
            sz_arr_of_blocks = other.sz_arr_of_blocks;
            begin_block = other.begin_block;
            end_block = other.end_block;
            begin_pos = other.begin_pos;
            return *this;
        }
        size_t size() const {
            return sz;
        }
        T& operator[](size_t ind) {
            return *(begin() + ind);
        }
        const T& operator[](size_t ind) const {
            return *(begin() + ind);
        }
        T& at(size_t ind) {
            if (begin_pos + ind < begin_pos || begin_pos + ind >= begin_pos + sz) throw std::out_of_range("");
            return *(begin() + ind);
        }
        const T& at(size_t ind) const {
            if (begin_pos + ind < begin_pos || begin_pos + ind >= begin_pos + sz) throw std::out_of_range("");
            return *(begin() + ind);
        }
        void clear_new_arr_of_blocks(T** arr, bool at_the_end, bool created_new_block = false) {
            if (created_new_block) {
                delete[] reinterpret_cast<char*>(arr[at_the_end ? (sz_arr_of_blocks + 1 + end_block) : (sz_arr_of_blocks + begin_block)]);
            }
            delete[] arr;
        }
        T** add_blocks(size_t new_sz_arr_of_blocks, bool at_the_end) {
            if (new_sz_arr_of_blocks == 0) return nullptr;
            T** new_arr_of_blocks = new T*[new_sz_arr_of_blocks];
            try {
                for (size_t i = sz_arr_of_blocks + 1 + begin_block; i < sz_arr_of_blocks + 1 + end_block; i++) {
                    new_arr_of_blocks[i] = arr_of_blocks[i - sz_arr_of_blocks - 1];
                }
                new_arr_of_blocks[at_the_end ? (sz_arr_of_blocks + 1 + end_block) : (sz_arr_of_blocks + begin_block)] = reinterpret_cast<T*>(new char*[size_of_block * sizeof(T)]);
            } catch(...) {
                clear_new_arr_of_blocks(new_arr_of_blocks, at_the_end, false);
                throw;
            }
            return new_arr_of_blocks;
        }
        void add_element_with_new_memory(bool at_the_end, const T& el) {
            size_t new_sz = sz + 1;
            size_t new_sz_arr_of_blocks = (sz_arr_of_blocks + 1) * 3;
            size_t new_begin_block = begin_block + sz_arr_of_blocks + (at_the_end ? 1 : 0);
            size_t new_end_block = end_block + sz_arr_of_blocks + (at_the_end ? 2 : 1);
            Position new_begin_pos = begin_pos + (sz_arr_of_blocks + 1) * size_of_block;
            T** new_arr_of_blocks = add_blocks(new_sz_arr_of_blocks, at_the_end);
            try {
                new (iterator(new_arr_of_blocks, at_the_end ? (new_begin_pos + sz) : (new_begin_pos - 1)).get_ptr()) T(el);
            } catch(...) {
                clear_new_arr_of_blocks(new_arr_of_blocks, at_the_end, true);
                throw;
            }
            delete[] arr_of_blocks;
            std::swap(new_sz, sz);
            std::swap(new_sz_arr_of_blocks, sz_arr_of_blocks);
            begin_pos = at_the_end ? new_begin_pos : (new_begin_pos - 1);
            std::swap(begin_block, new_begin_block);
            std::swap(end_block, new_end_block);
            std::swap(new_arr_of_blocks, arr_of_blocks);
        }
        void push_back(const T& el) {
            if (capacity_end() == 0) {
                add_element_with_new_memory(true, el);
                return;
            }
            bool created_block = false; 
            if (begin_pos + sz == Position(end_block - 1, size_of_block - 1)) {
                arr_of_blocks[end_block++] = reinterpret_cast<T*>(new char[size_of_block * sizeof(T)]);
                created_block = true;
            }
            try {
                new (end().get_ptr()) T(el);
                ++sz;
            } catch(...) {
                if (created_block) {
                    delete[] reinterpret_cast<char*>(arr_of_blocks[--end_block]);
                }
                throw;
            }
        }
        void pop_back() {
            (--end())->~T();
            --sz;
        }
        void push_front(const T& el) {
            if (capacity_beg() == 0) {
                add_element_with_new_memory(false, el);
                return;
            }
            bool created_block = false; 
            if (begin_pos == Position(begin_block, 0)) {
                arr_of_blocks[--begin_block] = reinterpret_cast<T*>(new char[size_of_block * sizeof(T)]);
                created_block = true;
            }
            try {
                new ((--begin()).get_ptr()) T(el);
                ++sz;
                --begin_pos;
            } catch(...) {
                if (created_block) {
                    delete[] reinterpret_cast<char*>(arr_of_blocks[begin_block++]);
                }
                throw;
            }
        }
        void pop_front() {
            (begin())->~T();
            --sz;
            ++begin_pos;
        }
        void move_elements(bool in_front, iterator it_beg, iterator it_end) {
            if (in_front) {
                if (it_end != it_beg) {
                    new (it_end.get_ptr()) T(*(it_end - 1));
                }
                --it_end;
                for (; it_end > it_beg; --it_end) {
                    *(it_end) = *(it_end - 1);
                }
            } else {
                for (; it_beg != it_end; ++it_beg) {
                    *(it_beg - 1) = *(it_beg);
                }
            }
        }
        void insert(iterator it, const T& el) {
            if (capacity_end() == 0) {
                ssize_t it_pos_add = it.get_position() - begin_pos;
                begin_pos = begin_pos + (sz_arr_of_blocks + 1) * size_of_block;
                T** new_arr_of_blocks = add_blocks((sz_arr_of_blocks + 1) * 3, true);
                begin_block = sz_arr_of_blocks + 1 + begin_block;
                end_block = sz_arr_of_blocks + 2 + end_block;
                sz_arr_of_blocks = (sz_arr_of_blocks + 1) * 3;
                delete[] arr_of_blocks;
                arr_of_blocks = new_arr_of_blocks;
                it = begin() + it_pos_add;
            } else if (begin_pos + sz == Position(end_block - 1, size_of_block - 1)) {
                arr_of_blocks[end_block++] = reinterpret_cast<T*>(new char[size_of_block * sizeof(T)]);
            }
            move_elements(true, it, end());
            if (it == end()) {
                new (it.get_ptr()) T(el);
            } else {
                *(it) = el;
            }
            ++sz;
        }
        void erase(iterator it) {
            move_elements(false, it + 1, end());
            (--end())->~T();
            --sz;
        }
        ~Deque() {
            clear_memory(arr_of_blocks, begin(), end(), begin_block, end_block);
        }
    };
