#include <iostream>
#include <cstdint>

template <size_t N>
class StackStorage {
    char memory[N];
    size_t first_free = 0;
public:
    StackStorage() {}
    StackStorage(const StackStorage&) = delete;
    char* allocate_memory(size_t n, size_t allignment) {
        size_t memory_as_int = reinterpret_cast<std::uintptr_t>(memory);
        size_t first_suitable = ((memory_as_int + first_free) % allignment == 0) ? (memory_as_int + first_free) : (memory_as_int + first_free + allignment - (memory_as_int + first_free) % allignment);
        first_suitable -= memory_as_int;
        if (N >= first_suitable + n) {
            first_free = first_suitable + n;
            return memory + first_suitable;
        }
        throw std::bad_alloc();
    }
};

template <typename T, size_t N> 
class StackAllocator {
    template <typename U, size_t M>
    friend class StackAllocator;
    StackStorage<N>* memory;
public:
    using value_type = T;
    using propagate_on_contatiner_copy_assignment = std::true_type;
    StackAllocator() = delete;
    StackAllocator(StackStorage<N>& memory): memory(&memory) {}
    template <typename U>
    StackAllocator(const StackAllocator<U, N>& other): memory(other.memory) {}
    template <typename U>
    StackAllocator& operator=(const StackAllocator<U, N>& other) {
        memory = other.memory;
    }
    T* allocate(size_t n) const {
        return reinterpret_cast<T*>(memory->allocate_memory(n * sizeof(T), sizeof(T)));
    } 
    void deallocate(__attribute__((unused)) T* ptr, __attribute__((unused)) size_t n) const {}
    template <typename U>
    struct rebind {
        using other = StackAllocator<U, N>;
    };
};

template <typename T, typename Alloc = std::allocator<T>>
class List: Alloc  {
public:
    struct BaseNode {
        BaseNode* prev = nullptr;
        BaseNode* next = nullptr;
        BaseNode() {
            prev = this;
            next = this;
        }
    };
    struct Node: BaseNode {
        T el;
        Node() = default;
        Node(const T& other): el(other) {}
    };
    using AllocTraits = std::allocator_traits<Alloc>;
    using NodeAlloc = typename AllocTraits::template rebind_alloc<Node>;
    using BaseNodeAlloc = typename AllocTraits::template rebind_alloc<BaseNode>;
    using NodeAllocTraits = typename std::allocator_traits<NodeAlloc>;
    using BaseNodeAllocTraits = typename std::allocator_traits<BaseNodeAlloc>;
    template <typename Value>
    class base_iterator {
        BaseNode* elem = nullptr;
    public:
        using difference_type = ptrdiff_t;
        using value_type = Value;
        using reference = value_type&;
        using iterator_category = std::bidirectional_iterator_tag;
        using pointer = Value*;
        base_iterator() = default;
        base_iterator(BaseNode* elem): elem(elem) {}
        reference operator*() const {
            return static_cast<Node*>(elem)->el;
        }
        value_type* operator->() const {
            return static_cast<Node*>(elem);
        }
        base_iterator& operator++() {
            elem = elem->next;
            return *this;
        }
        base_iterator operator++(int) {
            auto copy = *this;
            ++*this;
            return copy;   
        }
        base_iterator& operator--() {
            elem = elem->prev;
            return *this;
        }
        base_iterator operator--(int) {
            auto copy = *this;
            --*this;
            return copy;   
        }
        operator base_iterator<const Value>() const {
            base_iterator<const Value> res(elem);
            return res;
        }
        bool operator==(const base_iterator& other) const {
            return elem == other.elem;
        }
        Node* get_ptr() const {
            return static_cast<Node*>(elem);
        }
        BaseNode* get_base_ptr() const {
            return elem;
        }
    };
public:
    using iterator = base_iterator<T>;
    using const_iterator = base_iterator<const T>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    iterator begin() {
        return ++end();
    }
    const_iterator begin() const {
        return const_iterator(endNode.next);
    }
    const_iterator cbegin() const {
        return const_iterator(endNode.next);
    }
    iterator end() {
        return iterator(&endNode);
    }
    const_iterator end() const {
        return const_iterator(endNode.next->prev);
    }
    const_iterator cend() const {
        return const_iterator(endNode.next->prev);
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
    size_t sz;
    BaseNode endNode;
    void clear(BaseNode* fakeNode, NodeAlloc allocator, Node* last = nullptr, bool last_allocated = false) {
        if (last_allocated) {
            NodeAllocTraits::deallocate(allocator, last, 1);
        }
        BaseNode* ptr = fakeNode->next;
        while (ptr != fakeNode) {
            BaseNode* cptr = ptr;
            ptr = ptr->next;
            NodeAllocTraits::destroy(allocator, static_cast<Node*>(cptr));
            NodeAllocTraits::deallocate(allocator, static_cast<Node*>(cptr), 1);
        }
        fakeNode->next = fakeNode;
        fakeNode->prev = fakeNode;
    }
    void pop(bool at_the_end) {
        Node* ptr = (at_the_end) ? (--end()).get_ptr() : begin().get_ptr();
        if (at_the_end) {
            pop_from_the_end(&endNode);
        } else {
            pop_from_the_front(&endNode);
        }
        NodeAlloc newNodeAlloc = get_allocator(); 
        NodeAllocTraits::destroy(newNodeAlloc, ptr);
        NodeAllocTraits::deallocate(newNodeAlloc, ptr, 1);
        --sz;
    }
    void pop_from_the_front(BaseNode* fakeNode) {
        fakeNode->next->next->prev = fakeNode;
        fakeNode->next = fakeNode->next->next;
    }
    void pop_from_the_end(BaseNode* fakeNode) {
        fakeNode->prev->prev->next = fakeNode;
        fakeNode->prev = fakeNode->prev->prev;
    }
    void push_at_the_end(BaseNode* fakeNode, Node* newNode) {
        newNode->next = fakeNode;
        newNode->prev = fakeNode->prev;
        fakeNode->prev->next = newNode;
        fakeNode->prev = newNode;
    }
    void push_at_the_front(BaseNode* fakeNode, Node* newNode) {
        newNode->next = fakeNode->next;
        newNode->prev = fakeNode;
        fakeNode->next->prev = newNode;
        fakeNode->next = newNode;
    }
    Node* create_node(const T& el) {
        Node* newNode;
        bool allocated = false;
        NodeAlloc newNodeAlloc = get_allocator(); 
        try {
            newNode = NodeAllocTraits::allocate(newNodeAlloc, 1);
            allocated = true;
            NodeAllocTraits::construct(newNodeAlloc, newNode, el);
        } catch(...) {
            if (allocated) {
                NodeAllocTraits::deallocate(newNodeAlloc, newNode, 1);
            }
            throw;
        }
        return newNode;
    }
    void destroy_node(Node* node) {
        NodeAlloc newNodeAlloc = get_allocator(); 
        NodeAllocTraits::destroy(newNodeAlloc, node);
        NodeAllocTraits::deallocate(newNodeAlloc, node, 1);
    }
    void push(bool at_the_end, const T& el) {
        Node* newNode = create_node(el);
        if (at_the_end) {
            push_at_the_end(&endNode, newNode);
        } else {
            push_at_the_front(&endNode, newNode);
        }
        ++sz;
    }
    enum class type_of_construction {
        default_c, 
        copy_one_elem,
        copy_all_elems
    };
    template <type_of_construction type>
    void construct(BaseNode* fakeNode, NodeAlloc allocator, size_t n, const T* el = nullptr, const_iterator it_other = const_iterator(nullptr)) {
        if (n == 0) return;
        bool last_allocated = false;
        Node* ptr = nullptr;
        try {
            for (size_t i = 0; i < n; ++i) {
                last_allocated = false;
                ptr = NodeAllocTraits::allocate(allocator, 1);
                last_allocated = true;
                if constexpr (type == type_of_construction::default_c) {
                    NodeAllocTraits::construct(allocator, ptr);
                } else if constexpr (type == type_of_construction::copy_one_elem) {
                    NodeAllocTraits::construct(allocator, ptr, *el);
                } else if constexpr (type == type_of_construction::copy_all_elems) {
                    NodeAllocTraits::construct(allocator, ptr, *it_other);
                    ++it_other;
                }
                push_at_the_end(fakeNode, ptr);
            }
        } catch(...) {
            clear(fakeNode, allocator, ptr, last_allocated);
            throw;
        }
    }
public:
    NodeAlloc get_allocator() const {
        return static_cast<Alloc>(*this);
    }
    List(): sz(0) {}
    List(size_t n): sz(n) {
        construct<type_of_construction::default_c>(&endNode, get_allocator(), n);
    }
    List(size_t n, const T& el): sz(n) {
        construct<type_of_construction::copy_one_elem>(&endNode, get_allocator(), el);
    }
    List(const Alloc& allocator): Alloc(allocator), sz(0) {}
    List(size_t n, const Alloc& allocator): Alloc(allocator), sz(n) {
        construct<type_of_construction::default_c>(&endNode, get_allocator(), n);
    }
    List(size_t n, const T& el, const Alloc& allocator): Alloc(allocator), sz(n) {
        construct<type_of_construction::copy_one_elem>(&endNode, get_allocator(), n, el);
    }
    List(const List& other): Alloc(AllocTraits::select_on_container_copy_construction(other.get_allocator())), sz(other.sz) {
        construct<type_of_construction::copy_all_elems>(&endNode, get_allocator(), sz, nullptr, other.begin());
    }
    List& operator=(const List& other) {
        Alloc newAlloc = (std::is_base_of<std::true_type, typename AllocTraits::propagate_on_container_copy_assignment>::value) ? other.get_allocator() : get_allocator();
        BaseNode newEndNode;
        construct<type_of_construction::copy_all_elems>(&newEndNode, newAlloc, other.sz, nullptr, other.begin());
        clear(&endNode, get_allocator());
        newEndNode.prev->next = &endNode;
        newEndNode.next->prev = &endNode;
        endNode = newEndNode;
        if (std::is_base_of<std::true_type, typename AllocTraits::propagate_on_container_copy_assignment>::value) {
            static_cast<Alloc&>(*this) = newAlloc;
        }
        sz = other.sz;
        return *this;
    }
    void push_back(const T& el) {
        push(true, el);
    }
    void push_front(const T& el) {
        push(false, el);
    }
    void pop_back() {
        pop(true);
    }
    void pop_front() {
        pop(false);
    }
    void insert(const_iterator it, const T& el) {
        Node* newNode = create_node(el);
        push_at_the_end(it.get_base_ptr(), newNode);
        ++sz;
    }
    void erase(const_iterator it) {
        pop_from_the_end(it.get_base_ptr()->next);
        destroy_node(it.get_ptr());
        --sz;
    } 
    size_t size() const {
        return sz;
    }
    ~List() {
        clear(&endNode, get_allocator());
    }
};