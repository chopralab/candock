#ifndef IT_H
#define IT_H
#include <functional>
#include <map>
#include <memory>
#include <vector>
#include "candock/helper/debug.hpp"
#include "candock/helper/error.hpp"

namespace candock {
namespace molib {

template <class T, class U, class P>
class template_vector_iterator {  // Now define it
    T& __mols;
    U index;

   public:
    template_vector_iterator(T& mols) : __mols(mols), index(mols.__m.begin()) {}
    template_vector_iterator(T& mols, bool)
        : __mols(mols), index(mols.__m.end()) {}  // the end sentinel
    template_vector_iterator(T& mols, bool, bool)
        : __mols(mols), index(mols.__m.rbegin()) {}  // the rbegin sentinel
    template_vector_iterator(const template_vector_iterator& rv)
        : __mols(rv.__mols), index(rv.index) {}  // the copy constructor
    template_vector_iterator operator+(int offset) {
        template_vector_iterator other(*this);
        other.index += offset;
        return other;
    }
    template_vector_iterator& operator++() {
        ++index;
        return *this;
    }
    template_vector_iterator& operator++(int) { return operator++(); }
    template_vector_iterator& operator--() {
        --index;
        return *this;
    }
    template_vector_iterator& operator--(int) { return operator--(); }
    P& current() const { return **index; }
    P& operator*() const { return current(); }
    P* operator->() const { return &current(); }
    bool operator==(const template_vector_iterator& rv) const {
        return index == rv.index;
    }  // comparison tests
    bool operator!=(const template_vector_iterator& rv) const {
        return index != rv.index;
    }  // comparison tests
};

template <class T, class P>
class template_vector_container {
    typedef std::vector<T> __vector;
    typedef typename __vector::iterator __iterator;
    typedef typename __vector::const_iterator __const_iterator;
    typedef typename __vector::const_reverse_iterator __const_reverse_iterator;
    __vector __m;

   public:
    P& add(T t) {
        __m.push_back(std::move(t));
        return *__m.back();
    }
    void clear() { __m.clear(); }
    friend class template_vector_iterator<template_vector_container<T, P>,
                                          __iterator, P>;
    typedef template_vector_iterator<template_vector_container<T, P>,
                                     __iterator, P>
        iterator;
    friend class template_vector_iterator<const template_vector_container<T, P>,
                                          __const_iterator, P>;
    typedef template_vector_iterator<const template_vector_container<T, P>,
                                     __const_iterator, P>
        const_iterator;
    friend class template_vector_iterator<const template_vector_container<T, P>,
                                          __const_reverse_iterator, P>;
    typedef template_vector_iterator<const template_vector_container<T, P>,
                                     __const_reverse_iterator, P>
        const_reverse_iterator;
    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(*this, true); }
    const_iterator begin() const { return const_iterator(*this); }
    const_iterator end() const { return const_iterator(*this, true); }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(*this, true, true);
    }
    P& first() const { return begin().current(); }
    P& last() const {
        return rbegin().current();
    }  // expressions such as iterator it = end() uses copy constructor and not
       // assignment operator!!!
    P& operator[](int i) const { return *__m[i]; }
    P& element(int i) const { return operator[](i); }
    bool empty() const { return (begin() == end()); }
    size_t size() const { return __m.size(); }
    //~ void erase(int i) { __m.erase(i); }
    void erase(int i) { __m.erase(__m.begin() + i); }
    void remove_if(std::function<bool(const P&)> fptr) {
        for (__iterator it = __m.begin(), ite = __m.end(); it != ite;) {
            if (fptr(*it))
                it = __m.erase(it);
            else
                ++it;
        }
    }
};

template <class T, class U, class P>
class template_map_iterator {  // Now define it
    T& __mols;
    U index;

   public:
    template_map_iterator(T& mols) : __mols(mols), index(mols.__m.begin()) {}
    template_map_iterator(T& mols, bool)
        : __mols(mols), index(mols.__m.end()) {}  // the end sentinel
    template_map_iterator(T& mols, bool, bool)
        : __mols(mols), index(mols.__m.rbegin()) {}  // the rbegin sentinel
    template_map_iterator(const template_map_iterator& rv)
        : __mols(rv.__mols), index(rv.index) {}  // the copy constructor
    //~ template_map_iterator& operator=(const template_map_iterator &rv) { cout
    //<< "in operator="<<endl; __mols = &rv.__mols; index = rv.index; return
    //*this; }  // strange std::move ...
    template_map_iterator& operator++() {
        ++index;
        return *this;
    }
    template_map_iterator& operator++(int) { return operator++(); }
    template_map_iterator& operator--() {
        --index;
        return *this;
    }
    template_map_iterator& operator--(int) { return operator--(); }
    P& current() const { return *(index->second); }
    P& operator*() const { return current(); }
    P* operator->() const { return &current(); }
    bool operator==(const template_map_iterator& rv) const {
        return index == rv.index;
    }  // comparison tests
    bool operator!=(const template_map_iterator& rv) const {
        return index != rv.index;
    }  // comparison tests
};

template <class T, class U, class V, class Z = int>
class template_map_container {
    typedef std::map<Z, std::unique_ptr<T> > __map;
    typedef typename __map::iterator __iterator;
    typedef typename __map::const_iterator __const_iterator;
    typedef typename __map::const_reverse_iterator __const_reverse_iterator;
    __map __m;
    V* __br;

   public:
    T& aadd(T* t, U* u) {
        t->set_br(u);
        __m.insert(make_pair(__m.size(), std::unique_ptr<T>(t)));
        return *__m[__m.size() - 1];
    }
    T& aadd(Z p, T* t, U* u) {
        t->set_br(u);
        __m.insert(make_pair(p, std::unique_ptr<T>(t)));
        return *__m[p];
    }
    //~ T& aadd_no_br(Z p, T *t) { __m.insert(make_pair(p, unique_ptr<T>(t)));
    //return *__m[p]; }
    void clear() { __m.clear(); }
    friend class template_map_iterator<template_map_container<T, U, V, Z>,
                                       __iterator, T>;
    typedef template_map_iterator<template_map_container<T, U, V, Z>,
                                  __iterator, T>
        iterator;
    friend class template_map_iterator<const template_map_container<T, U, V, Z>,
                                       __const_iterator, T>;
    typedef template_map_iterator<const template_map_container<T, U, V, Z>,
                                  __const_iterator, T>
        const_iterator;
    friend class template_map_iterator<const template_map_container<T, U, V, Z>,
                                       __const_reverse_iterator, T>;
    typedef template_map_iterator<const template_map_container<T, U, V, Z>,
                                  __const_reverse_iterator, T>
        const_reverse_iterator;
    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(*this, true); }
    const_iterator begin() const { return const_iterator(*this); }
    const_iterator end() const { return const_iterator(*this, true); }
    const_reverse_iterator rbegin() const {
        return const_reverse_iterator(*this, true, true);
    }
    const V& br() const { return *__br; }
    void set_br(V* br) { __br = br; }
    T& first() const { return begin().current(); }
    T& last() const {
        return rbegin().current();
    }  // expressions such as iterator it = end() uses copy constructor and not
       // assignment operator!!!
    T& operator[](Z p) const {
        __const_iterator i = __m.find(p);
        if (i == __m.end()) {
            throw Error("die : wrong map index...");
        }
        return *(i->second);
    }
    T& element(Z p) const { return operator[](p); }
    bool has_element(Z p) const { return (__m.find(p) != __m.end()); }
    bool empty() const { return (begin() == end()); }
    size_t size() const { return __m.size(); }
    void erase(Z p) { __m.erase(p); }
    void erase_shrink(Z p) {
        __m[p].swap(__m.rbegin()->second);
        __m.erase(__m.rbegin()->first);
    }
    void remove_if(std::function<bool(const T&)> fptr) {
        for (__iterator it = __m.begin(), ite = __m.end(); it != ite;) {
            if (fptr(*(it->second)))
                it = __m.erase(it);
            else
                ++it;
        }
    }
};
}
}

#endif
