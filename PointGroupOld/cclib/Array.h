#ifndef _Array_h_
#define _Array_h_

#include <iostream.h>
#include <assert.h>

template <class Item>
class Array {
 public:

  const int ARRAY_PRINT_NUM_PER_LINE = 8;

  // Default constructor
  Array( int size = 1);

  // Destructor
  ~Array();

  // Copy constructor
  Array(const Array<Item>& right);

  // Assignment operator
  const Array<Item>& operator=( const Array<Item>& right);

  // Overloaded subscript operator for non-const Arrays
  // reference return creates an lvalue
  Item& operator[]( const int subscript ) {
      assert( 0 <= subscript && subscript < Size );
      return Data[ subscript ]; // reference return
  }

  // Overloaded subscript operator for const Arrays
  const Item& operator[]( const int subscript ) const {
      assert( 0 <= subscript && subscript < Size );
      return Data[ subscript ]; // reference return
  }

  // Get the array size.
  int GetSize() const {return Size; }
 private:
  int Size;
  Item *Data;
};

// Overloading the iostream operators for Array
template <class Item>
ostream& operator<<(ostream& s, const Array<Item>& array);

template <class Item>
istream& operator>>(istream& s, Array<Item>& array);

#include "Array.cpp"

#endif /* Array_h */
